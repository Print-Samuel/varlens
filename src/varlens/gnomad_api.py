"""Look up population allele frequency from gnomAD's public GraphQL API.

Use this instead of downloading a full gnomAD VCF when you only need
frequencies for a handful of variants (e.g. one patient's VCF): no download,
no API key. For more than a couple thousand variants, downloading a gnomAD
subset and using ``annotate.load_frequencies`` is faster and puts far less
load on gnomAD's free, shared, public service -- that is exactly why this
module enforces ``max_variants`` and a conservative pace between requests.

Uses only the standard library (``urllib``), deliberately, to match the rest
of the project: no extra dependency just to make one HTTP call.
"""
from __future__ import annotations

import json
import time
import urllib.error
import urllib.request
from typing import Callable, Dict, Iterable, MutableMapping, Optional

from . import __version__
from .vcf import Key

API_URL = "https://gnomad.broadinstitute.org/api"

# A descriptive User-Agent -- not a browser impersonation -- so gnomAD's
# operators can see what's calling; the default "Python-urllib/x.y" agent is
# what many APIs' anti-abuse rules flag first.
_USER_AGENT = f"varlens/{__version__} (educational bioinformatics tool; +https://github.com/)"

# Keep the query minimal: allele count and allele number from exomes and
# genomes are all we need to compute a combined allele frequency.
_QUERY = """
query VarlensLookup($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    exome { ac an }
    genome { ac an }
  }
}
"""

_RETRYABLE_STATUS = {429, 500, 502, 503, 504}


class GnomadAPIError(RuntimeError):
    """Raised when gnomAD's API can't be reached, or refuses the request outright."""


def variant_id(key: Key) -> str:
    """gnomAD's variant-id string for a key: 'chrom-pos-ref-alt' (no 'chr' prefix)."""
    chrom, pos, ref, alt = key
    return f"{chrom}-{pos}-{ref}-{alt}"


def _combined_af(payload: dict) -> Optional[float]:
    """Combine exome + genome allele counts into one frequency.

    gnomAD v4 keeps exome and genome cohorts separate rather than always
    publishing a pre-combined number, so we sum AC and AN ourselves across
    whichever cohorts the variant was actually observed in.
    """
    ac = an = 0
    seen = False
    for level in ("exome", "genome"):
        block = payload.get(level)
        if block and block.get("an"):
            ac += block["ac"]
            an += block["an"]
            seen = True
    if not seen or an == 0:
        return None
    return ac / an


def _post(query: str, variables: dict, timeout: float) -> dict:
    """POST one GraphQL request. Raises urllib.error.HTTPError / URLError as-is."""
    body = json.dumps({"query": query, "variables": variables}).encode("utf-8")
    request = urllib.request.Request(
        API_URL,
        data=body,
        headers={"Content-Type": "application/json", "User-Agent": _USER_AGENT},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def _http_error_detail(exc: urllib.error.HTTPError) -> str:
    """Best-effort read of the response body, which usually has the real reason
    (gnomAD's GraphQL errors normally come with a JSON body explaining what was
    wrong with the request -- 'Bad Request' alone doesn't tell us that)."""
    try:
        body = exc.read().decode("utf-8", errors="replace")
    except Exception:
        return ""
    return body[:300]


def _retry_after_seconds(exc: urllib.error.HTTPError) -> Optional[float]:
    header = exc.headers.get("Retry-After") if exc.headers else None
    if not header:
        return None
    try:
        return float(header)
    except ValueError:
        return None  # Retry-After can also be an HTTP date; not handled here


def _is_rate_limit_message(detail: str) -> bool:
    return "rate limit" in detail.lower()


def query_frequencies(
    keys: Iterable[Key],
    dataset: str = "gnomad_r4",
    timeout: float = 10.0,
    max_retries: int = 4,
    delay: float = 2.0,  # starting pace between *successful* requests; ramps up if rate-limited
    max_delay: float = 30.0,  # ceiling for the adaptive pace above
    backoff_base: float = 8.0,  # first per-attempt wait after a rate-limit/5xx, doubled each retry
    max_variants: int = 2000,
    on_progress: Optional[Callable[[int, int], None]] = None,
    on_skip: Optional[Callable[[Key, str], None]] = None,
    cache: Optional[MutableMapping[str, dict]] = None,
) -> Dict[Key, float]:
    """Query gnomAD's public API one variant at a time.

    Returns ``{key: allele_frequency}``. Variants gnomAD has never seen are
    left out of the result -- the same convention ``annotate.load_frequencies``
    uses for "no data found".

    gnomAD signals "you're going too fast" as an HTTP **400** with a JSON body
    containing the text "rate limit" -- not the HTTP 429 status the standard
    would suggest. This function checks the response body for that phrase (in
    addition to honouring standard 429/5xx codes), and when it sees it: (1)
    backs off and retries that one variant, and (2) permanently increases the
    pace used for every request afterwards (``delay`` ramps towards
    ``max_delay``), since hitting the limit once means the current pace is
    too fast for the rest of the batch too, not just this one variant.

    A single variant that gnomAD's API refuses for a genuine, non-rate-limit
    reason (e.g. a malformed variant ID) does **not** abort the whole batch:
    it is skipped (``on_skip`` is called with the key and a reason, if given)
    and the run continues. Only if *every* variant in the batch fails does
    this raise ``GnomadAPIError``, since that pattern means something
    systemic is wrong (no network, wrong dataset name, etc.).

    ``cache``, if given, is a plain ``{variant_id_string: {...}}`` dict
    (JSON-serialisable, so the caller can persist it to disk between runs).
    A variant already in the cache is never re-queried over the network --
    its cached frequency (or cached "confirmed absent") is used directly.
    Every *successfully resolved* variant this call makes (found or
    confirmed absent) is written into it as it goes, so a caller that saves
    ``cache`` after the call builds up a growing local copy of gnomAD over
    repeated runs. Skipped/failed lookups are deliberately *not* cached,
    so they get a fresh attempt next run instead of failing forever.

    Raises ``ValueError`` if ``keys`` has more than ``max_variants`` entries.
    """
    keys = list(keys)
    if len(keys) > max_variants:
        raise ValueError(
            f"{len(keys)} variants requested, over max_variants={max_variants}. "
            "Querying gnomAD's public API one variant at a time is meant for small "
            "batches; for more, download a gnomAD subset instead (see the README)."
        )

    freqs: Dict[Key, float] = {}
    skipped = 0
    attempted = 0  # keys that actually needed a network call (excludes cache hits)
    current_delay = delay
    for i, key in enumerate(keys):
        vid_str = variant_id(key)

        if cache is not None and vid_str in cache:
            cached = cache[vid_str]
            if "af" in cached:
                freqs[key] = cached["af"]
            if on_progress:
                on_progress(i + 1, len(keys))
            continue  # cache hit: no network call, no pacing delay

        attempted += 1
        skip_reason: Optional[str] = None
        found_af: Optional[float] = None
        confirmed_absent = False
        for attempt in range(max_retries):
            try:
                data = _post(_QUERY, {"variantId": vid_str, "dataset": dataset}, timeout)
            except urllib.error.HTTPError as exc:
                detail = _http_error_detail(exc)
                rate_limited = exc.code in _RETRYABLE_STATUS or _is_rate_limit_message(detail)
                if not rate_limited:
                    # A real, non-rate-limit error: won't fix itself by waiting, but
                    # not worth aborting the other variants over -- skip just this one.
                    skip_reason = f"HTTP {exc.code} ({exc.reason})" + (f": {detail}" if detail else "")
                    break
                current_delay = min(current_delay * 2, max_delay)  # slow down for the rest of the batch
                if attempt == max_retries - 1:
                    skip_reason = f"rate limited after {max_retries} attempts (HTTP {exc.code})"
                    break
                wait = _retry_after_seconds(exc) or backoff_base * (2 ** attempt)
                time.sleep(wait)
                continue
            except (urllib.error.URLError, TimeoutError) as exc:
                if attempt == max_retries - 1:
                    skip_reason = f"network error after {max_retries} attempts: {exc}"
                    break
                time.sleep(backoff_base * (2 ** attempt))
                continue

            if not data.get("errors"):
                variant = (data.get("data") or {}).get("variant")
                if variant:
                    found_af = _combined_af(variant)
                confirmed_absent = True  # a clean response, whether or not it had data, is definitive
            else:
                confirmed_absent = True  # gnomAD explicitly saying "not found" is also definitive
            break  # got a usable response

        if skip_reason:
            skipped += 1
            if on_skip:
                on_skip(key, skip_reason)
        else:
            if found_af is not None:
                freqs[key] = found_af
            if confirmed_absent and cache is not None:
                cache[vid_str] = {"af": found_af} if found_af is not None else {"absent": True}

        if on_progress:
            on_progress(i + 1, len(keys))
        if i + 1 < len(keys):
            time.sleep(current_delay)

    if attempted and skipped == attempted:
        raise GnomadAPIError(
            f"All {attempted} looked-up variant(s) failed -- this looks systemic "
            "(no network, wrong --dataset, or gnomAD is fully down), not a per-variant issue."
        )
    return freqs
