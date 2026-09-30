import io
import json
import urllib.error
from unittest.mock import patch

import pytest

from varlens.gnomad_api import GnomadAPIError, query_frequencies


class _FakeResponse:
    def __init__(self, payload: dict):
        self._body = json.dumps(payload).encode("utf-8")

    def read(self):
        return self._body

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _api_payload(exome=None, genome=None, errors=None):
    if errors:
        return {"errors": errors}
    return {"data": {"variant": {"exome": exome, "genome": genome}}}


def _http_error(code, headers=None, body=b""):
    return urllib.error.HTTPError(
        url="https://gnomad.broadinstitute.org/api",
        code=code,
        msg="error",
        hdrs=headers or {},
        fp=io.BytesIO(body),
    )


def test_combines_exome_and_genome_counts():
    payload = _api_payload(exome={"ac": 1, "an": 100}, genome={"ac": 3, "an": 300})
    with patch("varlens.gnomad_api.urllib.request.urlopen", return_value=_FakeResponse(payload)):
        freqs = query_frequencies([("1", 100, "A", "T")], delay=0, backoff_base=0)
    assert freqs[("1", 100, "A", "T")] == pytest.approx(4 / 400)


def test_missing_variant_is_absent_not_zero():
    payload = _api_payload(errors=[{"message": "Variant not found"}])
    with patch("varlens.gnomad_api.urllib.request.urlopen", return_value=_FakeResponse(payload)):
        freqs = query_frequencies([("1", 999, "A", "T")], delay=0, backoff_base=0)
    assert ("1", 999, "A", "T") not in freqs


def test_only_genome_data_is_used_when_no_exome():
    payload = _api_payload(exome=None, genome={"ac": 5, "an": 1000})
    with patch("varlens.gnomad_api.urllib.request.urlopen", return_value=_FakeResponse(payload)):
        freqs = query_frequencies([("1", 1, "A", "T")], delay=0, backoff_base=0)
    assert freqs[("1", 1, "A", "T")] == pytest.approx(0.005)


def test_transport_error_raises_after_retries():
    with patch(
        "varlens.gnomad_api.urllib.request.urlopen",
        side_effect=urllib.error.URLError("boom"),
    ):
        with pytest.raises(GnomadAPIError):
            query_frequencies([("1", 1, "A", "T")], max_retries=2, delay=0, backoff_base=0)


def test_429_retries_and_then_succeeds():
    ok_payload = _api_payload(exome={"ac": 1, "an": 10})
    responses = [_http_error(429), _FakeResponse(ok_payload)]

    def _fake_urlopen(*_a, **_k):
        result = responses.pop(0)
        if isinstance(result, Exception):
            raise result
        return result

    with patch("varlens.gnomad_api.urllib.request.urlopen", side_effect=_fake_urlopen):
        freqs = query_frequencies([("1", 1, "A", "T")], delay=0, backoff_base=0, max_retries=3)
    assert freqs[("1", 1, "A", "T")] == pytest.approx(0.1)


def test_429_honours_retry_after_header():
    with patch(
        "varlens.gnomad_api.urllib.request.urlopen",
        side_effect=_http_error(429, headers={"Retry-After": "7"}),
    ):
        with patch("varlens.gnomad_api.time.sleep") as fake_sleep:
            with pytest.raises(GnomadAPIError):
                query_frequencies([("1", 1, "A", "T")], max_retries=2, delay=0, backoff_base=99)
    # 7 (from Retry-After) should win over backoff_base=99
    waits = [call.args[0] for call in fake_sleep.call_args_list]
    assert 7.0 in waits


def test_400_with_rate_limit_message_is_treated_as_retryable():
    """This is the real behaviour observed from gnomAD: HTTP 400, not 429, with a
    JSON body saying 'rate limit'. It must be retried, not skipped as a bad request."""
    rate_limited = _http_error(400, body=b'{"errors":[{"message":"Query rate limit exceeded. Please try again in a few minutes."}]}')
    ok_payload = _api_payload(exome={"ac": 1, "an": 4})
    responses = [rate_limited, _FakeResponse(ok_payload)]

    def _fake_urlopen(*_a, **_k):
        result = responses.pop(0)
        if isinstance(result, Exception):
            raise result
        return result

    with patch("varlens.gnomad_api.urllib.request.urlopen", side_effect=_fake_urlopen):
        freqs = query_frequencies([("1", 1, "A", "T")], delay=0, backoff_base=0, max_retries=3)
    assert freqs[("1", 1, "A", "T")] == pytest.approx(0.25)


def test_rate_limit_permanently_slows_down_subsequent_requests():
    """Once gnomAD signals its rate limit, the pace for the REST of the batch should
    increase too, not just the retry of the one variant that got limited."""
    rate_limited = _http_error(400, body=b'{"errors":[{"message":"Query rate limit exceeded"}]}')
    ok = _api_payload(exome={"ac": 1, "an": 10})
    # variant 1: rate-limited once then succeeds; variant 2: succeeds first try
    responses = [rate_limited, _FakeResponse(ok), _FakeResponse(ok)]

    def _fake_urlopen(*_a, **_k):
        result = responses.pop(0)
        if isinstance(result, Exception):
            raise result
        return result

    sleeps = []
    with patch("varlens.gnomad_api.urllib.request.urlopen", side_effect=_fake_urlopen):
        with patch("varlens.gnomad_api.time.sleep", side_effect=lambda s: sleeps.append(s)):
            query_frequencies(
                [("1", 1, "A", "T"), ("1", 2, "A", "T")],
                delay=1.0,
                backoff_base=0,
                max_retries=3,
            )
    # the inter-variant pause after the rate-limited variant should be > the starting delay (1.0)
    assert max(sleeps) > 1.0


def test_cache_hit_skips_network_entirely():
    with patch("varlens.gnomad_api.urllib.request.urlopen", side_effect=AssertionError("should not be called")):
        freqs = query_frequencies(
            [("1", 1, "A", "T")],
            delay=0,
            backoff_base=0,
            cache={"1-1-A-T": {"af": 0.42}},
        )
    assert freqs[("1", 1, "A", "T")] == pytest.approx(0.42)


def test_cache_records_confirmed_absent_and_found():
    found_payload = _api_payload(exome={"ac": 1, "an": 4})
    absent_payload = _api_payload(errors=[{"message": "not found"}])
    responses = [_FakeResponse(found_payload), _FakeResponse(absent_payload)]

    def _fake_urlopen(*_a, **_k):
        return responses.pop(0)

    cache: dict = {}
    with patch("varlens.gnomad_api.urllib.request.urlopen", side_effect=_fake_urlopen):
        query_frequencies(
            [("1", 1, "A", "T"), ("1", 2, "A", "T")],
            delay=0,
            backoff_base=0,
            cache=cache,
        )
    assert cache["1-1-A-T"] == {"af": pytest.approx(0.25)}
    assert cache["1-2-A-T"] == {"absent": True}


def test_skipped_variant_is_not_cached():
    bad = _http_error(400, body=b'{"errors":[{"message":"Invalid variant id"}]}')
    ok_payload = _api_payload(exome={"ac": 1, "an": 4})
    responses = [bad, _FakeResponse(ok_payload)]

    def _fake_urlopen(*_a, **_k):
        result = responses.pop(0)
        if isinstance(result, Exception):
            raise result
        return result

    cache: dict = {}
    with patch("varlens.gnomad_api.urllib.request.urlopen", side_effect=_fake_urlopen):
        query_frequencies(
            [("1", 1, "A", "T"), ("1", 2, "A", "T")],
            delay=0,
            backoff_base=0,
            cache=cache,
        )
    assert "1-1-A-T" not in cache  # the failed one must stay out, so it's retried next run
    assert cache["1-2-A-T"] == {"af": pytest.approx(0.25)}


def test_non_retryable_http_error_is_skipped_not_fatal_when_other_variants_exist():
    """A single 400 shouldn't take down a whole batch: the other variant still resolves."""
    ok_payload = _api_payload(exome={"ac": 2, "an": 20})
    bad = _http_error(400, body=b'{"errors":[{"message":"Invalid variant id"}]}')
    responses = [bad, _FakeResponse(ok_payload)]

    def _fake_urlopen(*_a, **_k):
        result = responses.pop(0)
        if isinstance(result, Exception):
            raise result
        return result

    skipped = []
    with patch("varlens.gnomad_api.urllib.request.urlopen", side_effect=_fake_urlopen):
        freqs = query_frequencies(
            [("1", 1, "A", "T"), ("1", 2, "A", "T")],
            delay=0,
            backoff_base=0,
            on_skip=lambda key, reason: skipped.append((key, reason)),
        )
    assert ("1", 1, "A", "T") not in freqs
    assert freqs[("1", 2, "A", "T")] == pytest.approx(0.1)
    assert len(skipped) == 1
    assert skipped[0][0] == ("1", 1, "A", "T")
    assert "400" in skipped[0][1]
    assert "Invalid variant id" in skipped[0][1]  # body detail must be surfaced, not just "Bad Request"


def test_all_variants_failing_is_treated_as_systemic_and_raises():
    with patch("varlens.gnomad_api.urllib.request.urlopen", side_effect=_http_error(400)):
        with pytest.raises(GnomadAPIError):
            query_frequencies([("1", 1, "A", "T")], max_retries=2, delay=0, backoff_base=0)


def test_max_variants_guard_rejects_large_batches():
    keys = [("1", i, "A", "T") for i in range(5)]
    with pytest.raises(ValueError):
        query_frequencies(keys, max_variants=3, delay=0, backoff_base=0)


def test_progress_callback_is_called_once_per_variant():
    payload = _api_payload(exome={"ac": 0, "an": 100})
    calls = []
    with patch("varlens.gnomad_api.urllib.request.urlopen", return_value=_FakeResponse(payload)):
        query_frequencies(
            [("1", 1, "A", "T"), ("1", 2, "A", "T")],
            delay=0,
            backoff_base=0,
            on_progress=lambda done, total: calls.append((done, total)),
        )
    assert calls == [(1, 2), (2, 2)]
