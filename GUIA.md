# Guía de varlens (instalación, uso y comandos)

> Esta guía está en español, pensada como referencia rápida propia. El `README.md` (en inglés)
> es la documentación pública del repositorio para quien lo visite en GitHub; ambos describen
> la misma herramienta.

## 1. Qué es

**varlens** es una herramienta de línea de comandos escrita en Python puro (sin dependencias
externas) que:

1. Lee un archivo **VCF** (las variantes genéticas de una persona).
2. Las cruza con **ClinVar** (qué se sabe sobre si esas variantes causan enfermedad) y con la
   **frecuencia poblacional** (gnomAD): cuánta gente sana tiene esa misma variante.
3. Aplica unas reglas transparentes, inspiradas en los criterios ACMG/AMP que se usan en
   genética clínica real, para clasificar cada variante en uno de 4 niveles de prioridad.
4. Genera un informe HTML navegable (y opcionalmente una tabla TSV).

**Importante:** es una herramienta educativa/de portfolio, no un clasificador clínico. Así lo
dice el propio informe HTML que genera.

## 2. Cómo funciona (los módulos, por si te preguntan en una entrevista)

- `vcf.py` — lee el VCF línea a línea, separa filas multialélicas en una variante por alelo, y
  sabe filtrar por región o por un conjunto de variantes de interés sin cargar el archivo entero
  en memoria.
- `annotate.py` — cruza las variantes con ClinVar (clasifica su significado clínico en
  patogénica/benigna/VUS/conflictiva, y calcula las estrellas de revisión) y con la frecuencia
  poblacional.
- `gnomad_api.py` — alternativa a descargar gnomAD entero: consulta su API pública variante por
  variante, con reintentos, control del límite de peticiones (gnomAD lo señala como un HTTP 400
  con un mensaje de texto, no como el 429 estándar) y una caché en disco para no repetir
  consultas ya hechas.
- `rules.py` — el corazón de la clasificación: BA1/BS1 (demasiado común para ser la causa de una
  enfermedad rara), PM2 (rara o ausente), LoF (pérdida de función) y el propio ClinVar, combinados
  en 4 niveles: **Alta prioridad**, **Revisión manual**, **Incierta/sin anotar**, **Baja prioridad**.
- `report.py` — genera el HTML autocontenido (con modo oscuro automático) y la tabla TSV.
- `cli.py` — une todas las piezas y expone las opciones de línea de comandos.

## 3. Instalación desde cero (Windows, cmd)

### 3.1. Comprobar/instalar Python 3.9 o superior

```
py -0
```

Si no aparece ninguna versión 3.9+, descárgala de [python.org/downloads](https://www.python.org/downloads/)
marcando la casilla **"Add python.exe to PATH"** durante la instalación.

### 3.2. Crear el entorno virtual

Dentro de la carpeta `varlens` (donde está `pyproject.toml`):

```
py -3.11 -m venv .venv
```

(cambia `3.11` por la versión que tengas instalada, según lo que viste en `py -0`)

### 3.3. Activar el entorno

```
.venv\Scripts\activate.bat
```

Si usas PowerShell y da un error de política de ejecución, o usa `cmd` en su lugar, o ejecuta
antes: `Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass`.

### 3.4. Actualizar pip e instalar el proyecto

```
python -m pip install --upgrade pip
pip install -e ".[dev]"
```

### 3.5. Ejecutar los tests

```
pytest -q
```

Deberías ver `49 passed` (si no tienes Node.js instalado, el test que compara la demo de la landing con las reglas de Python no hace nada y también cuenta como pasado).

## 4. Cada vez que abres `cmd` de nuevo

Solo hace falta repetir esto (el entorno virtual se "olvida" al cerrar la ventana):

```
cd Desktop\varlens
.venv\Scripts\activate.bat
```

## 5. Probar con los datos de ejemplo (sintéticos, incluidos en el repo)

```
varlens data\example_patient.vcf --clinvar data\example_clinvar.vcf --freq data\example_freq.vcf --out report.html --tsv report.tsv
```

## 6. Descargar datos reales

| Qué | De dónde |
|---|---|
| Genoma real de referencia (HG001/NA12878, sano, GRCh38) | https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/ |
| ClinVar real (GRCh38) | https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz |

Guárdalos, por ejemplo, en una carpeta `data_real\` dentro de `varlens`.

## 7. Limitar qué parte del genoma escanear

**Opción recomendada — `--region`, directo sobre el genoma completo, sin crear archivos nuevos:**

```
varlens data_real\HG001_GRCh38.vcf.gz --region 21 --clinvar data_real\clinvar.vcf.gz --out report.html
varlens data_real\HG001_GRCh38.vcf.gz --region 21:20000000-26000000 --clinvar data_real\clinvar.vcf.gz --out report.html
```

`--region` acepta un cromosoma entero (`21`) o una ventana de coordenadas (`21:inicio-fin`).
Ojo: el principio de los cromosomas acrocéntricos (13, 14, 15, 21, 22) suele estar vacío de
variantes llamadas (es zona de heterocromatina cerca del centrómero) — si una región no
devuelve nada, prueba más adentro (por ejemplo, a partir de 20.000.000).

**Opción antigua — recortar un archivo nuevo con `subset_vcf.py`** (útil si quieres un archivo
fijo y pequeño que reutilizar varias veces sin volver a filtrar cada vez):

```
python scripts\subset_vcf.py data_real\HG001_GRCh38.vcf.gz data_real\HG001_chr21.vcf.gz --chrom 21 --max-variants 1500
```

## 8. Frecuencia poblacional: API de gnomAD vs archivo descargado

**Sin descargar nada (API pública, recomendado para pocas variantes o una región concreta):**

```
varlens data_real\HG001_GRCh38.vcf.gz --region 21 --clinvar data_real\clinvar.vcf.gz --freq-api --out report.html --tsv report.tsv
```

Parámetros útiles:
- `--max-api-variants` (por defecto 4000): tope de seguridad de variantes a consultar.
- `--freq-api-delay` (por defecto 2.0s) y `--freq-api-max-delay` (por defecto 30s): ritmo entre
  peticiones; sube solo si gnomAD avisa de que vas demasiado rápido.
- `--freq-api-cache` (por defecto `.varlens_gnomad_cache.json`) y `--no-cache`: las variantes ya
  consultadas se guardan, así que volver a correr el mismo VCF no repite peticiones.

**Con un archivo descargado** (mejor para más variantes o sin conexión durante la ejecución):

```
varlens sample.vcf --clinvar clinvar.vcf.gz --freq gnomad_subset.vcf.gz --out report.html
```

## 9. Crear un control positivo real (para comprobar que detecta patogénicas de verdad)

```
python scripts\make_positive_control.py --clinvar data_real\clinvar.vcf.gz --chrom 21 --list
python scripts\make_positive_control.py --clinvar data_real\clinvar.vcf.gz --chrom 21 --patient data_real\HG001_GRCh38.vcf.gz --out data_real\HG001_positive_control.vcf.gz
varlens data_real\HG001_positive_control.vcf.gz --clinvar data_real\clinvar.vcf.gz --region 21 --out report_pc.html
```

Este script saca una variante patogénica real y verificada (con sus estrellas de revisión) de tu
propio `clinvar.vcf.gz`, y la inyecta como heterocigota en una copia del VCF del paciente.

## 10. Ajustar los umbrales de clasificación

```
varlens sample.vcf --clinvar clinvar.vcf.gz --freq-api --ba1 0.05 --bs1 0.01 --pm2 0.0001 --out report.html
```

- `--ba1`: frecuencia a partir de la cual se considera benigna sin más (por defecto 0.05 = 5%).
- `--bs1`: frecuencia demasiado alta para una enfermedad rara (por defecto 0.01 = 1%).
- `--pm2`: frecuencia por debajo de la cual se considera "rara" (por defecto 0.0001 = 0.01%).

## 11. Tabla resumen de todos los comandos

| Acción | Comando |
|---|---|
| Crear entorno virtual | `py -3.11 -m venv .venv` |
| Activarlo | `.venv\Scripts\activate.bat` |
| Instalar el proyecto | `pip install -e ".[dev]"` |
| Ejecutar tests | `pytest -q` |
| Demo con datos sintéticos | `varlens data\example_patient.vcf --clinvar data\example_clinvar.vcf --freq data\example_freq.vcf --out report.html` |
| Recortar un VCF a un cromosoma/tamaño | `python scripts\subset_vcf.py in.vcf.gz out.vcf.gz --chrom 21 --max-variants 1500` |
| Escanear solo una región, sin recortar | `varlens in.vcf.gz --region 21:20000000-26000000 --clinvar clinvar.vcf.gz --out report.html` |
| Frecuencia vía API (sin descargar gnomAD) | `varlens in.vcf.gz --clinvar clinvar.vcf.gz --freq-api --out report.html` |
| Frecuencia vía archivo local | `varlens in.vcf.gz --clinvar clinvar.vcf.gz --freq gnomad.vcf.gz --out report.html` |
| Ver candidatas patogénicas reales de un cromosoma | `python scripts\make_positive_control.py --clinvar clinvar.vcf.gz --chrom 21 --list` |
| Inyectar una patogénica real como control positivo | `python scripts\make_positive_control.py --clinvar clinvar.vcf.gz --chrom 21 --patient in.vcf.gz --out pc.vcf.gz` |
| Benchmark ciego contra ClinVar | `varlens-benchmark clinvar.vcf.gz --markdown BENCHMARK.md --embed docs\index.html` |

## 12. Benchmark ciego contra ClinVar

El benchmark comprueba cuánto acierta el programa **sin mirar la etiqueta de ClinVar** (si la
mirara sería circular: la respuesta correcta sería también un dato de entrada). Usa la frecuencia
(`AF_EXAC`, `AF_TGP`, `AF_ESP`, que ya vienen en el propio VCF de ClinVar) y la consecuencia molecular
(`MC`), y usa la etiqueta patogénica/benigna solo al final, para puntuar.

```
varlens-benchmark data_real\clinvar.vcf.gz --markdown BENCHMARK.md --embed docs\index.html
```

- Lee el ClinVar completo una vez (unos minutos) y escribe `benchmark_results.json`, `BENCHMARK.md` y,
  con `--embed`, los números dentro de la landing page.
- Más rápido para probar: `--region 21`. Respuesta más estricta: `--min-stars 2`.
- Resultados (todos con intervalo de confianza del 95 %): % de patogénicas marcadas para revisión, %
  de benignas marcadas por error, **% de patogénicas descartadas como baja prioridad** (la falsa
  tranquilidad, el error más grave), % de benignas reconocidas como baja prioridad, y el valor
  predictivo positivo junto con la prevalencia de la que depende.
- La etiqueta `high` no puede salir en este modo (necesita ClinVar), así que lo que mide es cuánta señal
  recuperan por sí solas las reglas de frecuencia y de pérdida de función.

## 13. Landing page (`docs/index.html`)

Es una página única y autocontenida (HTML, CSS y JavaScript, sin cargar nada de fuera) con una demo
interactiva de las reglas, el esquema del pipeline, la tabla de reglas, el benchmark y los comandos.

Para publicarla con GitHub Pages: en el repositorio, **Settings > Pages > Deploy from a branch**,
rama `main`, carpeta `/docs`. Quedará en `https://Print-Samuel.github.io/varlens/`. Hasta que
ejecutes el benchmark con `--embed`, la sección de resultados muestra un aviso de "todavía no
generados" y **ningún número**.

## 14. Subir esto a GitHub

1. Crea un repositorio nuevo y vacío en GitHub (sin README, sin licencia — ya los trae el
   proyecto) y copia su URL (algo como `https://github.com/Print-Samuel/varlens.git`).
2. Dentro de la carpeta `varlens`, con el entorno activado:

```
git init
git add .
git commit -m "varlens: pipeline de anotación y priorización de variantes"
git branch -M main
git remote add origin https://github.com/Print-Samuel/varlens.git
git push -u origin main
```

3. Abre `README.md` y sustituye `Print-Samuel` por tu nombre de usuario real de GitHub (aparece
   en el badge de tests y en el ejemplo de `git clone`).

Si `git` no está instalado, bájalo de [git-scm.com](https://git-scm.com/downloads) — durante la
instalación, las opciones por defecto valen para esto.
