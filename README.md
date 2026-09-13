# MVA Hackathon 2026 — starter kit

Base reproducible y segura para trabajar primero en **Track 1: Variant Prediction** y reutilizar después los hallazgos en **Track 2: Drug Repurposing**.

## Decisión de estrategia

1. Empezar con el VCF ya generado (~315 MB) y el documento clínico.
2. Inspeccionar el contenido y confirmar GRCh38, muestras, anotaciones y términos HPO.
3. Construir una priorización de variantes y pares *compound heterozygous*.
4. Usar los FASTQ (~85 GB en total) únicamente para comprobaciones que el VCF no resuelva: calidad, cobertura, mosaicismo, CNV/SV o realineamiento.
5. No gastar ninguno de los seis envíos hasta validar localmente el formato y revisar científicamente los diez candidatos.

## Privacidad crítica

- `data/`, `artifacts/` y archivos genómicos están excluidos de Git.
- No subas FASTQ, VCF, el documento clínico ni resultados intermedios sensibles a GitHub, chats, unidades compartidas o servicios externos.
- El repositorio público debe contener código reproducible, documentación metodológica y únicamente los resultados derivados que las reglas permitan publicar.
- Conserva los datos solo durante el periodo autorizado y cumple la eliminación y confirmación exigidas por los organizadores.

## Requisitos

- Windows 10/11.
- Python 3.11 o 3.12 de 64 bits.
- Cuenta de Hugging Face con acceso aprobado al dataset.
- Al menos 2 GB libres para esta primera fase. Los FASTQ requerirán 100–150 GB más adelante.

## Inicio rápido en PowerShell

Abre PowerShell dentro de esta carpeta y ejecuta:

```powershell
Set-ExecutionPolicy -Scope Process Bypass
./scripts/setup.ps1
```

Inicia sesión en Hugging Face sin guardar el token dentro del proyecto:

```powershell
./.venv/Scripts/hf.exe auth login
```

Descarga solo el VCF, su índice y el documento clínico:

```powershell
./scripts/download-minimal.ps1
```

Genera el informe técnico de entrada:

```powershell
./scripts/run-intake.ps1
```

El resultado queda en `artifacts/intake_report.json`. Contiene estructura técnica agregada y términos HPO detectados; no copia el texto clínico ni variantes individuales.

## Crear y validar una predicción

Copia el ejemplo y reemplaza sus valores ficticios por candidatos reales:

```powershell
Copy-Item ./config/candidates.example.json ./config/candidates.local.json
./.venv/Scripts/python.exe -m mva_hackathon.submission build `
  ./config/candidates.local.json `
  ./submissions/predictions.csv
./.venv/Scripts/python.exe -m mva_hackathon.submission validate `
  ./submissions/predictions.csv
```

El validador comprueba las columnas oficiales, `PROBAND01`, GRCh38 por convención, pares completos, EPCR entre 0 y 1, orden descendente y el máximo de diez filas.

## Estructura

```text
mva-hackathon-starter/
├── config/                 # Ejemplos; candidates.local.json no se versiona
├── data/raw/               # Datos restringidos; nunca se versionan
├── artifacts/              # Informes locales; nunca se versionan
├── docs/                   # Estrategia, reglas y hoja de ruta
├── scripts/                # Comandos PowerShell reproducibles
├── src/mva_hackathon/      # Código Python
├── submissions/            # CSV locales; se revisan antes de enviar
└── tests/                  # Pruebas unitarias
```

## Siguiente hito

Ejecutar `run-intake.ps1`. Con su salida técnica podremos elegir correctamente la anotación (VEP/ANN/ClinVar/gnomAD), definir filtros y construir el primer ranking sin adivinar el esquema real del VCF.
