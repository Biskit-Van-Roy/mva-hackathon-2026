# Hoja de ruta competitiva — MVA Hackathon 2026

Fecha de inicio: 29 de agosto de 2026. Cierre oficial: 24 de octubre de 2026 a las 23:59 UTC.

## Qué estamos resolviendo

El Track 1 es un diagnóstico genómico N-of-1. Debemos ordenar hasta diez hipótesis. Cada hipótesis puede ser una sola variante o un par *compound heterozygous*. El evaluador oficial revela que la respuesta clínica contiene dos variantes en combinación; por eso el análisis por pares y por gen es central.

El Track 2 reutiliza el mecanismo biológico del Track 1 para proponer fármacos ya aprobados como hipótesis investigables. No es una recomendación clínica ni una afirmación de eficacia.

## Fases

### Fase 0 — Acceso y seguridad (29–30 agosto)

- Aceptar las reglas del dataset con la cuenta de Hugging Face.
- Descargar únicamente VCF, índice y fenotipo.
- Ejecutar el informe de entrada.
- Confirmar GRCh38, una sola muestra, campos de calidad y anotaciones disponibles.

Salida: `artifacts/intake_report.json`.

### Fase 1 — Fenotipo computable (31 agosto–2 septiembre)

- Revisar síntomas y términos HPO localmente.
- Normalizar HPO, edad de inicio y hallazgos ausentes explícitos.
- Crear un mapa fenotipo → enfermedad → gen con fuentes versionadas.

Salida: perfil HPO depurado y lista de genes priorizados con evidencia.

### Fase 2 — Variante y herencia (3–10 septiembre)

- Normalizar y descomponer el VCF contra GRCh38.
- Anotar consecuencia, frecuencia poblacional, ClinVar y predictores funcionales.
- Aplicar filtros de calidad y rareza sin eliminar de forma rígida variantes plausibles.
- Construir pares heterocigotos dentro del mismo gen y evaluar fase cuando sea posible.
- Mantener carriles separados para SNV/indel, splicing, CNV/SV y mosaicismo.

Salida: tabla auditable de candidatos y razones de inclusión/exclusión.

### Fase 3 — Ranking multimodal (11–20 septiembre)

- Combinar evidencia fenotípica, genética, herencia, calidad y literatura.
- Calibrar EPCR para que el orden sea explícito y reproducible.
- Ejecutar análisis de sensibilidad: con y sin conocimiento previo de MVA.
- Revisar manualmente lecturas/coverage de los principales candidatos si es necesario.

Salida: primer top 10 científico, todavía sin enviar.

### Fase 4 — Primer envío y aprendizaje (21–25 septiembre)

- Validar el CSV localmente.
- Congelar código, configuración y reporte del modelo 1.
- Usar el primer envío y registrar puntajes sin ajustar a ciegas.

Salida: baseline reproducible y resultado del leaderboard.

### Fase 5 — Profundización e innovación (26 septiembre–10 octubre)

- Revisar limitaciones del baseline.
- Usar FASTQ solo para preguntas concretas no resueltas.
- Añadir evidencia estructural, de splicing o mosaicismo si los datos lo sostienen.
- Preparar Track 2: mecanismo, red de vías, candidatos de fármacos y riesgos.

### Fase 6 — Entrega final (11–23 octubre)

- Reservar al menos dos intentos de Track 1 para mejoras justificadas.
- Finalizar repositorio público sin datos restringidos.
- Finalizar reporte con trazabilidad y limitaciones.
- Grabar pitch de tres minutos para Track 2.
- Auditar privacidad, reproducibilidad, enlaces y licencias.

## Estrategia de puntuación

- Pareja correcta en rango 1: 100 puntos.
- Pareja correcta en rangos 2–3: 50 puntos.
- Pareja correcta en rangos 4–5: 25 puntos.
- Pareja correcta en rangos 6–10: 10 puntos.
- Una sola variante correcta del par recibe la mitad del puntaje de ese rango.
- F-max se calcula a nivel de variantes individuales recorriendo umbrales EPCR.

Consecuencia práctica: el primer renglón debe ser la mejor pareja completa, no solo el mejor alelo aislado. Agregar ruido por encima del candidato correcto perjudica precisión y rango.

## Reglas de trabajo

- Ningún dato restringido se versiona o comparte.
- Cada decisión de filtrado queda registrada con versión de herramienta y base de datos.
- Los resultados negativos y limitaciones se documentan.
- Ningún fármaco se presenta como tratamiento; todas las propuestas son hipótesis para validación independiente.
