# Especificación oficial verificada

Verificada el 29 de agosto de 2026 contra la aplicación y el código público del Hackathon.

## Track 1

- Máximo seis envíos por participante; solo cuenta el mejor.
- Un CSV con 1–10 hipótesis ordenadas por EPCR.
- `proband_id` aceptado: `PROBAND01`.
- Coordenadas: GRCh38.
- Una fila puede representar una variante única o un par *compound heterozygous*.
- Entregables adicionales obligatorios: repositorio público de GitHub y reporte PDF o Markdown.
- El evaluador oficial indica que el *ground truth* es un par *compound heterozygous*.

Columnas exactas:

```text
proband_id,chrom_1,pos_1,ref_1,alt_1,chrom_2,pos_2,ref_2,alt_2,epcr,finding_type,notes
```

Puntaje de rango: 100 para rango 1; 50 para rangos 2–3; 25 para rangos 4–5; 10 para rangos 6–10. Recuperar solo una de las dos variantes obtiene la mitad.

## Track 2

- La interfaz activa acepta un solo envío por equipo, aunque una descripción anterior haya indicado tres.
- Reporte PDF o Markdown.
- Repositorio público de GitHub.
- Video de tres minutos en YouTube o Vimeo.
- Evaluación: rigor científico 35%, impacto 25%, innovación 25%, escalabilidad 15%.

## Datos y reglas

- Dataset restringido de aproximadamente 85 GB y once archivos.
- Primera descarga recomendada por este proyecto: documento clínico, VCF de 315 MB e índice.
- No redistribuir los datos ni intentar contactar al sujeto, su familia o contactos de MVA Society.
- Eliminar todos los datos e intermedios dentro de los 30 días posteriores al cierre y confirmar la eliminación por correo a los organizadores.
- Los resultados, el código y los reportes se publican con licencia CC BY.

## Fuentes

- Dataset: https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data
- Track 1: https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/submit_track1.py
- Evaluador: https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/evaluation.py
- Track 2: https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/submit_track2.py
- Reglas: https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/rules.py
