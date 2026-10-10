# Proyecto final — equipo C 

## Qué es este repositorio

El análisis del conjunto de datos que se le asignó a tu equipo, de principio a
fin, de forma que **cualquier persona pueda reproducirlo desde cero**. Esa es la
condición que se evalúa (criterio A2 de la rúbrica, 30 %): no que el resultado sea
bonito, sino que corra en otra máquina sin intervención manual.

## Integrantes

| Nombre | Qué hizo |
|---|---|
| Mónica Reyes Ramírez| |
| Natalie B. Pineda Morán | |


## Cómo reproducir este análisis

```bash
git clone <la-url-de-este-repositorio>
cd <el-repositorio>
quarto render reporte.qmd
```

Para reproducir este análisis se debe ***crear un archivo .venv***
```bash
python -m venv .venv
source .venv/bin/activate        # Mac/Linux
.venv\Scripts\activate           # Windows
pip install -r requirements.txt
```

**Los datos no están en el repositorio** y no deben estarlo: `.gitignore` los
excluye. Declara aquí de dónde salen.

| Dato | Dónde vive | Cómo se obtuvo |
|---|---|---|
| | `/ruta/en/el/cluster` | |

## Estructura

| Carpeta | Qué va aquí |
|---|---|
| `reporte.qmd` | El reporte ejecutable. Es el entregable |
| `R/` | Funciones que uses en más de un lugar |
| `datos/` | Vacía a propósito: los datos no se versionan |
| `figuras/` | Salida del reporte. Tampoco se versiona |

## Bitácora de decisiones

Cada vez que tomes una decisión de análisis —un umbral, un método, un filtro—
anótala aquí con su razón. **La tabla de decisiones es el criterio A3 de la
rúbrica** (15 %), y llenarla al final, de memoria, se nota.

| Fecha | Decisión | Por qué | Qué se probó antes |
|---|---|---|---|
|   | Umbral de cuentas mínimas | Casi todos los barcodes son gotas vacías con muy pocas cuentas. En el histograma se ven dos grupos separados (vacías y células), y 500 cae en el espacio entre ellos. La curva de rango también cae de golpe en ~4,000 barcodes| se probaron valores de 200,300, 400, 500, 1000, 2000, 5000, los resultados se mantienen estables de 200 a 2000 y 5000 recorta el pico de células |
|   | Umbral de genes detectados | Una célula dañada o un resto celular detecta pocos genes, y por eso se quita lo que está por debajo de 1,551, un grupo bajo separado del grupo principal. Arriba se dejó un límite amplio (6,249), porque los dobletes (dos células en una gota) se revisan mejor con otra herramienta| Tres niveles de exigencia (k = 3, 4 y 5). Con el más estricto se perdía el 18 % de los barcodes, demasiado, porque ya cortaba células normales|
|   | Fracción mitocondria máxima | Cuando una célula se rompe, pierde casi todo su RNA excepto el mitocondrial, así que un % mt alto indica célula dañada. La mayoría de tus células está en ~4.4 %, y el corte solo quita la cola larga. No usamos un valor de específico, porque no sabemos el tejido| Tres niveles de exigencia (cortes de 6.95 %, 7.80 % y 8.65 %). Se eligió el más permisivo porque no hay un grupo claramente anormal en ese rango.|





## Referencias usadas:
- Luecken MD, Theis FJ. Current best practices in single-cell RNA-seq analysis: a tutorial. Mol Syst Biol. 2019 Jun 19;15(6):e8746. doi: 10.15252/msb.20188746. PMID: 31217225; PMCID: PMC6582955.
