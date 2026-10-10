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
| 09 / 10 / 2026 | Umbral de cuentas mínimas | Hay un valle entre gotas vacías y células en el histograma, consistente con el acantilado de la curva de rango (~4000) | se probaron valores de 200,300, 400, 500, 1000, 2000, 5000, los resultados se mantienen estables de 200 a 2000 y 5000 recorta el pico de células |
