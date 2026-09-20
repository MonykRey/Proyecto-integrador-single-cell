# Proyecto final — equipo N

> Plantilla del módulo de análisis de células individuales, LCG 2027-1.
> **Cambia esta línea y el título por los de tu equipo.** Todo lo demás está
> pensado para que lo edites, no para que lo borres.

## Qué es este repositorio

El análisis del conjunto de datos que se le asignó a tu equipo, de principio a
fin, de forma que **cualquier persona pueda reproducirlo desde cero**. Esa es la
condición que se evalúa (criterio A2 de la rúbrica, 30 %): no que el resultado sea
bonito, sino que corra en otra máquina sin intervención manual.

## Integrantes

| Nombre | Qué hizo |
|---|---|
| | |
| | |

Esta tabla no es burocracia: el criterio B de la rúbrica es una respuesta
individual dirigida, y saber quién trabajó en qué parte es lo que la hace justa.

## Cómo reproducir este análisis

```bash
git clone <la-url-de-este-repositorio>
cd <el-repositorio>
quarto render reporte.qmd
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
| | | | |
