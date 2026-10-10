# Proyecto final — equipo C 

## Qué es este repositorio

El análisis del conjunto de datos que se le asignó a tu equipo, de principio a
fin, de forma que **cualquier persona pueda reproducirlo desde cero**. Esa es la
condición que se evalúa (criterio A2 de la rúbrica, 30 %): no que el resultado sea
bonito, sino que corra en otra máquina sin intervención manual.

## Integrantes

| Nombre | Qué hizo |
|---|---|
| Mónica Reyes Ramírez| Control de calidad y Normalización/selección de genes  |
| Natalie B. Pineda Morán | Preguntas iniciales y Reducción de dimensionalidad /agrupamiento |


## Cómo reproducir este análisis

```bash
git clone <url-del-repositorio>
cd Proyecto-integrador-single-cell
```
### Requisitos 
- Visual Studio Code
- Python 3.11
- conda (recomendado) o pip

### Datos
- Descarga los datos del equipo C a tu computadora, define la ruta donde se va a guardar esa carpeta que debe contener los archivos **barcodes.tsv, genes.tsv y matrix.mtx.**
- Luego de clonar el repositorio, abre el archivo `analisis_C.py` en Visual Studio Code y cambia la ruta (que se almacena en la variable RUTA_DATOS) a la ruta donde tienes los datos **barcodes.tsv, genes.tsv y matrix.mtx.**

### Librerias necesarias para este análisis
- numpy
- pandas
- matplotlib
- scanpy
- python-igraph
- leidenalg
- scikit-misc
- scrublet

### Comandos para instalar las librerías en un entorno (recomendable)
**En tu terminal de linux/Mac**

```bash
conda create -n sc-equipoC python=3.11 -y
conda activate sc-equipoC
conda install -c conda-forge numpy pandas matplotlib scanpy python-igraph leidenalg scikit-misc jupyter ipykernel -y
pip install scrublet
```
> Si `pip` no funciona, prueba con `pip3`.

### Correr el análisis
```bash
conda activate sc-equipoC
python analisis_C.py
```

### Salidas
- Las figuras se guardan automáticamente en `figuras/` al terminar el script.
- Los datos podrás visualizarlos en tu terminal de Visual Studio Code

### Otras formas de poder tener las dependencias

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

| Fecha | Valor decidido | Decisión | Por qué | Qué se probó antes |
|---|---|---|---| ---|
|   | 500| Umbral de cuentas mínimas | Casi todos los barcodes son gotas vacías con muy pocas cuentas. En el histograma se ven dos grupos separados (vacías y células), y 500 cae en el espacio entre ellos. La curva de rango también cae de golpe en ~4,000 barcodes| Se probaron valores de 200,300, 400, 500, 1000, 2000, 5000, los resultados se mantienen estables de 200 a 2000 y 5000 recorta el pico de células |
|   | Entre 1,551 y 6,249| Umbral de genes detectados | Una célula dañada o un resto celular detecta pocos genes, y por eso se quita lo que está por debajo de 1,551, un grupo bajo separado del grupo principal. Arriba se dejó un límite amplio (6,249), porque los dobletes (dos células en una gota) se revisan mejor con otra herramienta| Tres niveles de exigencia (k = 3, 4 y 5). Con el más estricto se perdía el 18 % de los barcodes, demasiado, porque ya cortaba células normales|
|   | 8.65% | Fracción mitocondrial máxima | Cuando una célula se rompe, pierde casi todo su RNA excepto el mitocondrial, así que un % mt alto indica célula dañada. La mayoría de tus células está en ~4.4 %, y el corte solo quita la cola larga. No usamos un valor de específico, porque no sabemos el tejido| Tres niveles de exigencia (cortes de 6.95 %, 7.80 % y 8.65 %). Se eligió el más permisivo porque no hay un grupo claramente anormal en ese rango.|
|   | 3 | Mínimo de células por gen | Quitamos genes detectados en 1-2 células de 3,749 células que pueden ser ruido. No encontramos un punto de quiebre donde pudiéramos definir claramente el umbral, entonces nos guiamos del estándar de la literatura| 1, 3, 5, 10 y 20, se ve una caída gradual sin puntos de quiebre, con esto nos quedan 23,801 genes en un mínimo de 3 células|
|   | No se filtró | Filtro de dobletes | El histograma de Scrublet no muestra dos grupos separados así que el umbral automático no tiene respaldo. Pero estos valores se anotaron | Scrublet tiene un umbral automático, encontramos 26 dobletes (0.7%) |
|   | Escalar cada célula a 10,000 cuentas totales y aplicar logaritmo | Normalización | Revisamos el histograma del tamaño de biblioteca y vimos que cada célula tenía una cantidad muy distinta de cuentas (de 2,679 a 29,371; mediana 7,928), lo cual refleja diferencias de secuenciación y no biología. Por eso escalamos cada célula a 10,000 cuentas totales. Elegimos 10,000 porque es un valor estándar y porque queda cerca de la mediana de nuestros datos, así que no distorsiona las cuentas de la mayoría de las células. Después de normalizar, todas suman exactamente 10,000, y el logaritmo reduce el peso de los valores extremos | Hicimos un histograma del tamaño de la biblioteca y comparamos contra 1e4 contra usar la mediana como alternativa|
|   | Quedarnos con los 2,000 genes que más cambian entre células  | Selección de genes variables | Elegimos 2,000 genes con seurat_v3 sobre las cuentas originales, porque es un número estándar y porque los datos lo respaldan: el gen en el lugar 2,000 todavía varía más de lo esperado por azar (varianza normalizada de 1.41, cuando sin variación extra sería ~1), y la mediana de los HVG es 1.78 contra 0.92 del resto de los genes. Además, casi no hay genes técnicos entre ellos (0 ribosomales y 2 mitocondriales), y los primeros son marcadores de células inmunes | Comparamos  1,000, 2,000 y 3,000 genes, pero como cada conjunto incluye al anterior esa comparación no aportó información, así que nos basamos en la variabilidad del último gen elegido|





## Referencias usadas:
- Luecken MD, Theis FJ. Current best practices in single-cell RNA-seq analysis: a tutorial. Mol Syst Biol. 2019 Jun 19;15(6):e8746. doi: 10.15252/msb.20188746. PMID: 31217225; PMCID: PMC6582955.
