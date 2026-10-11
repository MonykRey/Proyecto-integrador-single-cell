# Proyecto final — equipo C 🐱

## Qué es este repositorio ❓

El análisis del conjunto de datos que se le asignó a tu equipo, de principio a
fin, de forma que **cualquier persona pueda reproducirlo desde cero**. Esa es la
condición que se evalúa (criterio A2 de la rúbrica, 30 %): no que el resultado sea
bonito, sino que corra en otra máquina sin intervención manual.

## Integrantes 🧑‍🤝‍🧑

| Nombre | Qué hizo |
|---|---|
| Mónica Reyes Ramírez| Control de calidad y Normalización/selección de genes |
| Natalie B. Pineda Morán | Preguntas iniciales y Reducción de dimensionalidad /agrupamiento |


## Cómo reproducir este análisis 🧐

El análisis está pensado para correrse **de preferencia en una computadora personal con Visual Studio Code**, no en el clúster.

```bash
git clone https://github.com/MonykRey/Proyecto-integrador-single-cell/blob/master/README.md
cd Proyecto-integrador-single-cell
```

### Requisitos
- Visual Studio Code (con las extensiones de Python y, opcionalmente, Quarto)
- Python 3.11
- conda (recomendado) o pip
- [Quarto](https://quarto.org/docs/get-started/) instalado, para renderizar el reporte

### Datos
Los datos no se suben al repositorio. Para obtenerlos:

1. Descarga a tu computadora la carpeta de datos del equipo C. Si están en el clúster, puedes copiarla con:
```bash
   scp -r <usuario>@<servidor>: /mnt/data/bioinfo3/compartido/10x/datasets/proyecto-2026/proyecto-equipo-C ./data/
```
2. La carpeta debe contener **barcodes.tsv, genes.tsv y matrix.mtx**.
3. Cambia la ruta **en los dos archivos**. En ambos es la primera variable:
   - `analisis_C.py`: variable `RUTA_DATOS`, al inicio del script.
   - `reporte.qmd`: variable `RUTA_DATOS`, en el chunk `rutas`.
  

```python
   RUTA_DATOS = "/ruta/donde/guardaste/proyecto-equipo-C"
```
   Si guardas la carpeta en `data/proyecto-equipo-C/` dentro del repositorio, no necesitas cambiar nada. 

### Librerías necesarias
- numpy
- pandas
- matplotlib
- scikit-learn
- scanpy
- python-igraph
- leidenalg
- scikit-misc
- scrublet
- jupyter, pyyaml, nbformat y nbclient (para renderizar el reporte con Quarto)

### Instalar las librerías en un entorno (recomendado)
**En la terminal de Linux/Mac** (o en la terminal integrada de Visual Studio Code):

```bash
conda create -n sc-equipoC python=3.11 -y
conda activate sc-equipoC
conda install -c conda-forge numpy pandas matplotlib scikit-learn scanpy python-igraph leidenalg scikit-misc jupyter ipykernel pyyaml nbformat nbclient -y
pip install scrublet
```
> Si `pip` no funciona, prueba con `pip3`.

En Visual Studio Code, selecciona el intérprete del entorno: `Ctrl/Cmd + Shift + P` → **Python: Select Interpreter** → `sc-equipoC`.

### Correr el análisis
En la terminal de Visual Studio Code:

```bash
conda activate sc-equipoC
python analisis_C.py
```

### Renderizar el reporte con quarto
Con el **entorno activo**, en la **terminal de Rstudio** desde la carpeta del repositorio (donde esta el .qmd):

```bash
conda activate sc-equipoC
quarto render reporte.qmd
```

Esto genera `reporte.html`. Hazlo desde la terminal con `sc-equipoC` activo, porque Quarto usa el Python del entorno activo. Si falla con "Jupyter is not available", instala `jupyter` y `pyyaml` en ese entorno.

### Salidas
- Las figuras se guardan automáticamente en `figuras/`.
- Los resultados numéricos se imprimen en la terminal de Visual Studio Code.
- Los objetos intermedios (`.h5ad`) se guardan en `data/`.

### Otra forma de instalar las dependencias (sin conda)
Se puede crear un **entorno virtual** `.venv`:

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
| Matriz de cuentas 10x del equipo C (`matrix.mtx`, `genes.tsv`, `barcodes.tsv`) | `/mnt/data/bioinfo3/compartido/10x/datasets/proyecto-2026/proyecto-equipo-C` | Asignada al equipo C por el curso de Bioinformática y Estadística 3 |

## Estructura 🌳

| Carpeta o archivo | Qué va aquí |
|---|---|
| `reporte.qmd` | El reporte ejecutable (Quarto + Python). |
| `analisis_C_part3.py` | Script del análisis por celdas con todo el código completo del pipeline |
| `data/` | Vacía a propósito:  Aquí se coloca la matriz 10x (`matrix.mtx`, `genes.tsv`, `barcodes.tsv`) y se guardan los objetos intermedios `.h5ad` |
| `figuras/` | Salida del reporte y del script. |
| `requirements.txt` | Paquetes de Python necesarios para reproducir el análisis |
| `README.md` | Descripción del proyecto |

## Bitácora de decisiones 📖

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
|   | 20 | Número de PCs | El codo tiene un primer quiebre en el PC 5 y un segundo cerca del PC 20-22, después de lo cual la curva es una meseta (cada PC aporta ~0.0013-0.0018). Con 20 PCs el agrupamiento es estable | 10, 15, 20, 30 y 40 PCs. El ARI contra la referencia (20 PCs) fue de 0.975, 0.945 y 0.951 con 15, 30 y 40 PCs; con 10 PCs bajó a 0.888 y se perdió un grupo |
|   | 0.5 | Resolución de Leiden | Tiene la silueta más alta (0.355) y el grupo más chico tiene 28 células. Con resolución ≥ 0.8 la silueta cae a 0.235 y aparecen grupos de 15 células; con 1.5 hay 20 grupos (sobre-fragmentación). Los 11 grupos resultantes corresponden a tipos celulares reconocibles con los marcadores | 0.2, 0.5, 0.8, 1.0 y 1.5. Con 0.2 hay 9 grupos y silueta similar (0.351), pero se fusionan dos grupos que con 0.5 sí se separan |
|   | 15 | Número de vecinos | Es el valor por defecto de scanpy | No se probaron otros valores |
|   | Ninguna | Regresión de covariables técnicas | Ningún PC está dominado por técnica: la correlación más alta con cuentas totales, % mt o puntaje de doblete fue 0.39 (PC1 con cuentas totales). En los grupos finales las medianas de % mt (3.57 a 5.91 %) y de doblete (0.01 a 0.07) no indican grupos técnicos | Correlación de PC1 a PC10 con cuentas totales, % mt y puntaje de doblete, y medianas de esas covariables por grupo |



## Referencias usadas 📑 :
- Luecken MD, Theis FJ. Current best practices in single-cell RNA-seq analysis: a tutorial. Mol Syst Biol. 2019 Jun 19;15(6):e8746. doi: 10.15252/msb.20188746. PMID: 31217225; PMCID: PMC6582955.
