# %% Rutas (único lugar donde se declaran)
RUTA_DATOS = "data/proyecto-equipo-C"   # matrix.mtx, genes.tsv, barcodes.tsv

# %% Semilla y paquetes
# --- Librerías estándar ---
import random                  # generación de números aleatorios (útil para fijar semillas)

# --- Configuración del backend de matplotlib ---
import matplotlib
try:
    # get_ipython() solo existe dentro de entornos interactivos 
    # Si existe, se mantiene el backend por defecto y las gráficas salen en el panel.
    get_ipython()
except NameError:
    # Si falla, estamos ejecutando desde la terminal 
    # "Agg" es un backend sin ventanas: las figuras no se muestran, solo se guardan en archivo.
    matplotlib.use("Agg")

# --- Librerías de análisis de datos ---
import numpy as np             # operaciones numéricas y arreglos
import pandas as pd            # manejo de tablas (DataFrames)
import matplotlib.pyplot as plt  # creación de gráficas (importar DESPUÉS de fijar el backend)

# --- Análisis de scRNA-seq ---
import scanpy as sc            # análisis de datos de célula única (QC, normalización, clustering, UMAP)

# --- Utilidades de sistema ---
from pathlib import Path       # manejo de rutas de archivos de forma portable 

# --- Reproducibilidad ---
SEMILLA = 20260903             # semilla fija
random.seed(SEMILLA)           # fija la semilla del módulo random de Python
np.random.seed(SEMILLA)        # fija la semilla de NumPy 

# ============================ Carpeta de salida ======================================
# Crea la carpeta "figuras" si no existe; exist_ok=True evita error si ya está creada
Path("figuras").mkdir(exist_ok=True)

def guardar(fig, nombre):
    """Guarda la figura en figuras/ y la muestra (o la cierra si no hay pantalla)."""
    # Guarda la figura como PNG a 150 dpi; bbox_inches="tight" recorta márgenes en blanco
    fig.savefig(f"figuras/{nombre}.png", dpi=150, bbox_inches="tight")

    # Si el backend es "Agg" (terminal, sin ventanas), no hay nada que mostrar:
    # se cierra la figura para liberar memoria
    if matplotlib.get_backend().lower() == "agg":
        plt.close(fig)
    # Si hay entorno interactivo, se muestra en el panel
    else:
        plt.show()

# %% 1. Importar
# ====================================== Importar los datos ======================================
# La carpeta 10x guarda genes x células; AnnData guarda células x genes -> scanpy transpone al leer.
adata = sc.read_10x_mtx(RUTA_DATOS, var_names="gene_symbols",
                        make_unique=True, cache=False)
print(adata)

# %% P1 · ¿Cuántos genes y cuántas células?
# ========================== Pregunta 01 ===================================================
n_celulas, n_genes = adata.shape
print(f"{n_celulas} células (filas) x {n_genes} genes (columnas)")
print("Primeros barcodes:", list(adata.obs_names[:3]))
print("Primeros genes:   ", list(adata.var_names[:3]))

# %% P2 · ¿Cuentas crudas? (enteros no negativos)
# ========================== Pregunta 02 ===================================================

X = adata.X
print(pd.Series({
    "minimo": X.min(),
    "maximo": X.max(),
    "hay_decimales": bool(np.any(X.data % 1 != 0)),
    "hay_negativos": bool(np.any(X.data < 0)),
}, dtype=object).to_string())

# %% P3 · ¿Lote o condición? (solo hay el sufijo del barcode)
# ========================== Pregunta 03 ===================================================

barcodes = pd.Series(adata.obs_names)
print("Columnas de obs:", list(adata.obs.columns))
print("Longitud del barcode:\n",
      barcodes.str.split("-").str[0].str.len().value_counts().to_string())
print("Sufijos:\n", barcodes.str.split("-").str[-1].value_counts().to_string())

# %% P4 · ¿raw o filtered? (curva de rango de barcodes)
# ========================== Pregunta 04 ===================================================

totales = np.asarray(X.sum(axis=1)).ravel()
print(pd.Series({
    "barcodes": len(totales),
    "con_0_cuentas": int((totales == 0).sum()),
    "con_menos_de_100_cuentas": int((totales < 100).sum()),
    "mediana_cuentas": float(np.median(totales)),
}, dtype=object).to_string())

rango = np.sort(totales)[::-1]
fig, ax = plt.subplots(figsize=(5, 4))
ax.loglog(np.arange(1, len(rango) + 1), np.maximum(rango, 1))
ax.set(xlabel="Rango del barcode (log)", ylabel="Cuentas totales (log)",
       title="Curva de rango de barcodes")
guardar(fig, "01_curva_rango")

# %% Anotación de genes
ids = adata.var["gene_ids"].astype(str)
print("Prefijos de ID:\n", ids.str[:4].value_counts().head().to_string())
print("¿IDs con versión?:", bool(ids.str.contains(r"\.\d+$").any()))
print("Genes MT-:", int(adata.var_names.str.upper().str.startswith("MT-").sum()))

# %% 3. Forma de la matriz (práctica de la sesión 01)
# Con una matriz raw, esta fracción de ceros se acerca a 1 porque casi todos
# los barcodes son gotas vacías; hay que recalcularla después de llamar células.
prop_ceros = 1 - X.nnz / (X.shape[0] * X.shape[1])
dispersa_MB = (X.data.nbytes + X.indices.nbytes + X.indptr.nbytes) / 1e6
densa_MB = X.shape[0] * X.shape[1] * 8 / 1e6
print(f"Fracción de ceros: {prop_ceros:.6f}")
print(f"Tipo {type(X).__name__} | dispersa {dispersa_MB:.1f} MB | densa (no se crea) {densa_MB:.1f} MB")

# No se convierte toda la matriz a densa: solo esquinas pequeñas.
# Se toman las 5 células con más cuentas y los 5 genes más expresados.
suma_genes = np.asarray(X.sum(axis=0)).ravel()
top = np.argsort(suma_genes)[::-1][:5]
mejores = np.argsort(totales)[::-1][:5]
print(pd.DataFrame(X[mejores][:, top].toarray(),
                   index=adata.obs_names[mejores], columns=adata.var_names[top]).to_string())

# %% 4. Métricas de calidad por barcode
# Se calculan y se miran; todavía NO se filtra nada.
# Genes mitocondriales: "MT-" en humano, "mt-" en ratón (upper() cubre ambos).
# El % mitocondrial es un indicador de calidad, no un resultado biológico.
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                           log1p=False, inplace=True)
cols = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
print(adata.obs[cols].describe().round(2).to_string())
# pct_counts_mt queda vacío (NaN) en los barcodes con 0 cuentas: es esperado.

# %% Histogramas de TODOS los barcodes (escala log)
fig, axes = plt.subplots(1, 2, figsize=(10, 3.5))
axes[0].hist(np.log10(adata.obs["total_counts"] + 1), bins=80)
axes[0].set(xlabel="log10(cuentas totales + 1)", ylabel="barcodes")
axes[1].hist(np.log10(adata.obs["n_genes_by_counts"] + 1), bins=80)
axes[1].set(xlabel="log10(genes detectados + 1)")
guardar(fig, "02_hist_todos_los_barcodes")

# %% Solo gotas con señal, SOLO PARA VISUALIZAR (no es un filtro)
UMBRAL_VISUAL = 100   # arbitrario: sirve para ver mejor, no para decidir qué es célula
vis = adata.obs[adata.obs["total_counts"] >= UMBRAL_VISUAL]
print("Barcodes con >=", UMBRAL_VISUAL, "cuentas:", len(vis))
print(vis[cols].describe().round(2).to_string())

fig, axes = plt.subplots(1, 3, figsize=(13, 3.5))
for ax, col in zip(axes, cols):
    ax.hist(vis[col].dropna(), bins=60)
    ax.set(title=col, xlabel=col, ylabel="barcodes")
plt.tight_layout()
guardar(fig, "03_hist_qc_visual")

# %% Cuentas totales vs genes detectados, coloreado por % mitocondrial
fig, ax = plt.subplots(figsize=(5.5, 4.5))
pts = ax.scatter(vis["total_counts"], vis["n_genes_by_counts"],
                 c=vis["pct_counts_mt"], s=4, cmap="viridis")
ax.set(xscale="log", yscale="log", xlabel="Cuentas totales (log)",
       ylabel="Genes detectados (log)")
fig.colorbar(pts, label="% mitocondrial")
guardar(fig, "04_cuentas_vs_genes")