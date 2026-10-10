# %% Rutas (único lugar donde se declaran)
RUTA_DATOS = "/Users/monicareyes/Desktop/proyecto-equipo-C"   # matrix.mtx, genes.tsv, barcodes.tsv

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
import scrublet as scr          # detección de dobletes

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

# %% 5. Umbral de cuentas: histograma con candidatos
tc = adata.obs["total_counts"]
candidatos = [200, 300, 400, 500]

fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
axes[0].hist(np.log10(tc[tc > 0]), bins=120)
axes[0].set(xlabel="log10(cuentas totales)", ylabel="barcodes", yscale="log",
            title="Todos los barcodes con cuentas > 0")
axes[1].hist(np.log10(tc[tc >= 100]), bins=120)
axes[1].set(xlabel="log10(cuentas totales)", ylabel="barcodes",
            title="Zoom: barcodes con >= 100 cuentas")
for ax in axes:
    for u in candidatos:
        ax.axvline(np.log10(u), color="red", ls="--", lw=0.8)
plt.tight_layout()
guardar(fig, "05_hist_filtro_cuentas")

# %% Qué pasa con cada candidato
mt = adata.var["mt"].to_numpy()
filas = []
for u in candidatos:
    m = (tc >= u).to_numpy()         # máscara en numpy, no Serie
    sub = adata.X[m]                  # matriz dispersa solo de esos barcodes
    filas.append({
        "umbral_cuentas": u,
        "barcodes": int(m.sum()),
        "genes_detectados": int((np.asarray(sub.sum(axis=0)).ravel() > 0).sum()),
        "genes_mt_detectados": int((adata.X[m][:, adata.var["mt"].to_numpy()].sum(axis=0) > 0).sum()),
        "mediana_cuentas": float(np.median(tc.to_numpy()[m])),
        "mediana_genes_detectados": float(np.median(adata.obs["n_genes_by_counts"].to_numpy()[m])),
        "mediana_pct_mt": float(np.median(adata.obs["pct_counts_mt"].to_numpy()[m])),
    })
print(pd.DataFrame(filas).round(2).to_string(index=False))

# %% 6. aplicar el filtro de cuentas (solo barcodes con cuentas >= 500)
# aun no se filtran genes, solo barcodes

umbral_filtro = 500
cel = (adata.obs["total_counts"] >= umbral_filtro).to_numpy()  # máscara en numpy, no Serie
adata = adata[cel].copy()  # se hace copia para no modificar el original
print(f"Filtradas {cel.sum()} células con >= {umbral_filtro} cuentas")

# %% 7. Umbrales de genes y %mt 
def limites_mad (x, k):
    """Devuelve los límites inferior y superior de x según k veces la MAD."""
    mediana = np.median(x)
    mad = np.median(np.abs(x - mediana))
    return mediana - k * mad, mediana + k * mad

obs = adata.obs
metricas = {
    "log1p_total_counts": np.log1p(obs["total_counts"]).to_numpy(),
    "log1p_n_genes_by_counts": np.log1p(obs["n_genes_by_counts"]).to_numpy(),
    "pct_counts_mt": obs["pct_counts_mt"].to_numpy(),
}

limites = []
for nombre, x in metricas.items():
    for k in [3,4,5]:
        lo, hi = limites_mad(x, k)
        fuera= int (((x < lo) | (x > hi)).sum())
        limites.append({
            "metrica": nombre,
            "k": k,
            "limite_inferior": lo,
            "limite_superior": hi,
            "barcodes_fuera": fuera,
        })
print(pd.DataFrame(limites).round(2).to_string(index=False))

# %% graficas con limites candidatos (k = 3, 5)
fig, axes = plt.subplots(1, 3, figsize=(13, 3.5))
x_g = obs["n_genes_by_counts"].to_numpy()
axes[0].hist(x_g, bins=60)
axes[0].set(title="Genes detectados", xlabel="n_genes_by_counts", ylabel="barcodes")
for k in [3, 5]:
    lo, hi = limites_mad(x_g, k)
    axes[0].axvline(lo, color="red", ls="--", lw=0.8)
    axes[0].axvline(hi, color="red", ls="--", lw=0.8)
x_mt = obs["pct_counts_mt"].to_numpy()
axes[1].hist(x_mt, bins=60)
axes[1].set(title="% mitocondrial", xlabel="pct_counts_mt", ylabel="barcodes")
for k in [3, 5]:
    lo, hi = limites_mad(x_mt, k)
    axes[1].axvline(lo, color="red", ls="--", lw=0.8)
    axes[1].axvline(hi, color="red", ls="--", lw=0.8)
x_tc = obs["total_counts"].to_numpy()
axes[2].hist(x_tc, bins=60)
axes[2].set(title="Cuentas totales", xlabel="total_counts", ylabel="barcodes")
for k in [3, 5]:
    lo, hi = limites_mad(x_tc, k)
    axes[2].axvline(lo, color="red", ls="--", lw=0.8)
    axes[2].axvline(hi, color="red", ls="--", lw=0.8)
plt.tight_layout()
guardar(fig, "06_hist_qc_mad")

# %% 8 que tan sensibles son los resultados a los umbrales de QC
obs = adata.obs
g = np.log1p(obs["n_genes_by_counts"].to_numpy())
mt = obs["pct_counts_mt"].to_numpy()

def lim(x, k):
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    return med - k * mad, med + k * mad

limites_qc = []
for k in [3, 4, 5]:
    g_lo, g_hi = lim(g, k)
    mt_lo, mt_hi = lim(mt, k)
    fuera = ((g < g_lo) | (g > g_hi) | (mt < mt_lo) | (mt > mt_hi)).sum()
    limites_qc.append({
        "k": k,
        "barcodes_fuera": int(fuera),
        "porcentaje_fuera": float(fuera / len(obs) * 100),
    })
print(pd.DataFrame(limites_qc).round(2).to_string(index=False))

# %% 9 aplicar los filtros de calidad sobre adata (barcodes con >= 500 cuentas)
obs = adata.obs
g = np.log1p(obs["n_genes_by_counts"].to_numpy())
mt = obs["pct_counts_mt"].to_numpy()

g_min = lim(g, 4)[0]    # genes, inferior: k = 4
g_max = lim(g, 5)[1]    # genes, superior: k = 5
mt_max = lim(mt, 5)[1]  # % mt, superior: k = 5

f_g_bajo, f_g_alto, f_mt_alto = (g < g_min), (g > g_max), (mt > mt_max)
print(f"Genes: {np.expm1(g_min):.0f} a {np.expm1(g_max):.0f} | % mt máx: {mt_max:.2f}")
print(f"Fuera por pocos genes: {f_g_bajo.sum()} | por muchos genes: {f_g_alto.sum()} | por % mt: {f_mt_alto.sum()}")

mantener = ~(f_g_bajo | f_g_alto | f_mt_alto)
print(f"Antes del filtro: {adata.n_obs} barcodes")
adata_f = adata[mantener].copy()     # adata (4,088 barcodes) queda intacto
print(f"Después del filtro: {adata_f.n_obs} barcodes ({100 * adata_f.n_obs / adata.n_obs:.1f} %)")

# %% 10. filtrar genes: cuantos sobreviven con un mínimo de células
n_cel_por_gen = np.asarray((adata_f.X > 0).sum(axis=0)).ravel()
tabla_genes = pd.DataFrame({
    "min_celulas": [1, 3, 5, 10, 20],
    "genes_retenidos": [int((n_cel_por_gen >= m).sum()) for m in [1, 3, 5, 10, 20]],
})
tabla_genes["porcentaje_retenido"] = 100 * tabla_genes["genes_retenidos"]  / adata_f.n_vars
tabla_genes["genes_eliminados"] = adata_f.n_vars - tabla_genes["genes_retenidos"]
print(tabla_genes.to_string(index=False))

# %% grafica para saber cuantas celuas dectadas por gen
fig, ax = plt.subplots(figsize=(5, 4))
ax.hist(n_cel_por_gen, bins=60)
ax.set(xlabel="Células detectadas por gen", ylabel="genes", yscale="log", title="Genes detectados en >= 1 célula")
for m in [1 , 3, 5, 10, 20]:
    ax.axvline(m, color="red", ls="--", lw=0.8)
guardar(fig, "07_hist_genes_detectados")

# %% 11. filtrar genes: mínimo de 3 células
min_celulas = 3
genes_a_mantener = (n_cel_por_gen >= min_celulas)
print(f"Genes antes del filtro: {adata_f.n_vars}")
adata_f = adata_f[:, genes_a_mantener].copy()
print(f"Genes después del filtro: {adata_f.n_vars} ({100 * adata_f.n_vars / adata.n_vars:.1f} %)")  

# %%12 deteccion de dobletes con scrublet
# Un doblete es una gota que contiene dos células, lo que puede confundir el análisis. Scrublet estima la probabilidad de que cada célula sea un doblete.
# esto lo hace con un modelo de simulación y comparación de perfiles de expresión. Se recomienda filtrar los dobletes antes de hacer análisis downstream.

sc.pp.scrublet(adata_f, random_state=SEMILLA)

print(f"Dobletes predichos: {adata_f.obs['predicted_doublet'].sum()} "
      f"de {adata_f.n_obs} ({100 * adata_f.obs['predicted_doublet'].mean():.1f} %)")
print(f"Umbral: {adata_f.uns['scrublet']['threshold']:.3f}")

fig, ax = plt.subplots(figsize=(5, 3.5))
ax.hist(adata_f.obs["doublet_score"], bins=60, log=True)
ax.axvline(adata_f.uns["scrublet"]["threshold"], color="red", ls="--", lw=1)
ax.set(xlabel="Puntaje de doblete", ylabel="barcodes (log)",
       title="Scrublet (línea roja = umbral automático)")
guardar(fig, "08_scrublet")

# %% 12b. Decisión sobre dobletes
# El histograma no muestra dos poblaciones separadas, así que el umbral automático
# no es confiable. Se anotan los puntajes en adata_f.obs pero no se elimina ningún
# barcode.
print("Barcodes con puntaje > 0.366:", (adata_f.obs["doublet_score"] > 0.366).sum())
print("Se conservan todos:", adata_f.n_obs, "barcodes")