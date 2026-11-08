# Mataquito River Watershed Analysis

Code, data and outputs for a study of erosion rates and sediment sources in the Mataquito River watershed, central Chile (~35°S). The study combines basin-averaged in situ cosmogenic ¹⁰Be erosion rates with detrital zircon U-Pb geochronology from 11 nested river-sand samples (CT-1 through CT-11). The samples run along an Andes-to-coast transect from the high Principal Cordillera, through the Central Depression, to the Coastal Cordillera and the Pacific outlet.

The analyses cover:

- watershed delineation
- ¹⁰Be production rates and subwatershed erosion rates
- detrital zircon age distributions and top-down mixture modeling (detritalPy-mix)
- quartz and zircon fertility and mineral-flux ratios
- downstream mixing below the Teno–Lontué confluence
- river-profile (χ, k_sn, erodibility) and catchment-geometry analysis
- subwatershed geology
- a sediment mass balance at the outlet

---

## Associated manuscript

This repository accompanies the submitted manuscript:

> **Combining in situ ¹⁰Be with Detrital Zircon U-Pb Geochronology to Quantify Erosion Rates and Sediment Sources: Mataquito River, Chile**
>
> Grant P. Long¹\*, Felipe Aron², George E. Hilley¹
>
> ¹ Department of Earth and Planetary Sciences, Stanford University, Stanford, California 94305, USA
> ² Departamento de Geología, Universidad de Chile, Santiago, RM 8370450, Chile
>
> \* Corresponding author: glong1@stanford.edu · ORCID [0009-0008-8330-9316](https://orcid.org/0009-0008-8330-9316)

Citation metadata for this code is in [`CITATION.cff`](CITATION.cff). The [Supplemental material](#supplemental-material) section lists the manuscript's supplemental figures and tables.

---

All analysis code and data live in the [`Mataquito/`](Mataquito/) folder. **Unless stated otherwise, file paths in this README are relative to `Mataquito/`.**

## Directory structure

```
Mataquito/
├── README.md
├── data/
│   ├── sample_metadata/        MataquitoSampleData.xlsx (canonical ¹⁰Be sample table), Quartz mass fractions.xlsx, mataquito_area_km2.xlsx
│   ├── dem/                    fd_utm30m, area_utm30m, CT_{1,5,6,7}_mask (committed); dem_utm30m.tif and failure-index GeoTIFFs (external)
│   ├── geochronology/
│   │   └── zircon_upb/         Per-sample U-Pb lab data tables + combined CT-1-11_ZrUPb_datasheets.xlsx
│   ├── geology/                subwatershed_geologic_units_data.xlsx, combined_subwatershed_summary.xlsx,
│   │   ├── Geologic_Map/       colored geologic-map shapefile and PNG (NB09)
│   │   └── landslides/         2010 Maule earthquake landslide inventory (reference data)
│   ├── hydrology/              DGA discharge/precipitation .xls, CHELSA precipitation GeoTIFF (external)
│   │   └── CAMELS-CL/          camels_cl_<gauge>/ daily discharge and catchment attributes
│   └── shapefiles/
│       ├── sample_points/      flow-corrected sample points, subwatersheds30m.shp
│       ├── flowlines/          corrected_flowlines.shp (NB07)
│       ├── geology/            1:1,000,000 geologic map shapefiles (units, faults, volcanic centers)
│       ├── slope/              CT-1 / CT-7 low-slope (<5°, <10°) rasters and polygons
│       └── *.shp, *_Slope.tif  Central Depression polygon, adjusted CT-5/CT-6 areas, thrusts, slope rasters (external)
├── notebooks/                  01–10 analysis notebooks (+ LaTeX tables written by NB09)
├── mataquito/                  Shared Python package and batch mixing scripts
├── tests/                      pytest unit tests for the package
├── figures/
│   ├── manuscript/             Manuscript figures
│   ├── zircon/                 DZ age distributions (NB04)
│   │   └── mixing/             detritalPy-mix figures (NB05 / batch script)
│   ├── fertility/              Fertility/flux violin and distribution plots (NB06)
│   ├── mixing/                 Downstream confluence mixing plots (NB07)
│   └── mass_balance/           Mass-balance figures (NB10, optional)
├── results/
│   ├── detritalpy_output/      detritalPy CSV/KML exports (NB04)
│   ├── mixing_coefficients/    detritalPy-mix workbooks, summaries and bootstrap CSVs (NB05 / batch script)
│   ├── fertility/              Fertility and flux tables (NB06)
│   ├── mass_balance/           Area accounting and budget summary (NB10)
│   └── geology_tables/         Geologic summary tables (NB09)
└── supplemental/               Manuscript supplemental figures S1–S4 and tables S1–S9
```

Some folders are created only when the corresponding save lines are enabled: `figures/river_profiles/` and `results/river_profiles/` (NB08), and `figures/geology/` (NB09).

---

## Contents

1. [Directory structure](#directory-structure)
2. [Supplemental material](#supplemental-material)
3. [Study area and drainage network](#study-area-and-drainage-network)
4. [Setup](#setup)
5. [Workflow and run order](#workflow-and-run-order)
6. [Notebook reference](#notebook-reference)
7. [Batch mixing scripts](#batch-mixing-scripts)
8. [The `mataquito` package](#the-mataquito-package)
9. [Data](#data)
10. [Tests](#tests)
11. [License](#license)

---

## Supplemental material

The supplemental figures and tables for the manuscript are in [`Mataquito/supplemental/`](Mataquito/supplemental/). Three tables (S6, S7 and S8) are written directly by a notebook. The others are copies of a notebook output or of a data file in this repository, as listed below.

| Item | Contents | Produced by / source |
|------|----------|----------------------|
| `Figure_S1.pdf` | Geologic map of the Mataquito watershed with the unit legend | NB09, geologic-map section (1:1,000,000 geology) |
| `Figure_S2.pdf` | Quartz and zircon fertility distributions (log₁₀) for the CT-8.2 & CT-3\* → CT-2 and CT-8.2 & CT-3\* → CT-9 scenarios | NB06, Section 13 (`dist_*` plots) |
| `Figure_S3.pdf` | Quartz and zircon fertility distributions (log₁₀) for the CT-5 & CT-6 → CT-10 / CT-11 / CT-8.2 scenarios | NB06, Section 13 (`dist_*` plots) |
| `Figure_S4.pdf` | detritalPy-mix results for CT-8.2 & CT-3 → CT-2: parent and child age distributions and mixing-proportion violin plot | NB05 / batch mixing script (`figures/zircon/mixing/`) |
| `Table_S1.xlsx` | Geologic-unit areas (km²) and percentages for each subwatershed | NB09 (`combined_subwatershed_summary.xlsx`) |
| `Table_S2.xlsx` | Quartz mass fractions of the 250–500 µm fraction used for ¹⁰Be | `data/sample_metadata/Quartz mass fractions.xlsx` |
| `Table_S3.csv` | Inter-sample comparison matrices (similarity, likeness, K-S, Kuiper, cross-correlation, …) | NB04, "Export sample comparison matrices" (`results/detritalpy_output/interSampleComparisonKDE.csv`) |
| `Table_S4.csv` | Zircon fertility and zircon flux ratios for every scenario and both mixing metrics (median, IQR, linear and log₁₀) | NB06, Section 15 (`results/fertility/zircon_fertility_and_flux_by_metric.csv`) |
| `Table_S5.xlsx` | Age-group definitions (min/max age and geologic justification) used to bin zircon ages | Compiled for the manuscript; used for the NB04 age-proportion bins |
| `Table_S6.xlsx` | Best-fit mixing coefficients (V_max and r²-KDE) for each parent pair and child | **Written by NB05**, last section |
| `Table_S7.xlsx` | Part A: sediment budget at CT-9 (PC = CT-8.2). Part B: quartz and zircon mineral-flux ratios | **Written by NB10**, Section 8 |
| `Table_S8.xlsx` | k_sn, erosion rate and erodibility K for the Teno/Lontué outlets, with θ, n and lowland sensitivity | **Written by NB08**, Section 4.9 |
| `Table_S9.xlsx` | Sample locations, dates, drainage areas, grain size, ¹⁰Be concentrations and erosion rates | `data/sample_metadata/MataquitoSampleData.xlsx` |

---

## Study area and drainage network

The 11 sample sites form a nested hierarchy, so downstream sites integrate sediment from their upstream tributaries:

```
CT-7 ──→ CT-5 ──┐
                 ├──→ CT-4 ──→ CT-10 ──→ CT-11 ──→ CT-8.2 ──→ CT-2 ──┐
CT-1 ──→ CT-6 ──┘                                            ├──→ CT-9
                                                CT-3 ─────────┘
```

- **Teno branch:** CT-7 (upper Teno) → CT-5 (lower Teno).
- **Lontué branch:** CT-1 (Río Claro / upper Lontué) → CT-6 (lower Lontué).
- **Principal Cordillera children** below the Teno–Lontué confluence: CT-10, CT-11 and CT-8.2.
- **Central Depression / Coastal Cordillera:** CT-2. CT-3 is a small Coastal Cordillera tributary that also serves as the lithologic proxy "CT-3\*" for the coastal part of the basin.
- **Outlet:** CT-9, the Mataquito near the coast.
- **CT-4** lies just below the confluence. It is used in several notebooks as an exploratory scenario.
- **Naming:** CT-8.2 is stored as `CT-8` in `MataquitoSampleData.xlsx` and in `mataquito.sample_data`. The zircon data and figures use `CT-8.2`.

---

## Setup

### 1. Python environment

All analyses were run with Python 3.11. A conda environment is the easiest way to get the geospatial stack:

```bash
conda create -n mataquito -c conda-forge python=3.11 numpy pandas scipy matplotlib openpyxl xlrd xlsxwriter seaborn scikit-learn statsmodels geopandas shapely rasterio pyproj geopy matplotlib-scalebar numba pysheds folium simplekml peakutils kdepy jupyterlab pytest
```

```bash
conda activate mataquito
```

```bash
pip install detritalpy
```

| Package(s) | Needed by |
|------------|-----------|
| numpy, pandas, scipy, matplotlib, openpyxl | Everything |
| detritalpy (≥ 1.4.6; developed with 1.4.14), KDEpy, scikit-learn, folium, simplekml, peakutils, statsmodels | NB04, NB05, batch mixing scripts |
| geopandas, shapely, rasterio, pyproj | NB01, NB07, NB08, NB09 |
| geopy | NB07 |
| pysheds, numba, scikit-learn, xlrd | NB08 |
| seaborn, matplotlib-scalebar, xlsxwriter | NB09 |
| **TopoAnalysis** (not on PyPI) | NB01, NB08, `mataquito.river_profiles` |
| pytest | `tests/` |

**TopoAnalysis** is the Stanford Tectonic Geomorphology Lab DEM toolkit: <https://github.com/STGL/TopoAnalysis>. Clone it and install it into the same environment:

```bash
git clone https://github.com/STGL/TopoAnalysis.git
```

```bash
pip install -e TopoAnalysis
```

### 2. Get the code and external data

```bash
git clone <this repository>
```

```bash
cd MataquitoWatershedAnalysis/Mataquito
```

GeoTIFFs (`*.tif`) are gitignored, so the large rasters are not in the repository. Put them where [Data → External files](#external-files-not-in-the-repository) says before running NB01, NB07, NB08 or NB09. Every other notebook runs from files already committed.

### 3. Launch the notebooks from `notebooks/`

Every notebook finds the package and the data relative to its own folder. They use `sys.path.insert(0, '..')`, `Path.cwd().parent` or `'../data/...'`, so start Jupyter from the `notebooks/` directory or open the notebooks with that as the working directory:

```bash
cd notebooks
```

```bash
jupyter lab
```

The batch mixing scripts and the tests run from `Mataquito/` instead (see their sections).

---

## Workflow and run order

Each notebook is self-contained. The diagram shows which notebooks feed others through files on disk. Every intermediate file is already committed, so you can open any notebook on its own. You only need to follow the arrows if you want to regenerate everything from scratch.

```
NB01 DEM delineation ──► data/dem/ (fd, area, CT_*_mask) ──► NB07 confluence mixing
        (external DEM)                                    └─► NB08 river profiles / geometry

NB02 production rates      (standalone)
NB03 subwatershed rates    (standalone)

NB04 zircon age distributions ──► results/detritalpy_output/   (Table S3)

NB05 zircon mixing  ─┐ (interactive, one scenario at a time)
batch mixing script ─┴► results/mixing_coefficients/*.xlsx, *.txt, *.csv   (Table S6, Fig. S4)
                              │
                              ├─► NB06 fertility & flux ──► results/fertility/   (Table S4, Figs. S2–S3)
                              │                                   │
                              │                                   └─► NB10 mass balance ──► results/mass_balance/  (Table S7)
                              └─► NB07 (mixing length scales)

NB09 subwatershed geology  ──► results/geology_tables/   (Table S1, Fig. S1)
```

A full regeneration runs in this order: **01 → 02 → 03 → 04 → batch mixing script (or 05) → 06 → 07 → 08 → 09 → 10.**

**Saving outputs.** Several notebooks are set up not to overwrite committed outputs:

- **NB06** has `SAVE_FIGS = False` in Section 1. Set it to `True` to write every figure PDF. Each PDF embeds a timestamp, so rewriting unchanged figures would otherwise touch every file.
- **NB08 and NB10** have every figure and CSV save line commented out. Uncomment the ones you need. The supplemental tables (S7, S8) are always written.
- **NB04** detritalPy plot calls use `savePlot=False`. Set `savePlot=True` to write them.

---

## Notebook reference

All paths below are relative to `Mataquito/`.

### 01 — DEM watershed delineation
`notebooks/01_dem_watershed_delineation.ipynb`

- **Purpose:** Fill the 30 m DEM, compute D8 flow direction and drainage area with TopoAnalysis, and delineate a basin mask and area for each of the 11 sample points. The masks are also converted to polygon shapefiles for QGIS.
- **Inputs:** `dem_utm30m.tif` (external, UTM 19S). Sample outlet coordinates are defined in the notebook (EPSG:32719).
- **Configure:** set `DEM_DATA_DIR` (the folder holding the DEM) and `MASK_DIR` (where masks are written) in the two `EXTERNAL DATA PATHS` / `EXTERNAL OUTPUT PATH` cells. The notebook `os.chdir`s into these folders.
- **Outputs** (written into `DEM_DATA_DIR` / `MASK_DIR`):
  - `dem_utm_filled30m`, `fd_utm30m`, `area_utm30m` and `logarea_utm30m` grids
  - one `CT_<n>_mask` per sample
  - `*_polygon.shp` basin polygons

  The `fd_utm30m`, `area_utm30m` and `CT_*_mask` grids that later notebooks use are committed in `data/dem/`.
- **Runtime:** filling and routing the full DEM takes several minutes.

### 02 — Production rates
`notebooks/02_production_rates.ipynb`

- **Purpose:** ¹⁰Be surface production rates at each sample site from Stone (2000) latitude/altitude scaling. It uses the 30–40° latitude coefficients, elevation → pressure, and an SLHL reference of 5.1 atoms g⁻¹ yr⁻¹. It also computes apparent production rates back-calculated from each sample's measured concentration and erosion rate (P = N·E·ρ/Λ).
- **Inputs:** `data/sample_metadata/MataquitoSampleData.xlsx`, via `mataquito.production_rates` and `mataquito.sample_data`.
- **Outputs:** printed tables only.

### 03 — Subwatershed erosion rates
`notebooks/03_erosion_rate_uncertainties.ipynb`

- **Purpose:** Deterministic mass balance through the drainage network, E_sub = (E_down·A_down − Σ E_up·A_up) / A_sub. It walks through CT-5, CT-6, CT-10, CT-11, CT-8, CT-2 and CT-9 and shows the flux breakdown at each node. Negative subwatershed rates indicate net storage or deposition.
- **Inputs:** `MataquitoSampleData.xlsx` and the network topology in `mataquito.flow_network`.
- **Outputs:** printed tables and an in-notebook figure.

### 04 — Zircon age distributions
`notebooks/04_zircon_age_distributions.ipynb`

- **Purpose:** detritalPy plots and statistics for all 11 samples:
  - CDFs and KDEs (bandwidth 2 Myr) over 0–25 Ma and 0–400 Ma
  - age-proportion bar graph
  - interactive sample map
  - maximum depositional ages
  - multidimensional scaling (MDS, stress and Shepard plots, QQ and heat-map matrices)
  - inter-sample comparison matrices
  - distribution and age/error exports
  - a critical-D test on the K-S statistics
- **Inputs:** `data/geochronology/zircon_upb/CT-1-11_ZrUPb_datasheets.xlsx`, which combines the per-sample lab data tables in the same folder.
- **Local kernels:** the notebook loads `mataquito/detritalFuncs_local.py` and swaps detritalPy's adaptive-KDE kernels for `mataquito/adaptiveKDE_local.py` (see [package](#the-mataquito-package)).
- **Outputs:**
  - `figures/zircon/DZageDistributions0-25Ma.pdf` and `figures/zircon/DZageDistributions0-400Ma.pdf`.
  - CSV exports (`ageProportionsCT-1-11.csv`, `MDA.csv`, `interSampleComparison*.csv`, `Dist_to_CSV.csv`, `AgeError_toCSV.csv`) are written by detritalPy to `notebooks/Output/`.
  - The copies used downstream and in Table S3 are in `results/detritalpy_output/`. The critical-D cell reads `results/detritalpy_output/interSampleComparisonKDE.csv`.

### 05 — Zircon mixing (detritalPy-mix)
`notebooks/05_zircon_mixing.ipynb`

- **Purpose:** Top-down (forward) mixture modeling of detrital zircon ages with bootstrapped uncertainty (Sharman & Johnstone, detritalPy-mix). Each child is modeled as a mixture of two parents. The notebook steps through:
  - best-fit mixture
  - bootstrapped mixing coefficients
  - age-distribution, objective-function, distribution-comparison and spatial plots
  - Excel export
- **Choosing a scenario:** edit `parent_list` / `child_list` in Section II. Every scenario used in the study is listed there, commented out. Set `objective_metric` (`'vmax'` or `'r2-kde'`), `bw` and `nBootstrapIterations` in the cells that follow.
- **Iterations:** the notebook defaults to `nBootstrapIterations = 100` (about a minute). The workbooks used downstream were produced with 10,000 iterations by the [batch mixing script](#batch-mixing-scripts).
- **Outputs:**
  - `results/mixing_coefficients/mix_coeffs_all_<parents>_to_<child>_<metric>.csv`
  - `results/mixing_coefficients/dPy-mix_<parents>_to_<child>_results_<metric>.xlsx`
  - figures in `figures/zircon/mixing/`
- **Table S6:** the last section reads the `dPy-mix_*_summary_<metric>.txt` files in `results/mixing_coefficients/` and writes `supplemental/Table_S6.xlsx`.

### 06 — Quartz and zircon fertility and flux
`notebooks/06_fertility_analysis.ipynb`

- **Purpose:** Monte Carlo propagation that combines flux-ordered ¹⁰Be erosion rates with the detritalPy-mix bootstrap coefficients. It reports four ratios per parent/child scenario for both metrics (`vmax`, `r2kde`):
  - quartz fertility Q_A/Q_B
  - quartz flux ratio F_QA/F_QB
  - zircon fertility Z_A/Z_B
  - zircon flux ratio F_ZA/F_ZB

  Sample A is the Principal Cordillera parent (CT-5 or CT-8.2) and B is CT-6 or CT-3\*. Erosion-rate uncertainty is external throughout.
- **Sections:**
  1. Imports and paths
  2. Sample data
  3. Monte Carlo draws
  4. Flux ordering
  5. Flux-ordered rate plots
  6. Equations
  7. Scenario registry (CT-3\* areas)
  8. Mixing results (both metrics)
  9. Coefficient violin plots
  10. Quartz fertility and flux
  11. Zircon fertility and flux
  12. Fertility and flux violin plots
  13. Distribution plots
  14. Flux-ordering sensitivity (rerun on unordered rates)
  15. Summary tables

  The last cell is a preflight check that lists any missing input files.
- **Inputs:**
  - `data/sample_metadata/MataquitoSampleData.xlsx`
  - `results/mixing_coefficients/dPy-mix_*_results_{vmax,r2kde}.xlsx`, both the bootstrap sheet and the `Best-fit` sheet
- **Outputs:**
  - `results/fertility/erosion_rates_flux_ordered.csv`
  - `results/fertility/zircon_fertility_and_flux_by_metric.csv` (Table S4)
  - `results/fertility/mixing_coefficients_both_metrics.csv`
  - `quartz_fertility_and_flux.csv` and `flux_ordering_sensitivity.csv`, whose save lines are commented out in Section 15
  - figure PDFs in `figures/fertility/` when `SAVE_FIGS = True`

### 07 — Downstream confluence mixing
`notebooks/07_downstream_confluence_mixing.ipynb`

- **Purpose:** Trace D8 flowlines from the sample points to the outlet, measure stream distance from the Teno–Lontué confluence and from the outlet, and plot how ¹⁰Be concentration, erosion rate and zircon mixing proportions change downstream (CT-4 → CT-10 → CT-11 → CT-8.2). CT-5 and CT-6 reference bands show whether the downstream samples are well mixed.
- **Inputs:**
  - `data/dem/fd_utm30m` and `data/dem/area_utm30m`
  - `data/shapefiles/sample_points/mataquitosamplepoints_corrected_for_flow_30mDEM.shp`, the sample points snapped to the 30 m flow network
  - `results/mixing_coefficients/mix_coeffs_all_ct-5_ct-6_to_{ct-4,ct-10,ct-11,ct-8.2}_vmax.csv`
  - `MataquitoSampleData.xlsx`
- **Outputs:**
  - `data/shapefiles/flowlines/corrected_flowlines.shp`
  - in `figures/mixing/`: `corrected_confluence_distances_visualization.pdf`, `downstream_mixing_length_scales_filtered.pdf`, `downstream_mixing_Be-10.pdf`, `downstream_mixing_erosionrate.pdf`, `outlet_distances_visualization.pdf`, `downstream_mixing_erosionrate_from_outlet.pdf` and `downstream_mixing_be10concentration_from_outlet.pdf`

### 08 — River profiles and catchment geometry
`notebooks/08_river_profile_analysis.ipynb`

- **Purpose:**
  1. χ–elevation profiles for six Coastal Cordillera esteros near the outlet and for the Teno (CT-7, CT-5) and Lontué (CT-1, CT-6) outlets, giving k_sn and erodibility K from the ¹⁰Be erosion rates.
  2. Sensitivity of k_sn and K to θ, n, Central Depression lowland cells and a block bootstrap.
  3. Teno/Lontué catchment geometry.
  4. Upper-catchment topography of CT-7 versus CT-1.
  5. Teno/Lontué runoff.
- **Sections:**
  1. **Setup.** Paths and every analysis parameter: θ = 0.4, A₀ = 1 km², a 1 km² channel threshold, a 500 m lowland elevation, `N_BOOT = 1000` with seed 42, and others.
  2. **χ–z scatter and concavity scan** over θ = 0.1–1.0.
  3. **Sample-outlet χ–z analysis** at θ = 0.4: watershed extraction, distance–elevation profiles, regression fits, k_sn and K solution. The coastal esteros are assigned the CT-3 erosion rate.
  4. **Teno/Lontué sensitivity:**
     - 4.1 drainage-area check against Table 1
     - 4.2 k_sn in proper units
     - 4.3 θ = 0.35–0.50
     - 4.4 cells above 500 m only
     - 4.5 block bootstrap of k_sn into K
     - 4.6 K for n = 1
     - 4.7 Teno/Lontué ratios
     - 4.8 reproducibility check
     - 4.9 **writes `supplemental/Table_S8.xlsx`**
  5. **Catchment geometry.** The DEM is cropped and re-routed with pysheds, with epsilon depression filling and HAND. This section covers trunk channels, valley-floor width transects, lowland area, the nested CT-1 and CT-7 catchments, hillslope-angle distributions, width functions and planform shape.
  6. **Upper-catchment topography (CT-7 vs CT-1):**
     - hypsometry
     - 5 and 10 km local relief
     - slope distributions, exceedance curves and threshold-slope fraction
     - Gaussian-mixture hillslope population
     - area below 5° and below 10°
     - the Moon et al. (2011) failure index across a transmissivity sweep, using CHELSA and uniform precipitation
     - optional GeoTIFF export of the failure index for QGIS
  7. **Stream discharge and runoff.** CAMELS-CL daily discharge (Alvarez-Garreton et al., 2018) through 2020, extended to 2024 with DGA records. Teno is gauge 7104002; Lontué is the same-day sum of 7112001 and 7115001. This section computes annual and daily runoff and tests the Lontué − Teno difference.
- **Inputs:**
  - `data/dem/`: `dem_utm30m.tif` (external), plus the committed `fd_utm30m`, `area_utm30m` and `CT_{1,5,6,7}_mask`
  - `data/shapefiles/CT_7_Slope.tif` and `CT_1_Slope.tif` (external)
  - `data/hydrology/CHELSA_pr_annual_mm_1981-2010_Mataquito_utm19s.tif` (external)
  - `data/hydrology/CAMELS-CL/camels_cl_{7104002,7112001,7115001}/`
  - `data/hydrology/Caudales Medios Diarios 13_10_2026 15_21.xls` (DGA)
  - `MataquitoSampleData.xlsx`
- **Outputs:**
  - `supplemental/Table_S8.xlsx` is always written.
  - All other figures, CSV tables and failure-index GeoTIFFs have their save lines commented out. They write to `figures/river_profiles/`, `results/river_profiles/` and `data/dem/` respectively when enabled. Create those folders first.
- **Runtime:** this is the longest notebook. The concavity scan prints a time estimate before it runs, and the 1,000-iteration block bootstrap (4.5) and the full-grid relief and failure-index calculations each take several minutes.

### 09 — Subwatershed geology
`notebooks/09_subwatershed_geology.ipynb`

- **Purpose:** Areal proportions of geologic units in each CT subwatershed, chi-squared and Kruskal–Wallis tests across subwatersheds, pie charts, and a geologic map of the watershed (1:1,000,000) with faults, volcanic centers, hillshade and scale bar.
- **Inputs:**
  - `data/geology/subwatershed_geologic_units_data.xlsx`: unit areas per subwatershed, already extracted
  - `data/shapefiles/geology/` (`geo_centro_u19.shp`, `fallas_centro_u19.shp`, `fallas-simb_centro_u19.shp`, `centros-line_centro_u19.shp`)
  - `data/shapefiles/sample_points/subwatersheds30m.shp`
  - `data/shapefiles/hillshademauleregion.tif` (external)
- **Configure:** paths are set in the `EXTERNAL DATA PATHS` cell (cell 3).
- **Outputs:**
  - `results/geology_tables/geologic_area-percents.xlsx` and `results/geology_tables/mata_geologic_units.xlsx`
  - `combined_subwatershed_summary.xlsx` (Table S1), `geologic_summary_complete.tex` and `summary_table_only.tex`, written to the working directory (`notebooks/`). Copies are kept in `results/geology_tables/` and `data/geology/`.
  - `figures/geology/mataquito_geology_1-1000000.{png,pdf}`
  - `data/geology/Geologic_Map/complete_geologic_map.png` and `geo_centro_colored.shp`

### 10 — Mass balance at the outlet
`notebooks/10_mass_balance.ipynb`

- **Purpose:** Volumetric erosional sediment flux (m³/yr) leaving the watershed at CT-9, from the domains that supply it. Note that 1 m/Myr × 1 km² = 1 m³/yr.
  - **Route A (primary):** sum of the Principal Cordillera parents (CT-5 + CT-6), a Central Depression transfer zone (0), and the Coastal Cordillera represented by CT-3\*.
  - **Route B (check):** un-weights the measured downstream rates by quartz fertility. It reduces algebraically to Route A.
  - **Single-sample budgets (comparison):** S1 is the outlet only; S2 is one downstream PC sample plus CT-3\*.
  - **Configurations:** everything is run twice, with CT-10 (`PC_CT10`) or CT-8.2 (`PC_CT82`) as the downstream Principal Cordillera reference sample.
- **Monte Carlo settings:** N = 100,000 draws, seed 17, with the same flux ordering as NB06 (`mataquito.erosion`).
- **Sections:**
  1. Draws and flux ordering
  2. Area accounting
  3–5. Route A, Route B and single-sample budgets
  6. Cross-configuration check
  7. Summary table and key findings
  8. Table S7
  9. Figures (commented out)
  10. Assumptions
- **Inputs:**
  - `MataquitoSampleData.xlsx`
  - for Table S7 Part B, `results/fertility/quartz_fertility_and_flux.csv` and `results/fertility/zircon_fertility_and_flux_by_metric.csv` from NB06
- **Outputs:**
  - `results/mass_balance/area_accounting.csv`
  - `results/mass_balance/mass_balance_summary.csv` (median, IQR, 2.5–97.5 percentile per quantity per configuration)
  - **`supplemental/Table_S7.xlsx`**
  - optional figures in `figures/mass_balance/`

---

## Batch mixing scripts

`mataquito/all_mixing_scenarios_seeded.py` runs detritalPy-mix for every two-parent scenario in the study, without editing the notebook. Each scenario is run with both `vmax` and `r2-kde`, using 10,000 bootstrap iterations by default. This script produced the mixing results used in NB06, NB07, NB10, Table S6 and Figure S4.

The scenarios are:

- **CT-7 & CT-1 →** each of CT-4, CT-10, CT-11, CT-8.2, CT-2 and CT-9
- **CT-5 & CT-6 →** each of CT-4, CT-10, CT-11, CT-8.2, CT-2 and CT-9, plus a combined CT-4/10/11/8.2 run
- **{CT-4, CT-10, CT-11, CT-8.2} & CT-3 →** CT-2
- **{CT-4, CT-10, CT-11, CT-8.2, CT-2} & CT-3 →** CT-9

Run it from `Mataquito/`:

```bash
python mataquito/all_mixing_scenarios_seeded.py
```

```bash
python mataquito/all_mixing_scenarios_seeded.py ct5ct6_to_CT10 --metric vmax -n 100
```

The first command runs every scenario with both metrics. The second is a quick test: one scenario, one metric, 100 iterations.

- **Options:**
  - `--iterations/-n` sets the bootstrap iterations.
  - `--metric {vmax,r2-kde,both}` chooses the metric.
  - `--seed` sets the base seed.
  - Scenario names are positional arguments; see `SCENARIOS` in the script.
- **Run plan:** the script prints the run plan and asks for confirmation before it starts.
- **Reproducibility:** each (scenario, metric) run gets its own seed, derived from an MD5 hash of the base seed, scenario name and metric. Re-running one scenario on its own therefore reproduces its draws from a full batch.
- **Settings:** the age range for the mixing calculation is 0–4500 Myr in 1 Myr steps, with a KDE bandwidth of 2 Myr and 2σ input uncertainties. Perturbed ages are clipped to the evaluation grid so that young grains cannot go negative.
- **Outputs:**
  - in `results/mixing_coefficients/`: `dPy-mix_*_results_<metric>.xlsx` (sheets: bootstrap ensemble per child, `Best-fit`, `Model_fit`), `dPy-mix_*_summary_<metric>.txt`, `mix_coeffs_all_*_<metric>.csv` and `run_log.txt`
  - in `figures/zircon/mixing/`: age-distribution and distribution-comparison PDFs, each at 0–25 Ma and 0–400 Ma, plus objective-function and spatial-comparison plots
- **Runtime:** a full run takes several hours. Progress is appended to `run_log.txt`.

`mataquito/all_mixing_scenarios.py` is an unseeded variant. It uses the installed detritalPy mixer and a split [1.5, 5] Myr bandwidth, and it is not used for the manuscript results.

---

## The `mataquito` package

These are shared functions imported by the notebooks with `sys.path.insert(0, '..')`. Each module can be imported on its own.

| Module | Description | Key functions / contents | Used by |
|--------|-------------|--------------------------|---------|
| `sample_data` | Loads the canonical sample spreadsheet | `load_samples()`, `get_erosion_rates()`, `get_areas()`, `get_uncertainties()`, `get_production_rates()` | NB03, NB06, NB07, NB08, NB10 |
| `production_rates` | Stone (2000) latitude/altitude scaling for ¹⁰Be | `elevation_to_pressure()`, `stone2000_scaling()`, `stone2000_production_rate()` | NB02 |
| `flow_network` | Drainage topology and deterministic subwatershed mass balance | `FLOW_NETWORK`, `TRAVERSAL_ORDER`, `subwatershed_area()`, `all_subwatershed_erosion_rates()` | NB03 |
| `erosion` | Monte Carlo erosion-rate draws and the downstream-increasing-flux (flux-ordering) constraint; returns index-aligned parent/child triples | `generate_mc_samples()`, `flux_order_samples()`, `raw_triples()` | NB06, NB10 |
| `fertility` | Quartz and zircon fertility and mineral-flux ratios, the CT-3\* area convention, and the mixing-results loader | `load_mixing_results()`, `quartz_fertility()`, `zircon_fertility()`, `ct3_star_area()`, `ratio_stats()`, `THRESHOLD` | NB06, NB10 |
| `flowlines` | D8 flowline tracing and along-stream distances | `create_flowlines_from_raster()`, `calculate_stream_distances_from_confluence()`, `calculate_stream_distances_from_outlet()` | NB07 |
| `river_profiles` | χ-elevation, k_sn and erodibility helpers (needs TopoAnalysis) | `extract_watershed()`, `get_ksn()`, `compute_r2()`, `solve_K()`, `fit_and_draw()` | NB08 |
| `detritalFuncs_local` | Local copy of detritalPy's `detritalFuncs` (G. R. Sharman et al.) with the plotting and KDE bandwidth handling used here | — | NB04, NB05, batch scripts |
| `detritalMixer_local` | Local copy of detritalPy-mix's `detritalMixer` (S. A. Johnstone & G. R. Sharman) | — | NB05, batch scripts |
| `adaptiveKDE_local` | Local copy of the Shimazaki & Shinomoto adaptive KDE kernels, patched into `detritalpy.adaptiveKDE` | `sskernel()`, `ssvkernel()` | NB04, NB05, batch scripts |
| `all_mixing_scenarios_seeded`, `all_mixing_scenarios` | Batch runners for detritalPy-mix | See [Batch mixing scripts](#batch-mixing-scripts) | — |

The modules in `tests/` also serve as usage examples.

---

## Data

### Sample metadata

`data/sample_metadata/MataquitoSampleData.xlsx` (Sheet1, 11 samples) is the single source of ¹⁰Be data for every notebook, read through `mataquito.sample_data`. Its columns are:

| Column | Description |
|--------|-------------|
| `Sample_ID` | CT-1 … CT-11 (CT-8.2 stored as `CT-8`) |
| `Sample_Type`, `Latitude`, `Longitude` | Sample description and location (decimal degrees) |
| `Source_Area` | Upstream drainage area (km²), via `get_areas()` |
| `Elevation`, `Virtual_Elevation` | Sample and basin-averaged effective elevation (m) |
| `Surface_Production_Rate` | Basin-averaged ¹⁰Be production rate (atoms g⁻¹ yr⁻¹), via `get_production_rates()` |
| `Shielding_Factor` | Topographic shielding |
| `Quartz`, `9Be carrier (g)`, `10Be/9Be`, `10Be/9Be_uncertainty` | Quartz mass, carrier and AMS ratio |
| `10Be_concentration`, `10Be_concentration_uncertainty` | ¹⁰Be concentration (atoms/g) |
| `Erosion_rate` | Basin-averaged ¹⁰Be erosion rate (m/Myr), via `get_erosion_rates()` |
| `Erosion_rate_uncertainty_internal`, `Erosion_rate_uncertainty_external` | 1σ uncertainties (m/Myr). The external value is used for all Monte Carlo draws, via `get_uncertainties()` |

The full published sample table, with coordinates, dates, grain size and concentrations, is `supplemental/Table_S9.xlsx`.

### Committed data

The repository includes:

- the zircon U-Pb data tables
- the sample spreadsheets
- the D8 flow-direction and drainage-area grids and catchment masks (`data/dem/`)
- all shapefiles (sample points, subwatersheds, flowlines, geologic map)
- the CAMELS-CL gauge folders and DGA discharge files
- every intermediate result in `results/`

### External files (not in the repository)

`*.tif` is gitignored. Obtain or regenerate these files and place them at the paths shown:

| File | Path | Source | Needed by |
|------|------|--------|-----------|
| 30 m DEM, UTM 19S | `data/dem/dem_utm30m.tif` (and the `DEM_DATA_DIR` set in NB01) | Copernicus WorldDEM-30 (see attribution below), projected to WGS 84 / UTM zone 19S (EPSG:32719) | NB01, NB08 |
| Slope rasters (degrees) | `data/shapefiles/CT_7_Slope.tif`, `data/shapefiles/CT_1_Slope.tif` | Derived from the DEM (in part with TopoAnalysis) and clipped to the CT-7 and CT-1 masks | NB08 §6 |
| Mean annual precipitation 1981–2010 | `data/hydrology/CHELSA_pr_annual_mm_1981-2010_Mataquito_utm19s.tif` | CHELSA climatology, regridded to the DEM grid | NB08 §6.8 |
| Hillshade | `data/shapefiles/hillshademauleregion.tif` | Derived from the DEM | NB09 map |
| Low-slope rasters | `data/shapefiles/slope/CT_{1,7}_slope<{5,10}.tif` | Derived from the slope rasters (the `.gpkg` polygons are committed) | Map figures |

**DEM attribution:** This map was produced using Copernicus WorldDEM-30 © DLR e.V. 2010-2014 and © Airbus Defence and Space GmbH 2014-2018 provided under COPERNICUS by the European Union and ESA; all rights reserved. Plotted in coordinate reference system EPSG:32719 - WGS 84 / UTM zone 19S. The hillshade, slope rasters and all other topographic derivatives in this repository were produced from this DEM.

Discharge is from CAMELS-CL (Alvarez-Garreton et al., 2018, <https://camels.cr2.cl>) and the Dirección General de Aguas (DGA) online records. The geologic map shapefiles (central-zone sheets, `*_centro_u19`) are modified from Servicio Nacional de Geología y Minería (2004), 1:1,000,000 map of Chile: Servicio Nacional de Geología y Minería, 2004, Mapa Geológico de Chile: Versión digital: SERNAGEOMIN Publicación Geológica Digital, No. 7 CD-ROM, <https://tiendadigital.sernageomin.cl/es/geologia-basica/3273-mapa-geologico-de-chile-version-digital.html>. `data/geology/landslides/` holds the comprehensive earthquake-induced landslide inventory for the 2010 Chile (Maule) megathrust earthquake, used as reference data.

---

## Tests

From `Mataquito/`:

```bash
pytest tests/ -v
```

| Test file | Covers |
|-----------|--------|
| `test_sample_data.py` | Loading `MataquitoSampleData.xlsx` |
| `test_production_rates.py` | Stone (2000) scaling |
| `test_flow_network.py` | Network topology and subwatershed mass balance |
| `test_fertility.py` | Mixing-results loader, CT-3\* areas, quartz/zircon fertility and flux |
| `test_flowlines.py` | Line distance and substring helpers (needs shapely, rasterio) |
| `test_river_profiles.py` | k_sn, R² and K solution (needs TopoAnalysis) |

If shapely/rasterio or TopoAnalysis are not installed, skip those files:

```bash
pytest tests/ -v --ignore=tests/test_flowlines.py --ignore=tests/test_river_profiles.py
```

---

## License

MIT; see [`LICENSE`](LICENSE). For questions about the code or data, contact Grant P. Long (glong1@stanford.edu).
