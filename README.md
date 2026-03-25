# Mataquito River Watershed Analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19212420.svg)](https://doi.org/10.5281/zenodo.19212420)

Erosion rates, sediment provenance, and mineral fertility in the Mataquito River watershed, central Chile (~35S). Uses cosmogenic 10Be, detrital zircon U-Pb geochronology, quartz/zircon fertility analysis, DEM analysis, downstream confluence mixing, and subwatershed geology.

This repository supports an Andes-to-coast transect study of 11 nested catchment samples (CT-1 through CT-11) across the Mataquito watershed, which drains from the high Andes through the Central Depression to the Pacific coast. The analysis links cosmogenic nuclide-derived erosion rates to sediment provenance using detrital zircon mixing models and mineral fertility.

## Drainage Network

The 11 sample sites form a nested hierarchy where downstream sites integrate sediment from upstream tributaries:

```
CT-7 ──→ CT-5 ──┐
                 ├──→ CT-10 ──→ CT-11 ──→ CT-8 ──→ CT-2 ──┐
CT-1 ──→ CT-6 ──┘                                          ├──→ CT-9
                                              CT-3 ─────────┘
```

- **Headwaters**: CT-1 (Rio Claro), CT-3 (small tributary), CT-7 (upper Teno)
- **Outlet**: CT-9 (Mataquito at coast)
- **CT-4**: Excluded from the hierarchy (problematic data); its area is lumped into the CT-10 subwatershed calculation

## Directory Structure

```
Mataquito/
├── data/
│   ├── sample_metadata/      # Sample locations, areas, erosion rates
│   ├── dem/                   # Flow direction and drainage area rasters
│   ├── shapefiles/            # Sample points and flowlines
│   ├── geochronology/         # Zircon U-Pb data sheets
│   └── geology/               # Subwatershed geologic unit data
├── notebooks/
│   ├── 01_dem_watershed_delineation.ipynb
│   ├── 02_production_rates.ipynb
│   ├── 03_erosion_rate_uncertainties.ipynb
│   ├── 04_zircon_age_distributions.ipynb
│   ├── 05_zircon_mixing.ipynb
│   ├── 06_fertility_analysis.ipynb
│   ├── 07_downstream_confluence_mixing.ipynb
│   ├── 08_river_profile_analysis.ipynb
│   └── 09_subwatershed_geology.ipynb
├── mataquito/                 # Shared Python package
├── tests/                     # Unit tests for mataquito package
├── figures/
│   ├── manuscript/            # Publication figures
│   ├── fertility/             # Quartz and zircon fertility plots
│   ├── mixing/                # Downstream confluence mixing plots
│   ├── zircon/                # Detrital zircon age distributions
│   ├── erosion/               # Erosion rate figures
│   └── geology/               # Geologic map figures
└── results/
    ├── detritalpy_output/     # detritalPy CSV/KML outputs
    ├── mixing_coefficients/   # Bootstrap mixing results
    └── geology_tables/        # Geologic summary tables
```

## Notebooks

| # | Notebook | Purpose | Key Inputs | Key Outputs | External Data? |
|---|----------|---------|------------|-------------|----------------|
| 01 | DEM watershed delineation | Delineate drainage basins from 30m DEM using D8 flow routing | `dem_utm30m.tif`, sample coordinates | Basin masks, area calculations, polygon shapefiles | Yes |
| 02 | Production rates | Calculate cosmogenic 10Be production rates via Stone (2000) scaling | `MataquitoSampleData.xlsx` | Production rates (atoms/g/yr) per sample | No |
| 03 | Erosion rate uncertainties | Compute subwatershed erosion rates by mass balance, propagate uncertainty | `MataquitoSampleData.xlsx`, flow network | Subwatershed erosion rates with MC uncertainty | No |
| 04 | Zircon age distributions | Plot detrital zircon U-Pb age distributions using detritalPy | Zircon U-Pb spreadsheets | KDE/PDP plots, age proportion CSVs | No |
| 05 | Zircon mixing | Model downstream zircon provenance via detritalPy-mix | Zircon data + parent/child pairs | Mixing coefficient CSVs, distribution plots | No |
| 06 | Fertility analysis | Quantify quartz and zircon fertility from mixing + cosmogenic data | Mixing coefficients, erosion rates | Fertility violin plots, erosion rate comparisons | No |
| 07 | Downstream confluence mixing | Trace D8 flowlines and analyze 10Be/erosion rate evolution downstream | Flow direction/area rasters, sample points | Distance vs. erosion rate and 10Be plots | Yes |
| 08 | River profile analysis | Chi-elevation profiles for coastal streams near Mataquito outlet | `dem_utm30m.tif`, `fd_utm30m`, `area_utm30m` | Chi vs. elevation plots for 6 coastal streams | Yes |
| 09 | Subwatershed geology | Calculate geologic unit proportions and statistical tests per subwatershed | Geologic map shapefiles, geology data | Unit percentages, chi-squared tests, summary tables | Yes |

## `mataquito` Package

Shared Python package extracted from the notebooks. Each module can be imported independently.

| Module | Description | Key Functions |
|--------|-------------|---------------|
| `sample_data` | Load sample metadata from Excel | `load_samples()`, `get_erosion_rates()`, `get_areas()` |
| `flow_network` | Drainage network topology and mass balance | `all_subwatershed_erosion_rates()`, `subwatershed_area()` |
| `erosion` | Monte Carlo uncertainty propagation for erosion rates | `generate_mc_samples()`, `flux_order_samples()` |
| `production_rates` | Stone (2000) latitude/altitude scaling for 10Be | `stone2000_production_rate()`, `elevation_to_pressure()` |
| `fertility` | Quartz and zircon fertility calculations | `load_wct()`, `quartz_fertility()`, `zircon_fertility()` |
| `flowlines` | D8 flowline tracing and stream distance calculations | `create_flowlines_from_raster()`, `calculate_stream_distances_from_confluence()` |

See `tests/` for usage examples.

## Data

**`data/sample_metadata/MataquitoSampleData.xlsx`** — Canonical sample data (Sheet1, 11 samples):

| Column | Description |
|--------|-------------|
| `Sample_ID` | Sample identifier (CT-1 through CT-11) |
| `Source_Area` | Upstream drainage area (km²) |
| `Erosion_rate` | Basin-averaged erosion rate (m/Myr) |
| `Ext_uncert` | External uncertainty on erosion rate (m/Myr) |
| `P_rate` | Cosmogenic 10Be production rate (atoms/g/yr) |

**`results/mixing_coefficients/`** — Bootstrap mixing proportion CSVs from detritalPy-mix (NB05). Each row is a bootstrap replicate; columns are parent contributions.

## External Data

Some notebooks require external data files not included in this repository (large DEMs, geologic map shapefiles, hillshade rasters). These notebooks have clearly marked `EXTERNAL DATA PATHS` sections at the top that must be configured for your local environment.

Affected notebooks: 01, 07, 08, 09.

## Running Tests

```bash
pip install pytest
cd Mataquito
pytest tests/ -v
```

Note: `test_flowlines.py` requires shapely and rasterio, which may have binary compatibility issues in some environments. To skip it: `pytest tests/ -v --ignore=tests/test_flowlines.py`

## Dependencies

**Required**: numpy, pandas, matplotlib, openpyxl, scipy

**Optional** (needed by specific notebooks):
- shapely, rasterio, geopandas — NB01, NB07, NB08 (DEM/spatial analysis)
- detritalPy — NB04, NB05 (zircon age distributions and mixing)
- TopoAnalysis — NB01, NB08 (DEM processing and chi analysis)
- seaborn — NB09 (geology plotting)
