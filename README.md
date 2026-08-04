# TerraClimate Agroclim

R pipeline to compute global, gridded agroclimatic indices for viticulture
from [TerraClimate](https://www.climatologylab.org/terraclimate.html) data,
and to relate them to the world's wine-growing regions.

It implements the Tonietto & Carbonneau (2004) multicriteria geoviticultural
classification (Huglin Heliothermal Index, Cool Night Index, Riou Dryness
Index) together with a Koppen-Geiger climate classification and a set of
complementary seasonal/extreme indices, at ~4 km global resolution, for
observed climate (2001-2020, 1961-1990, 1981-2010) and for +2°C / +4°C
warming-level climatologies.

**Code only.** The input/output data (NetCDF downloads, per-tile CSV
extractions, index rasters, figures) are not included in this repository -
see [Data sources](#data-sources) below.

## Pipeline

| # | Script | Purpose |
|---|--------|---------|
| 01 | `scripts/01_download_terraclimate.R` | Download yearly TerraClimate NetCDF files (tmin, tmax, vpd, pet, ppt) |
| 02 | `scripts/02_build_climatologies.R` | Average yearly files into a multi-year monthly climatology |
| 03 | `scripts/03_extract_annual_series_and_indices.R` | Extract every land pixel's monthly series and compute agroclimatic indices, year by year |
| 04 | `scripts/04_extract_climatology_indices.R` | Same computation, applied directly to a climatology (one averaged "year") |
| 05 | `scripts/05_build_agroclim_rasters_annual.R` | Compile the per-tile annual-series CSVs into summary NetCDF rasters (one per index) |
| 06 | `scripts/06_build_agroclim_rasters_climatology.R` | Compile the per-tile climatology CSVs into NetCDF rasters |
| 07 | `scripts/07_extract_vgdb_point_data.R` | Extract index rasters (all periods/scenarios) and terrain (elevation, slope) at Vineyard Geodatabase points |
| 08 | `scripts/08_describe_vgdb_points.R` | Descriptive statistics of the indices at vineyard points |
| 09 | `scripts/09_explore_vgdb_agroclim_features.R` | Exploratory analysis of agroclimatic features / Koppen classes at vineyard points |
| 10 | `scripts/10_filter_vineyard_maps_general.R` | Map areas climatically suitable for viticulture, general criteria |
| 11 | `scripts/11_filter_vineyard_maps_specific.R` | Same, for a specific set of criteria (e.g. a given style/region) |
| 12 | `scripts/12_statistics_general.R` | Surface-area statistics of the potentially suitable areas (general criteria) |
| 13 | `scripts/13_statistics_specific.R` | Surface-area statistics of the potentially suitable areas (specific criteria) |

Shared, reusable code lives in `R/`:

- `R/agroclim_indices.R` - `calc_agroclim_indices()`: the agroclimatic index
  calculations, shared by scripts 03 and 04.
- `R/koppen.R` - `koppen.fast()`: Koppen-Geiger climate classification.
- `R/tiling.R` - `make_pixel_tiles()` / `write_tile_list()`: splitting a
  global raster into batches for extraction.

`archive/` keeps older superseded versions of some scripts (see their header
comments), for traceability only - they are not maintained.

## Setup

1. R (>= 4.3) with the packages used across the scripts: `tidyverse`, `sf`,
   `terra`, `tidyterra`, `geodata`, `data.table`, `Rfast`, `ncdf4`,
   `downloader`, `elevatr`, `pacman`.
2. Copy `R/config.R.example` to `R/config.R` and edit the two paths for your
   own machine (`R/config.R` is git-ignored).
3. Run the scripts in `scripts/` in numeric order (each stage reads the
   previous stage's output - see the table above).

## Data sources

- **TerraClimate**: Abatzoglou, J.T., S.Z. Dobrowski, S.A. Parks, K.C.
  Hegewisch (2018). TerraClimate, a high-resolution global dataset of
  monthly climate and climatic water balance from 1958-2015. *Scientific
  Data*, 5, 170191. https://doi.org/10.1038/sdata.2017.191
- **+2°C / +4°C warming-level years**: Wang, X., Jiang, D., & Lang, X.
  (2018). Climate Change of 4°C Global Warming above Pre-industrial Levels.
  *Advances in Atmospheric Sciences*, 35(7), 757-770.
  https://doi.org/10.1007/s00376-018-7160-4
- **Winegrape Vineyard Geodatabase (VGDB)**: not publicly redistributed by
  this repository; contact the author for access conditions.
- **World administrative boundaries / DEM / Koppen reference table**: third-
  party reference GIS layers, see `R/config.R.example` for the expected
  directory layout.

## Method references

- Tonietto, J., & Carbonneau, A. (2004). A multicriteria climatic
  classification system for grape-growing regions worldwide. *Agricultural
  and Forest Meteorology*, 124(1-2), 81-97.
  https://doi.org/10.1016/j.agrformet.2003.06.001
- Hall, A., & Jones, G. V. (2010). Spatial analysis of climate in
  winegrape-growing regions in Australia. *Australian Journal of Grape and
  Wine Research*, 16(3), 389-404.
- Peel, M. C., Finlayson, B. L., & McMahon, T. A. (2007). Updated world map
  of the Koppen-Geiger climate classification. *Hydrology and Earth System
  Sciences*, 11(5), 1633-1644. https://doi.org/10.5194/hess-11-1633-2007

## Known limitations / planned work

- The Gladstones Biologically Effective Degree Days (BEDD) index is not yet
  implemented (see the `TODO` placeholder in `R/agroclim_indices.R`).
- Two small bugs present in the original scripts (incorrect days-in-month
  used in the soil evaporation term of the Dryness Index) were fixed while
  factoring the code into `R/agroclim_indices.R`; see `tests/` for a
  regression check that quantifies and documents the fix.

## License

Code released under the [MIT License](LICENSE). This does not cover the
input datasets, which have their own terms (see Data sources above).
