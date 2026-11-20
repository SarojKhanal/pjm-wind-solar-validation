# 2025 PJM wind and solar validation

How closely does generation estimated from weather match the wind and solar power that PJM reports each hour? This EN.530.637 Energy Meteorology project compares ERA5 and HRRR-MET wind estimates and ERA5 solar PV estimates with PJM observations for 2025. The comparison covers the PJM system (RTO) and its western, mid-Atlantic and southern regions.

[Dataset downloads](#dataset-downloads) · [Weather-to-power code](src/validate.py) · [Run the analysis](#run-the-analysis) · [Hourly metrics](results/tables/metrics.csv) · [Example week](results/figures/hourly_example.png)

## Project layout

| Folder | Contents |
| --- | --- |
| [src/](src/) | Python scripts for fleet preparation, weather extraction, hourly validation, energy reconciliation and plotting. |
| [data/raw/](data/raw/) | Original ERA5 NetCDFs, PJM Data Miner exports, USGS inventories and shapefiles, EIA archives and the SAM turbine library. |
| [data/inputs/](data/inputs/) | The prepared fleet, PJM-specific EIA and HRRR-MET subsets, zone geometry and turbine curves used by the scripts. |
| [results/tables/](results/tables/) | Hourly series, monthly energy totals, metrics and audit tables. |
| [results/figures/](results/figures/) | Plant locations, example hourly generation, and seasonal and diurnal patterns. |

## Dataset downloads

The full hourly run reads three raw data sources. If a raw file is missing from your checkout, download the matching version and place it at the path below. The prepared inputs for the other sources are listed separately.

| Dataset | What to download and where to put it |
| --- | --- |
| PJM Data Miner 2 | Export hourly wind and solar generation as two CSVs. Select **January 1 through December 31, 2025 in Eastern prevailing time**, all six area labels and both UTC and EPT hour-beginning columns. Save as `data/raw/pjm/wind_gen.csv` and `data/raw/pjm/solar_gen.csv`. Each export has 52,560 data rows (8,760 hours × 6 areas). [Wind generation](https://dataminer2.pjm.com/feed/wind_gen) · [Solar generation](https://dataminer2.pjm.com/feed/solar_gen) |
| Copernicus ERA5 | Choose **ERA5 hourly data on single levels**, product type *reanalysis*, all UTC hours, January–December 2025, 0.25° grid, NetCDF. Set the area to north **43**, west **−91**, south **34.25**, east **−71.25**. Download one calendar month per variable group using the names below (36 files). [Climate Data Store download form](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=download) |
| USGS wind and solar inventories | Save the USWTDB **v9.0 (June 26, 2026)** shapefile ZIP as `data/raw/usgs/uswtdb_V9_0_20260626.zip`. Both wind scripts read turbine attributes directly from this archive. For solar, extract the USPVDB **v4.0 (April 14, 2026)** facility CSV into `data/raw/usgs/uspvdbCSV/`. [Wind turbine database downloads](https://energy.usgs.gov/uswtdb/data/) · [Solar PV v4.0 release](https://www.usgs.gov/data/united-states-large-scale-solar-photovoltaic-database-ver-40-april-2026) |

For ERA5, replace `YYYYMM` with the six-digit month (for example, `202501`). Keep the NetCDF coordinates `valid_time`, `latitude` and `longitude`.

| ERA5 fields | Expected file under `data/raw/era5/` |
| --- | --- |
| 100 m eastward and northward wind components (`u100`, `v100`) | `winds/2025/era5_winds_YYYYMM.nc` |
| 2 m air temperature (`t2m`) | `temps/2025/era5_temps_YYYYMM.nc` |
| Surface solar radiation downwards (`ssrd`) and total sky direct solar radiation at surface (`fdir`) | `solar/2025/era5_solar_YYYYMM.nc` |

The remaining inputs in [data/inputs/](data/inputs/) are the 2025 HRRR-MET wind series at the modeled plants, PJM subsets of EIA-860 and EIA-923, zone polygons, and SAM turbine curves. The [source and refresh notes](#source-files-and-optional-refresh) link to the original datasets and give the commands for rebuilding the EIA and HRRR-MET subsets.

## Run the analysis

From the project root, create a Python 3.12 environment. These commands use Conda to create the environment and pip to install the package versions used for the analysis.

```sh
conda create -n pjm-validation python=3.12 pip
conda activate pjm-validation
python -m pip install numpy==1.26.4 pandas==2.2.2 xarray==2023.6.0 netCDF4==1.7.4 pvlib==0.16.1 geopandas==1.1.1 matplotlib==3.9.2
```

With the raw files and prepared inputs in place, run the four scripts in order. This sequence runs offline.

```sh
python src/prepare_fleet.py
python src/validate.py
python src/reconcile_energy.py
python src/figures.py
```

1. [Prepare the fleet](src/prepare_fleet.py) joins the EIA inventory to USGS plant locations and assigns PJM regions.
2. [Validate hourly generation](src/validate.py) converts wind speed and solar irradiance to plant output, sums by region and matches modeled MW to PJM observations by UTC hour.
3. [Check monthly energy](src/reconcile_energy.py) compares modeled plant totals with EIA-923 and PJM. Use January–November for a comparison across complete months.
4. [Generate the figures](src/figures.py) plots plant locations, an example week, and seasonal and diurnal patterns. It also calculates sensitivity and low-output statistics.

## Reading the outputs

To see the system-wide result, open [metrics.csv](results/tables/metrics.csv) and filter to `region=RTO` and `split=all`. Under the baseline 15% wind-loss assumption, hourly wind RMSE is **923 MW for ERA5** and **730 MW for HRRR-MET**. ERA5 solar has a **+39.8% energy bias**. The final six Eastern hours of December 31 lack weather coverage, leaving 8,754 matched hours. Select `split=season` or `split=hour_ept` to see when errors are largest. The [seasonal and diurnal figure](results/figures/seasonal_diurnal.png) shows the corresponding generation profiles.

| File in `results/tables/` | Contents |
| --- | --- |
| `hourly_validation.csv.gz` | UTC hour beginning, observed MW, modeled MW, model and region. Each of three simulations has 8,754 paired rows for each of four areas. The two wind simulations share observations. |
| `metrics.csv` | Matched hours, means, bias, MAE, RMSE, energy bias as a percentage of observed energy and Pearson correlation. Filter `split=all` for full-period metrics. Other splits identify month, meteorological season, local hour and a daylight subset. |
| `observed_area_audit.csv` | Full 8,760-hour PJM coverage, MWh, extremes and negative-hour counts for every reported area. |
| `fleet_audit.csv` and `hourly_fleet_diagnostics.csv.gz` | Generator-level modeled hours and share with directly matched turbine curves, plus hourly active capacity and modeled output for fleet sensitivities. |
| `modeled_plant_monthly.csv.gz` and `monthly_energy_check.csv` | Simulated plant-month MWh and comparisons with PJM and EIA-923. December EIA totals cover the full month. December paired weather totals omit the final six local hours. |
| `sensitivity.csv`, `low_output.csv`, `paired_wind_rmse.csv` | Wind loss and PV fleet scenarios, lowest observed output deciles and paired 168-hour-block bootstrap for the difference in system wind RMSE. |

Hourly bias is simulated power minus observed power. Each MW value represents a one-hour interval, so its numerical contribution to energy is in MWh. An undefined correlation or a relative bias with a nonpositive reference total is left blank. Daylight uses midpoint solar elevation above 10° at 39.5° N, 79° W. For the low-output check, each technology uses its own observed tenth percentile, restricting PV to daylight first. The bootstrap samples paired squared errors with a fixed seed and circular blocks. It conditions on this year, its fleet and the baseline model.

## Source files and optional refresh

To rebuild a prepared EIA subset from the spreadsheets, run `python -m pip install openpyxl`. For HRRR-MET extraction, run `python -m pip install scipy "fsspec[http]" h5py requests`.

- **PJM and ERA5.** The 2025 wind and solar generation feeds are in `data/raw/pjm/`. The 36 monthly ERA5 files for 2025 are in `data/raw/era5/`, with 12 files each for wind, temperature and solar radiation. Compare each regional series directly. WEST, MIDATL and SOUTH account for the RTO total. Adding RFC or OTHER would double count generation. [PJM Data Miner](https://dataminer2.pjm.com/feed/wind_gen/definition) · [ERA5 single levels](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels)
- **Fleet.** The original USWTDB v9.0 and USPVDB v4.0 files in `data/raw/usgs/` supply locations and equipment attributes. The included `data/inputs/eia_pjm_2025.csv` is a selected-column subset of EIA-860 operating generators and 2025 retirements at plants with balancing authority PJM. Year-end out-of-service and standby units are excluded. The status snapshot can miss units that operated earlier in the year. `data/inputs/fleet.csv` preserves source coordinates, commissioning months, equipment fields and regional assignments. Refresh the subset with `python src/prepare_fleet.py --eia-zip data/raw/eia/eia8602025.zip`. [EIA-860 archive](https://www.eia.gov/electricity/data/eia860/xls/eia8602025.zip)
- **Geography.** The included `data/inputs/pjm_zone.geojson` is the 22-polygon PJM-21-bus-system layer used in the ICARUS-PJM-Dataset. It includes a separate MECK polygon grouped into SOUTH. The regional groupings follow PJM's broad regions. WEST includes AEP, APS, ATSI, COMED, DAY, DEOK, DUQ, EKPC and OVEC. MIDATL includes AECO, BGE, DPL, JCPL, METED, PECO, PENELEC, PEPCO, PPL, PSEG and RECO. The zone polygons and six transmission-owner fallbacks approximate a settlement-level crosswalk. Unresolved assets retain the label UNASSIGNED. [Polygon source](https://github.com/zhentong-shao/PJM-21-bus-system) · [PJM Manual 14B, regional definitions](https://www.pjm.com/-/media/DotCom/committees-groups/committees/mrc/2024/20240821/20240821-item-07b---4-manual-14b-revisions---clean.pdf)
- **HRRR-MET.** The included `data/inputs/hrrr_2025.npz` holds all 8,760 source timestamps, the 114 plant identifiers and selected grid point identifiers, grid coordinates, distance in km, and raw 100 m wind speed hundredths of m/s. Divide by the included scale factor to obtain m/s. The native f02 product uses the archive's valid timestamps. To refresh, run `python src/fetch_hrrr.py` after preparing the fleet. It reads coordinate and HDF5 chunk metadata and about 1.26 GB of selected wind chunks from the 3.28 TB public annual object. [Official OEDI archive](https://data.openei.org/submissions/8636)
- **Wind and PV conversion.** `data/inputs/wind_curves.csv` contains speed and normalized power for 19 selected SAM curves paired to USWTDB model strings. Curves are linearly interpolated. The fixed wind-loss fraction stands in for combined turbine and plant losses. The wind baseline uses plant centroids, a static turbine mix and a common vertical-shear exponent. The PV model gives priority to the EIA generator's tracking flags and geometry. USPVDB plant attributes fill gaps, followed by fixed tilt 25°, azimuth 180° and DC/AC ratio 1.25 as applicable. Single-axis geometry is standardized across sites. [SAM library](https://github.com/NatLabRockies/SAM/blob/develop/deploy/libraries/Wind%20Turbines.csv) · [pvlib](https://pvlib-python.readthedocs.io/en/stable/)
- **Energy check.** `data/inputs/eia923_pjm_2025.csv` contains plant/technology monthly net generation from the final 2025 EIA-923 generation and fuel sheet. The retained rows match the weather-covered plant IDs, fuel SUN or WND and prime movers PV, WT or WS. Refresh after the hourly validation with `python src/reconcile_energy.py --eia-zip data/raw/eia/f923_2025.zip`. [EIA-923 archive](https://www.eia.gov/electricity/data/eia923/xls/f923_2025.zip)

SHA-256 fingerprints for the source files used in this run are listed below, with paths relative to the project root.

```
2b27929d26cc1da9ad6a530da1cac6718cd8496098ee560c4d47cfffd96d08ba  data/raw/eia/eia8602025.zip
1bff7092a86f6678c069c986829da5a23104048fc1aafb65b0776e3c904284de  data/raw/eia/f923_2025.zip
dcfbf9149a82dce3eaa1de07002ede19c73e4eea31456109eceacaf27fc0b2c1  data/raw/sam/Wind Turbines.csv
6f5935061eeeb5a66eb3dafa22bff2358b285eea2cccd7bd94efc771b51a4e26  data/inputs/pjm_zone.geojson
```
