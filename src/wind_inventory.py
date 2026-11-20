"""Read the USWTDB attributes used by fleet preparation and wind validation."""
import geopandas as gpd
import numpy as np

from project_paths import RAW_DATA


def read_wind_turbines():
    archive = RAW_DATA / "usgs" / "uswtdb_V9_0_20260626.zip"
    columns = ["eia_id", "p_year", "t_cap", "t_model", "t_rd", "t_hh", "ylat", "xlong"]
    turbines = gpd.read_file(archive, columns=columns, ignore_geometry=True)
    # The shapefile marks missing numbers with -9999 and missing model names with "missing".
    return turbines.replace([-9999, "missing"], np.nan)
