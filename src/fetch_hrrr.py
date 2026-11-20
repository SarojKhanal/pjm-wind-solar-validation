"""Subset the uncompressed public HRRR-MET HDF5 object using read-only S3 range requests."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import local
import time

import fsspec
import h5py
import numpy as np
import pandas as pd
import requests
from scipy.spatial import cKDTree

from project_paths import INPUTS

URL = "https://nrel-pds-wtk.s3.us-west-2.amazonaws.com/hrrr_met_toolkit/v1.0.0/hrrr_nat_f02_conus_2025.h5"
SESSION = local()


def xyz(latlon):
    lat, lon = np.radians(latlon).T
    return np.column_stack((np.cos(lat)*np.cos(lon), np.cos(lat)*np.sin(lon), np.sin(lat)))


def read_range(start, size):
    if not hasattr(SESSION, "http"):
        SESSION.http = requests.Session()
    for attempt in range(4):
        try:
            end = start+size-1
            r = SESSION.http.get(URL, headers={"Range": f"bytes={start}-{end}"}, timeout=(10,60))
            r.raise_for_status()
            if r.status_code != 206 or not r.headers.get("Content-Range","").startswith(f"bytes {start}-{end}/") or len(r.content)!=size:
                raise ValueError("S3 returned an unexpected byte range")
            return r.content
        except (requests.RequestException, ValueError):
            if attempt==3:
                raise
            time.sleep(2**attempt)


def select_chunk(job, shape, dtype):
    row, offset, size, columns, targets = job
    raw = np.frombuffer(read_range(offset, size), dtype=dtype).reshape(shape)
    return row, targets, raw[:, columns].copy()


def main():
    fleet = pd.read_csv(INPUTS/"fleet.csv")
    plants = fleet.loc[fleet.resource.eq("wind"), ["eia_id","lat","lon"]].drop_duplicates().sort_values("eia_id")
    assert plants.eia_id.is_unique
    with fsspec.open(URL, block_size=262144, cache_type="blockcache") as file, h5py.File(file) as h:
        points = h["coordinates"][:]
        distances, sites = cKDTree(xyz(points)).query(xyz(plants[["lat","lon"]].to_numpy()))
        time_utc = h["time_index"][:].astype(str)
        field = h["windspeed_100m"]
        chunks, dtype = field.chunks, field.dtype
        assert not field.id.get_create_plist().get_nfilters() and dtype == np.dtype("<u2")
        fill, scale = field.attrs["fill_value"], field.attrs["scale_factor"]
        data = np.full((len(time_utc),len(sites)),fill,dtype=dtype)
        jobs = []
        for column in np.unique(sites // chunks[1] * chunks[1]):
            idx = np.where((sites >= column) & (sites < column+chunks[1]))[0]
            for row in range(0, len(time_utc),chunks[0]):
                chunk = field.id.get_chunk_info_by_coord((row,int(column)))
                assert chunk.size == np.prod(chunks)*dtype.itemsize
                jobs.append((row,chunk.byte_offset,chunk.size,sites[idx]-column,idx))
    print(f"Reading {len(jobs)} chunks ({sum(j[2] for j in jobs)/1e9:.2f} GB) for {len(plants)} plants",flush=True)
    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = [executor.submit(select_chunk,j,chunks,dtype) for j in jobs]
        for i, future in enumerate(as_completed(futures),1):
            row, targets, values = future.result()
            count = min(len(data)-row, len(values))
            data[row:row+count,targets] = values[:count]
            if i%100==0:
                print(f"{i}/{len(jobs)} chunks",flush=True)
    assert not np.any(data==fill)
    with TemporaryDirectory(dir=INPUTS) as tmp:
        out = Path(tmp)/"hrrr_2025.npz"
        np.savez_compressed(out, plant_id=plants.eia_id.to_numpy(), time_utc=time_utc, speed_100m=data,
                            scale=scale,grid_id=sites, site_latlon=points[sites],
                            distance_km=6371*2*np.arcsin(distances/2))
        out.replace(INPUTS/out.name)


if __name__ == "__main__":
    main()
