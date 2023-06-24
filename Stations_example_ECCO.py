"""
Example demonstrating how to plot center, corner and velocity
points from ECCO dataset when smapling isolated, discrete coordinates
"""


import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import oceanspy as ospy
from matplotlib import gridspec

from oceanspy.llc_rearrange import mates, rotate_vars

import copy as _copy


# from dask.distributed import Client
# client = Client()
# client

# Directory
ECCO_url = "catalog_ECCO.yaml"
ECCOod = ospy.open_oceandataset.from_catalog("LLC", ECCO_url)
od = ECCOod

print(od)

# stations from Poseidon Viewer

p=[{"type":"Point","coordinates":[-39.05481151032802,49.148182344281935]},{"type":"Point","coordinates":[-41.36769261284706,48.556348747481366]},{"type":"Point","coordinates":[-42.65262655869097,46.65061953371659]},{"type":"Point","coordinates":[-42.65262655869097,41.86937828918002]},{"type":"Point","coordinates":[-41.62467940201584,34.8253018629576]},{"type":"Point","coordinates":[-35.45699646196506,38.73449125129926]},{"type":"Point","coordinates":[-30.959727651511365,33.65689272953311]},{"type":"Point","coordinates":[-26.33396544647326,38.634189174007446]},{"type":"Point","coordinates":[-22.09368342518836,36.288699983098]},{"type":"Point","coordinates":[-21.32272305768201,41.0024199325517]},{"type":"Point","coordinates":[-30.31726067858941,42.44090828314131]},{"type":"Point","coordinates":[-35.58548985654945,45.49180909068974]},{"type":"Point","coordinates":[-33.4011021486148,49.315993610524]},{"type":"Point","coordinates":[-26.84793902481084,45.94038509073192]},{"type":"Point","coordinates":[-25.30601828979814,48.556348747481366]},{"type":"Point","coordinates":[-21.451216452266404,44.67529627981841]}]

lons, lats = ospy.utils.viewer_to_range(p)

# manual inspection of corner / velocity /scalar grid points


ds = _copy.deepcopy(od._ds)
co_list = [var for var in ds.coords]
ds = mates(ds.reset_coords())
ds = ds.set_coords(co_list)


if not ds.xoak.index:
    ds.xoak.set_index(["XC", "YC"], "scipy_kdtree")
cdata = {"XC": ("station", lons), "YC": ("station", lats)}
ds_data = xr.Dataset(cdata)
nds = ds.xoak.sel(XC=ds_data["XC"], YC=ds_data["YC"])

stations=nds.station.values


xg = []
yg = []
iX, iY, iface = (nds[f"{i}"].data for i in ("X", "Y", "face"))

for i in range(len(stations)):
    xrange = slice(iX[i], iX[i]+2)
    yrange = slice(iY[i], iY[i]+2)
    args = {'face': iface[i], 
            'Xp1': xrange,
            'Yp1': yrange,
           }
    xg.append(ds['XG'].isel(**args))
    yg.append(ds['YG'].isel(**args))
xg = np.array(xg)
yg = np.array(yg)


fig = plt.figure(figsize=(20, 10))

colors=[['#000099', '#0080FF'], ['darkred', 'red']]
for i in range(len(iface)):
    XG = xg[i]
    YG = yg[i]
    XU0 = XG[0, 0]
    XU1 = XG[0, 1]
    YV0 = YG[0, 0]
    YV1 = YG[1, 0]
    if iface[i] in np.arange(6):
        ii=0
        color=colors[ii][-1]
    else:
        ii=1
        color=colors[ii][-1]
    plt.plot(XG, YG, color, ls='', marker='x', markersize=2);
    plt.plot(XU0, nds.YC.values[i], colors[ii][0], marker='>', markersize=8)
    plt.plot(XU1, nds.YC.values[i], colors[ii][1], marker='>', markersize=8)
    plt.plot(nds.XC.values[i], YV0, colors[ii][0], marker='^', markersize=8)
    plt.plot(nds.XC.values[i], YV1, colors[ii][1], marker='^', markersize=8)
    plt.plot(nds.XC.values[i], nds.YC.values[i], color, marker='o', markersize=5, alpha=0.5)
# plt.savefig('faces_10_2_center_corner_points_original.png')
plt.show()


# corrected version

# note the tranposed behavior at rotated facets
fig = plt.figure(figsize=(20, 10))

colors=[['#000099', '#0080FF'], ['darkred', 'red']]
for i in range(len(iface)):
    if iface[i] in np.arange(6):
        ii=0
        color=colors[ii][-1]
        XG = xg[i]
        YG = yg[i]
        XU0 = XG[0, 0]
        XU1 = XG[0, 1]
        YV0 = YG[0, 0]
        YV1 = YG[1, 0]
    else:
        ii=1
        color=colors[ii][-1]
        XG = xg[i].T
        YG = yg[i].T[::-1, :]
        XU0 = XG[0, 0]
        XU1 = XG[0, 1]
        YV0 = YG[0, 0]
        YV1 = YG[1, 0]
    plt.plot(XG, YG, color, ls='', marker='x', markersize=2);
    plt.plot(XU0, nds.YC.values[i], colors[ii][0], marker='>', markersize=8)
    plt.plot(XU1, nds.YC.values[i], colors[ii][1], marker='>', markersize=8)
    plt.plot(nds.XC.values[i], YV0, colors[ii][0], marker='^', markersize=8)
    plt.plot(nds.XC.values[i], YV1, colors[ii][1], marker='^', markersize=8)
    plt.plot(nds.XC.values[i], nds.YC.values[i], color, marker='o', markersize=5, alpha=0.5)
# plt.savefig('faces_10_2_center_corner_points_corrected.png')
plt.show()

print('done')











