#!/usr/bin/env python
# coding: utf-8

import sys
import local_pcangsd as lp
import os
import lostruct
from skbio.stats.ordination import pcoa
import matplotlib.pyplot as plt
import seaborn as sns
import dask
import pandas as pd
import numpy as np


input = sys.argv[1]
store = sys.argv[2]


if not os.path.exists(store):
	lp.beagle_to_zarr(input, store, chunksize=100000)


print('finished beagle to zarr!')


ds = lp.load_dataset(store, chunks=100000) # open the Dataset

# Display the first 10 rows
#print(ds.head(100))

window_size = int(sys.argv[9])
var_num = int(sys.argv[10])

# For large chromosomes
# window size of 50000
# min variants of 500

# For small chromosomes
# window size of 20000
# min variants of 100

ds = lp.window(ds, type='position', size=window_size, min_variant_number = var_num)
print('loaded data')

# record the window loci
windf = pd.DataFrame({'window_start' : ds.window_start.to_numpy(),
		      'window_stop' : ds.window_stop.to_numpy(),
		      })
#windf
print('starting pca')

#%%time

pca_zarr_store = lp.pca_window(
    ds,
    store=sys.argv[3], # where to store the result
    tmp_folder=sys.argv[4], # need a tmp folder, /tmp/tmp_local_pcangsd is default
    k=5, # number of PCs to retain
    min_maf=0.05, # applied on each window
)

ds_pca = lp.load_dataset(sys.argv[3])
#ds_pca

results = lp.to_lostruct(ds_pca)

print(f"Results on {results.shape[0]} windows")

pc_dists = lostruct.get_pc_dists(results, jax=False)

mds = pcoa(pc_dists)

window_center = lp.get_window_center(ds_pca)

plt.figure()
plt.scatter(x=window_center, y=mds.samples["PC1"])
_ = plt.title("MDS Coordinate 1 (y-axis) compared to Window (x-axis)")
_ = plt.xlabel("Position of window (bp)")
_ = plt.ylabel("MDS 1")


output_plot_name = sys.argv[5]
output_plot_dir = sys.argv[6]

# Save the plot to a file
plt.savefig(output_plot_dir + output_plot_name + "_mds_plot.png", dpi=300, bbox_inches="tight")

# Close the plot
plt.close()


mds_12 = mds.samples.loc[:, ["PC1", "PC2"]].copy()
xy = mds_12.to_numpy()

corners = lostruct.corners(xy, prop=0.05) # no k parameter implemented yet, default=3, prop=0.05, default=1?
print(corners[:10])

mds_12['corner'] = 'other'
mds_12['window'] = range(pc_dists.shape[0])
mds_12['window_center'] = window_center

for i in range(3):
	mds_12.iloc[corners[:,i], 2] = f'corner {i+1}'

_ = sns.scatterplot(
	data=mds_12[mds_12.corner=='other'], x="PC1", y="PC2", color='gray',
)
_ = sns.scatterplot(
	data=mds_12[mds_12.corner!='other'], x="PC1", y="PC2", hue="corner",
	hue_order=["corner 1", "corner 2", "corner 3"],
	palette='colorblind',
)

# Save the plot to a file
plt.savefig(output_plot_dir + output_plot_name + "_pc1_pc2_scatterplot.png", dpi=300, bbox_inches="tight")

# Close the plot
plt.close()

_ = sns.scatterplot(
	data=mds_12[mds_12.corner=='other'],
	x='window_center', y='PC1', color='gray',
)
_ = sns.scatterplot(
	data=mds_12[mds_12.corner!='other'],
	x="window_center", y="PC1", hue="corner",
	hue_order=["corner 1", "corner 2", "corner 3"],
	palette='colorblind',
)
_ = plt.xlabel("Position of window (bp)")


# Save the plot to a figure
plt.savefig(output_plot_dir + output_plot_name + "_pc1_position_scatterplot.png", dpi=300, bbox_inches="tight")

# Close the plot
plt.close()

_ = sns.scatterplot(
	data=mds_12[mds_12.corner=='other'],
	x='window_center', y='PC2', color='gray',
)
_ = sns.scatterplot(
	data=mds_12[mds_12.corner!='other'],
	x="window_center", y="PC2", hue="corner",
	hue_order=["corner 1", "corner 2", "corner 3"],
	palette='colorblind',
)
_ = plt.xlabel("Position of window center (bp)")

# Save the plot to a figure
plt.savefig(output_plot_dir + output_plot_name + "_pc2_position_scatterplot.png", dpi=300, bbox_inches="tight")

# Close the plot
plt.close()

column_index = int(sys.argv[7])

corner_pca = lp.pcangsd_merged_windows(ds, corners[:, column_index], k=5)

input_file = sys.argv[8]
meta = pd.read_csv(input_file, sep='\t')

#meta = pd.read_csv('/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/16-local_pcangsd/male_samples_pops.tsv', sep='\t')
#ids = pd.read_csv('data/hilo/other_files/HILO_MAIZE55_PARV50_ids.list', names=['ID'])
#meta = meta.set_index(meta.ID).reindex(ids.ID)
#meta = meta.set_index(meta.Sample)
#pop_group = sys.argv[8]

# Custom palette
custom_palette = sns.color_palette([
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
    "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
    "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854"
])

# palette="colorblind"

_ = sns.scatterplot(
	x=corner_pca[3][0], y=corner_pca[3][1],
	hue=meta['Pop'],
	palette=custom_palette,
)
_ = plt.xlabel("PC1")
_ = plt.ylabel("PC2")
_ = plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

# Save the plot to a figure
plt.savefig(output_plot_dir + output_plot_name + "_merged_intervals_scatterplot.png", dpi=300, bbox_inches="tight")

# Close the plot
plt.close()
