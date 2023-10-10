#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to query and generate specific maps using neuromaps toolbox
Maps to be downloaded are specified in the categorized.json file

Created on 2023
"""
import sys
import os
import nibabel as nib
import numpy as np
import json
import shutil
from neuromaps.datasets import fetch_annotation, fetch_atlas
from neuromaps import transforms, images

# Change working dir to script dir
script_path = os.path.realpath(sys.argv[0])
if sys.argv[0]:
    # Call from terminal
    script_path = os.path.dirname(script_path)
os.chdir(script_path)

# Directories
tmp_dir  = os.path.join(script_path, r'../tmp')
maps_dir = os.path.join(script_path, r'../maps')

# Ensure these directories are empty
for dir_path in [tmp_dir, maps_dir]:
    if os.path.isdir(dir_path):
        shutil.rmtree(dir_path)
    os.mkdir(dir_path)

# Surface and Volume maps
for map_type in ['surface', 'volume']:
    os.mkdir(os.path.join(maps_dir, map_type))

# Surface maps will be in the FsAveage space (164k vertices per hemisphere)
atlases_fsaverage = fetch_atlas("fsaverage", "164k")

# Load neurotransmitters info for neuromaps from JSON file
# NeuroTX > NeuroTXSubtypes > TracerMaps
with open(r'../categorized.json') as selected_maps:
    data = json.load(selected_maps)
neuro_tx_maps = data["maps"]

# Process each neurotransmitter
for neuro_tx, neuro_tx_subtypes in neuro_tx_maps.items():
    # Make forlder of neuro_tx
    os.mkdir(os.path.join(maps_dir, 'surface', f'{neuro_tx}'))
    for neuro_tx_subtype, tracer_maps in neuro_tx_subtypes.items():
        for map_info in tracer_maps:
            source = map_info["source"]
            desc   = map_info["description"]
            space  = map_info["space"]
            den    = map_info["density"]
            # Fetch original map
            original_map = fetch_annotation(source=source, desc=desc, space=space, den=den, data_dir=tmp_dir)
            # Not transforming these specific maps because neuromaps warns that they are best used in the provided fsaverage space
            # https://github.com/netneurolab/neuromaps/blob/abc085a/neuromaps/datasets/annotations.py#L238
            if source == 'norgaard2021' or (source == 'beliveau2017' and desc == 'cimbi36'):
                # Make a copy of the surface files and set the name in the same format as the other ones
                if space == 'fsaverage':
                    for gii, hemi in zip(original_map, ['l', 'r']):
                        output = os.path.join(maps_dir, f"./surface/{neuro_tx}/source-{source}_desc-{desc}_space-fsaverage_den-164k_{hemi}h.shape.gii")
                        shutil.copyfile(gii, output)
                continue
            # Generate surface images in FsAverage 164k space
            gii_images = transforms.mni152_to_fsaverage(original_map, '164k')
            # Set values in medial wall to NaN
            surface_images = []
            for mask_gii, hemi_gii in zip(atlases_fsaverage['medial'], gii_images):
                mwall_mask = nib.load(mask_gii).agg_data()
                receptor_hemi_data = hemi_gii.agg_data()
                receptor_hemi_data[mwall_mask == 0] = np.nan
                surface_images.append(images.construct_shape_gii(receptor_hemi_data))
            # Save maps
            for gii, hemi in zip(surface_images, ['l', 'r']):
                output_filename = os.path.join(maps_dir, f'./surface/{neuro_tx}/source-{source}_desc-{desc}_space-fsaverage_den-164k_{hemi}h.shape.gii')
                nib.save(gii, output_filename)
# Delete tmp folder
shutil.rmtree(tmp_dir)


