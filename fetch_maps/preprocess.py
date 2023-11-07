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

# Directory for full resolution maps
tmp_dir  = os.path.join(script_path, r'../tmp')

# Ensure this directory is empty
if os.path.isdir(tmp_dir):
    shutil.rmtree(tmp_dir)
os.mkdir(tmp_dir)

# Surface and Volume maps
for map_type in ['surface', 'volume']:
    os.mkdir(os.path.join(tmp_dir, map_type))

# Surface maps will be in the FsAveage space (164k vertices per hemisphere)
atlases_fsaverage = fetch_atlas("fsaverage", "164k")

# Load info for brain maps (annotations) from JSON file
# BrainMap > BrainMap Subtype > Description for BrainMap Subtype
with open(r'../categorized.json') as selected_maps:
    data = json.load(selected_maps)
brain_maps = data["maps"]

# Process each brain map
for brain_map, brain_map_subtypes in brain_maps.items():
    # Make folder for brain map
    os.mkdir(os.path.join(tmp_dir, 'surface', f'{brain_map}'))
    os.mkdir(os.path.join(tmp_dir, 'volume',  f'{brain_map}'))
    for brain_map_subtype, brain_map_descs in brain_map_subtypes.items():
        for brain_map_desc in brain_map_descs:
            source = brain_map_desc["source"]
            desc   = brain_map_desc["description"]
            space  = brain_map_desc["space"]
            den    = brain_map_desc["density"]
            N      = brain_map_desc["N"]
            age    = round(brain_map_desc["Age"])
            # Fetch original map
            original_map = fetch_annotation(source=source, desc=desc, space=space, den=den, data_dir=tmp_dir)
            if space == 'MNI152':
                output_filename = os.path.join(tmp_dir, f"./volume/{brain_map}/source-{source}_desc-{desc}_N-{N}_Age-{age}_space-mni152_den-{den}.nii.gz")
                shutil.copyfile(original_map, output_filename)
            # Not transforming these specific maps because neuromaps warns that they are best used in the provided fsaverage space
            # https://github.com/netneurolab/neuromaps/blob/abc085a/neuromaps/datasets/annotations.py#L238
            if source == 'norgaard2021' or (source == 'beliveau2017' and desc == 'cimbi36'):
                # Make a copy of the surface files and set the name in the same format as the other ones
                if space == 'fsaverage':
                    for gii, hemi in zip(original_map, ['l', 'r']):
                        output = os.path.join(tmp_dir, f"./surface/{brain_map}/source-{source}_desc-{desc}_N-{N}_Age-{age}_space-fsaverage_den-164k_{hemi}h.shape.gii")
                        shutil.copyfile(gii, output)
                        os.remove(output)
                continue
            # Generate surface images in FsAverage 164k space
            gii_images = transforms.mni152_to_fsaverage(original_map, '164k')
            os.remove(original_map)
            # Set values in medial wall to NaN
            surface_images = []
            for mask_gii, hemi_gii in zip(atlases_fsaverage['medial'], gii_images):
                mwall_mask = nib.load(mask_gii).agg_data()
                receptor_hemi_data = hemi_gii.agg_data()
                receptor_hemi_data[mwall_mask == 0] = np.nan
                surface_images.append(images.construct_shape_gii(receptor_hemi_data))
            # Save maps
            for gii, hemi in zip(surface_images, ['l', 'r']):
                output_filename = os.path.join(tmp_dir, f"./surface/{brain_map}/source-{source}_desc-{desc}_N-{N}_Age-{age}_space-fsaverage_den-164k_{hemi}h.shape.gii")
                nib.save(gii, output_filename)
shutil.rmtree(os.path.join(tmp_dir, 'annotations'))


