#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script is part of the bst-neuromaps plugin
https://neuroimage.usc.edu/brainstorm/Tutorials/Neuromaps

The plugin integrates curated brain annotations (aka "brain maps" or "mas") and
tools to further expand the accessibility and inclusivity of brain-mapping into
Brainstorm.

This script fetchs the maps listed in "categorized.json" from the neuromaps.
Each entry in "categorized.json" follows the format:

    "Category": {                     # to organize map files
      "Subtype": [                    # to organize map files
        {
          "source":      "savli2012", # to fetch and rename map file
          "description": "way100635", # to fetch and rename map file
          "space":       "MNI152",    # to fetch and rename map file
          "density":     "3mm",       # to fetch and rename map file
          "N":           35,          # to rename map file
          "Age":         26.3         # to rename map file
        }

Where, a Category can have multiple Subtypes, a Subtype can habe multiple Maps,
each Map has the keys "source" "description" "space" "density", "N" and "age"

Maps are fetched in their original space and density. Then they are transformed
to the surface-based coordinate system used throughout FreeSurfer i.e.,
'fsaverage' with 164k vertices per hemisphere

The space 'fsaverage' '164k' is selected as target as this is a common surface
between neuromaps and Brainstorm

- Volume maps: saved as they are, projected to 'fsaverage' 164k and renamed
- Surface maps: projected to 'fsaverage' 164k and renamed

Finally all maps are saved in the "../tmp" dir following this structure:

../temp
├── surface
│   ├── Category1
│   │   ├── subtype-SUBTYPE_source-SOURCE_desc-DESCRIPTION_N-N_Age-AGE_space-fsaverage_den-164k_lh.shape.gii
│   │   └── subtype-SUBTYPE_source-SOURCE_desc-DESCRIPTION_N-N_Age-AGE_space-fsaverage_den-164k_rh.shape.gii
│   └── Category2
│       ├── subtype-SUBTYPE_source-SOURCE_desc-DESCRIPTION_N-N_Age-AGE_space-fsaverage_den-164k_lh.shape.gii
│       └── subtype-SUBTYPE_source-SOURCE_desc-DESCRIPTION_N-N_Age-AGE_space-fsaverage_den-164k_rh.shape.gii
└── volume
    └── Category1
        └── subtype-SUBTYPE_source-SOURCE_desc-DESCRIPTION_N-N_Age-AGE_space-mni152_den-1mm.nii.gz

Authors: Le Thuy Duong Nguyen, Raymundo Cassani 2023
"""

import sys
import os
import nibabel as nib
import numpy as np
import json
import shutil
from neuromaps.datasets import fetch_annotation, fetch_atlas
from neuromaps import transforms, images
from decimal import Decimal, ROUND_HALF_UP

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

# Directories for Surface and Volume maps
for map_type in ['surface', 'volume']:
    os.mkdir(os.path.join(tmp_dir, map_type))

# Surface maps will be in the 'fsaverage'space with 164k vertices per hemisphere
atlases_fsaverage = fetch_atlas("fsaverage", "164k")

# Load info for brain maps (annotations) from JSON file
# Category > Subtype > MapFields
with open(r'../categorized.json') as selected_maps:
    data = json.load(selected_maps)
maps = data["maps"]

# Process each map
for map_category, map_subtypes in maps.items():
    # Make folder for category
    os.mkdir(os.path.join(tmp_dir, 'surface', f'{map_category}'))
    os.mkdir(os.path.join(tmp_dir, 'volume',  f'{map_category}'))
    # Process all map_subtypes
    for map_subtype, map_objs in map_subtypes.items():
        for map_obj in map_objs:
            # Repace / with --- (3 dashes), to have have valid filepath
            map_subtype = map_subtype.replace('/', '---')
            source = map_obj["source"]
            desc   = map_obj["description"]
            space  = map_obj["space"]
            den    = map_obj["density"]
            N      = map_obj["N"]
            age    = int(Decimal(map_obj["Age"]).to_integral_value(rounding=ROUND_HALF_UP))
            # Fetch original map, saved in ../tmp/annotations/
            original_map = fetch_annotation(source=source, desc=desc, space=space, den=den, data_dir=tmp_dir)
            # Clear to avoid saving last iteration surfaces
            gii_maps = None
            
            # Convert the original map to surface maps: space='fsaverage' with den='164k'
            if space == 'MNI152':
                # Copy to ../tmp/volume/
                output_filename = os.path.join(tmp_dir, f"./volume/{map_category}/subtype-{map_subtype}_source-{source}_desc-{desc}_N-{N}_Age-{age}_space-mni152_den-{den}.nii.gz")
                shutil.copyfile(original_map, output_filename)
                # Not transforming these specific maps because neuromaps warns that they are best used in the provided fsaverage 164k space
                # https://github.com/netneurolab/neuromaps/blob/abc085a/neuromaps/datasets/annotations.py#L238
                if source == 'norgaard2021' or source == 'beliveau2017':
                    original_map = fetch_annotation(source=source, desc=desc, space='fsaverage', den='164k', data_dir=tmp_dir)
                    # Copy to ../tmp/surface/
                    for gii, hemi in zip(original_map, ['l', 'r']):
                        output = os.path.join(tmp_dir, f"./surface/{map_category}/subtype-{map_subtype}_source-{source}_desc-{desc}_N-{N}_Age-{age}_space-fsaverage_den-164k_{hemi}h.shape.gii")
                        shutil.copyfile(gii, output)
                        os.remove(gii)
                    # Nothing else to do for these surface maps
                    continue
                else:
                    # Generate surface maps in 'fsaverage' 164k space from MNI152 volumes
                    gii_maps = transforms.mni152_to_fsaverage(original_map, '164k')
                    os.remove(original_map)

            elif space == 'fsLR':
                # Generate surface maps in 'fsaverage' 164k space from 'fsLR'
                gii_maps = transforms.fslr_to_fsaverage(original_map, '164k')

            elif space == 'civet':
                # Generate surface maps in 'fsaverage' 164k space from 'civet'
                gii_maps = transforms.civet_to_fsaverage(original_map, '164k')

            elif space == 'fsaverage' and den != '164k'
                # Generate surface maps in 'fsaverage' 164k space from 'fsaverage' with different density
                gii_maps = transforms.fsaverage_to_fsaverage(original_map, '164k')
            
            # For maps transformed to 'fsaverage' 164k, set values in medial wall to NaN
            surface_maps = []
            # Process each hemisphere (one gii file)
            for gii_hemi_mask, gii_hemi in zip(atlases_fsaverage['medial'], gii_maps):
                mwall_hemi_mask = nib.load(gii_hemi_mask).agg_data()
                hemi_data = gii_hemi.agg_data()
                hemi_data[mwall_hemi_mask == 0] = np.nan
                surface_maps.append(images.construct_shape_gii(hemi_data))
            # Save maps to ../tmp/surface/
            for gii, hemi in zip(surface_maps, ['l', 'r']):
                output_filename = os.path.join(tmp_dir, f"./surface/{map_category}/subtype-{map_subtype}_source-{source}_desc-{desc}_N-{N}_Age-{age}_space-fsaverage_den-164k_{hemi}h.shape.gii")
                nib.save(gii, output_filename)
# Delete ../tmp/annotations which was created by fetch_annotation()
shutil.rmtree(os.path.join(tmp_dir, 'annotations'))


