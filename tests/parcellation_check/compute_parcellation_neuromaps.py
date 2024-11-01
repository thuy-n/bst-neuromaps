#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compute_dk_parcellation_neuromaps.py

Parcellate surface and volume annotations using the Desikan-Killiany atlas
and save results in two .csv files:

    - 'parcellation_dk_surface_neuromaps.csv'
    - 'parcellation_dk_volume_neuromaps.csv'

Surface annotations are in the fsaverage 164k (per hemisphere)
Volume annotations are in the MNI152 space

Both annotations have to be previously fetched with 'preprocess.py', and
are located at ../../tmp/surface and ../../tmp/volume respectively.

Created on 2024
"""

import sys
import os
import re
import numpy as np
import pandas as pd
import nibabel as nib
from neuromaps.parcellate import Parcellater
from neuromaps.images import annot_to_gifti

#%% Find surface and volume files
# Change working dir to script dir
script_path = os.path.realpath(sys.argv[0])
if sys.argv[0]:
    # Call from terminal
    script_path = os.path.dirname(script_path)
os.chdir(script_path)
# Directory for full resolution maps
tmp_dir = os.path.join(script_path, r'../../tmp')
srf_dir = os.path.join(tmp_dir, r'surface')
vol_dir = os.path.join(tmp_dir, r'volume')

# Surface, find all left-hemisphere GIfTI files
gii_files = []
surf_shortfiles = []
for root, _, files in os.walk(srf_dir):
    for file in files:
        if file.endswith('lh.shape.gii'):
            gii_files.append(os.path.join(root, file))
            # Create name Ntx__subtype_desc_source
            match = re.search(r'^.*surface/(.*?)/subtype-(.*?)_source-(.*?)_desc-(.*?)_', os.path.join(root, file))
            surf_shortfiles.append(f'{match.group(1)}__{match.group(2)}_{match.group(4)}_{match.group(3)}')

# Volume, find all gunzipped NIfTI files
nii_files = []
vol_shortfiles = []
for root, _, files in os.walk(vol_dir):
    for file in files:
        if file.endswith('.nii.gz'):
            nii_files.append(os.path.join(root, file))
            # Create name Ntx__subtype_desc_source
            match = re.search(r'^.*volume/(.*?)/subtype-(.*?)_source-(.*?)_desc-(.*?)_', os.path.join(root, file))
            vol_shortfiles.append(f'{match.group(1)}__{match.group(2)}_{match.group(4)}_{match.group(3)}')

# Checks fetched brain annotations
if len(gii_files) == 0 or len(nii_files) == 0 or len(gii_files)!= len(nii_files):
    sys.exit("Error: Fetch brain annotations with 'preprocess.py' before running this scrip")


#%% Surface parcellation

# Get Desikan-Killiany parcellation as .gii files
# The lh.aparc and rh.aparc .annot files come from freesurfer7.3/freesurfer/subjects/fsaverage/label
dk_hemisphere_filebases = ['lh.aparc', 'rh.aparc']
dk_labels = []
dk_giis = []
ix_offset = 0
for dk_hemisphere_filebase in dk_hemisphere_filebases:
    # Read .annot file to get labels
    annot_file = './' + dk_hemisphere_filebase + '.annot'
    annot_data = nib.freesurfer.io.read_annot(annot_file, orig_ids=True)
    # Find index for 'corpuscallosum' and 'unknown' in labels and remove them
    ids_labels = [bytes(byte_list).decode('utf-8')for byte_list in annot_data[2]]
    ids_labels.pop(ids_labels.index('corpuscallosum'))
    ids_labels.pop(ids_labels.index('unknown'))
    # Add hemisphere tag
    ids_labels = [ids_label + '_' + dk_hemisphere_filebase[0].upper() for ids_label in ids_labels]
    dk_labels = dk_labels + ids_labels
    # Convert .annot files to .GIfTI
    # gifti_file = './' + dk_hemisphere_filebase + '.gii'
    gifti_data = annot_to_gifti(annot_file)[0]
    # Add offset, so labels are different to non-zero labels
    ix_nz = np.nonzero(gifti_data.darrays[0].data)[0]
    gifti_data.darrays[0].data[ix_nz] = gifti_data.darrays[0].data[ix_nz] + ix_offset
    ix_offset = ix_offset + (len(ids_labels))
    # nib.save(gifti_data, gifti_file)
    dk_giis.append(gifti_data)

# Parcellater
parc = Parcellater(tuple(dk_giis), 'fsLR')

# Get GIfTI files per map
parcel_srf_results = np.zeros([len(dk_labels), len(gii_files)])
for ix, l_file in enumerate(gii_files):
    r_file = l_file.replace('_lh.', '_rh.')
    map_imgs = tuple([l_file, r_file])
    # Apply parcellation, background label is ignored
    parcellated_srf = parc.fit_transform(map_imgs, 'fsLR')
    parcel_srf_results[:,ix] = parcellated_srf

# Create CSV with data
df = pd.DataFrame(parcel_srf_results, columns=surf_shortfiles)
df.index = dk_labels
df.index.name = 'Row'
df.to_csv('parcellation_dk_surface_neuromaps.csv')

#%% Volume parcellation

# 'MNI152NLin2009aSym' was pre-processed with FreeSurfer 7.3,
# then 'aparc+aseg.mgz' was saved as 'Desikan-Killiany_volatlas.nii'
# already in the MNI512 space and volume size [193, 229, 193]

# The cortical labels are the same as in 'colortable_desikan_killiany.txt',
# except that left hemisphere has 1000 added to the index and the right has 2000 added.
# as indicated in the 'FreeSurferColorLUT.txt' file

# The 'aparc+aseg.mgz' contains 110 labels:
#  1 : Unknown
# 41 : ASEG labels
# 68 : Cortical labels ordered as in the desikan_killiany.txt

parc = Parcellater('Desikan-Killiany_volatlas.nii', 'mni152', resampling_target='parcellation')
# Use 'parcellation' as resampling_target to upsample 3mm brain annotations to 1mm
# if parcellation is downsampled, there is the risk of removing labels
#
# Parcellater, ignores the Unknown label, thus it has [109] parcels:
# [ 0 : 41] 41 Subcortical labels
# [41 :   ] 68 Cortical labels

# Get NIfTI file per map
parcel_vol_results = np.zeros([len(dk_labels), len(nii_files)])
for ix, nii_file in enumerate(nii_files):
    # Apply parcellation, background label is ignored
    parcellated_vol = parc.fit_transform(nii_file, 'mni152')
    # Keep only cortical parcells
    parcel_vol_results[:,ix] = parcellated_vol[0, 41: ]

# Create CSV with data
df = pd.DataFrame(parcel_vol_results, columns=vol_shortfiles)
df.index = dk_labels
df.index.name = 'Row'
df.to_csv('parcellation_dk_volume_neuromaps.csv')