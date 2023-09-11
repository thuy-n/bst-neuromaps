import nibabel as nib
import numpy as np
import json
import shutil
from neuromaps.datasets import fetch_annotation, fetch_atlas
from neuromaps import transforms, images

# Directories
data_dir = r'../tmp'
maps_dir = r'../maps'

# FsAveage atlases (164k vertices per hemisphere)
atlases_fsaverage = fetch_atlas("fsaverage", "164k")

# Load info for neuromaps for neurotransmitters from JSON file
with open(r'../categorized.json') as selected_maps:
    data = json.load(selected_maps)
neuro_tx_maps = data["maps"]
# NeuroTX > NeuroTXSubtypes > TracerMaps

# Process each neurotransmitter
for neuro_tx, neuro_tx_subtypes in neuro_tx_maps.items():
    for neuro_tx_subtype, tracer_maps in neuro_tx_subtypes.items():
        for map_info in tracer_maps:
            source = map_info["source"]
            desc   = map_info["description"]
            space  = map_info["space"]
            den    = map_info["density"]

            # Fetch original map
            original_map = fetch_annotation(source=source, desc=desc, space=space, den=den, data_dir=data_dir)
            
            # Not transforming these specific maps because neuromaps warns that they are best used in the provided fsaverage space 
            if source == 'norgaard2021' or (source == 'beliveau2017' and desc == 'cimbi36'):
                
                # Make a copy of the surface files and set the name in the same format as the other ones
                if space == 'fsaverage':
                    for gii, hemi in zip(original_map, ['l', 'r']):
                        output = f"maps/transformed_fsav164k/source-{source}_desc-{desc}_space-fsaverage_den-164k_{hemi}h.shape.gii"
                        shutil.copyfile(gii, output)
                continue

            # Generate surface images in a format supported by Brainstorm
            gii_images = transforms.mni152_to_fsaverage(original_map, '164k')

            # Set values in medial wall to NaN
            surface_images = []
            for mask_gii, hemi_gii in zip(atlases_fsaverage['medial'], gii_images):
                mwall_mask = nib.load(mask_gii).agg_data()
                receptor_hemi_data = hemi_gii.agg_data()
                receptor_hemi_data[mwall_mask == 0] = np.nan
                surface_images.append(images.construct_shape_gii(receptor_hemi_data))

            for gii, hemi in zip(surface_images, ['l', 'r']):
                output_filename = f"../maps/transformed_fsav164k/source-{source}_desc-{desc}_space-fsaverage_den-164k_{hemi}h.shape.gii"
                nib.save(gii, output_filename)

