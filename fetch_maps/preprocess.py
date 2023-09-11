import neuromaps
import nibabel as nib
import numpy as np
import json
from neuromaps.datasets import fetch_annotation, fetch_atlas
from neuromaps import transforms, images


# Load maps from JSON file
with open('categorized.json') as selected_maps:
    data = json.load(selected_maps)

data_dir = '/thuy-n/neuromaps/tree/main/maps' 
maps = data["maps"]

atlases_fsaverage = fetch_atlas("fsaverage", "164k")


# Process each neurotransmitter
for neurotransmitter, subtypes in maps.items():
    for subtype, maps in subtypes.items():
        for map_info in maps:
            source = map_info["source"]
            desc = map_info["desc"]
            space = map_info["space"]
            den = map_info["den"]

            original_map = fetch_annotation(source=source, desc=desc, space=space, den=den, data_dir=data_dir)
            
            # Not transforming these specific maps because neuromaps warns that they are best used in the provided fsaverage space 
            if source == 'norgaard2021' or (source == 'beliveau2017' and desc == 'cimbi36'):
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
                output_filename = f"maps/transformed_fsav164k/source-{source}_desc-{desc}_space-fsaverage_den-164k_{hemi}h.shape.gii"
                nib.save(gii, output_filename)