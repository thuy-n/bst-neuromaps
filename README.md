![bst-neuromaps](https://github.com/thuy-n/testing/assets/130115390/0b14b8e0-8b94-4928-b83f-3202a3b6bbba)

# neuromaps as a Brainstorm plugin

This repository serves to host the plugin to use [neuromaps](https://github.com/netneurolab/neuromaps) directly in [Brainstorm](https://neuroimage.usc.edu/brainstorm/), as part of an [Open Science](https://www.mcgill.ca/neuro/open-science) initiative. Extending the toolbox from Python to the Brainstorm MATLAB environment provides researchers worldwide with the newest cutting-edge neuroimaging tools in an intuitive point-and-click user environment.

For this first iteration, we focused on the neurotransmitter receptors and transporters. Thirty different maps from the neuromaps toolbox were selected, covering nine different neurotransmitter systems: dopamine, norepinephrine, serotonin, acetylcholine, glutamate, GABA, histamine, cannabinoid, and opioid. These maps are sourced from open-access repositories, addressing the need for a comprehensive tool that integrates standardized analytic workflows for both surface and volumetric data. 

The plugin consists of three main parts:

1. Python code to query and generate specific curated maps using the neuromaps toolbox
   
3. Files for neuromaps in a format supported by Brainstorm

4. Brainstorm processes to perform statistical analyses with the brain maps in this plugin


## Key Features

- Two different coordinate systems:
  - MNI152 for volumetric maps- 152 normative MRI scans developed by the Montreal Neurological Institute.
  - FsAverage for surface maps- the default system used by the FreeSurfer software (164k vertices per hemisphere).
- Flexible framework designed to accommodate future expansions and updates.
- Repository of precomputed maps sourced from the published literature for both volumetric and surface systems, ensuring accessibility and ease of use. All the information for these maps, including the appropriate citations to use can be found in [this Excel sheet](https://docs.google.com/spreadsheets/d/1R0usElQw1HCYaIGMpgJk-u3HcL6N1nQ1/edit?usp=sharing&ouid=114237437498686296895&rtpof=true&sd=true). 


## Contributing
We welcome contributions from the community to help improve and expand the functionality of the neuromaps plugin. Feel free to submit pull requests, report issues, or provide suggestions! Your feedback is invaluable in ensuring a user-friendly experience for researchers worldwide. We believe that open science is most impactful when the countless everyone is provided with equal access to the newest and greatest resources in the field.

## Citations
Please cite the original conception of the neuromaps toolbox [Markello, Hansen et al., 2022](https://www.nature.com/articles/s41592-022-01625-w) when using this plugin. If you used any of the included maps, please also cite the original papers that publish the data. References for each map can be found in [the Excel sheet](https://docs.google.com/spreadsheets/d/1R0usElQw1HCYaIGMpgJk-u3HcL6N1nQ1/edit?usp=sharing&ouid=114237437498686296895&rtpof=true&sd=true).

## License
The neuromaps plugin for Brainstorm is distributed under the [GNU General Public License Version 3](https://opensource.org/licenses/GPL-3.0), dated 29 June 2007. See the [LICENSE](LICENSE) file for more details.
