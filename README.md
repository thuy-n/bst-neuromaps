![bst-neuromaps](https://github.com/thuy-n/testing/assets/130115390/0b14b8e0-8b94-4928-b83f-3202a3b6bbba)

# Neuromaps Brainstorm plugin

This repository serves to host the plugin to use [neuromaps](https://github.com/netneurolab/neuromaps) directly in [Brainstorm](https://neuroimage.usc.edu/brainstorm/), as part of an [Open Science](https://www.mcgill.ca/neuro/open-science) initiative. Neuromaps is a toolbox for assessing, transforming and analyzing structural and functional brain annotations (Markello, Hansen et al., 2022). 

For this first iteration, we focused on the neurotransmitter receptors and transporters. Thirty different brain annotations from the neuromaps toolbox were selected, covering nine different neurotransmitter systems: dopamine, norepinephrine, serotonin, acetylcholine, glutamate, GABA, histamine, cannabinoid, and opioid. These maps are sourced from open-access repositories, addressing the need for a comprehensive tool that integrates standardized analytic workflows for both surface and volumetric data. 


The plugin consists of three main parts:
1. Python code to query and generate specific curated maps using the neuromaps toolbox
2. Files for neuromaps in a format supported by Brainstorm
3. Brainstorm processes to perform statistical analyses with the brain maps in this plugin

See the [Brainstorm website](https://neuroimage.usc.edu/brainstorm/Tutorials/Neuromaps) for a full tutorial on how to run the plugin. 


## Motivation

The growing field of neuroimaging has allowed researchers in the past decades to capture thousands of ‘snapshots’ of brain activity, expanding our understanding of the functional organization of the brain. With the increasing number of these images, the neuroscientific field needs a systematic way to compare them and establish connections between brain features. While the neuromaps toolbox has begun to transform how neuroscientists approach neuroimaging, its requirement for computer programming experience poses a significant barrier for many researchers seeking to use it. To address this, we implemented curated brain annotations from neuromaps, along with correlation processes that account for spatial autocorrelation within Brainstorm, a collaborative open-source software dedicated to advanced analyses of brain recordings. Extending these pioneering tools to the MATLAB environment provides researchers worldwide with cutting-edge neuroimaging tools in an intuitive point-and-click user environment. As technological and data sharing advances have increasingly directed neuroscience research towards integrative questions rooted in data science, we believe that Open Science is most impactful when everyone is provided with equal access to the newest resources in the field. This implementation holds great promise to advance integrative research and bridge the gap between neurophysiological and neurochemical systems, allowing researchers to identify potential neurochemical targets for future clinical treatments. The reasons behind this implementation can be summarized in three main points:
1. **Increased Accessibility**: Provide researchers access to neuromaps without any prior computer programming experience.
2. **Effortless Integration**: Make it readily available to Brainstorm’s 41,000+ users in an intuitive point-and-click interface for streamlined exploration.
3. **Enhanced Capabillities**: Brainstorm's advanced functionalities enhances the overall user experience by providing greater flexibility and control over the visualization and analysis of brain annotations.


## Key Features

- Two different coordinate systems:
  - FsAverage for surface maps- the default system used by the FreeSurfer software (164k vertices per hemisphere).
  - MNI152 for volumetric maps- 152 normative MRI scans developed by the Montreal Neurological Institute.
- Flexible framework designed to accommodate future expansions and updates.
- Repository of precomputed annotations sourced from the published literature for both volumetric and surface systems, ensuring accessibility and ease of use. All the information for these maps, including the appropriate citations to use can be found in [this spreadsheet](https://docs.google.com/spreadsheets/d/1R0usElQw1HCYaIGMpgJk-u3HcL6N1nQ1/edit?usp=sharing&ouid=114237437498686296895&rtpof=true&sd=true).
   - To obtain the surfaces, the original maps offered in MNI152 space were transformed to FsAverage using the registration fusion framework proposed in neuromaps ([Buckner et al., 2011](https://journals.physiology.org/doi/full/10.1152/jn.00339.2011); [Wu et al., 2018](https://onlinelibrary.wiley.com/doi/10.1002/hbm.24213)).
- Correlation processes that account for spatial autocorrelation.


## Citations
- Please cite the original conception of the neuromaps toolbox when using this plugin:
> Markello, R.D., Hansen, J.Y., Liu, ZQ. et al. neuromaps: structural and functional interpretation of brain maps. Nat Methods 19, 1472–1479 (2022). [DOI: 10.1038/s41592-022-01625-w](https://doi.org/10.1038/s41592-022-01625-w).
    
- If you used any of the included maps, please also cite the original papers that publish the data. References for each map can be found in [the spreadsheet](https://docs.google.com/spreadsheets/d/1R0usElQw1HCYaIGMpgJk-u3HcL6N1nQ1/edit?usp=sharing&ouid=114237437498686296895&rtpof=true&sd=true).

- If you used the surface maps which were transformed using the registration fusion framework, please also cite:
> Wu, J. et al. Accurate nonlinear mapping between MNI volumetric and FreeSurfer surface coordinate systems. Hum. Brain Mapp. 39, 3793–3808 (2018). [DOI: 10.1002/hbm.24213](https://doi.org/10.1002/hbm.24213).


## Contributing
We welcome contributions from the community to help improve and expand the functionality of the neuromaps plugin. Feel free to submit pull requests, report issues, or provide suggestions! Your feedback is invaluable in ensuring a user-friendly experience for researchers worldwide!

## License
The neuromaps plugin for Brainstorm is distributed under the [GNU General Public License Version 3](https://opensource.org/licenses/GPL-3.0), dated 29 June 2007. See the [LICENSE](LICENSE) file for more details.
