![bst-neuromaps](https://github.com/thuy-n/testing/assets/130115390/0b14b8e0-8b94-4928-b83f-3202a3b6bbba)

# Brainstorm neuromaps plugin 

This repository serves to host the plugin to use [`neuromaps`](https://github.com/netneurolab/neuromaps) directly in [`Brainstorm`](https://neuroimage.usc.edu/brainstorm/), as part of an [Open Science](https://www.mcgill.ca/neuro/open-science) initiative. `Neuromaps` is a toolbox for assessing, transforming and analyzing structural and functional brain annotations (Markello, Hansen et al., 2022). 

For this first iteration, we focused on the neurotransmitter receptors and transporters. Twenty six different brain annotations from the `neuromaps` toolbox were selected, covering nine different neurotransmitter systems: dopamine, norepinephrine, serotonin, acetylcholine, glutamate, GABA, histamine, cannabinoid, and opioid. These maps are sourced from open-access repositories, addressing the need for a comprehensive tool that integrates standardized analytic workflows for both surface and volumetric data. 


The plugin consists of three main parts:
1. Python code to query and generate specific curated maps using the `neuromaps` toolbox.
2. Ready-to-use brain annotation files from `neuromaps` in a format supported by `Brainstorm`.
3. `Brainstorm` processes to perform statistical analyses with the brain maps in this plugin.

See the [Brainstorm website](https://neuroimage.usc.edu/brainstorm/Tutorials/Neuromaps) for a full tutorial on how to run the plugin. 


# Table of Contents
* [Motivation](#motivation)
* [Key Features](#key-features)
* [Citations](#citations)
* [Contributing](#contributing)
* [License](#license)

## <a name="motivation"></a>Motivation

The growing field of neuroimaging has allowed researchers in the past decades to capture thousands of ‘snapshots’ of brain activity, expanding our understanding of the functional organization of the brain. With the increasing number of neuroimaging data, the neuroscientific field needed a systematic way to compare these and establish connections between brain features. <BR><BR> While the[`neuromaps`](https://github.com/netneurolab/neuromaps) toolbox has begun to transform how neuroscientists approach neuroimaging, its requirement for computer programming experience poses a significant barrier for many researchers seeking to use it. To address this, we implemented curated brain annotations from `neuromaps`, along with correlation processes that account for spatial autocorrelation within `Brainstorm`, a collaborative open-source software dedicated to advanced analyses of brain recordings. Extending these pioneering tools to the MATLAB environment provides researchers worldwide with cutting-edge neuroimaging tools in an intuitive point-and-click user environment. 

As technological and data sharing advances have increasingly directed neuroscience research towards integrative questions rooted in data science, we believe that Open Science is most impactful when everyone is provided with equal access to the newest resources in the field. This implementation holds great promise to advance integrative research and bridge the gap between neurophysiological and neurochemical systems, allowing researchers to identify potential neurochemical targets for future clinical treatments. 

The reasons behind this implementation can be summarized in three main points:

1. **Increased Accessibility**: Provide researchers access to `neuromaps` without any prior computer programming experience.
2. **Effortless Integration**: Make it readily available to `Brainstorm`’s 41,000+ registered users in an intuitive point-and-click interface for streamlined exploration.
3. **Enhanced Capabillities**: `Brainstorm`'s advanced functionalities enhances the overall user experience by providing greater flexibility and control over the visualization and analysis of brain annotations.


## <a name="key-features"></a>Key Features

- Flexible framework designed to accommodate future expansions and updates.
- Repository of precomputed annotations sourced from the published literature in the FsAverage coordinate system (the default system used by the FreeSurfer software), ensuring accessibility and ease of use.
   - All the information for these maps, including the appropriate citations to use can be found in [this spreadsheet](https://mcgill-my.sharepoint.com/:x:/g/personal/le_thuy_nguyen_mail_mcgill_ca/EThdovlDN1tIiK85qmzBbZsBmON9MKwWSsBtgURTrDO5tg?e=xDEeSX).
   - To obtain the surfaces, the original maps offered in MNI152 space were transformed to FsAverage using the registration fusion framework proposed in neuromaps ([Buckner et al., 2011](https://journals.physiology.org/doi/full/10.1152/jn.00339.2011); [Wu et al., 2018](https://onlinelibrary.wiley.com/doi/10.1002/hbm.24213)).
- One-click installation within `Brainstorm`, enabling seamless visualization and facilitating integrative multimodal analysis.
- Correlation processes featuring spatially informed null models that account for spatial autocorrelation.

![features](https://github.com/user-attachments/assets/e00ddd86-802d-43fc-972a-0a2429fa21e2)

## <a name="citations"></a>Citations
- Please cite the original conception of the neuromaps toolbox when using this plugin:
  > Markello, R.D., Hansen, J.Y., Liu, ZQ. et al. neuromaps: structural and functionalinterpretation of brain maps. _Nat Methods_ **19**, 1472–1479 (2022). [DOI: 10.1038/s41592-022-01625-w](https://doi.org/10.1038/s41592-022-01625-w).

- If you used any of the included maps, please cite the original papers that publish the data. References for each map can be found in [the spreadsheet](https://mcgill-my.sharepoint.com/:x:/g/personal/le_thuy_nguyen_mail_mcgill_ca/EThdovlDN1tIiK85qmzBbZsBmON9MKwWSsBtgURTrDO5tg?e=8Uh9wE).

- If you used the surface maps which were transformed using the registration fusion framework, please also cite:
  > Wu, J. et al. Accurate nonlinear mapping between MNI volumetric and FreeSurfer surface coordinate systems. _Hum. Brain Mapp._ **39**, 3793–3808 (2018). [DOI: 10.1002/hbm.24213](https://doi.org/10.1002/hbm.24213).


## <a name="contributing"></a>Contributing
We welcome contributions from the community to help improve and expand the functionality of the neuromaps plugin. Feel free to submit pull requests, report issues, or provide suggestions! Your feedback is invaluable in ensuring a user-friendly experience for researchers worldwide! <BR><BR> For any questions, comments, or concerns, please contact Thuy at le.thuy.nguyen@mail.mcgill.ca or use the [feedback box](https://neuroimage.usc.edu/brainstorm/Tutorials/Neuromaps#Contributing) from the tutorial page.

## <a name="license"></a>License
The `Brainstorm-neuromaps plugin` is distributed under the [GNU General Public License Version 3](https://opensource.org/licenses/GPL-3.0), dated 29 June 2007. See the [LICENSE](LICENSE) file for more details.
