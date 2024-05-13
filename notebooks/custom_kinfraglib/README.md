# Custom KinFragLib notebooks
## 1. Custom filtering steps
### `1_1_custom_filters_unwanted_substructures.ipynb`
This notebook filters out unwanted substructures which are defined by PAINS
substructures and/or by the list from Brenk et al.
### `1_2_custom_filters_drug_likeness.ipynb`
This notebook filters out fragments not fulfilling the Rule of Three and the Quantitative Estimate
of Druglikeness (QED), which both reflect the molecular properties of the fragments.
### `1_3_custom_filters_synthesizability.ipynb`
This notebook filters the fragments for synthesizability using a buyable building block
filter and the SYnthetic Bayesian Accessibility (SYBA).
### `1_4_custom_filters_pairwise_retrosynthesizability.ipynb`
This notebook builds fragment pairs using only those fragments that passed all custom filtering steps.
Next, it uses the ASKCOS API to check if a one-step retrosynthetic route for this pair can be found and children, building this fragment pair, are returned from ASKCOS.
Finally, it will compare the children retrieved from ASKCOS with the fragments building the pair. If the fragments are
substructures of the children, then they pass this filter.
## 2. Custom filtering pipeline
### `2_1_custom_filters_pipeline.ipynb`
This notebook applies all the filtering steps explained before and gives the possibility to easily (de-)activate
single filtering steps and modify the parameters.
### `2_2_custom_filters_analysis.ipynb`
This notebook analyzes the custom-filtered fragment library and compares it with the pre-filtered and the reduced fragment set from the previous study.
### `2_3_custom_filters_paper.ipynb` 
This notebook contains all the code used to generate figures for the CustomKinFragLib paper. 
## Requirements
To succesfully apply the pairwise retrosynthesizability filters in `1_4_custom_filters_pairwise_retrosynthesizability.ipynb` and `2_1_custom_filters_pipeline.ipynb`, the ASCKOS API might need to be installed following the [Documenation](https://askcos-docs.mit.edu/guide/1-Introduction/1.1-Introduction.html). Note that for the data used in the notebooks, neccessary results from ASKCOS are already precomputed. Thus, ASKCOS does not need to be installed in that case. However, if the data is modified ASCKOS needs to be installed and running (using `make start` within the `askcos2_core` directory which is obtained following the [installion](https://askcos-docs.mit.edu/guide/1-Introduction/1.1-Introduction.html)).