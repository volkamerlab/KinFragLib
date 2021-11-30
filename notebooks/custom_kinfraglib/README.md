# Custom KinFragLib notebooks
## 1. Custom filtering steps
### `1_1_custom_filters_unwanted_substructures.ipynb`
This notebook filters out unwanted substructures which are defined by PAINS
substructures and/or by the list from Brenk et al.
### `1_2_custom_filters_drug_likeness.ipynb`
This notebook filters out fragments not fulfilling the Rule of Three and the Quantitative Estimate
of Druglikeness (QED), which both reflect the molecular properties of the fragments.
### `1_3_custom_filters_synthesizability.ipynb`
This notebook is filtering the fragments for synthesizability using a buyable building block
filter and the SYnthetic Bayesian Accessibility (SYBA).
### `1_4_custom_filters_pairwise_retrosynthesizability.ipynb`
This notebook is building fragment pairs using only those fragments that passed all custom filtering step.
Next, it requests retrosynthesis information using the ASKCOS API to check if a one-step retrosynthetic route for this pair can be found.
Finally, it will compare the children retrieved from ASKCOS with the fragments building the pair. If the fragments are
substructures of the children, then they pass this filter.
## 2. Custom filtering pipeline
### `2_1_custom_filters_pipeline.ipynb`
This notebook compares all the filtering steps explained before and gives the possibility to easily (de-)activate
single filtering steps and modify the parameters.
### `2_2_custom_filters_analysis.ipynb`
This notebook is analyzing the filtered fragment library and comparing it with the pre-filtered and the reduced set.