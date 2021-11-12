# custom KinFragLib
## 1. custom filters
Custom filters that can be (de-)activated by the user to retrieve a custom filtered fragment
library.
### `1_1_custom_filters_unwanted_substructures.ipynb`
This notebook filters out unwanted substructures which are defined by PAINS
substructures and/or by the list from Brenk et al.
### `1_2_custom_filters_drug_likeness.ipynb`
This notebook filters out fragments not fulfilling the Rule of Three and the Quantitative Estimate
of Druglikeness (QED), which both reflect the molecular properties of the fragments.
### `1_3_custom_filters_synthesizability.ipynb`
This notebook is filtering the fragments for syntheszibalility using a buyable building block
filter and the SYnthetic Bayesian Accessibility (SYBA).
