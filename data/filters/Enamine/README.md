# Enamine Building Blocks Comparison
We use an RDKit-based substructure search to identify Enamine Building Blocks matching our fragments. 
The created building block file is used in 
`notebooks/custom_kinfraglib/1_3_custom_filters_synthesizability.ipynb`. 

## Enamine Building Blocks
To perform a substructure match, we need to first download the building blocks from the [Enamine website](https://enamine.net/building-blocks/building-blocks-catalog). Please deposit the building blocks SDF file in the `../data/filters/Enamine` folder. 


## Substructure search 
To save computational runtime, we have implemented the substructure search of fragments in our library with the Enamine building blocks separately. It is generated by executing the `kinfraglib/filters/enamine_substructures.py` file which performs a substructure search of a given fragment library against the Enamine building blocks. If a substructure is found, the smallest matching Enamine building block is written to `data/filters/Enamine/Enamine_Building_Blocks.sdf`:  
```bash
python3 kinfraglib/filters/enamine_substructures.py -e data/filters/Enamine/enamine_building_blocks_original.sdf -f data/fragment_library -o data/filters/Enamine/Enamine_Building_Blocks.sdf
```
