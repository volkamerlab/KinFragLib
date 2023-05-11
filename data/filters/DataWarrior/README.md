# Use DataWarrior to create an Enamine Building Block file
The created DataWarrior file is used in 
`notebooks/custom_kinfraglib/1_3_custom_filters_synthesizability.ipynb`.

## Table of contents
1. Prepare KinFragLib fragment library for DataWarrior
2. Install DataWarrior
3. Get the Enamine Building Blocks


## 1. Prepare KinFragLib fragment library for DataWarrior

### a. Open a new notebook in the `notebooks/custom-kinfraglib` folder

### b. Load the fragment library

```python
%load_ext autoreload
%autoreload 2
```

```python
import pandas as pd
from pathlib import Path
from kinfraglib import utils, filters
```

```python
# define Path to data
HERE = Path(_dh[-1])
PATH_DATA = HERE / '../../data'
PATH_DataWarrior = HERE / '../../data/filters/DataWarrior'
```

```python
fragment_library = utils.read_fragment_library(PATH_DATA / 'fragment_library')
```

```python
pd.concat(fragment_library).reset_index(drop=True).shape
```

Apply pre-filters
- removing fragments in pool X
- removing duplicates
- removing fragments without dummy atoms (unfragmented ligands)
- removing fragments only connecting to pool X

For more information to the pre-filtering steps please have a look at 
`notebooks/kinfraglib/3_1_fragment_library_reduced.ipynb`.
```python
fragment_library = filters.prefilters.pre_filters(
    fragment_library)
```

### c. Save the fragment library to be processed with DataWarrior
```python
# save smiles without dummy atoms in sdf files for DataWarrior
filters.utils.save_fragments_wo_dummy(fragment_library, PATH_DataWarrior)
```

## 2. Install DataWarrior
Go to the [DataWarrior](https://openmolecules.org/datawarrior/download.html) web page and install 
DataWarrior for your OS

## 3. Get the Enamine Building Blocks
- open DataWarrior
- open file `fragments_wo_dummy.sdf` in the `data/filters/DataWarrior` directory 
  * choose "yes" when DataWarrior asks if it should interpret stereo centers as absolute
- mark all structures 
  * go to Database 
  * Search Enamine Building Blocks
  * choose superstructures at the top dropdown menu
  * set maximum row count to 100000 
  * set minimun package size to 1 
  * choose "retrieve cheapest..." at the bottom 
  * run query
- go to file dropdown menu 
  * save special 
  * choose SDfile 
  * name it `Enamine_Building_Blocks.sdf` and save it in a folder called DataWarrior in same 
    directory as your fragment library

Note: if you only want to use exact matches choose "exact structures" in the dropdown menu
