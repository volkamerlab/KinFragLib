# Use DataWarrior to create Enamine Building Block file

## Table of contents
1. Open new notebook in notebooks/custom-kinfraglib folder
2. Load fragment library
3. save fragment library for DataWarrior
4. Install DataWarrior
5. Get the Enamine Building Blocks


## 1. Prepare KinFragLib fragment library for DataWarrior

### a. Open a new notebook in the `notebooks/custom-kinfraglib` folder

## b. Load the fragment library

```python
%load_ext autoreload
%autoreload 2
```

```python
from kinfraglib import utils
```

```python
# define Path to data
HERE = Path(_dh[-1])
PATH_DATA = HERE / '../../data'
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
```python
fragment_library_pre_filtered = filters.prefilters.pre_filters(
    fragment_library)
```

## c. Save the fragment library to be processed with DataWarrior
```python
# save smiles without dummy atoms in sdf files for DataWarrior
filters.utils.save_smiles_wo_dummy(fragment_library, PATH_DATA)
```


## 4. Install DataWarrior
Go to the [DataWarrior](https://openmolecules.org/datawarrior/download.html) web page and install 
DataWarrior for your OS

## 5. Get Enamine Building Blocks
- open DataWarrior
- open file `smiles_wo_dummy.sdf` in the `data/fragment_library` directory 
  * choose 'yes' when DataWarrior asks if it should interpret stereo centers as absolute
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
  * name it "Enamine_Building_Blocks.sdf" and save it in a folder called DataWarrior in same 
    directory as your fragment library

Note: if you only want to use exact matches choose "exact structures" in the dropdown menu
