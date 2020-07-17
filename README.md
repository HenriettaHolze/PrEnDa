# Local Dashboard

## Installation

- clone this repo
```
git clone https://github.com/HenriettaHolze/<repo_name>.git
```
- make conda env from yaml file
```
conda env create -f <repo_name>/env.yaml
source activate PrEnDa
```

## Quick guide

To start the dashboard using one of the provided example data sets, run:
```
source activate PrEnDa
python <repo_name>/dashboard_local.py
```
When you click the link in the terminal, the dashboard should open in your browser. 

To use the other example dataset:
```sh
DATA_PATH="example_data/gssm" python dashboard_local.py
```


## Requirements

To explore your own data, you must provide a directory containing the following files:

- variant data csv `data.csv`
    - variant column: unique variant ID
    - parent column: variant ID of parent(s), separated by space if multiple (merge of 2 variants) or empty for the origin/wild-type
    - mutation column: 1-indexed mutations that are new in relation to parent variant, space separated
    - attribute column: one numeric attribute of variants

- settings file `settings.json`
    - json file that maps column headers of the data.csv file to the different column types
    - columns that are not specified, will be ignored

```
{
    "name": "4WOR_A_GSSM",
    "column_definitions": {
        "variant": "variant",
        "parent": "parent",
        "mutation": "mutation",
        "attribute": "ddG"
    }
}
```

- protein sequence fasta file `protein.fasta`
    - protein sequence of wild type or origin of directed engineering experiment

- PDB file `protein.pdb`

- MSA fasta file with homologous sequences `msa.fasta`
    - can be made with [HHblits](https://toolkit.tuebingen.mpg.de/tools/hhblits)
    - resulting a3m file can be converted to fasta file with `reformat.pl` [script](https://github.com/soedinglab/hh-suite/blob/master/scripts/reformat.pl) from the HHsuite  
    `./reformat.pl a3m fas <HHBLITS_A3M_FILE> <HHBLITS_FASTA_FILE> -r`

In order to run the dashboard with your own data, you have to set the environment variable `DATA_PATH`. You can simply run the dashboard with `DATA_PATH="<data_path>" python dashboard_local.py`. This will set the environment variable.  

## Data

2 example datasets are provided: one from a gene site saturation mutagenesis (GSSM) experiment and one simulating a directed protein engineering approach (directed).
The original data was retrived from the FireProt database. 

Give proper reference to FireProt!
Describe data preparation

## docker
test reproducibility with Docker
docker build -t prenda_image <path>
docker run -it -rm -p "8050:8050" -name prenda_container prenda_image

