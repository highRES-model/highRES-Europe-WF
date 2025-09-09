# highRES-model

[![Documentation Status](https://readthedocs.org/projects/highres-europe-wf/badge/?version=latest)](https://highres-europe-wf.readthedocs.io/en/latest/?badge=latest)

Welcome to the repository for the European version of the high temporal and spatial resolution electricity system model (highRES). The model is used to plan least-cost electricity systems for Europe and specifically designed to analyse the effects of high shares of variable renewables and explore integration/flexibility options. It does this by comparing and trading off potential options to integrate renewables into the system including the extension of the transmission grid, interconnection with other countries, building flexible generation (e.g. gas power stations), renewable curtailment and energy storage.

highRES is written in GAMS and its objective is to minimise power system investment and operational costs to meet hourly demand, subject to a number of system constraints. The transmission grid is represented using a linear transport model. To realistically model variable renewable supply, the model uses spatially and temporally-detailed renewable generation time series that are based on weather data.

## Getting started

To run the full workflow, two datapackages are needed they can be downloaded from [Zenodo](https://zenodo.org/records/15536457):

1. Resources (3.1 MB compressed, 6.3 MB uncompressed) <https://zenodo.org/records/15536457/files/resources.zip?download=1>.
2. Weatherdata (4.3 GB) <https://zenodo.org/records/15536457/files/weatherdata.zip?download=1>
3. Geodata (260 MB compressed, 517 MB uncompressed) <https://zenodo.org/records/15536457/files/geodata.zip?download=1>

## Windows
1. Clone the repository
    - To also get the submodule, go into the cloned foder and run these two commands
        ```
        git submodule init
        git submodule update
        ```
2. Install snakemake
    - Download miniforge windows exe <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe>
    - Install Miniforge
    - Open Miniforge Prompt from the start menu
    - Install the environment from the provided yaml file `mamba env create -f workflow/envs/highres_environment.yaml`
3. Activate the snakemake environment `mamba activate highres`
4. Navigate to the repository in your snakemake conda environment shell
5. Get the required input files
    ```
   curl -L -b cookies.txt "https://zenodo.org/records/15735203/files/resources.zip?download=1" --output resources.zip
   curl -L -b cookies.txt "https://zenodo.org/records/15735203/files/weatherdata.zip?download=1" --output weatherdata.zip
   curl -L -b cookies.txt "https://zenodo.org/records/15735203/files/geodata.zip?download=1" --output geodata.zip
   ```
6. Extract the required input files
    ```
    unzip resources.zip
    unzip weatherdata.zip
    unzip geodata.zip
    ```
7. Create a folder for shared input and move the geodata and weatherdata to that folder.
    ```
    mkdir shared_input
    mv geodata shared_input
    mv weatherdata shared_input
    ```
8. Make sure GAMS is installed and licensed and that gamspath is set correctly in the config file
9. Run `snakemake -c all --configfile config/config_ci.yaml`


## Documentation 
[Documentation](https://highres-europe-wf.readthedocs.io/en/latest/)
