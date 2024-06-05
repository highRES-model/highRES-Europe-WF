# highRES-model

[![Documentation Status](https://readthedocs.org/projects/highres-europe-wf/badge/?version=latest)](https://highres-europe-wf.readthedocs.io/en/latest/?badge=latest)

Welcome to the repository for the European version of the high temporal and spatial resolution electricity system model (highRES). The model is used to plan least-cost electricity systems for Europe and specifically designed to analyse the effects of high shares of variable renewables and explore integration/flexibility options. It does this by comparing and trading off potential options to integrate renewables into the system including the extension of the transmission grid, interconnection with other countries, building flexible generation (e.g. gas power stations), renewable curtailment and energy storage.

highRES is written in GAMS and its objective is to minimise power system investment and operational costs to meet hourly demand, subject to a number of system constraints. The transmission grid is represented using a linear transport model. To realistically model variable renewable supply, the model uses spatially and temporally-detailed renewable generation time series that are based on weather data.

## Getting started

To run the full workflow, two datapackages are needed they can be downloaded from:

1. (~80MB compressed, ~300MB uncompressed) <https://uio-my.sharepoint.com/:u:/g/personal/tobiasvh_uio_no/Eftsg10mEK9Mpi4TSN8aS9kBWlooGJ_99YDDaYcGiQvrYQ?e=xeu9Lk&download=1>.
2. (~10GB) <https://uio-my.sharepoint.com/:u:/g/personal/tobiasvh_uio_no/EdEmFkUQoL5Imy3-OumK_o0BcFqilpjB3CQOCbUwi_1T8g?e=O0kq50&download=1>

## Windows
1. Clone the repository
2. Install snakemake
    - Download miniforge windows exe <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe>
    - Install Miniforge
    - Open Miniforge Prompt from the start menu
    - Run the minimal install of the snakemake environment `mamba create -c bioconda -c conda-forge -n snakemake snakemake-minimal pandas zstd`
3. Activate the snakemake environment `mamba activate snakemake`
4. Navigate to the repository in your snakemake conda environment shell
5. Get the required input files
    ```
   curl -o shared_input.tar -L -b cookies.txt "https://uio-my.sharepoint.com/:u:/g/personal/tobiasvh_uio_no/EdEmFkUQoL5Imy3-OumK_o0BcFqilpjB3CQOCbUwi_1T8g?e=O0kq50&download=1" -o resources.tar.zst -L -b cookies.txt "https://uio-my.sharepoint.com/:u:/g/personal/tobiasvh_uio_no/Eftsg10mEK9Mpi4TSN8aS9kBWlooGJ_99YDDaYcGiQvrYQ?e=xeu9Lk&download=1"
   ```
6. Extract the required input files
    ```
    zstd -d resources.tar.zst
    mkdir resources
    tar xf resources.tar -C resources
    mkdir shared_input
    tar xf shared_input.tar -C shared_input
    ```
7. Make sure GAMS is installed and licensed and that gamspath is set correctly in the config file
8. Run `snakemake -c all --use-conda`


## Documentation 
[Documentation](https://highres-europe-wf.readthedocs.io/en/latest/)
