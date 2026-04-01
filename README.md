# highRES-model

[![Documentation Status](https://readthedocs.org/projects/highres-europe-wf/badge/?version=latest)](https://highres-europe-wf.readthedocs.io/en/latest/?badge=latest)

Welcome to the repository for the WIMBY project version of the high temporal and spatial resolution electricity system model (highRES). This version was used to model the core scenarios in [The electricity system value of the local acceptance of onshore wind in Europe](https://arxiv.org/abs/2603.28279).

highRES-Europe is used to plan least-cost electricity systems for Europe and specifically designed to analyse the effects of high shares of variable renewables and explore integration/flexibility options. It does this by comparing and trading off potential options to integrate renewables into the system including the extension of interconnection between countries, building flexible generation (e.g. gas power stations), renewable curtailment and energy storage.

highRES is a combination of a optimisation module, written in GAMS, with inputs prepared using a snakemake workflow. Its objective is to minimise power system investment and operational costs to meet hourly demand, subject to a number of system constraints. The transmission grid is represented using a linear transport model. To realistically model variable renewable supply, the model uses spatially and temporally-detailed renewable generation time series that are based on weather data.

## Getting started

To run the full workflow, three datapackages are needed:

1. Resources (1.7 MB)
2. Weatherdata (4.8 GB)
3. Geodata (350 MB)

The datapackages can be downloaded from [Zenodo](https://zenodo.org/records/19087632).

## Windows
1. Clone the repository
    - To also get the submodule, go into the cloned folder and run these two commands
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
    zenodo_get 10.5281/zenodo.19087631
   ```
6. Extract the required input files
    ```
    unzip resources.zip
    unzip weatherdata.zip
    unzip geodata.zip
    ```
7. Create a folder for shared input and move the geodata and weatherdata to that folder
    ```
    mkdir shared_input
    mv geodata shared_input
    mv weatherdata shared_input
    ```
8. Make sure GAMS is installed and licensed and that gamspath is set correctly in the config file
9. To run the nine core scenarios from the paper `snakemake -c all --configfile config/config_james.yaml`. The scenario combinations across the social and environmental dimensions are created via running all combinations of the config file parameter:
	```
	environmental: [
	   "low",
	   "medium",
	   "high",
	]
	social: [
	   "low",
	   "medium",
	   "high",
	]
	```
	
**Please note** the full workflow run for each of the scenarios requires (mostly for the GAMS optimisation part):

- ~90GB of RAM
- 10 cores
- ~24 hours of wall clock time per scenario

## Demo version

A demo of the model is available by running:  

`snakemake -c all --configfile config/config_demo.yaml`

This runs one scenario (low-low) for one month (~730 hours) as opposed to the usual ~8760 hours which should complete in ~20 minutes. Expected output is a result.gdx file with a total system cost (the costs variable) value of 27288.3.

## System requirements

This version of highRES has been tested using:

- Operating systems: Microsoft Windows 11 Pro Version 10.0.26100 Build 26100 and Red Hat Enterprise Linux Server release 7.9 (Maipo)
- GAMS 43.4.1

## Documentation 
[General highRES Documentation](https://highres-europe-wf.readthedocs.io/en/latest/)
