.. _getting-started-label:

Getting started
================

Requirements
-------------

To use highRES, a GAMS license is required. See the `GAMS website <https://www.gams.com/latest/docs/UG_License.html>`_ for more information.

Data bundle
------------

To run the full workflow, three datapackages are needed:

1. Resources (3.2 MB compressed, 6.4 MB uncompressed)
2. Weatherdata (4.3 GB)
3. Geodata (260 MB compressed, 518 MB uncompressed)

The datapackages can be downloaded from Zenodo https://doi.org/10.5281/zenodo.14223617.

Windows
----------------
1. Clone the repository
   - To also get the submodule, go into the cloned folder and run these two commands
   .. code-block:: console
      git submodule init
      git submodule update

2. Install snakemake
   - Download miniforge windows exe https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe
   - Install Miniforge
   - Open Miniforge Prompt from the start menu
   - Install the environment from the provided yaml file
   .. code-block:: console
      mamba env create -f workflow/envs/highres_environment.yaml

3. Activate the snakemake environment
.. code-block:: console
   mamba activate highres

4. Navigate to the repository in your snakemake conda environment shell

5. Get the required input files
.. code-block:: console
   zenodo_get 10.5281/zenodo.14223617

6. Extract the required input files
.. code-block:: console
   unzip resources.zip
   unzip weatherdata.zip
   unzip geodata.zip

7. Create a folder for shared input and move the geodata and weatherdata to that folder
.. code-block:: console
   mkdir shared_input
   mv geodata shared_input
   mv weatherdata shared_input

8. Make sure GAMS is installed and licensed and that gamspath is set correctly in the config file

9. Run
   .. code-block:: console
      snakemake -c all --configfile config/config_ci.yaml
   - Specify your own config file by changing the file name
