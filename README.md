# PyPCN: Protein Contact Networks in PyMOL

## A user-friendly PyMOL plugin for computation, visualization and analysis of Protein Contact Networks

Protein Contact Networks (PCNs) are a way to represent the tridimensional structure of a protein, allowing at the same time to simplify the description of protein complexity and apply the typical network formalism in the description of the structure-function relationship in proteins. Inter-residue contacts are described as binary adjacency matrices, which are derived from the graph representation of the α-carbons and distances according to defined thresholds. Algorithms for functional characterization, i.e. clustering techniques, centrality measures and community extractions metrics, are computed on binary adjacency matrices to unveil allosteric, dynamic and interaction mechanisms in proteins. Such strategies are commonly applied in a combinatorial way, albeit rarely found in seamless and user-friendly implementations. In this context, PCN-Miner is a Python module for integrating different algorithms and metrics dedicated to the analyses of PCNs. We have now integrated PCN-Miner in PyPCN, **an open-source PyMOL plugin to provide a GUI for assisting PCNs analyses**. A dedicated GUI, together with the visual support provided by PyMOL, makes the analysis more intuitive and simple, in a way that broadens the applicability of the analysis of proteins as PCNs.

Some general features:

- Handling of either PDBs or pre-computed adjacency matrices as inputs.
- Mapping of the results onto 3D protein structures providing an intelligible visualization. 
- Support for more than 24 algorithms for PCN analyses. 
- Visualization of contact matrices as interactive plots. 

... and so more, see the Quick Guide for further details

## Requirements

**Minimal requirement**: a recent version of PyMOL installed on your computer. 

PyPCN is compatible with incentive PyMOL builds distributed by [Schrodinger](https://pymol.org/2/ "Schrodinger website") (required PyMOL version >= 2.3.4) and open source builds (required PyMOL version >= 2.3.0).

PyPCN is distributed freely to the public and it has been tested and runs on Windows, macOS and Linux versions of PyMOL.

(Some incompatibilities may arise with the usage of PyMOL version 2.5.x if ‘undo’ function is enabled, which in PyMOL 2.5.2 still shows some shortcomings. Therefore, when the plugin is opened, the ‘undo’ function is automatically disabled and it is strongly suggested to keep it disabled when using the plugin.)


## Download

PyPCN plugin ZIP file: [download from here]( "PyPCN plugin ZIP file direct download") 


## Plugin Installation 
 
PyPCN is installed via the PyMOL plugin manager:

* First download the latest version of the plugin ZIP file [here](  "PyPCN plugin ZIP file direct download") 

* Launch PyMOL and use the *Plugin* → *Plugin Manager* command from the main menu of PyMOL. The plugin manager window of PyMOL will open.

* Click on *Install New Plugin* and press the *Choose File…* button. Select the **PyPCN ZIP file** which you have downloaded before. 
You will be asked to give the path of the directory in which to install the plugin files. Just select the default option if you are unsure about what to do (the location of the plugin files does not make any difference when running the plugin).


## Dependencies 

Some dependencies, which are not distributed along with PyMOL, must be installed. PyPCN is equipped with an **automatic installation** process, but only available on Incentive PyMOL version 2.5.

Please follow the instructions reported below for those cases in which the automatic installation is not supported (i.e. Open-Source PyMOL, Incentive PyMOL versions lower than 2.5, PyMOL setups in which the automatic installation fails for not widely known reasons). 

**Before continuing**: the protocol reported next is one of the multiple ways to do that, any alternative is plausible, as long as it allows for the dependencies to be correctly loaded from PyMOL. To check it, just try to import the module of interest in the PyMOL command-line. 
PyMOL in Ubuntu Linux, macOS and windows operating systems can be installed in several ways. However, we would suggest installing it in a dedicated Conda environment to ensure full compatibility with the dependencies. Please, follow the steps as reported.

**Install The “Conda package manager”.**

The “Conda package manager” has two versions, Anaconda and Miniconda, which are both functional for our purposes. If the “Conda package manager” is not installed, it can be downloaded at the following links:

    Anaconda: https://www.anaconda.com/products/distribution 
    
    Miniconda: https://docs.conda.io/en/latest/miniconda.html 

**Create and activate a dedicated environment** 

Environments in “Conda” can help in avoiding conflicts between installed packages. 
A documentation about that can be found here:

    https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

**To create a dedicated environment, the conda create command should be used:**

    conda create --name environment_name

**Then, the environment must be activated:**

    conda activate environment_name

**Install PyMOL in the environment**

When the “Conda package manager” is installed, and the dedicated environment created and activated:

   * open source PyMOL can be installed as reported at the next link: https://anaconda.org/conda-forge/pymol-open-source 
   * Incentive PyMOL can be installed as reported on the official website at the next link: https://pymol.org/2/

**Install PyPCN** 

PyPCN’s installation in open source PyMOL is the same as explained in “Section 2”.

**Install the dependencies** 

For each dependency, type the preferred command reported in Table 1. 

Usually, in a newly created environment, Biopython and matplotlib are not installed either. This lack would block the opening of PyPCN. Both are Python modules needed by some PyMOL’s features and its plugins. 

    Biopython can be installed by making use of the conda install command as reported at the link: https://anaconda.org/conda-forge/biopython  

    matplotlib can be installed by making use of the conda install command as reported at the link: https://anaconda.org/conda-forge/matplotlib


|      Name     |       Installation with pip        |            Other                 |      
|:-------------:|:----------------------------------:|:--------------------------------:|
|     2.5.2     | Linux (Ubuntu 20.04.3 LTS), 64-bit |          Incentive               |
|     2.5.1     | Linux (Ubuntu 21.04), 64-bit       |   Incentive (Conda package)      |
|     2.5.2     | Linux (Ubuntu 18.04.2 LTS), 64-bit |          Incentive               |
|     2.5.0     | Linux (Ubuntu 21.04), 64-bit       |   Open source (Conda package)    |
|     2.5.0     | Linux (Ubuntu 20.04.3 LTS), 64-bit |   Open source (Conda package)    |
|     2.4.1     | Windows (v.10 Home), 64-bit        |          Incentive               |
|     2.5.2     | Windows (v.10 Pro), 64-bit         |          Incentive               |
|     2.5.0     | Windows (v.10 Pro), 64-bit         |  Open source (Conda package)     |
|     2.5.2     | MacOS (High Sierra v.10.13.6)      |          Incentive               |
|     2.4.0     | MacOS (High Sierra v.10.13.6)      |  Open source (Conda package)     |
|     2.5.2     | MacOS (Monterey v.12.2.1)          |          Incentive               |
|     2.5.2     | MacOS (Big Sur v.11.6.5)           |          Incentive               |
|     2.4.0     | MacOS (Big Sur v.11.6.5)           |  Open source (Conda package)     |


## Threads 



## Tested platforms 

| PyMOL version |          Operating system          |         PyMOL source             |      
|:-------------:|:----------------------------------:|:--------------------------------:|
|     2.5.2     | Linux (Ubuntu 20.04.3 LTS), 64-bit |          Incentive               |
|     2.5.1     | Linux (Ubuntu 21.04), 64-bit       |   Incentive (Conda package)      |
|     2.5.2     | Linux (Ubuntu 18.04.2 LTS), 64-bit |          Incentive               |
|     2.5.0     | Linux (Ubuntu 21.04), 64-bit       |   Open source (Conda package)    |
|     2.5.0     | Linux (Ubuntu 20.04.3 LTS), 64-bit |   Open source (Conda package)    |
|     2.4.1     | Windows (v.10 Home), 64-bit        |          Incentive               |
|     2.5.2     | Windows (v.10 Pro), 64-bit         |          Incentive               |
|     2.5.0     | Windows (v.10 Pro), 64-bit         |  Open source (Conda package)     |
|     2.5.2     | MacOS (High Sierra v.10.13.6)      |          Incentive               |
|     2.4.0     | MacOS (High Sierra v.10.13.6)      |  Open source (Conda package)     |
|     2.5.2     | MacOS (Monterey v.12.2.1)          |          Incentive               |
|     2.5.2     | MacOS (Big Sur v.11.6.5)           |          Incentive               |
|     2.4.0     | MacOS (Big Sur v.11.6.5)           |  Open source (Conda package)     |

