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

The automatic installation of the dependencies is only supported by PyMOL version >= 2.5


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

