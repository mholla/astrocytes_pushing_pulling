# Astrocytes based white matter growth during cortical folding

Code to reproduce the results of pre-print *"Astrocytes in white matter respond to tensile cues during cortical folding: a numerical study"* by Karan Taneja, Kengo Saito, Hiroshi Kawasaki and Maria Holland.

### Folder defintions
- `simulation-files`:  
  - `Abaqus-scripts`: Contains the python scripts and VUMAT files to run the simulations
  - `odb-data-files`: Contains the coordinates of the top of the cortex (pial surface) from the simulations in a .csv file, and the python script to extract them from the odb file.
  - `sim-analysis-GI`: Contains the python script that postprocesses the .csv files of the simulations from sub-folder `odb-data-files`, to generate the gyrification index (GI) file at each timepoint in that simulation.
- `image-processing`:
  - `Raw_image`: Contains raw histological images at P10 (n=2) and P16 (n=3)
  - `Raw_extracted_points`: Contains extracted cortical plate points from the images using ImageJ.
  - `Image_analysis`: Contains python script to calculate the GI values from extracted points and create figures for convex hulls on the raw images.
- `manuscript-figures`: Contains python scripts used to generate subplots in Figures 1 and 2 in the manuscript.
 - `GI_compare_model_images.py`: Python script that takes the .csv files from sub-folder `sim-analysis-GI` and values from the python code in `Image_analysis` to generate Figure 8 in the manuscript, which is stored in `manuscript-figures`. 
 - `python_requirements.txt`: A list of python libraries used to postprocess the results.




### Prerequisites

- [Abaqus 2024](https://www.3ds.com/products/simulia/abaqus)
- [Python](https://www.anaconda.com/download)
- [ImageJ](https://imagej.net/ij/)



## Running the files

### Abaqus python scripts

The python scripts in subfolders `Abaqus-scripts` and `odb-data-files` can be run on the terminal as 

    abaqus cae noGUI=<script_name>.py

### GI Postprocessing scripts

The python script to postprocess the .csv files, can be run on your preferred python IDE or on the terminal as

    ./<script_name>.py 


### Note: Install the python dependencies if needed

The python libraries mentioned in `python_requirements.txt` can be installed on your machine by going to the terminal and running the following command
    python -m pip install -r python_requirements.txt


## Authors

  - [Karan Taneja](https://scholar.google.com/citations?hl=en&user=j2vT-84AAAAJ)
  - [Kengo Saito](https://scholar.google.com/citations?hl=en&user=PL0U1YQAAAAJ)
  - [Hiroshi Kawasaki](https://scholar.google.com/citations?hl=en&user=mJ4WHW0AAAAJ)
  - [Maria Holland](https://scholar.google.com/citations?hl=en&user=dUTauN0AAAAJ&view_op=list_works&sortby=pubdate)

  
## Acknowledgments
The following folks contributed by helping in formulating the equations, discussing the validity of the results and working with ImageJ (fun).  
  - [Johannes Krotz](https://johanneskrotz.org/) 
  - [CoMMaND Lab at University of Notre Dame](https://commandlab.nd.edu/people/)
  - [Julian Najera](https://timelab.nd.edu/people/)

