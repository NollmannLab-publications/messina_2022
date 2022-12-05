# Messina, et al., 2022

This repository contains the code and data needed to generate the figures for the following manuscript:

**3D chromatin interactions involving Drosophila insulators are infrequent and arise before TADs and transcription**

*Olivier Messina [1], Flavien Raynal [2], Julian Gurgo [1], Jean-Bernard Fiche [1], Vera Pancaldi [2,3], Marcelo Nollmann [1]*

[1] Centre de Biologie Structurale, Univ Montpellier, CNRS UMR 5048, INSERM U1054, Montpellier, France.
[2] Université de Toulouse, Inserm, CNRS, Université Toulouse III-Paul Sabatier, Centre de Recherches en Cancérologie de Toulouse, Toulouse, France.
[3] Barcelona Supercomputing Center, Barcelona, Spain.

## File contents

The repository contains the following scripts:

- HiC_Aggregation_Analysis.m 

This is a Matlab script used to compute Aggregation Peak Analysis for chromatin binding factors. 
The functions used in this script are encoded in an additional file called functionsContainer.m

- HiM_matrices_plots.m

This is a Matlab script used to compute pair-wise distance and proximity maps.

### Data

The latest version of the data can be found at: https://osf.io/aqtxj/.

You will be able to download a Numpy array with the pair-wise distance maps for each of the datasets available: 

1. *Drosophila melanogaster*, dpp locus, nuclear cycle 12
2. *Drosophila melanogaster*, dpp locus, nuclear cycle 14
