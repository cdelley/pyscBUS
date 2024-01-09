# BUS format module for python

Functions to read and write uncompressed BUS files with python. The main use for this module
is as submodules in code realted to single-cell experiment preprocessing.
    
The file format was introduced by Páll Melsted, Sina Booeshaghi, Lior Pachtr
and coworkers for kallisto.

*Melsted, Páll, Booeshaghi, A. Sina et al. Modular and efficient pre-processing
of single-cell RNA-seq. BioRxiv (2019): 673285, doi.org/10.1101/673285.*   

See file specifications at:

https://github.com/BUStools/BUS-format
https://github.com/BUStools/bustools

Compressed BUS format is currently not supported. These need to be created (or decompressed)
with the BUStools. Check for details:

*Einarsson, Pétur Helgi, and Páll Melsted. “BUSZ: compressed BUS files.” 
Bioinformatics (Oxford, England) vol. 39,5 (2023): btad295. doi:10.1093/bioinformatics/btad295*

## Deprecated Compressed File Format From My Earlier Projects

In adition, contains functions to write a condensed file format that I designed and used
in my earlier projects. However, I recommend to use the BUS format in future projects.
