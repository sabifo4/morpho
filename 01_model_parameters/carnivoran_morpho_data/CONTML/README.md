# Details of CONTML input and output files

* **`19s.tree` file**:   
   Input file for `CONTML` which contains the tree topology fixed to estimate the morphological tree. It can be found inside [`with_correlation`](https://github.com/sabifo4/morpho/blob/master/01_model_parameters/carnivoran_morpho_data/CONTML/with_correlation)
   and [`without_correlation`](https://github.com/sabifo4/morpho/blob/master/01_model_parameters/carnivoran_morpho_data/CONTML/without_correlation) directories.
* **`Z_19sp_Morpho.aln` file** and :  
   Input file for `CONTML` which contains the alignment with the morphological continuous characters.
   This alignment was used when estimating the ML tree topology for quantitative morphological data when accounting for within-lineage character correlation. 
   This file can be found [here](https://github.com/sabifo4/morpho/blob/master/01_model_parameters/carnivoran_morpho_data/CONTML/with_correlation) .
* **`Ms_19sp_Morpho.aln` file** and :  
   Input file for `CONTML` which contains the alignment with the morphological continuous characters.
   This alignment was used when estimating the ML tree topology for quantitative morphological data when ignoring within-lineage character correlation. 
   This file can be found [here](https://github.com/sabifo4/morpho/blob/master/01_model_parameters/carnivoran_morpho_data/CONTML/without_correlation) .
* **`outtree` file**:  
   Estimated morphological tree output by `CONTML` when fixing the tree topology detailed in the input `19s.tree` file.
* **`outfile` file**:  
   Output log file by `CONTML` after estimating the morphological tree `outtree`.

# Options used in the CONTML software for this analysis

```
Continuous character Maximum Likelihood method version 3.695

Settings for this run:
  U                       Search for best tree?  No, use user trees in input
  L                Use lengths from user trees?  No
  C  Gene frequencies or continuous characters?  Continuous characters
  O                              Outgroup root?  Yes, at species number 15
  M                 Analyze multiple data sets?  No
  0         Terminal type (IBM PC, ANSI, none)?  IBM PC
  1          Print out the data at start of run  No
  2        Print indications of progress of run  Yes
  3                              Print out tree  Yes
  4             Write out trees onto tree file?  Yes

  Y to accept these or type the letter for one to change
Y
contml.exe: can't find input tree file "intree"
Please enter a new file name> 19s.tree

Output written to file "outfile"

Tree also written onto file "outtree"

Done.

Press enter to quit.
```
  
# Estimated morphological tree 

The estimated tree with morphological branches can be represented as it follows:

<p align="center">
  <img width="700" height="500" src="https://github.com/sabifo4/morpho/blob/master/figs/Fig6_MorphoTree.png">
</p>
