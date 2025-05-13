# Phylo-pPCD

This is a lite version of [PhyloProfile](https://github.com/BIONF/PhyloProfile), which is pre-configured only for its use with the data of 235 [potential plant cell wall degrading (pPCD) enzymes across the tree of life](https://applbio.biologie.uni-frankfurt.de/download/fDOG_pPCD/cellulase_pp_data.tar.gz).

Users can display the full data set, customized subsets of the pPCD enzymes and taxa, zoom in on individual ortholog pairs, and ultimately display the domain architectures of the orthologs. Interactive links connect the information about taxon, protein sequence and Pfam (El-Gebali, et al. 2019) or SMART domains (Letunic, et al. 2009) with the corresponding public databases.

# Table of Contents
* [Usage](#usage)
* [Bugs](#bugs)
* [License](#license)
* [Contact](#contact)

# Usage
Feature-aware phylogenetic profile of the 235 pPCD enzymes across ~18,000 taxa from three domains of the tree of life.

* Taxa (x-axis) are ordered with increasing taxonomic distance to *Rhizoctonia solani*.
* Rows indicate the pPCD proteins.
* Dots indicate orthologs of each protein in the respective taxon. To inspect the domain architecture comparison between the protein in reference species (e.g. *R. solani*) and the one from the datapoint, click on the dot and select Domain plot.
* The color scheme represents the similarity of annotated features between the seed proteins and the respective orthologs using the FAS scores ([Dosch, et al. 2023](https://doi.org/10.1093/bioinformatics/btad226)) by assigning values between 0 and 1. The score is penalized for missing features (FAS_F; dot color) and added features (FAS_B; background color) in the orthologous protein. For instance, a protein with an identical feature architecture to the reference protein in reference species would have FAS_F = 1 and FAS_B = 1, and would appear as a blue dot in a white background according to the color scale.
* For further details on the biological interpretation of the data, please refer to the article (link will be soon available).

# Standalone installation

From an R terminal, PhyloCellulase can be installed locally using *devtools*:

```r
if (!requireNamespace("devtools"))
    install.packages("devtools")
devtools::install_github("BIONF/PhylopPCD", INSTALL_opts = c('--no-lock'), dependencies = TRUE)
```

Then, to run it, enter:

```r
library(PhylopPCD)
runPhylopPCD()
```

Check your web browser, *Phylo-pPCD* will be displayed there ;-) 

For the first time running, the tool will download some pre-calculated data. Please be patient until you see a the profile plot.

# Bugs
Any bug reports or comments, suggestions are highly appreciated. Please [open an issue on GitHub](https://github.com/BIONF/PhylopPCD/issues/new) or be in touch via email.

# License
This tool is released under [MIT license](https://github.com/BIONF/PhylopPCD/blob/master/LICENSE).

# Contact
Ingo Ebersberger
ebersberger@bio.uni-frankfurt.de
