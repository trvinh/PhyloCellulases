# PhyloCellulase

This is a lite version of [PhyloProfile](https://github.com/BIONF/PhyloProfile), which is pre-configured only for its use with the data of ~~303 [eukaryotic ribosome biogenesis factors](https://applbio.biologie.uni-frankfurt.de/download/RibosomeBiogenesis/PhyloRBF/)~~. The phylogenetic profiles for proteins in the ~~RBFeuk set across 982 taxa provide the first unifying resource for tracing eukaryotic ribosome biogenesis factors across the organismal diversity~~.

Users can display the full data set, customized subsets of cellulases and taxa, zoom in on individual ortholog pairs, and ultimately display the domain architectures of the orthologous cellulases. Interactive links connect the information about taxon, protein sequence and Pfam (El-Gebali, et al. 2019) or SMART domains (Letunic, et al. 2009) with the corresponding public databases.

# Table of Contents
* [Usage](#usage)
* [Bugs](#bugs)
* [License](#license)
* [Contact](#contact)

<!-- ![](https://github.com/BIONF/PhyloProfileCorona/blob/gh-pages/www/phyloprofile_interface.png) -->
![](phyloprofile_interface.png)

# Usage
Feature-aware phylogenetic profile of the 280 yeast and 23 human ribosome biogenesis factor genes across 982 taxa (727 Archaea, 23 Bacteria, 232 Eukaryotes).
* Taxa (x-axis) are ordered with increasing taxonomic distance to yeast *Saccharomyces cerevisiae*.
* Rows indicate the ribosome biogenesis factors.
* Dots indicate orthologs of each protein in the respective taxon. To inspect the domain architecture comparison between the protein in yeast and the one from the datapoint, click on the dot, select detailed plot, select the ortholog from the barplot, and click on show domain architecture.
* The color scheme represents the similarity of annotated features between yeast proteins and the respective orthologs using the FAS scores ([Koestler, et al. 2010](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-417)) by assigning values between 0 and 1. The score is penalized for missing features (FAS_F; dot color) and added features (FAS_B; background color) in the orthologous protein. For instance, a protein with an identical feature architecture to the reference protein in yeast would have FAS_F = 1 and FAS_B = 1, and would appear as a blue dot in a white background according to the color scale.
* For further details on the biological interpretation of the data, please refer to the article (link will be soon available).

# Standalone installation

From an R terminal, PhyloCellulase can be installed locally using *devtools*:

```r
if (!requireNamespace("devtools"))
    install.packages("devtools")
devtools::install_github("BIONF/PhyloCellulase", INSTALL_opts = c('--no-lock'), dependencies = TRUE, build_vignettes = TRUE)
```

Then, to run it, enter:

```r
library(PhyloCellulase)
runPhyloCellulase()
```

Check your web browser, *PhyloCellulase* will be displayed there ;-) For the first time running, the tool will download some pre-calculated data. Please be patient until you see a the profile plot.

# Bugs
Any bug reports or comments, suggestions are highly appreciated. Please [open an issue on GitHub](https://github.com/BIONF/PhyloCellulase/issues/new) or be in touch via email.

# License
This tool is released under [MIT license](https://github.com/BIONF/PhyloCellulase/blob/master/LICENSE).

# Contact
Ingo Ebersberger
ebersberger@bio.uni-frankfurt.de
