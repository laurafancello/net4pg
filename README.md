---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

`​``{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
`​``

# Handle Ambiguity of Protein Identifications from Shotgun Proteomics  

 <p style="text-align: center;">[ [GitLab](https://gitlab.com/laura_fancello/ccs4prot) ] 
</p>

**Analyze ambiguous protein identifications using graph connected components (CCs)** 
--- 

  Ambiguity of protein identifications is an important issue in shotgun proteomics and it is due to the
presence of shared peptides. Shared peptides are very frequent in higher eukaryotes and can generate quite complex peptide-to-protein
mapping structures. These structures can be efficiently represented using bipartite graphs, with peptides
and proteins as vertices and with edges featuring peptide to protein membership. The graph-based representation also facilitates the 
assessment and quantification of ambiguity in protein identifications, by means of graph theory and, in particular, using graph connected 
components (CCs). Graph CCs are the largest subgraphs in which any two vertices are connected to each 
other by a path and not connected to any other of the vertices in the supergraph. Proteins sharing one or more peptides are
thus gathered in the same CC (multi-protein CCs), while unambiguous protein identifications are represented by CCs with a single
protein vertex (single-protein CCs). The proportion of multi-protein versus single-protein CCs and the size (i.e. number of protein members) of
multi-protein CCs can be used as a measure of the ambiguity of protein identifications. 
CCs represent a peptide-centric strategy to group proteins and it is independent from the variety of protein_centric strategies of protein
grouping and protein inference.  As such, it does not require protein inference and it is widely applicable, reproducible and transparent. 


**Reduce ambiguity of protein identifications by transcriptome-informed post-hoc filtering** 
--- 
  The availability of an increasing number of matched proteomic and transcriptomic datasets can be exploited to reduce ambiguity of protein quantification. 
Indeed, according to the central dogma of biology, there can be no protein without the corresponding transcript. Following this, proteins identifications
for which the corresponding transcript is identified in the sample-matched transcriptome are more likely to be correct than protein identifications with no 
expressed transcript.
Based on this rationale we propose a  a transcriptome-informed post-hoc filtering strategy to reduce ambiguity of protein identifications. This strategy 
consists in the removal of proteins whose transcript is not identified in the sample-matched transcriptome. This can be operated in two ways: 
1. remove all proteins with no expressed transcript and all those peptides which only map on these proteins;
2. remove only those proteins with no expressed transcript which are exclusively identified by shared peptides. 
This latter option is more cautious as it does not cause the removal of any peptide identification and it filters out only ambiguous protein identifications (i.e. 
proteins with shared peptides). However, it can be interesting to investigate whether peptide identifications removed using the first option are of lower quality.
---

## Install ccs4prot as an R package

ccs4prot is as an **R package**. 

Download the package with the git clone command:

```bash
git clone https://gitlab.com/laura_fancello/ccs4prot.git
```

Initiate R and install the R package using devtools (devtools needs to be installed as well)

```{r}
library("devtools")
devtools::install("ccs4prot")
```


## Usage

To learn how to use ccs4prot, please refer to the introductory vignette posted at this link: 

* [https://gitlab.com/laura_fancello/ccs4prot/vignettes/IntroToCCs4prot.Rmd]https://gitlab.com/laura_fancello/ccs4prot/vignettes/IntroToCCs4prot.Rmd)

