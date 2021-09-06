---
output: github_document
---

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

<style>
body {
text-align: justify}
</style>

# Handle Ambiguity of Protein Identifications from Shotgun Proteomics  


**Analyze ambiguous protein identifications using graph connected components (CCs)** 
--- 


Protein inference is a central issue in proteomics, given the presence of shared peptides (*i.e.*, peptides that might originate from different proteins sharing homology, from different proteoforms due to alternative mRNA splicing, post-translational modifications, proteolytic cleavages, and/or allelic variants). Indeed, in bottom-up mass spectrometry-based proteomics, the most widely used proteomic approach, peptide-protein connectivity is lost for experimental reasons and protein identifications are to be inferred from peptide identifications. Shared peptides can generate quite complex peptide-to-protein mapping structures but these can be efficiently represented using bipartite graphs, with peptides
and proteins as vertices and with edges featuring peptide to protein membership.
Graph connected components (CCs) (*i.e.*, the largest subgraphs in which any two vertices are connected to each other by a path and not connected to any other of the vertices in the supergraph) can be used as a mesure of the level of ambgiuty in protein identifications.  
CCs represent a peptide-centric strategy to group proteins and it is independent from the variety of protein_centric strategies of protein
grouping and protein inference.  As such, it does not require protein inference and it is widely applicable, reproducible and transparent.   
The CCs4prot package allows to build graph from shotgun proteomic identifications and calculate its connected component.  


**Reduce ambiguity of protein identifications by transcriptome-informed filtering** 
--- 
  The availability of an increasing number of matched proteomic and transcriptomic datasets can be exploited to reduce ambiguity of protein quantification. 
Indeed, according to the central dogma of biology, there can be no protein without the corresponding transcript. Following this, proteins identifications
for which the corresponding transcript is identified in the sample-matched transcriptome are more likely to be correct than protein identifications with no 
expressed transcript.  
The CCs4prot package implements a transcriptome-informed filtering strategy and allows to measure the impact of the filtering on ambiguity of protein identifications.  


## Install the CCs4prot R package

Download the package with the git clone command:

```bash
git clone https://github.com/laurafancello/CCs4prot.git
```

Initiate R and install the R package using devtools (devtools needs to be installed as well)

```{r}
library("devtools")
devtools::install("CCs4prot")
```


## Usage

To learn how to use CCs4prot, please refer to the introductory vignette posted at this link: 

* [https://github.com/laurafancello/CCs4prot/blob/main/vignettes/IntroToCCs4prot.Rmd)

