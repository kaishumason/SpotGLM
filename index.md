# SpotGLM
SpotGLM is an R package for performing **cell type-specific differential testing** in **spatial omics** data. 
It adapts to different data modalities, such as **gene expression**, **chromatin accessibility**, or **isoform usage**, 
and identifies cell type specific changes associated with local tissue context. SpotGLM accounts for the 
mixed-cell composition inherent to spatial barcoding technologies.

SpotGLM is integrated with the [SPARROW package](https://kaishumason.github.io/SPARROW/), making it **fast and ultra-scalable** for data sets with millions of spots.

![Schematic of SpotGLM Analysis](man/figures/schematic.png)

## Installation

To install SpotGLM from GitHub:

```
install.packages("devtools")
devtools::install_github("kaishumason/SpotGLM") # install
```

We strongly recommend using SpotGLM with SPARROW, which is a power-preserving data reduction method which can speed up analyses.  SPARROW can also be installed from GitHUB:
```
install.packages("devtools")
devtools::install_github("kaishumason/SPARROW") # install
```


## Getting Started

Please follow these tutorials to get started on your data:

[Using SpotGLM with SPARROW on Visium HD](articles/Vignette_VisiumHD_Mouse_Kidney_analysis.html)

[Using SpotGLM on Visium data](articles/Visium_analysis.html)

[SpotGLM for isoform switching analysis on spatial long read data](articles/Spatial_Long_Read_analysis.html)

[SpotGLM for spatial ATAC analysis](articles/Spatial_ATAC_analysis.html)

Here is a more technical tutorial that applies spotGLM to simulation data.  By working with simulation data, you can see the generative model used in spotGLM:

[Demonstration of SpotGLM on simulation data](articles/Intro_to_SpotGLM.html)
