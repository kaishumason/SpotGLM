# SpotGLM
SpotGLM is an R package for performing **cell type-specific differential testing** in **spatial omics** data. 
It adapts to different data modalities, such as **gene expression**, **chromatin accessibility**, or **isoform usage**, 
and identifies cell type specific changes associated with local tissue context. SpotGLM accounts for the 
mixed-cell composition inherent to spatial barcoding technologies.

SpotGLM is integrated with the [SPARROW package](https://kaishumason.github.io/SPARROW/), making it **fast and ultra-scalable** for data sets with millions of spots.

![](man/figures/schematic.png)

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

[Using SpotGLM with SPARROW on Visium HD (kidney example)](articles/Vignette_VisiumHD_Mouse_Kidney_analysis.html)

[Using SpotGLM on Visium data (colorectal cancer example)](articles/Visium_analysis.html)

[SpotGLM for isoform switching analysis on spatial long read data (olfactory bulb example)](articles/Spatial_Long_Read_analysis.html)

[SpotGLM for spatial ATAC analysis (mouse brain example)](articles/Spatial_ATAC_analysis.html)

Here is a more technical tutorial that applies spotGLM to simulation data.  By simulating the data yourself, you can gain understanding of the underlying model:
[Demonstration of SpotGLM on simulation data](articles/Intro_to_SpotGLM.html)
