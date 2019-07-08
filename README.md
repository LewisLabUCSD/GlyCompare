# GlyCompare

Here, we present GlyCompare, a novel method wherein glycans from glycomic data are decomposed to a minimal set of intermediate substructures, thus incorporating shared intermediate glycan substructures into all comparisons of glycans. 

## Citation
Bokan Bao+, Benjamin P. Kellman+, Austin W. T. Chiang, Austin K. York, Mahmoud A. Mohammad, Morey W. Haymond, Lars Bode, and Nathan E. Lewis. 2019. “**Correcting for Sparsity and Non-Independence in Glycomic Data through a System Biology Framework.**” bioRxiv. https://doi.org/10.1101/693507.

### Demonstrations & Manuscript Analyses
- [All demos](https://github.com/LewisLabUCSD/GlyCompare/tree/master/example_notebook)
- [Figure 2 & 3](https://github.com/LewisLabUCSD/GlyCompare/blob/master/example_notebook/Fig2_Fig3_epo_analysis.ipynb)
- [Figure 4 & 5](https://github.com/LewisLabUCSD/GlyCompare/blob/master/example_notebook/Fig4_Fig5_hmo_analysis.ipynb)

![workflow](GlyCompare_flow.png)

## Disclaimer:

The GlyCompare framework provides several tools that account for the influence of the glycan substructure network in the analysis of glycomic data. However, for its effective use, there are two primary requirements for the data it processes. First, the tools require that the glycan is stored in a tree-like structure. Thus, neither cyclic nor glycan with undefined topology in glycoCT format can be processed. Second, during the substructure matching, in terms of the linkage specificity, the code can only handle two types of substructure analysis. One has the exact linkage specification, and one ignores all linkage specification and only accounts for topology. Currently, our tool cannot handle partial ambiguity in linkages of a glycan. The code and the manual are freely available and will be continually developed to enable its accessibility to all types of scientists. 


'pip glycompare'
