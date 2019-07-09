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


Pipeline guide.

### Create a new session and all data are required in ,/source_data/
- An abundance table (1.1) is necessary to run the full GlyCompare pipeline, if the glytoucan_id for each glycan is specified. For example, [iscience_data](https://github.com/LewisLabUCSD/GlyCompare/blob/master/example_data/test_iscience/source_data/abundace_table.csv). Complex data need more, [epo_data](https://github.com/LewisLabUCSD/GlyCompare/blob/master/example_data/paper_epo/source_data/). 
- 1.1. abundance_table.xls. An abundance_table has column: sample (glycoprofile), row: glycan (the glycan identifier can be glytoucan_id or custimized: m/z, hplc). 
    
    

- 1.2. (optional) external_profile_naming.json. Glycompare defaults to index names for each sample (e.g. 1,2,3). If you want to specify a name for each sample use. i.e. [paper_epo/source_data/external_profile_naming.json](https://github.com/LewisLabUCSD/GlyCompare/blob/master/example_data/paper_epo/source_data/external_profile_naming.json)
    
    

- 1.3. (optional) glycan_identifier_to_glytoucan_id.json. If the glycan identifiers in abundance_table.xls are not glytoucan_id (if uses m/z or costumized id). See 2.4. [paper_epo/source_data/glycan_identifier_to_glytoucan_id.json](https://github.com/LewisLabUCSD/GlyCompare/blob/master/example_data/paper_epo/source_data/glycan_identifier_to_glytoucan_id.json)
    
    
- 1.4. (optinal) source_data/glycoct/. If part of glycan structures are manually curated, we need a source_data/glycoct/ directory to store all of them.
     i.e. [paper_epo/source_data/glycoct](https://github.com/LewisLabUCSD/GlyCompare/tree/master/example_data/paper_epo/source_data/glycoct)


### There are several parameter should be set up first
    2.1. working_addr : root working dir

    2.2. project_name: usually same as the folder of root

    2.3. __init__.num_processors: number of processes needed

    2.4. complex_naming: 
        True if they use the mz, hplc naming or self-created naming rather than glytoucan_id in abundance table.
        
    2.5. complex_profile: 
        False, if glycan identifiers in abundance_table.xls is not glytoucan_id.
        
        True, if a glycan identifier, for example m/z, in abundance_table.xls has isomers, the m/z-structure matching 
        must be specified in glycan_identifier_to_glytoucan_id.json (1.3)
        
    #    glycan_profile['5'] = {
    #         "2244": "G04483SK",
    #         "2605": "G30460NZ",
    #         "2967": "G17689DH",
    #         "3416": "3416.1",
    #     }
    #    glycan_profile['6'] = {
    #         "1591": "G07483YN",
    #         "1836": "G80858MF",
    #         "2081": "G80393PG",
    #         "2326": "G00176HZ",
    #         "2530": "G79457WN",
    #         "2892": "G79412GP",
    #         "3096": "G40242TG", 
    #         "3457": "3457.1",
    #         "3777": "G76812VG",
    #         "4022": "G80264ZA",
    #         "4226": "G86696LV",
    #         "4587": "4587.1",
    #         "5037": "G49604DB",
    #     }
        False, if there is no ambiguity between abundance table index and glycan structure ID.

    
    2.6. external_profile_naming: as 1.2.
        True, 
        False

    2.7. structure_loader: for loading glytoucan id, could be:
        list of glycan_id
        glycan_dict

    2.8. glytoucan_db_addr: the addr for glytoucan database if needed

    # meta_name: a gff type table which includes glycan's naming information
