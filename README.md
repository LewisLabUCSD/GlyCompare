# GlyCompare
Comparing glycans, leveraging motif substructure arithmetic 

## Paper link
https://www.overleaf.com/16221947rjwhxdqwjbmq


## Definition and nomenclature for variables in Glycompare:
  - Each keyword represents a defined type of variable. 

### Glycan Part:
- **glycan_**
  - type: a glycan object. The object of gly.structure.glycan.glycan
  
- **glycan_dict**
  - type: dict. key is **glycan_id**. glycan_dict[glycan_id]=**glycan_**.
  
- **glycan_list**
  - type: list (of **glycan_**). 
  
- **glycan_id**
  - type: str. String ID used for a **glycan_**.
  
- **glycan_str**
  - type: str. glycoct structure of a **glycan_**. Generated by str(**glycan_**).
  
- **glycan_str_dict**
  - type: str. key is **glycan_id**. glycan_str_dict[glycan_id]=str(**glycan_**).
  
- **mz**
  - type: str. Annotated m/z for a **glycan_**

- **motif_(substructure)**
  - type: a glycan object. The object of gly.structure.glycan.glycan. Just a substructure that breaking down from **glycan_**.
  
- **motif_id** 
  - type: str. String index ID used for a **glycan_**.

- **motif_vec**
  - type: list (of **motif_**). A list of **glycan_** type obj. A list of substructures that breaking down from **glycan_**.

- **motif_dict**
  - type: dict. Key is the degree of **motif_**. motif_dict\[key\] is a **motif_vec**.

- **glycan_motif_dict**
  - type: dict. Key is the glycan_id. glycan_motif_dict\[glycan_id\] is a **motif_dict**.
  
- **match_vec**
  - type: list (of float). length = len(motif_vec), The existency(0/1/1+) within a **glycan_** for each **motif_** in **motif_vec**
  
- **match_mtrix**
  - type: 2D list (of float). List size: number of **glycan_**; Sublist size: len(**motif_vec**); simply a list of **match_vec**
  
### Distinguish the variables that have same type.
    - glycan_1, glycan_2 mean two different **glycan_**
    - motif_1, motif_2 mean two different **motif_**
  
### Glycan IO Part:
    - All glycan_ obj will be stored with glycoct format.
    - GlycoCT: 
    ```
    from glypy.io import glycoct
    glycan_ = glycoct.loads(json_load(target_address)) 
    json_store(target_address, str(glycan_))
    ```
  
### Glycoprofile Part

- **relative_abundance**
  - type: list (of float). Sum to 1. Relative abundance for each annotated glycan in a glycoprofile.

- **weighted_vec**
  - type: list (of float). = match_mtrix.T * relative_abundance
  
- **_class_** **glycan_profile_obj()**:
  - glycan_id_list: 
    - type: list. A list of glycan_id.
  
  - mz_list: 
    - type: list. A list of m/z for each glycan_.
  
  - weighted_vec: 
    - as mentioned
  
  - .relative_abundance: 
    - as mentioned
  
  


## Basic Workflow

### Conversion of glycans to motif space

Given different glycan representations, we will use glypy to translate these into motif space

- GlycoCT: 
```
from glypy.io import glycoct
glycan_ = glycoct.loads(glycan_glycoct) 
```


- IUPAC 
```
from glypy.io import iupac
glycan_ = iupac.loads(glycan_iupac) 
```




#### Demo: Drawing and exporting a glycan in glycoCT format
- Access glycan with GlyTouCan ID in localized database `glytoucan[ID]`
- If no ID, draw your glycan in GlyTouCan.org and export the glycoCT format to a local directory 
- Run ```load_glycan_local()``` with given directory name 
- A processed `NBT_for_motif_extraction.json` json file which contains all the substructures of a glycan will be saved in targeted directory 


### Existing Comparisons in GlyPy
Using glypy we can compare two glycans based on ____?
```
glypy.compare(glycan_1, glycan_2) 
```

### Compare glycans in motif-space

With minimal loss of depth information, simple arithmetic comparisons are now possible in motif space:

- Common Core Extraction: ```commonMotifs = glycan_motif_1*glycan_motif_2```
- Unique elements of glycan1: ```uniqueMotifs = glycan_motif_1-glycan_motif_2```
- Frequency of Motifs: ```motifFq = glycan_motif_1+glycan_motif_2```
- Unions of Motifs: ```motifUnion = (numpy.array(glycan_motif_1+glycan_motif_1)>0).aslist()```

### Visualize Results

Each element of the glycan motif corresponds to a glypy readable glycan. Unique motifs, and common cores can easily be visualized in glypy

```plot_glycan_list( motif_vector , index[commonMotif==1])```
```plot_glycan( motif_vector[ 0 ]```

## Advanced Functionality

### Extract motif to `a_glycan_motif_dict`
Once glycans are converted to `glypy.Glycan` objects we can project them into motif space

```
a_glycan_motif_dict = extract_motif.extract_motif(glycan_obj)
```

For the batch work, please refer:
```
extract_motif.get_motif_pip(gly_len=23, prior=True)
```
### Customize `motif_vector` and match `a_glycan_motif_dict` 
Generally, a customized substructre vector can be generated and all glycans are matched through pipeline
```
customize_motif_vec.customizing_motif_vec_pip()
motif_lib = motif_class.GlycanMotifLib(json.load(output_motif_dic_degree_list_addr)) # unicarbkb_motifs_12259.json
motif_lib.motif_vec
```
#### Customize `motif_vector`
- customize a motif_vector, please refer `customize_motif_vec.customizing_motif_vec_pip()`
```
customize_motif_vec.get_motif_dict_degree_list_pipe(NBT_glycan_dict_degree_list_glycoct_for_motif,
                                                                output_motif_dic_degree_list_addr)
```
#### Match `a_glycan_motif_dict` to `motif_vector`
- Glycan with extracted motifs can be matched to a provided motif_vector using 
```
match_dict[glycan_name]={}
customize_motif_vec.match_motif(motif_vector, a_glycan_motif_dict, glycan_name, match_dict, index=0)
```
- For the batch work, please refer `customize_motif_vec.customizing_motif_vec_pip()`
```
customize_motif_vec.motif_matching_wrapper(NBT_motif_dic_degree_list,
                                                        NBT_glycan_dict_degree_list_glycoct_for_motif,
                                                        output_matched_glycan_addr)
```
### Motif-based differential glycomics workflow

We have implimented functions for reading glycoprofiles, glycan profile data include 
- glycan structural data
A Glycan structural profile is like this:
`{'name1':{'m/z':'GlyTouCanID or customized ID'}}`
For example:
```
{'Gly01': {'2244': 'G04483SK',
           '4587': '4587.1',
           '5037': 'G49604DB',
           '5486': '5486.1'}}
```
- glycan abundance data table 
```
{'m/z_1':abundance_list1,
 'm/z_2':abundance_list2,
 'm/z_3':abundance_list3,
}
```
- example for a 5-profile abundance table
```
{'2244': [0, 0, 0, 0, 0],
 '4587': [0, 0, 0, 0, 0],
 '5037': [0, 0, 0, 0, 0],
 '5486': [0, 0, 0, 0, 0],}
```
### the example code
```
abundance_data_table = load_json(__init__.json_address + "NBT_dict_name_abundance_cross_profile.json")
merged_glycan_profile, _ = glycan_profile.load_glycan_profile_dic()
glycan_profiles_obj = glycan_profile.combine_profile_mz_with_motif_existance(merged_glycan_profile,
                                                                                  abundance_data_table)


```
## comparison across enzyme isoform
# extract motifs
glycanMotifMatrix1 = np.matrix( [mapMotifs( extractMotifs( g1 ) , motifsVector = motifVector ) for g1 in glycanList1 ] )
glycanMotifMatrix2 = np.matrix( [mapMotifs( extractMotifs( g2 ) , motifsVector = motifVector ) for g2 in glycanList2 ] )
# abunance weighted motif vector
glycanMotifVectorAbundance1 = glycanMotifMatrix1%*%abundance1
glycanMotifVectorAbundance2 = glycanMotifMatrix1%*%abundance2
# all motif vector
glycanMotifVector1 = numpy.sum( glycanMotifMatrix1 , axis=1)
glycanMotifVector2 = numpy.sum( glycanMotifMatrix2 , axis=1)
```
Now that we have mapped each glycan in the glycoprofile to our motif vector, we can search for the common core WITHIN our profiles individually. 
```
commonMotifs_glycoprofile1 = numpy.prod( glycanMotifMatrix1 , axis=1 )
commonMotifs_glycoprofile2 = numpy.prod( glycanMotifMatrix2 , axis=1 )
```
We can also find the common core BETWEEN our profiles
```
numpy.prod( commonMotifs_glycoprofile1 , commonMotifs_glycoprofile2 )
```

### Motif-based glycan clustering
We can also cluster glycoprofiles based on their motif contents:

- Cluster glycoprofiles on motif-abundance vector
```g = sns.clustermap( [glycanMotifVectorAbundance1 , glycanMotifVectorAbundance2] , metric='braycurtis')```
- Cluster glycoprofiles on common motif vector
```g = sns.clustermap( [commonMotifs_glycoprofile1 , commonMotifs_glycoprofile2] , metric='braycurtis')```

### Extact common motifs within a cluster
under construction

### Construct Dependency Structure for New Mofif Vector
```
glycan_distance = np.matrix(1,len(motifVector),len(motifVector))
for g1 in mofitVector:
  for g2 in motifVector:
    glypy.compare(g1,g2)
# contruct hierarchy from distance matrix
```

### Motif Statistics 

- Basic enrichment
```
ben glycoprofile test....whatever
```
- Motif Enrichment with heirarchical dependency 
```
glycoprofile compare with hierarchy
```
