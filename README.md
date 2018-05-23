# GlyCompare
Comparing glycans, leveraging motif substructure arithmetic 

## Paper link
https://www.overleaf.com/16221947rjwhxdqwjbmq

## Basic Workflow

### Conversion of glycans to motif space

Given different glycan representations, we will use glypy to translate these into motif space

- GlycoCT: 
```
from glypy.io import glycoct
glycoct.loads(glycan_glycoct) 
```


- IUPAC 
```
from glypy.io import iupac
iupac.loads(glycan_iupac) 
```

#### Demo: Drawing and exporting a glycan in glycoCT format
- Access glycan with GlyTouCan ID in localized database `glytoucan[ID]`
- If no ID, draw your glycan in GlyTouCan.org and export the glycoCT format to a local directory 
- Run ```load_glycan_local``` with given directory name 
- A processed json file which contains all the substructures of a glycan will be saved in targeted directory 

Once glycans are converted to glypy objects we can project them into motif space

```extract_motif(glycan_obj)```

Extracted motifs can be mapped to a provided motif vector using ```match_motif()```

```
motif_lib = motif_class.GlycanMotifLib(json.load('motifs.json')) # unicarbkb_motifs_12259
motif_lib.motif_vec
```
or 
```
customizing_motif_vec_pip()
    get_motif_dict_degree_list_pipe()
    motif_matching_wrapper()
```
### Existing Comparisons in GlyPy
Using glypy we can compare two glycans based on ____?
```
glypy.compare( glycan1,glycan2) 
```

### Compare glycans in motif-space

With minimal loss of depth information, simple arithmetic comparisons are now possible in motif space:

- Common Core Extraction: ```commonMotifs = glycanMotifs1*glycanMotifs2```
- Unique elements of glycan1: ```uniqueMotifs = glycanMotifs1-glycanMotifs2```
- Frequency of Motifs: ```motifFq = glycanMotifs1+glycanMotifs2```
- Unions of Motifs: ```motifUnion = (numpy.array(glycanMotifs1+glycanMotifs2)>0).aslist()```

### Visualize Results

Each element of the glycan motif corresponds to a glypy readable glycan. Unique motifs, and common cores can easily be visualized in glypy

```plot_glycan_list( motifVector , index[commonMotif==1])```
```plot_glycan( motifVector[ 0 ]```

## Advanced Functionality

### Construct a custom motif vector
- For every glycan in your database, run ```extractMotifs```
- Combine all extracted motifs into a global vector using ```mergeMotifs```
```
motifVector = mergeMotifs( [ extractMotifs(g) for g in glycans ] )
```
### Motif-based differential glycomics workflow

We have implimented functions for reading glycoprofiles
A Glycan profile is like this:
`{'name':{'m/z':'GlyTouCanID'}}`
For example:
```
{'Gly01': {'2244': 'G04483SK',
                                '4587': '4587.1',
                                '5037': 'G49604DB',
                                '5486': '5486.1'}}
```
glycanList1,abundance1 = read_glycoprofile( 'glycoprofile1.dic' )
glycanList2,abundance2 = read_glycoprofile( 'glycoprofile2.dic' )
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
