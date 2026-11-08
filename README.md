# StereoDist
## This is a fork

> [!NOTE] 
> Original code at https://bitbucket.org/oesperlab/stereodist/src/master/
> 
>This is an adaptation to allow reading Newick tree files originating from the HDC-Phylo code database ( https://github.com/cthadeusantos/hdc-phylo ).
>
>This adaptation also modified a "gap" in the original code, so that when there are fewer than two comparable mutations, CASet and DISC now return NaN instead of 0. Thus, these cases do not appear as distances computed as zero.
>In the .npy matrices, undefined values remain as NaN. In the plots, they appear in gray, with the legend "no comparable pairs"; the color scale remains fixed between 0 and 1. The sorted plot uses only the defined values and reports how many undefined comparisons were excluded. The minimum and maximum calculations also ignore undefined values.

### Original text from  https://bitbucket.org/oesperlab/stereodist/src/master/

CASet and DISC distance implementations for cancer evolutionary trees.

Trees should be input in Newick format:
```
      GPR158
      /    \
  MAP2K1  PLA2G16    -->    (MAP2K1,(LRRC16A)PLA2G16)GPR158;
             \
            LRRC16A            
```

Multi-labeled nodes use `{Label1,Label2}`:
```
   GPR158,SLC12A1
      /       \
  MAP2K1    PLA2G16,EXOC6B     -->     (MAP2K1,(LRRC16A){PLA2G16,EXOC6B}){GPR158,SLC12A1};
                \
               LRRC16A            
``` 

Trees should be in a single text file with one tree per line, ending in a semicolon. Any line beginnning with a '#' will be treated as a comment. Both distance measures will run using the intersection distance by default, only averaging over pairs of mutations appearing in both trees. 

To run CASet distance: 
```
python3 CASet.py SampleTrees.txt [-ouptm]
```

To run DISC distance:
```
python3 DISC.py SampleTrees.txt [-ouptm]
```

positional arguments:
```
  inputFile
```

optional arguments:
```
  -h, --help            show this help message and exit
  -o OUTPUTFILE, --outputFile OUTPUTFILE
  -u, --union
  -p PICKLE, --pickle PICKLEFILE
  -t, --treePrint
  -m, --minmax
```

Using the pickle option will write a pickled matrix to the specified file (see https://docs.python.org/3/library/pickle.html). This n-by-n matrix lists the pairwise distance between each tree where the [i,j]th entry is the distance between trees i and j in the Newick input file. 

Sample run on the trees above to produce a pickled matrix: 
```
cd stereodist
python3 CASet.py SampleTrees.txt -u -p CASetSample.pickle
```
Unpickling the file from above gives a matrix showing that the distance between the two trees is as follows: 
``` 
0.0 0.8
0.8 0.0
``` 
