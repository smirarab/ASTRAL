DESCRIPTION:
-----------
ASTRAL is a Java program for estimating a species tree given a set of unrooted gene trees. ASTRAL is statistically consistent under multi-species coalescent model (and thus is useful for handling ILS). It finds the tree that maximizes the number of induced quartet trees in the set of gene trees that are shared by the species tree. The algorithm has an exact version that can run for small datasets (less than 18 taxa) and a more useful version that can handle large datasets (103 taxa an 800 genes were analyzed in few minutes).


Warnings:
-----------

ASTRAL is under active development (see the [development branch](https://github.com/smirarab/ASTRAL/tree/development)). There are many features that we are adding to ASTRAL (multi-alleles, unresolved gene trees, etc.). However, the most important one that needs mention is the definition of set X (see our paper), which restricts the search space. 

The significance of set X is that all bipartitions in our species tree are restricted to those that are in X. By default, X is the set of bipartitions in the input gene set. We can build scenarios in simulations where our default strategy for setting X fails to produce good results. We need to enlarge this set in such situations (lots of taxa, too few genes, lots of ILS, too much missing data). There are two ways to do this: automatically using heuristic algorithms, and by using additional trees that might include reasonable hypotheses of the species tree. We are actively working on the first approach and our development branch has new methods that greatly alleviate this potential issue. As a user, you can try finding other trees, anything that might have a chance of containing bipartitions from the species tree, and add those to the set X (using `-e` option). However, if your dataset has one of the following criteria, we **strongly recommend that you test the latest version from the development branch** on your dataset:

* Very high levels of ILS
* Few genes. What is few depends on the number of taxa. 100 genes is quite a lot for 10 taxa, but not much at all for 100 taxa. 
* Lots of missing data. 

For more about this, see [Step 10 of our tutorial on the development branch](https://github.com/smirarab/ASTRAL/blob/development/astral-tutorial.md#step-10-automatic-addition-of-bipartitions-to-x).

Reference and contact
----------

The algorithm used is described in:

S. Mirarab, R. Reaz, Md. S. Bayzid, T. Zimmermann, M.S. Swenson, and T. Warnow1
"ASTRAL: Genome-Scale Coalescent-Based Species Tree Estimation", accepted in ECCB 2014 and to appear in Bioinformatics.

See our [tutorial](astral-tutorial.md) in addition to the rest of this README file. 

Email: `astral-users@googlegroups.com` for questions.



INSTALLATION:
-----------
There is no installation required to run ASTRAL. 
You simply need to download the [zip file](https://github.com/smirarab/ASTRAL/raw/master/Astral.4.4.4.zip) 
and extract the contents to a folder of your choice. Alternatively, you can clone the [github repository](https://github.com/smirarab/ASTRAL/). You can run make.sh to build the project or simply use the jar file that is included with the repository. 

ASTRAL is a java-based application, and should run in any environment (Windows, Linux, Mac, etc.) as long as java is installed. Java 1.5 or later is required. We have tested ASTRAL only on Linux and MAC.

To test your installation, go to the place where you uncompressed ASTRAL, and run:

```
java -jar astral.4.4.4.jar -in test_data/song_primates.50.gene.tre
```

This should quickly finish. There are also other sample input files under `test_data/` that can be used.

ASTRAL can be run from any directories. You just need to run `java -jar /path/to/astral/astral.4.4.4.jar`.
Also, you can move `astral.4.4.4.jar` to any location you like and run it from there, but note that you need
to move the `lib` directory as well. 

EXECUTION:
-----------
ASTRAL currently has no GUI. You need to run it through command-line. In a terminal, go the location where you have downloaded the software, and issue the following command:

```
  java -jar astral.4.4.4.jar
```

This will give you a list of options available in ASTRAL.

To find the species tree given a set of gene trees in a file called `in.tree`, use:

```
java -jar astral.4.4.4.jar -i in.tree
```

The results will be outputted to the standard output. To save the results in a file use the `-o` option:

```
java -jar astral.4.4.4.jar -i in.tree -o out.tre
```

Note that, currently, the input gene trees need to be fully resolved, and each taxon name should appear at most once in each gene tree (missing data are fine). 

### Bootstrapping:

To perform 100 replicates of multi-locus bootstrapping ([Seo 2008](http://www.ncbi.nlm.nih.gov/pubmed/18281270)), use:

```
java -jar astral.4.4.4.jar -i best_ml -b bs_paths 
```

In this command, `bs_paths` is a file that gives the location of gene tree bootstrap files, one line per gene. 
`best_ml` has all the "main" trees (e.g. best ML trees) in one file. 
ASTRAL outputs 100 bootstrapped replicates, then outputs a greedy consensus of the 100 bootstrap replicates, and then outputs one main tree (estimated based on the `best_ml` input) and 
draws support on this main tree using the 100 bootstrap replicates.

Also related to bootstrapping are `-r` (to set the number of replicates), `-g` (to enable gene/site resampling), `-s` (to set the seed number).

### Taxon names:
Leaves of a gene trees need to have distinct labels; these labels will appear as the leaf names in the species tree by default. If the gene trees contain multiple copies of the gene for the same taxa, the current version of ASTRAL is not able to handle the input (something we plant to fix in future), but you can always give different individuals different
names and it should run.  


### Memory:
For big datasets (say more than 100 taxon) increasing the memory available to Java can result in speed up. Note that you should give Java only as much as free memory you have available on your machine. So, for example, if you have 3GB of free memory, you can invoke ASTRAL using the following command to make all the 3GB available to Java:

```
java -Xmx3000M -jar astral.4.4.4.jar -i in.tree
```

Acknowledgment
-----------
ASTRAL code uses bytecode and some reverse engineered code from PhyloNet package (with permission from the authors).


Bug Reports:
-----------
contact ``astral-users@googlegroups.com``
