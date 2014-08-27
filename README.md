DESCRIPTION:
-----------
ASTRAL is a Java program for estimating a species tree given a set of unrooted gene trees. ASTRAL is statistically consistent under multi-species coalescent model (and thus is useful for handling ILS). It finds the tree that maximizes the number of induced quartet trees in the set of gene trees that are shared by the species tree. The algorithm has an exact version that can run for small datasets (less than 18 taxa) and a more useful version that can handle large datasets (103 taxa an 800 genes were analyzed in few minutes).

The algorithm used is described in:

S. Mirarab, R. Reaz, Md. S. Bayzid, T. Zimmermann, M.S. Swenson, and T. Warnow1
"ASTRAL: Genome-Scale Coalescent-Based Species Tree Estimation", accepted in ECCB 2014 and to appear in Bioinformatics.

See our [tutorial](astral-tutorial.md) in addition to the rest of this README file. 

Email: `astral-users@googlegroups.com` for questions.


INSTALLATION:
-----------
There is no installation required to run ASTRAL. 
You simply need to download the [zip file](https://github.com/smirarab/ASTRAL/raw/master/Astral.4.6.3test.2.zip) 
and extract the contents to a folder of your choice. Alternatively, you can clone the [github repository](https://github.com/smirarab/ASTRAL/). You can run make.sh to build the project or simply use the jar file that is included with the repository. 

ASTRAL is a java-based application, and should run in any environment (Windows, Linux, Mac, etc.) as long as java is installed. Java 1.5 or later is required. We have tested ASTRAL only on Linux and MAC.

To test your installation, go to the place where you uncompressed ASTRAL, and run:

```
java -jar astral.4.6.3test.2.jar -in test_data/song_primates.50.gene.tre
```

This should quickly finish. There are also other sample input files under `test_data/` that can be used.

ASTRAL can be run from any directories. You just need to run `java -jar /path/to/astral/astral.4.6.3test.2.jar`.
Also, you can move `astral.4.6.3test.2.jar` to any location you like and run it from there, but note that you need
to move the `lib` directory as well. 

EXECUTION:
-----------
ASTRAL currently has no GUI. You need to run it through command-line. In a terminal, go the location where you have downloaded the software, and issue the following command:

```
  java -jar astral.4.6.3test.2.jar
```

This will give you a list of options available in ASTRAL.

To find the species tree given a set of gene trees in a file called `in.tree`, use:

```
java -jar astral.4.6.3test.2.jar -i in.tree
```

The results will be outputted to the standard output. To save the results in a file use the `-o` option:

```
java -jar astral.4.6.3test.2.jar -i in.tree -o out.tre
```

The input gene trees can have missing taxa, polytommies (unresolved branches), and also multiple individuals per species. When multiple individuals from the same species are available, a mapping file needs to be provided using a `-a` option. This mapping file should have one line per species, and each line needs to be in one of two formats:

```
species_name [number of individuals] individual_1 individual_2 ...

species_name:individual_1,individual_2,...
```

### Bootstrapping:

To perform 100 replicates of multi-locus bootstrapping ([Seo 2008](http://www.ncbi.nlm.nih.gov/pubmed/18281270)), use:

```
java -jar astral.4.6.3test.2.jar -i best_ml -b bs_paths 
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
java -Xmx3000M -jar astral.4.6.3test.2.jar -i in.tree
```

Acknowledgment
-----------
ASTRAL code uses bytecode and some reverse engineered code from PhyloNet package (with permission from the authors).


Bug Reports:
-----------
contact ``astral-users@googlegroups.com``
