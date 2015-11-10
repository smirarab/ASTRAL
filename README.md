DESCRIPTION:
-----------
ASTRAL is a Java program for estimating a species tree given a set of unrooted gene trees. ASTRAL is statistically consistent under multi-species coalescent model (and thus is useful for handling ILS). It finds the tree that maximizes the number of induced quartet trees in the set of gene trees that are shared by the species tree. The algorithm has an exact version that can run for small datasets (less than 18 taxa) and a more useful version (its default) that can handle large datasets (tested for up to 1000 taxa and 1000 genes).

The algorithm used is described in:

* Mirarab, Siavash, Rezwana Reaz, Md. Shamsuzzoha Bayzid, Theo Zimmermann, M Shel Swenson, and Tandy Warnow. “ASTRAL: Genome-Scale Coalescent-Based Species Tree.” Bioinformatics (ECCB special issue) 30, no. 17 (2014): i541–i548. [doi:10.1093/bioinformatics/btu462](doi.org/10.1093/bioinformatics/btu462).
* Mirarab, Siavash, Tandy Warnow. “ASTRAL-II: Coalescent-Based Species Tree Estimation with Many Hundreds of Taxa and Thousands of Genes.”. Bioinformatics (ISMB special issue) 31, no. 12 (2015): i44–i52. [doi:10.1093/bioinformatics/btv234](http://bioinformatics.oxfordjournals.org/content/31/12/i44)

The code given here corresponds to ASTRAL-II.

See our [tutorial](astral-tutorial.md) in addition to the rest of this README file. Also, the chapter of my dissertation that describes ASTRAL in detail is provided [here](thesis-astral.pdf).

Email: `astral-users@googlegroups.com` for questions.

INSTALLATION:
-----------
There is no installation required to run ASTRAL. 
You simply need to download the [zip file](https://github.com/smirarab/ASTRAL/raw/master/Astral.4.8.0.zip) 
and extract the contents to a folder of your choice. Alternatively, you can clone the [github repository](https://github.com/smirarab/ASTRAL/). You can run `make.sh` to build the project or simply use the jar file that is included with the repository. 

ASTRAL is a java-based application, and should run in any environment (Windows, Linux, Mac, etc.) as long as java is installed. Java 1.5 or later is required. We have tested ASTRAL only on Linux and MAC.

To test your installation, go to the place where you uncompressed ASTRAL, and run:

```
java -jar astral.4.8.0.jar -i test_data/song_primates.424.gene.tre
```

This should quickly finish. There are also other sample input files under `test_data/` that can be used.

ASTRAL can be run from any directories. You just need to run `java -jar /path/to/astral/astral.4.8.0.jar`.
Also, you can move `astral.4.8.0.jar` to any location you like and run it from there, but note that you need
to move the `lib` directory as well. 

EXECUTION:
-----------
ASTRAL currently has no GUI. You need to run it through command-line. In a terminal, go the location where you have downloaded the software, and issue the following command:

```
  java -jar astral.4.8.0.jar
```

This will give you a list of options available in ASTRAL.

To find the species tree given a set of gene trees in a file called `in.tree`, use:

```
java -jar astral.4.8.0.jar -i in.tree
```

The results will be outputted to the standard output. To save the results in a file use the `-o` option (**Strongly recommended, unless you are using a pipeline**):

```
java -jar astral.4.8.0.jar -i in.tree -o out.tre
```

The input gene trees can have missing taxa, polytommies (unresolved branches), and also multiple individuals per species. When multiple individuals from the same species are available, a mapping file needs to be provided using a `-a` option. This mapping file should have one line per species, and each line needs to be in one of two formats:

```
species_name [number of individuals] individual_1 individual_2 ...

species_name:individual_1,individual_2,...
```

The code for handling multiple individuals is in its infancy and might not work well yet. Keep posted for improvements to this feature. As of July, 2015, we *strongly* recommend
that you use the evolving code in the [multiind](https://github.com/smirarab/ASTRAL/tree/multiind) branch for multi individuals.

### Bootstrapping:

To perform 100 replicates of multi-locus bootstrapping ([Seo 2008](http://www.ncbi.nlm.nih.gov/pubmed/18281270)), use:

```
java -jar astral.4.8.0.jar -i best_ml -b bs_paths -r 100
```

In this command, `bs_paths` is a file that gives the location of gene tree bootstrap files, one line per gene. 
`best_ml` has all the "main" trees (e.g. best ML trees) in one file. 

##### Bootstrap Output: 

The output file generated when using the bootstrapping feature with 100 replicates (`-r 100`) contains the following trees, in this order:

* 100 bootstrapped replicate trees; each tree is the result of running ASTRAL on a set of bootstrap gene trees (one per gene).
* A greedy consensus of the 100 bootstrapped replicate trees; this tree has support values drawn on branches based on the bootstrap replicate trees. Support values show the percentage of bootstrap replicates that contain a branch. 
* The “main” ASTRAL tree; this is the results of running ASTRAL on the `best_ml` input gene trees. This main tree also includes support values, which are again drawn based on the 100 bootstrap replicate trees.

If `-r` option is set to anything other than 100, the number of replicates would be accordingly adjusted.    
**Note** that by default (i.e., when no `-r` is given), ASTRAL only performs 100 replicates regardless of the number of replicates in your bootstrapped gene trees.
If you want to bootstrap with a different number of replicates, you must use `-r`. 

Also related to bootstrapping are `-g` (to enable gene/site resampling) and `-s` (to set the seed number) options. 


### Memory:
For big datasets (say more than 100 taxon) increasing the memory available to Java can result in speed ups. Note that you should give Java only as much free memory as you have available on your machine. So, for example, if you have 3GB of free memory, you can invoke ASTRAL using the following command to make all the 3GB available to Java:

```
java -Xmx3000M -jar astral.4.8.0.jar -i in.tree
```

Acknowledgment
-----------
ASTRAL code uses bytecode and some reverse engineered code from PhyloNet package (with permission from the authors).


Bug Reports:
-----------
contact ``astral-users@googlegroups.com``
