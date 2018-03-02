DESCRIPTION:
-----------
ASTRAL is a tool for estimating an unrooted species tree given a set of unrooted gene trees.
ASTRAL is statistically consistent under the multi-species coalescent model (and thus is useful for handling incomplete lineage sorting, i.e., ILS).
ASTRAL finds the species tree that has the maximum number of shared induced quartet trees with the set of gene trees, subject to the constraint that the set of bipartitions in the species tree comes from a predefined set of bipartitions. This predefined set is empirically decided by ASTRAL (but see tutorial on how to expand it). 



The current code corresponds to **ASTRAL-III** (see below for the publication).
The algorithm was designed by Tandy Warnow and Siavash Mirarab originally. ASTRAL-III incorporates many ideas by Chao Zhang and Maryam Rabiee.
[Code developers](https://github.com/smirarab/ASTRAL/graphs/contributors) are mainly Siavash Mirarab, Chao Zhang, Maryam Rabiee, and Erfan Sayyari.

Email: `astral-users@googlegroups.com` for questions.



#### Documentations

1. The rest of this README file
- **Our [tutorial](astral-tutorial.md)**.
- The chapter of Siavash Mirarab's dissertation that describes ASTRAL in detail is provided [here](thesis-astral.pdf).
- Publications below have scientific details
- A [developer guide](developer-guide.md).

##### Publications:

- The original algorithm (ASTRAL-I) is described in:
    - Mirarab, Siavash, Rezwana Reaz, Md. Shamsuzzoha Bayzid, Theo Zimmermann, M Shel Swenson, and Tandy Warnow. “ASTRAL: Genome-Scale Coalescent-Based Species Tree.” Bioinformatics (ECCB special issue) 30 (17): i541–i548. 2014. [doi:10.1093/bioinformatics/btu462](doi.org/10.1093/bioinformatics/btu462).
- All the versions between 4.7.4  and 5.1.0 corresponds to ASTRAL-II, described in:
    * Mirarab, Siavash, and Tandy Warnow. “ASTRAL-II: Coalescent-Based Species Tree Estimation with Many Hundreds of Taxa and Thousands of Genes.”. Bioinformatics (ISMB special issue) 31 (12): i44–i52. 2015. [doi:10.1093/bioinformatics/btv234](http://bioinformatics.oxfordjournals.org/content/31/12/i44)
- Since version 5.1.1, the code corresponds to ASTRAL-III, described in:
    * Zhang, Chao, Erfan Sayyari, and Siavash Mirarab. “ASTRAL-III: Increased Scalability and Impacts of Contracting Low Support Branches.” In Comparative Genomics: 15th International Workshop, RECOMB CG, 2017. [doi:10.1007/978-3-319-67979-2_4](https://doi.org/10.1007/978-3-319-67979-2_4). 
- Since version 4.10.0, ASTRAL can also compute branch length (in coalescent units) and a measure of support called “local posterior probability”, described here:
    * Sayyari, Erfan, and Siavash Mirarab. “Fast Coalescent-Based Computation of Local Branch Support from Quartet Frequencies.” Molecular Biology and Evolution 33 (7): 1654–68. 2016. [doi:10.1093/molbev/msw079](http://mbe.oxfordjournals.org/content/early/2016/05/12/molbev.msw079.short?rss=1)



INSTALLATION:
-----------
There is no installation required to run ASTRAL.
You simply need to download the [zip file](https://github.com/smirarab/ASTRAL/raw/master/Astral.5.11.1b.zip)
and extract the contents to a folder of your choice. Alternatively, you can clone the [github repository](https://github.com/smirarab/ASTRAL/). You can run `make.sh` to build the project or simply use the jar file that is included with the repository.

ASTRAL is a java-based application, and should run in any environment (Windows, Linux, Mac, etc.) as long as java is installed. Java 1.5 or later is required. We have tested ASTRAL only on Linux and MAC.

To test your installation, go to the place where you put the uncompressed ASTRAL, and run:

```
java -Djava.library.path=. -jar astral.5.11.1b.jar -i test_data/song_primates.424.gene.tre
```

This should quickly finish. There are also other sample input files under `test_data/` that can be used.

ASTRAL can be run from any directory. You just need to run `java -jar /path/to/astral/astral.5.11.1b.jar`.
Also, you can move `astral.5.11.1b.jar` to any location you like and run it from there, but note that you need
to move the `lib` directory as well.

EXECUTION:
-----------
ASTRAL currently has no GUI. You need to run it through the command-line. In a terminal, go the location where you have downloaded the software, and issue the following command:

```
  java -Djava.library.path=. -jar astral.5.11.1b.jar
```

This will give you a list of options available in ASTRAL.

To find the species tree given a set of gene trees in a file called `in.tree`, use:

```
java -Djava.library.path=. -jar astral.5.11.1b.jar -i in.tree
```

The results will be outputted to the standard output. To save the results in a file use the `-o` option (**Strongly recommended**):

```
java -Djava.library.path=. -jar astral.5.11.1b.jar -i in.tree -o out.tre
```
To save the logs (**also recommended**), run:

```
java -Djava.library.path=. -jar astral.5.11.1b.jar -i in.tree -o out.tre 2>out.log
```

###### Input: 
* The input gene trees are in the Newick format
* The input trees can have missing taxa, polytomies (unresolved branches), and also multiple individuals per species.
*  Taxon names cannot have quotation marks in their names (sorry!). This means you also cannot have weird characters like ? in the name (underscore is fine).
* When multiple individuals from the same species are available, you can ask ASTRAL to force them to be together in the species tree. To do this, a mapping file needs to be provided using the `-a` option. This mapping file should have one line per species, and each line needs to be in one of two formats:

```
species_name [number of individuals] individual_1 individual_2 ...

species_name:individual_1,individual_2,...
```
   Note that when multiple individuals exist for the same species, your species name should be different from the individual names.
   
###### Output: 
The output in is Newick format and gives: 

* the species tree topology, 
* branch lengths in coalescent units (only for internal branches or for terminal branches if that species has multiple individuals),
* branch supports measured as [local posterior probabilities](). 
* It can also annotate branches with other quantities, such as quartet support, as described in the [tutorial](astral-tutorial.md).



### Bootstrapping:

To perform 100 replicates of multi-locus bootstrapping ([Seo 2008](http://www.ncbi.nlm.nih.gov/pubmed/18281270)), use:

```
java -Djava.library.path=. -jar astral.5.11.1b.jar -i best_ml -b bs_paths -r 100
```

In this command, `bs_paths` is a file that gives the location (file path) of gene tree bootstrap files, one line per gene. See the [tutorial](astral-tutorial.md)
for more details.
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
For big datasets (say more than 200 taxa), increasing the memory available to Java can result in speedups. Note that you should give Java only as much free memory as you have available on your machine. So, for example, if you have 3GB of free memory, you can invoke ASTRAL using the following command to make all the 3GB available to Java:

```
java -Xmx3000M -Djava.library.path=. -jar astral.5.11.1b.jar -i in.tree
```

Acknowledgment
-----------
ASTRAL code uses bytecode and some reverse engineered code from PhyloNet package (with permission from the authors).


Bug Reports:
-----------
contact ``astral-users@googlegroups.com``
