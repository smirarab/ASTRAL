-------------------------
DESCRIPTION:
---------------
ASTRAL is a Java program for estimating a species tree given a set of unrooted gene trees. ASTRAL is statistically consistent under multi-species coalescent model (and thus is useful for handling ILS). It finds the tree that maximizes the number of induced quartet trees in the set of gene trees that are shared by the species tree. The algorithm has an exact version that can run for small datasets (less than 18 taxa) and a more useful version (default) that can handle large datasets (tested for up to 1000 taxa and 1000 genes).

The original ASTRAL algorithm is described in:

Mirarab, Siavash, Rezwana Reaz, Md. Shamsuzzoha Bayzid, Theo Zimmermann, M Shel Swenson, and Tandy Warnow. “ASTRAL: Genome-Scale Coalescent-Based Species Tree.” Bioinformatics 30, no. 17 (2014): i541–i548. [doi:10.1093/bioinformatics/btu462](http://bioinformatics.oxfordjournals.org/cgi/content/full/btu462?ijkey=U8qRWnpC7BM4Syn&keytype=ref).

The current ASTRAL code is what we call ASTRAL-II and is described in detail in an upcoming publication. 

Email: `astral-users@googlegroups.com` for questions


-------------------------
Tutorial Steps:
---------------

### Step 1: INSTALLATION:

There is no installation required to run ASTRAL. You simply need to download the [zip file](https://github.com/smirarab/ASTRAL/raw/master/Astral.4.7.8.zip) and extract the contents to a folder of your choice. Alternatively, you can clone the [github repository](https://github.com/smirarab/ASTRAL/) if you are familiar with git. If you clone the git repository, you can run `make.sh` to build the project, or simply use the jar file that is included with the repository. 

ASTRAL is a Java-based application, and should run in any environment (Windows, Linux, Mac, etc.) as long as Java is installed. Java 1.5 or later is required. We have tested ASTRAL only on Linux and MAC, but it should work on Windows too.

In the remaining of the tutorial, we will assume you have extracted the ASTRAL zip file into a directory with the path `~/astral-home/`. In the commands given below, substitute `~/astral-home/` with the directory you have chosen for ASTRAL. 


### Step 2: Running ASTRAL from command-line to see the help:

ASTRAL currently has no GUI. You need to run it through command-line. Open a terminal (in Windows, look for a program called `Command Prompt` and run that; In Linux you should know how to do this; In MAC, search for an application called `Terminal`). Once the terminal is opened up, go the location where you have downloaded the software (e.g. using `cd ~/astral-home/`), and issue the following command:

```
  Java -jar astral.4.7.8.jar
```

This will print the list of options available in ASTRAL. If no errors are printed, your ASTRAL installation is fine and you can proceed to the next steps. 

### Step 3: Running ASTRAL on a sample input dataset 

We will next run ASTRAL on an input dataset. From the ASTRAL directory, run:

```
Java -jar astral.4.7.8.jar -i test_data/song_mammals.424.gene.tre
```

The results will be outputted to the standard output. To save the results in an output file use the `-o` option:

```
Java -jar astral.4.7.8.jar -i test_data/song_mammals.424.gene.tre -o test_data/song_mammals.tre
```

Here, the main input is just a file that contains all the input gene trees in newick format. The input gene trees are treated as unrooted, whether or not they have a root. Note that the output of ASTRAL should also be treated as an unrooted tree. 

The test file that we are providing here is based on the [Song et. al.](http://www.pnas.org/content/109/37/14942.short) dataset of 37 mammalian species and 442 genes. We have removed 23 problematic genes (21 mislabeled genes and 2 genes we classified as outliers) and we have also re-estimated gene trees using RAxML on the alignments that authors of that paper kindly provided to us. 

Some points regarding the input file are worth mentioning upfront:

**Polytomies**: The input gene trees can have polytomies (unresolved branches) since [version 4.6.0](CHANGELOG.md). 

**Multiple Individuals**: ASTRAL can handling multiple individuals from the same species since [version 4.5.0](CHANGELOG.md). If `s1i1`,`s1i2`,`s1i3` all belong to a species called `S1`, and similarly, `s2i1`,`s2i2` all belong to a species called `S2`, a mapping file with one of the following formats needs to be provided to ASTRAL:

```
S1 3 s1i1 s1i2 s1i3
S2 2 s2i1 s2i2
```
or

```
S1:s1i1,s1i2,s1i3
S2:s2i1,s2i2
```

**However,** it should be emphasized that the current handling of multiple taxa is not well tested and needs more work. Please stay tuned for upcoming improvements that fix/improve this feature. 



### Step 4: Viewing results of ASTRAL:

The output of ASTRAL is a tree in newick format. These trees can be viewed in many existing tools. Here are few that are used by many people:

1. [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) is probably the most widely used tool. It produces nice looking figures and works under Linux, Windows and MAC. 
2. [Archaeopteryx](https://sites.google.com/site/cmzmasek/home/software/archaeopteryx) is very good for viewing large trees and also works under all three operating systems. 
3. [EvolView](http://www.evolgenius.info/evolview): online application. You don't even need to download. 

There are [many more](http://en.wikipedia.org/wiki/List_of_phylogenetic_tree_visualization_software).

For this tutorial, let's use the online viewer (EvolView) or any other tool you can manage to download and install. Using either of these applications open the `test_data/song_mammals.tre` file. We will explore various tree viewing options. Importantly, we will reroot the tree at the correct node, which is always necessary, since the rooting of the ASTRAL trees is arbitrary and meaningless.  


### Step 5: The exact version of ASTRAL

ASTRAL has an exact and a heuristic version. The heuristic version solves the optimization problem exactly subject to the constraint that all the bipartitions in the species tree should be present in at least one of the input gene trees, or in a set of few other bipartitions that (since [version 4.5.1](CHANGELOG.md)) ASTRAL automatically infers from the set of input gene trees as likely bipartitions in the specie tree. The constrained version of ASTRAL is the default. However, when you have small number of taxa (typically 17 taxa or less), the exact version of ASTRAL can also run in reasonable time (see the figure below). 

![Running time Figure](runningtime.png "ASTRAL running times for exact version as function of the number of taxa")

Since the mammalian dataset we have used so far has 37 taxa, the exact version cannot run on it. However, we have created a subset of this dataset that has all 9 primates, tree shrew, rat, rabbit, horse, and the sloth (a total of 14 taxa). We can run the exact version of ASTRAL on this reduced dataset. Run:

```
Java -jar astral.4.7.8.jar -i test_data/song_primates.424.gene.tre -o test_data/song_primates.424.exact.tre -x
```

Using the `-x` option results in running the exact version of the ASTRAL algorithm. This run should finish in about 30 seconds. Now, we will run ASTRAL on the same input using the default heuristic algorithm:

```
Java -jar astral.4.7.8.jar -i test_data/song_primates.424.gene.tre -o test_data/song_primates.424.default.tre
```
This time, ASTRAL finished in under a second. So, is there a difference between the output of the exact and the heuristic version? Open up the two trees in your tree viewer tool and compare them. You will notice they are identical. You could also compare the scores outputted by ASTRAL and notice that they are identical. 


### Step 6: Example where exact helps


The default primate dataset we used in the previous step had 424 genes and 14 taxa. Since we have a relatively large number of gene trees, we could reasonably expect the exact and heuristic versions to generate identical output. The key point here is that as the number of genes increases, the probability that each bipartition of the species tree appears in at least one input gene tree increases. Thus, with 424 genes all bipartitions from the species tree are in at least one input gene tree, and therefore the exact and the heuristic versions are identical. 

We tried hard to find a subset of genes in the biological primates dataset where the exact and the heuristic versions did not match. We couldn't! So we had to resort to simulations. We simulated a 14-taxon dataset with extreme levels of ILS (average 87% RF between gene trees and the species tree). Now, with this simulated dataset, if you take only 10 genes, something interesting happens.


Run

```
Java -jar astral.4.7.8.jar -i test_data/simulated_14taxon.gene.tre -o test_data/simulated_14taxon.default.tre
```

and then

```
Java -jar astral.4.7.8.jar -i test_data/simulated_14taxon.gene.tre -o test_data/simulated_14taxon.exact.tre -x
```

Now you see that the tree outputted by the exact version has a slightly higher score (4812=48.07% versus 4803=47.98%), and a slightly different topology compared to the heuristic version. Thus, in extreme cases (i.e. lots of ILS and/or gene tree estimation error and few available gene trees compared to the number of taxa), one could observe differences between exact and heuristic versions. Note that how many genes should be considered few depends on the number of taxa you have, and also how much missing data there is. 

The main point of our ASTRAL-II work is to make sure the heuristic and exact versions essentially work identically. We have tested the heuristic version and a condition similar to the exact condition in our ASTRAL-II paper and have observed no differences. Thus, we believe the ASTRAL-II heuristics do a very good job of searching the tree space. 


### Step 7: Providing ASTRAL with extra trees:
**Note:** Since ASTRAL-II (version 4.7.2), we do not see any particular need for adding extra trees, unless the number of input trees is extremely small. We provide the following information, but note that ASTRAL-II should be able to define a thorough search space without a need for extra trees. Nevertheless, if extra trees are available, adding them never hurts.  

What if we had few very discordant gene trees, but the number of taxa was sufficiently large that we couldn't run the exact version?
We always have a second option. Imagine that you are able to create a set of hypothetical trees using various methods. For example, maybe you have prior hypothesis of what the species tee should be. Or, maybe you have run concatenation and have potential species trees. Most realistically, maybe you have a collection of bootstrapped gene tress that can be used. ASTRAL allows you to provide these sets of alternative trees to expand the space of that ASTRAL considers. Thus, ASTRAL will solve the optimization problem subject to the constraint that each bipartition should come either from one of the input gene trees, or the ones it infers automatically, or these "extra" gene trees. The extra gene trees, however, do not contribute to the score. They just control the space bing searched. 

So, now run:

```
Java -jar astral.4.7.8.jar -i test_data/simulated_primates_5X.10.gene.tre -o test_data/simulated_primates_5X.10.species.tre -e test_data/simulated_primates_5X.10.bootstrap.gene.tre
```
Here, the `-e` option is used to input a set of extra trees that ASTRAL uses to expand its search space. The file provided simply has 200 bootstrap replicates for each of the these 10 simulated genes. 

### Step 8: Running on large datasets:
We will now run ASTRAL on a relatively large dataset. Run:

```
Java -jar astral.4.7.8.jar -i test_data/100-simulated-boot
```

The input file here is a simulated dataset with 100 sequences and 100 replicates of bootstrapped gene trees for 25 loci (thus 2,500 input trees). Note that ATRAL finishes on this dataset in a matter of seconds. 


### Step 9: Multi-locus Bootstrapping:
Astral can perform multi-locus bootstrapping ([Seo, 2008](http://www.ncbi.nlm.nih.gov/pubmed/18281270)). To be able to perform multi-locus bootstrapping, ASTRAL needs to have access to bootstrap replicates for each gene. To start multi-locus bootstrapping using ASTRAL, you need to provide the location of all gene tree bootstrap replicates. To run bootstrapping on our test input files, go to `test_data` directory, and decompress the file called `song_mammals.424genes.bs-trees.zip`. Now run

```
Java -jar ../astral.4.7.8.jar -i song_mammals.424.gene.tre -b bs-files
```

This will run 100 replicates of bootstrapping. The argument after `-i` (here `song_mammals.424.gene.tre`) contains all the maximum likelihood gene trees (just like the case where bootstrapping was not used). The `-b` option tells ASTRAL that bootstrapping needs to be performed. Following `-b` is the name of a file (here `bs-files`) that contains the location of gene tree bootstrap files, one line per gene. For example, the first line is `424genes/100/raxmlboot.gtrgamma/RAxML_bootstrap.allbs`, which tells ASTRAL that the gene tree bootstrap replicates of the first gene can be found in a file called `424genes/100/raxmlboot.gtrgamma/RAxML_bootstrap.allbs`.

By default, ASTRAL performs 100 bootstrap replicates, but the `-r` option can be used to perform any number of replicates. For example, 

```
Java -jar ../astral.4.7.8.jar -i song_mammals.424.gene.tre -b bs-files -r 150
```

will do 150 replicates. Note that your input gene tree bootstrap files need to have enough bootstrap replicates for the number of replicates requested using `-r`. For example, if you have `-r 150`, each file listed in `bs-files` should contain at least 150 bootstrap replicates.

**Gene resampling:** 
ASTRAL performs site-only resampling by default (see [Seo, 2008](http://www.ncbi.nlm.nih.gov/pubmed/18281270)). ASTRAL can also perform gene/site resampling, which can be activated with the `-g` option:

```
Java -jar ../astral.4.7.8.jar -i song_mammals.424.gene.tre -b bs-files -g -r 100
```

Note that when you perform gene/site resampling, you need more gene tree replicates than the number of multi-locus bootstrapping replicates you requested using `-r`. For example, if you have `-g -r 100`, you might need 150 replicates for some genes (and less than 100 replicates for other genes). This is because when genes are resampled, some genes will be sampled more often than others by chance.

Finally, since bootstrapping involves a random process, a seed number can be provided to ASTRAL to ensure reproducibility. The seed number can be set using the `-s` option (by default 692 is used). 

**Output:** 
ASTRAL first performs the bootstrapping, and for each replicate, it outputs the bootstrapped ASTRAL tree. So, if number of replicates is set to 100, it outputs 100 trees. Then, it outputs a greedy consensus of all the 100 bootstrapped trees (with support drawn on branches). Finally, it performs the main analysis (i.e. on trees provided using `-i` option) and draws branch support on this main tree using the bootstrap replicates. The tree outputted at the end therefore is the ASTRAL tree on main input trees, with support values drawn based on bootstrap replicates. Support values are shown as branch length (i.e. after a colon sign) and
are percentages. 

Step 10: Scoring a species tree
--------------
You can use the `-q` option in ASTRAL to score an existing species tree. The ASTRAL score is the fraction of the induced quartet trees in the input set that are in the species tree. So, a score of `0.9` would mean that 90% of the quartet trees induced by your gene trees are in your species tree. 

To score a tree using ASTRAL, run:

```
Java -jar astral.4.7.8.jar -q test_data/simulated_14taxon.default.tre -i test_data/simulated_14taxon.gene.tre
```

This will score the species tree given in `test_data/simulated_14taxon.default.tre` compared to the gene trees given in `test_data/simulated_14taxon.gene.tre`. It will output:

```
Quartet score is: 4803
Normalized quartet score is: 0.4798201798201798
```

This means 4803 induced quartet trees from the gene trees are in the species tree, and these 4803 quartets are 47.98% of all the quartet trees that could be found in the species tree. As mentioned before, this dataset is one with a *very high* ILS level. 

Miscellaneous :
---------------

### Automatic addition of bipartitions to X.
The search space explored by ASTRAL is limited to bipartitions in a given set X. This means
that the accuracy of ASTRAL can be impacted by what is in the set X. 
ASTRAL-I (the algorithm described in our 2014 paper) sets the set X to all bipartitions from the gene trees. This is a very good start and sufficient in many cases. However, in ASTRAL-II, we have added few extra mechanisms to add extra bipartitions to set X. In our upcoming paper on ASTRAL-II, we show that for the most challenging datasets, these additions dramatically increase the accuracy of ASTRAL, to the point that the search space seems sufficiently large in all the analyses we ran.

Some of the mechanisms used in ASTRAL-II are the following:

1. **Completing gene trees:** If the input tree has missing taxa, we use a built-in mechanisms to first complete those gene trees and then add bipartitions to X from the completed gene trees (note that we use the original incomplete gene tree for score calculations). Our current implemented heuristic for doing this is based on calculating the similarity between all pairs of taxa (measured as how often they appear together on the same side of a quartet tree in induced gene trees) and using the four point condition to figure out where the taxon belongs in the gene tree. 

- **Greedy-based addition**: By default, we use an approach based on greedy consensus to find out additional bipartitions and add them to the set X, regardless of whether input gene trees are complete or not. In short, we first find the greedy consensus of the input gene trees at various thresholds. For each polytomy in these consensus trees, we find resolutions by randomly sampling taxa from its pending branches and calculating greedy consensus for gene trees restricted to these random samples. We do this many times for each polytomy and add all the induced new bipartitions. We also calculate a UPGMA resolution of the polytomy using our similarity matrix and add the resulting bipartitions to set X. We also use the similarity matrix to add all possible caterpillar resolutions of the polytomies according to the similarity matrix. 

- **Similarity-based addition**: By default, we also use our similarity matrix (also used for completing trees) to add some more bipartitions using a UPGMA tree. 
 
- **Resolving unresolved gene trees**: When input gene trees are unresolved, we use a strategy similar to our Greedy-based addition strategy for inferring and adding extra bipartitions that effectively resolve the polytomies. Once again, the score calculation is based on the unresolved tree, and these additions only affect the set X. 

- **Auto-expansion**: If a cluster of taxa in the dynamic programming phase cannot be further decomposed into two parts, we add all possible resolutions of it, if it is small enough (less than 6 taxa), or some new resolutions for it if it is larger.  


- **Similarity-based extensive addition**: this feature is by default off, because it did not help in our extensive simulations and results in a very large number of additions. You can turn it on if you want. 

Mechanisms 1, 4, and 5 are always in place. Mechanisms 2 and 3 are by default activated, but you can turn them off to increase speed (not recommended). To turn off these two features, use `-p 0`. To enable mechanism 6, use `-p 2`. 

We do not recommend `-p 0` option, because greedy-based additions could be quite important, and they typically have a relatively small impact on the running time. 


### Memory:
For big datasets (say more than 100 taxon) increasing the memory available to Java might be necessary. Note that you should never give Java more memory than what you have available on your machine. So, for example, if you have 4GB of free memory, you can invoke ASTRAL using the following command to make 3GB available to Java:

```
Java -Xmx3000M -jar astral.4.7.8.jar -i in.tree
```

### Other options:

* `-m [a number]`: removes genes with less that the specified number of leaves in them. Thus, this is useful for requiring a certain level of taxon occupancy. 
* `-k completed`: To build the set X (and *not* to score the species tree), ASTRAL internally completes the gene trees. To see these completed gene trees, run this option. This option is usable only when you also have `-o`.
* `-k bootstraps` and `-k bootstraps_norun`: these output the bootstrap replicate inputs to ASTRAL. These are useful if you want to run ASTRAL separately on each bootstrap replicate on a cluster. 
* `-k searchspace_norun`: outputs the search space (constraint set X) and exits. 

### Acknowledgment

ASTRAL code uses bytecode and some reverse engineered code from PhyloNet package (with permission from the authors).


### Bug Reports:
contact: ``astral-users@googlegroups.com``
