
ASTRAL is a java program for estimating a species tree given a set of unrooted gene trees. ASTRAL is statistically consistent under multi-species coalescent model (and thus is useful for handling ILS). The optimization problem solved by ASTRAL seeks to find the tree that maximizes the number of induced quartet trees in gene trees that are shared by the species tree. The optimization problem is solved exactly for a constrained version of the problem that restricts the search space. An exact solution to the unconstrained version is also implemented and can run on small datasets (less than 18 taxa). The current repository (master branch) includes the ASTRAL-III algorithm.  

Read the [README](README.md) file in addition to this tutorial. 

Email: `astral-users@googlegroups.com` for questions. Please subscribe to the mailing list for infrequent updates. 



* [Installation](#installation)
* [Running ASTRAL](#running-astral)
  * [ASTRAL Help](#astral-help)
  * [Running on the sample mammalian dataset](#running-on-the-sample-mammalian-dataset)
  * [Running on larger datasets:](#running-on-larger-datasets)
  * [Running with unresolved gene trees](#running-with-unresolved-gene-trees)
  * [Running on a multi\-individual datasets](#running-on-a-multi-individual-datasets)
* [Interpreting output](#interpreting-output)
  * [Viewing results of ASTRAL:](#viewing-results-of-astral)
  * [Branch length and support](#branch-length-and-support)
  * [The ASTRAL Log information](#the-astral-log-information)
* [Scoring existing trees](#scoring-existing-trees)
* [Extensive branch annotations](#extensive-branch-annotations)
  * [Newick annotations](#newick-annotations)
  * [File export of branch annotations](#file-export-of-branch-annotations)
  * [Prior hyper\-parameter](#prior-hyper-parameter)
* [Multi-locus Bootstrapping:](#multi-locus-bootstrapping)
* [The Search space of ASTRAL](#the-search-space-of-astral)
  * [Exact version](#exact-version)
    * [Example where exact helps](#example-where-exact-helps)
  * [Providing ASTRAL with extra trees](#providing-astral-with-extra-trees)
* [Miscellaneous](#miscellaneous)
  * [Memory](#memory)
  * [Other options](#other-options)
  * [Acknowledgment](#acknowledgment)
  * [Bug Reports](#bug-reports)


Installation
---------------

Refer to the [README](README.md#installation)

Running ASTRAL
---------------

ASTRAL currently has no GUI. You need to run it through command-line. 

* Open a terminal (on Windows, look for a program called `Command Prompt` and run that; on Linux you should know how to do this; on MAC, search for an application called `Terminal`).
* Once the terminal opens, go the location where you have downloaded the software (e.g. using `cd ~/astral-home/`),

###  ASTRAL Help

 To see the help, issue the following command:

```
  java -jar __astral.jar__
```

This will print the list of options available in ASTRAL. If no errors are printed, your ASTRAL installation is fine and you can proceed to the next sections. 

###  Running on the sample mammalian dataset 

We will next run ASTRAL on an input dataset. From the ASTRAL directory, run:

```
java -jar __astral.jar__ -i test_data/song_mammals.424.gene.tre
```

The results will be outputted to the standard output. To save the results in an output file use the `-o` option:

```
java -jar __astral.jar__ -i test_data/song_mammals.424.gene.tre -o test_data/song_mammals.tre
```

Here, the main input is just a file that contains all the input gene trees in Newick format. The input gene trees are treated as unrooted, whether or not they have a root. Note that the **output of ASTRAL should also be treated as an unrooted tree**. 

The test file that we are providing here is based on the [Song et. al.](http://www.pnas.org/content/109/37/14942.short) dataset of 37 mammalian species and 442 genes. We have removed 23 problematic genes (21 mislabeled genes and 2 genes we classified as outliers) and we have also re-estimated gene trees using RAxML on the alignments that the authors of that paper kindly provided to us. 

The input gene trees can have polytomies (unresolved branches) since [version 4.6.0](CHANGELOG.md). 

### Running on larger datasets:
We will now run ASTRAL on a larger dataset. Run:

```
java -jar __astral.jar__ -i test_data/100-simulated-boot
```

The input file here is a simulated dataset with 100 sequences and 100 replicates of bootstrapped gene trees for 25 loci (thus 2,500 input trees). Note that ASTRAL finishes on this dataset in a matter of seconds. 

A larger real dataset from the [1kp](http://www.pnas.org/content/early/2014/10/28/1323926111) dataset is also included. This dataset includes
424 genes from 103 species. Run:

```
java -jar __astral.jar__ -i test_data/1KP-genetrees.tre -o test_data/1kp.tre 2> test_data/1kp.log
```

This takes about a minute to run on a laptop. On this dataset, notice in the ASTRAL log information that it originally starts with 11043 clusters in its search space, and using heuristics implemented in ASTRAL-II, it increases the search space slightly to 11085 clusters. For more challenging datasets (i.e., more discordance or fewer genes) this number might increase a lot. 

### Running with unresolved gene trees

In our [ASTRAL-III paper](https://doi.org/10.1007/978-3-319-67979-2_4) we showed that contracting very low support branches (e.g., below 10% bootstrap support) from gene trees can improve accuracy somewhat. 
Thus, we recommend removing very low support branches. 

To contract low support branches, you can use many tools,  including the [newick utilities](htpp://cegg.unige.ch/newick_utils). If you have newick utilities installed, you can use

```
nw_ed  1KP-genetrees.tre 'i & b<=10' o > 1KP-genetrees-BS10.tre
```

To create a file `1KP-genetrees-BS10.tre` that includes the 1KP dataset with branches of 10% support or lower contracted. If you don't have newick utilities, don't worry. The contracted file is part of the ASTRAL distribution. 

```
java -jar __astral.jar__ -i test_data/1KP-genetrees-BS10.tre -o test_data/1kp-BS10.tre 2> test_data/1kp-bs10.log
```

Compare the species tree generated here with that generated with the fully resolved gene trees. You can confirm that the tree topology has not changed in this case, but the branch lengths and the branch support have all changed (and that they tend to both increase). By comparing the log files you can also see that after contracting low support branches, the normalized quartet score increases to 0.92321 (from 0.89467 with no contraction). This is expected as low support branches tend to increase not decrease discordance. 

### Running on a multi-individual datasets

When multiple individuals from the same species are available, to force the species to be monophyletic, a mapping file needs to be provided using the `-a` option. This mapping file should have one line per species, and each line needs to be in one of two formats:

```
species_name [number of individuals] individual_1 individual_2 ...

species_name:individual_1,individual_2,...
```
Some rules about the mapping file:

* When multiple individuals exist for the same species, your species names should be different from the individual names.
* You cannot have empty names (e.g., `,,`)
* The same individual name should not be mapped to multiple species
* Each individual name should appear in at least one gene name

For example, you can run

```
java -jar __astral.jar__ -i test_data/1KP-genetrees.tre -o test_data/1KP-speciestrees-multiind.tre -a test_data/namemap-1kp.txt 2> test_data/1kp-speciestree-multiind.log
```

Interpreting output
-----

### Viewing results of ASTRAL:

The output of ASTRAL is a tree in Newick format. These trees can be viewed in many existing tools. Here are some that are used by many people:

1. [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) is probably the most widely used tool. It produces nice looking figures and works under Linux, Windows, and MAC. 
2. [Archaeopteryx](https://sites.google.com/site/cmzmasek/home/software/archaeopteryx) is very good for viewing large trees and also works under all three operating systems. 
3. [EvolView](http://www.evolgenius.info/evolview): online application. You don't even need to download. 

There are [many more tools](http://en.wikipedia.org/wiki/List_of_phylogenetic_tree_visualization_software).

For this tutorial, let's use the online viewer (EvolView) or any other tool you can manage to download and install. Using either of these applications open the `test_data/song_mammals.tre` file. We will explore various tree viewing options. Importantly, we will reroot the tree at the correct node, which is always necessary, since the rooting of the ASTRAL trees is arbitrary and meaningless. 

There have been some reports that FigTree and some other tools sometimes have difficulty opening ASTRAL trees. This is likely because ASTRAL does not generate terminal branch lengths. In the case of FigTree, opening the three two or three times in a row often works (who knows why!). In other tools, you may have to add an arbitrary branch length to each branch. [This script](https://github.com/smirarab/global/blob/master/src/mirphyl/utils/add-bl.py) may be of help.  

### Branch length and support

ASTRAL measures branch length in coalescent units and also has a fast way of measuring support without a need for bootstrapping. 
The algorithms to compute branch lengths and support and the meaning of support outputted is further described in [this paper](http://mbe.oxfordjournals.org/content/early/2016/05/12/molbev.msw079.short?rss=1).
We will return to these in later sections. Some points have to be emphasized:

* ASTRAL only estimates branch lengths for internal branches and those terminal branches that correspond to species with more than one individuals sampled. 
* Branch lengths are in coalescent units and are a direct measure of the amount of discordance in the gene trees. As such, they are prone to underestimation because of statistical noise in gene tree estimation.   
* Branch support values measure the support for a quadripartition (the four clusters around a branch) and not the bipartition, as is commonly done. 



### The ASTRAL Log information

ASTRAL outputs lots of useful information to your screen ([stderr](https://en.wikipedia.org/wiki/Standard_streams), really). You can capture this information
by directing your stderr to a file. Capturing the log is highly recommended. Here is how you would capture stderr:

```
java -jar __astral.jar__ -i test_data/song_mammals.424.gene.tre -o test_data/song_mammals.tre 2> song_mammals.log
```

Here are some of the important information captured in the log:

* Number of taxa, and their names. Double check these to make sure they are correct.
* Number of genes. 
* Version of ASTRAL used in your analysis.
* The normalized quartet score (proportion of input gene tree quartet trees satisfied by the species tree). This is a number between zero and one; the higher this number, the *less* discordant your gene trees are. 
* The final optimization score is similar to the above number, but is not normalized (the number of gene tree quartets satisfied by the species tree).
* Running time.
* More advanced info: the size of the search space in terms of the number of clusters and number of tripartitions (i.e., elements weighted).  

Scoring existing trees
---

You can use the `-q` option in ASTRAL to score an existing species tree to produce its quartet score, compute its branch lengths, and compute its ASTRAL branch support values. The ASTRAL score is the fraction of the induced quartet trees in the input set that are in the species tree. So, a score of `0.9` would mean that 90% of the quartet trees induced by your gene trees are in your species tree. 

To score a tree using ASTRAL, run:

```
java -jar __astral.jar__ -q test_data/simulated_14taxon.default.tre -i test_data/simulated_14taxon.gene.tre -o test_data/simulated_scored.tre 2> test_data/simulated_scored.log
```

This will score the species tree given in `test_data/simulated_14taxon.default.tre` with respect to the gene trees given in `test_data/simulated_14taxon.gene.tre`. It will output the following in the log:

```
Quartet score is: 4803
Normalized quartet score is: 0.4798201798201798
```

This means 4803 induced quartet trees from the gene trees are in the species tree, and these 4803 quartets are 47.98% of all the quartet trees that could be found in the species tree. As mentioned before, this dataset is one with a *very high* ILS level. 

In addition to giving an overall score, when you score a tree, branch lengths
and branch support are also computed and outputted. In the next section, we will
introduce ways to output even more information per branch. 

* When scoring a tree, you probably want to capture the stderr using `2>name_of_a_file` redirection, as described before. 

* Just like the use of ASTRAL for inferring species trees, you can combine `-a` and `-q` to compute branch support values for multi-individual dataset.  

Extensive branch annotations
--

Whether you are inferring a species tree or scoring one using the `-q` option, you will always get estimates of the branch length and local posterior support for each branch. 
In addition to these default annotations for each branch, you can ask ASTRAL to output other per branch information. 

**Three Possible Rearrangements:**
Around each branch in an unrooted tree, there are four groups. If you think about a rooted tree, the four groups defined by a branch are the first child (`L`), the second child (`R`), the sister group (`S`), and everything else (`O`). With these four groups, if we keep all the groups intact, we can have three unrooted topologies: `RL|SO`, `RS|LO`, and `RO|LS`. The first topology is what the
current tree has, and we refer it to as the main topology. The two others
are alternative topologies, and we refer to `RS|LO` and `RO|LS` as the first and the second alternatives, respectively. ASTRAL can output branch annotations for the main tree as well as the two alternatives.

* When the output includes annotations for alternative topologies, we always first show the main topology, followed by `RS|LO` and `RO|LS`, in that order.
 
* Note: a common question I get is how do I know which child is `L` and which child is `R`. Unfortunately, the answer depends on your tree viewing software. There is no general rule. 
	* You can always look at the newick file directly, and the first child is `L` and the second is `R`; thus, in the output newick format, we always have
`(LEFT,RIGHT)`. 
	* But most people would find this cumbersome. Instead, you can look at options `-t 16` and `-t 32` [below](#file-export-of-branch-annotations) that output a `.csv` file.

**`-t` option**: 
To enable extra per branch information, you need to use the `-t` option, and various types of information that can be turned on. Most `-t` options just annotate the newick file. Some of the will also turn on the output of the `.csv` file mentioned earlier.


### Newick annotations


* *no annotations* (`-t 0`): This turns off calculation and reporting of posterior probabilities and branch lengths.  
* *Quartet support* (`-t 1`): The percentage of quartets in your gene trees that agree with a branch (normalized quartet support) give use a nice way of measuring the amount of gene tree conflict around a branch. Note that the local posterior probabilities are computed based on a transformation of normalized quartet scores (see [Figure 2 of this paper](https://doi.org/10.1093/molbev/msw079)).
* *Alternative quartet topologies* (`-t 8`): Outputs `q1`, `q2`, `q3`; these three values show quartet support (as defined above) for the main topology (`LR|SO`), first alternative (`RS|LO`) and second alternative (`RO|LS`), respectively.
* *Local posterior* (`-t 3`): is the default where we show [local posterior probability](https://doi.org/10.1093/molbev/msw079) for the main topology.
* *Alternative posteriors* (`-t 4`): The output includes three local posterior probabilities: one for the main topology, and one for each of the two alternatives (`RS|LO` and `RO|LS`, in that order). The posterior of the three topologies adds up to 1. This is because of our locality assumption, which basically asserts that we assume the four groups around the branch (`L`, `R`, `S`, and  `O`) are each correct and therefore, there are only three possible alternatives. 
* *Full annotation* (`-t 2`): When you use this option, for each branch you get a lot of different measurements:
   * `q1`,`q2`,`q3`: these three values show quartet support (as defined in the description of `-t 1`) for the main topology, the first alternative, and the second alternative, respectively. 
   * `f1`, `f2`, `f3`: these three values show the total number of quartet trees in all the gene trees that support the main topology, the first alternative, and the second alternative, respectively.
   * `pp1`, `pp2`, `pp3`: these three show the local posterior probabilities (as defined in the description of `-t 4`) for the main topology, the first alternative, and the second alternative, respectively.
   * `QC`: this shows the total number of quartets defined around each branch (this is what our paper calls `m`).
   * `EN`: this is the effective number of genes for the branch. 
      * If you don't have any missing data, this would be the number of genes in your dataset. 
      * When there are missing data, some gene trees might have nothing to say about a branch and these would not count towards EN. Thus, the effective number of genes might be smaller than the total number of genes. 
      * Things become more complex if your gene trees include polytomies. A fully resolved gene tree would contribute 1 to `f1+f2+f3` if it has at least one quartet around that branch and would contribute 0 if it has no quartet around it. If a gene tree has polytomies, it can contribute less than 1 to `f1+f2+f3`. Let `M` be the total number of fully resolved quartets of species around a branch of interest in the species tree (`M` already takes into account missing data). Let's say that out of `M` quartets, only `R` are resolved in the gene tree due to polytomies. The contribution of this gene tree to `f1+f2+f3` is `R/M<1` in this scenario. Thus, `EN=f1+f2+f3` can be a non-integer when we have multifurcations in the input gene trees. 
 
* *Polytomy test* (`-t 10`): runs an experimental test to see if a null hypothesis that the branch is a polytomy could be rejected. See [this paper: doi:10.3390/genes9030132](http://www.mdpi.com/2073-4425/9/3/132). 

Run:

```
java -jar __astral.jar__ -q test_data/1kp.tre -i test_data/1KP-genetrees.tre -t 2 -o test_data/1kp-scored-t2.tre
```
```
java -jar __astral.jar__ -q test_data/1kp.tre -i test_data/1KP-genetrees.tre -t 4 -o test_data/1kp-scored-t4.tre
```
```
java -jar __astral.jar__ -q test_data/1kp.tre -i test_data/1KP-genetrees.tre -t 8 -o test_data/1kp-scored-t8.tre
```
```
java -jar __astral.jar__ -q test_data/1kp.tre -i test_data/1KP-genetrees.tre -t 10 -o test_data/1kp-scored-t10.tre
```
read all the values given for a couple of branches and try to make sense of them. 


###  File export of branch annotations

Since it is often hard to know which branch is `L` and which branch is `R`, understanding branch annotations is a bit hard for most users. To help, we have added a feature for outputting some of the branch annotations into a `.csv` file. 

* Note that we strongly suggest using [DiscoVista]() to visualize the quartet frequencies. If you find DiscoVista hard to install and use, you can instead use these .csv files. 

To get the `.csv` outputs, you can use `-t 16` and `-t 32`.

* *.csv*  output has the following format. 
	* The output file is always called `freqQuad.csv` and is written to the same directory as the input file (sorry for the ugliness!)
	* The file is tab-delimited. 
	* 1st column: a dummy name for the node. Note that each three lines in a row have the same node number. Around each node, we have three possible unrooted toplogies (NNI rearrangements). We show stats for these three rearrangements. 
	* 2nd column: the topology name for which we are giving the scores. Here, `t1` is always the main topology (observed in your species tree) and `t2` and `t3` are the two alternatives. 
	* 3rd column: Gives the actual topology with the format: `{A}|{B}#{C}|{D}`. This means that the quartet topology being scored is putting groups A and B together on one side, and groups C and D on the other side. Please remember that quartets are unrooted trees. Each of the groups is a comma-separate list of species. 
	* 5th column: number of gene trees that match the the topology in this line. 
		* You will note that this number is not always an integer number. The reason is that in each gene tree, groups A, B, C, and D may not be together. Those gene trees still count as 1 unit, but they can contribute a fraction of that total of 1 to each of the tree topologies. So a gene tree may count as 0.7 for one topology, 0.2 for another, and 0.1 for the third.
	* 6th column: This is the total number of gene trees that had any useful information about this branch.
		* If you have no missing data, this should equal the total number of gene trees. 
		* If you have missing data, some genes may be missing one of groups A, B, C, or D entirely. Those genes will be agnostic about this branch. This column gives the number of genes that have at least one species from each of A, B, C, and D.  
	* 4th column is likely you are interested in and it depends on whether you used `-t 16` or `-t 32`. 
		* `-t 16`: This column is the local posterior for the topology given in this line. Note that the local posterior probability is different from normalized quartet score. See [Figure 2 of this paper](https://doi.org/10.1093/molbev/msw079).
		* `-t 32`: This column is simply 5th column divided by 6 column. Thus, it gives the **normalized quartet score** for that topology. Note that the three lines with the same node name (1st column) will add up to one in their 4th column. 

Hints: 

* since `freqQuad.csv` is tab-delimited, if you open it with Excel, you may have better luck if you rename it to `freqQuad.txt` or `freqQuad.tsv`

* Please ignore the R script `freqQuadVisualization.R`. That script is not to be used stand-alone. It can only be used within DiscoVista. 

#### Example:

Run:
```
java -jar __astral.jar__ -q test_data/1kp.tre -i test_data/1KP-genetrees.tre -t 16 -o test_data/1kp-scored-t16.tre
```

and examine the `freqQuad.csv` output file. 

Then run:
```
java -jar __astral.jar__ -q test_data/1kp.tre -i test_data/1KP-genetrees.tre -t 32 -o test_data/1kp-scored-t32.tre
```

and examine `freqQuad.csv`.

Note how all columns are the same except for the 4th column. 


### Prior hyper-parameter

Our calculations of the local posterior probabilities and branch lengths use a Yule prior model for the branch lengths of the species tree. The speciation rate (in coalescent units) of the Yule process (lambda) is by default set to 0.5, which results in a flat prior for the quartet frequencies in the `[1/3,1]` range. Using `-c` option one can adjust the hyper-parameter for the prior. For example, you might want to estimate lambda from the data after one run and plug the estimate prior in a subsequent run. We have not yet fully explored the impact of lambda on the posterior. For branch lengths, lambda acts as a pseudocount and can have a substantial impact on the estimated branch length for very long branches. More specifically, if there is no, or very little discordance around a branch, the MAP lengths of the branch (which is what we report) is almost fully determined by the prior. 

Run the following two commands and compare the lengths of the longest branches:

```
java -jar __astral.jar__ -q test_data/1kp.tre -i test_data/1KP-genetrees.tre -c 2 -o test_data/1kp-scored-c2.tre
```

```
java -jar __astral.jar__ -q test_data/1kp.tre -i test_data/1KP-genetrees.tre -c 0.001 -o test_data/1kp-scored-cs.tre
``` 

Note that setting lambda to 0 results in reporting ML estimates of the branch lengths instead of MAP. However, for branches with no discordance, we cannot compute a branch lengths. For these, we currently arbitrarily set ML to 10 coalescent units (we might change this in future versions).


Multi-locus Bootstrapping:
---
ASTRAL outputs a branch support value even without bootstrapping. Our [analyses](http://mbe.oxfordjournals.org/content/early/2016/05/12/molbev.msw079.short?rss=1) show that this form of support is more reliable than bootstrapping (under the conditions we explored). Nevertheless, you may want to run bootstrapping. ASTRAL can perform multi-locus bootstrapping ([Seo, 2008](http://www.ncbi.nlm.nih.gov/pubmed/18281270)). 

To be able to perform multi-locus bootstrapping, ASTRAL needs to have access to bootstrap replicate trees for each gene. 

#### Example

To start multi-locus bootstrapping using ASTRAL, you need to provide the location of all gene tree bootstrap replicates. To run bootstrapping on our test input files, 


* go to `test_data` directory, and 
* decompress the file called `song_mammals.424genes.bs-trees.zip`. 
* Now run:

```
java -jar ../__astral.jar__ -i song_mammals.424.gene.tre -b bs-files -o song_mammals.bootstrapped.astral.tre
```

This will run 100 replicates of bootstrapping in addition to one run of ASTRAL on the main trees. 

#### Input

* The argument after `-i` (here `song_mammals.424.gene.tre`) contains all the maximum likelihood gene trees (just like the case where bootstrapping was not used). 
* The `-b` option tells ASTRAL that bootstrapping needs to be performed. Following `-b` is the name of a file (here `bs-files`) that contains the **file path** of gene tree bootstrap files, one line per gene. Thus, the input file is **not** a file full of trees. It's a file full of paths of files full of trees.
    * For example, the first line is `424genes/100/raxmlboot.gtrgamma/RAxML_bootstrap.allbs`, which tells ASTRAL that the gene tree bootstrap replicates of the first gene can be found in a file called `424genes/100/raxmlboot.gtrgamma/RAxML_bootstrap.allbs`. In this case, `bs-files` has 424 lines (the number of genes) and each file named in `bs-files` has 100 trees (the number of bootstrap replicates). 


#### Output

The output file (here, `song_mammals.bootstrapped.astral.tre`) includes:

1. 100 bootstrapped replicate trees; each tree is the result of running ASTRAL on a set of bootstrap gene trees (one per gene).
2. A greedy consensus of the 100 bootstrapped replicate trees; this tree has support values drawn on branches based on the bootstrap replicate trees. Support values show the percentage of bootstrap replicates that contain a branch.
3. The “main” ASTRAL tree; this is the results of running ASTRAL on the `best_ml` input gene trees. This main tree also includes support values, which are again drawn based on the 100 bootstrap replicate trees.

**Note**: Support values are shown as percentages, as opposed to local posterior probabilities that are shown as a number between 0 and 1. When you perform bootstrapping, local posterior probabilities are not computed. If you want those as well, use the output of bootstrapping as input to astral with `-q` to annotate branches with posterior probabilities (see [branch annotations](#extensive-branch-annotations)).

As ASTRAL performs bootstrapping, it continually outputs the bootstrapped ASTRAL tree for each replicate. So, if the number of replicates is set to 100, it first outputs 100 trees. Then, it outputs a greedy consensus of all the 100 bootstrapped trees (with support drawn on branches). Finally, it performs the main analysis (i.e., on trees provided using `-i` option) and draws branch support on this main tree using the bootstrap replicates. Therefore, in this example, the output file will include 102 trees. 

**What to use:** The most important tree is the tree outputted at the end; this is the ASTRAL tree on main input trees, with support values drawn based on bootstrap replicates. Our [analyses](https://academic.oup.com/sysbio/article-lookup/doi/10.1093/sysbio/syu063) have shown this tree to be better than the consensus tree in most cases.  


#### Number of replicates 

By default, ASTRAL performs 100 bootstrap replicates, but the `-r` option can be used to perform any number of replicates.
**Note:**  when no `-r` is given, ASTRAL only performs 100 replicates regardless of the number of replicates in your bootstrapped gene trees.
For example, 

```
java -jar ../__astral.jar__ -i song_mammals.424.gene.tre -b bs-files -r 150 -o song_mammals.bootstrapped.150.astral.tre
```

will do 150 replicates. Note that your input gene tree bootstrap files need to have enough bootstrap replicates for the number of replicates requested using `-r`. For example, if you have `-r 150`, each file listed in `bs-files` should contain at least 150 bootstrap replicates.



#### Gene+site resampling:
ASTRAL performs site-only resampling by default (see [Seo, 2008](http://www.ncbi.nlm.nih.gov/pubmed/18281270)). ASTRAL can also perform gene+site resampling, which can be activated with the `-g` option:

```
java -jar ../__astral.jar__ -i song_mammals.424.gene.tre -b bs-files -g -r 100 -o song_mammals.bootstrapped.gs.astral.tre
```

Note that when you perform gene/site resampling, you need more gene tree replicates than the number of multi-locus bootstrapping replicates you requested using `-r`. For example, if you have `-g -r 100`, you might need 150 replicates for some genes (and less than 100 replicates for other genes). This is because when genes are resampled, some genes will be sampled more often than others by chance.

#### Gene-only resampling:
ASTRAL can also perform gene-only bootstrapping using the `--gene-only` option. This form of bootstrapping requires only one input file, which is given using `-i`. Thus, for this, you don't need to use `-b`. The following performs bootstrapping by resampling genes in the input file:

```
java -jar ../__astral.jar__ -i song_mammals.424.gene.tre --gene-only -o song_mammals.bootstrapped.go.astral.tre
```


Finally, since bootstrapping involves a random process, a seed number can be provided to ASTRAL to ensure reproducibility. The seed number can be set using the `-s` option (by default 692 is used). 


The Search space of ASTRAL
---

### Exact version
ASTRAL has an exact and a heuristic version. The heuristic version solves the optimization problem exactly subject to the constraint that all the bipartitions in the species tree should be present in at least one of the input gene trees, or in a set of some other bipartitions that (since [version 4.5.1](CHANGELOG.md)) ASTRAL automatically infers from the set of input gene trees as likely bipartitions in the species tree. The constrained version of ASTRAL is the default. However, when you have a small number of taxa (e.g., 18 taxa or less), the exact version of ASTRAL can also run in reasonable time.

Since the mammalian dataset we have used so far has 37 taxa, the exact version cannot run on it. However, we have created a subset of this dataset that has all 9 primates, tree shrew, rat, rabbit, horse, and the sloth (a total of 14 taxa). We can run the exact version of ASTRAL on this reduced dataset. Run:

```
java -jar __astral.jar__ -i test_data/song_primates.424.gene.tre -o test_data/song_primates.424.exact.tre -x
```

Using the `-x` option results in running the exact version of the ASTRAL algorithm. This run should finish in about 30 seconds. Now, we will run ASTRAL on the same input using the default heuristic algorithm:

```
java -jar __astral.jar__ -i test_data/song_primates.424.gene.tre -o test_data/song_primates.424.default.tre
```
This time, ASTRAL finished in under a second. So, is there a difference between the output of the exact and the heuristic version? Open up the two trees in your tree viewer tool and compare them. You will notice they are identical. You could also compare the scores outputted by ASTRAL and notice that they are identical. 


####  Example where exact helps


The default primate dataset we used in the previous step had 424 genes and 14 taxa. Since we have a relatively large number of gene trees, we could reasonably expect the exact and heuristic versions to generate identical output. The key point here is that as the number of genes increases, the probability that each bipartition of the species tree appears in at least one input gene tree increases. Thus, with 424 genes all bipartitions from the species tree are in at least one input gene tree, and therefore, the exact and the heuristic versions are identical. 

We tried hard to find a subset of genes in the biological primates dataset where the exact and the heuristic versions did not match. We couldn't! So we had to resort to simulations. We simulated a 14-taxon dataset with extreme levels of ILS (average 87% RF between gene trees and the species tree). Now, with this simulated dataset, if you take only 10 genes, something interesting happens.


Run

```
java -jar __astral.jar__ -i test_data/simulated_14taxon.gene.tre -o test_data/simulated_14taxon.default.tre
```

and then

```
java -jar __astral.jar__ -i test_data/simulated_14taxon.gene.tre -o test_data/simulated_14taxon.exact.tre -x
```

Now you see that the tree outputted by the exact version has a slightly higher score (4812=48.07% versus 4803=47.98%), and a slightly different topology compared to the heuristic version. Thus, in extreme cases (i.e., lots of ILS and/or gene tree estimation error and few available gene trees compared to the number of taxa), one could observe differences between the exact and heuristic versions. Note that how many genes should be considered few depends on the number of taxa you have, and also how much missing data there is. 

The main point of our ASTRAL-II work is to make  the heuristic version as close to the exact version as possible. We have tested the heuristic version and a condition similar to the exact condition in our ASTRAL-II paper and have observed no real differences. Thus, we believe the ASTRAL-II heuristics do a very good job of searching the tree space. But when the exact version can be run, there is no reason not to. 


### Providing ASTRAL with extra trees


We always have another option for increasing the search space. Imagine that you are able to create a set of hypothetical trees using various methods. For example, maybe you have a prior hypothesis of what the species tree could be. Or, maybe you have run concatenation and have potential species trees. Most realistically, maybe you have a collection of bootstrapped gene trees that can be used. ASTRAL allows you to provide these sets of alternative trees to expand the space of that ASTRAL considers. Thus, ASTRAL will solve the optimization problem subject to the constraint that each bipartition should come either from one of the input gene trees, or the ones it infers automatically, or these "extra" gene trees. The extra gene trees, however, do not contribute to the calculation of the score, which is always computed against input gene trees. They just add to the space being searched, and thus may be beneficial for ASTRAL. 

Since ASTRAL-II (version 4.7.2), we have reduced the need for adding extra trees, unless the number of input trees is extremely small. Nevertheless, if extra trees are available, adding them never hurts accuracy and typically has minimal
impact on running time.  

To expand the search space, you can run:

```
java -jar __astral.jar__ -i test_data/simulated_primates_5X.10.gene.tre -o test_data/simulated_primates_5X.10.species.tre -e test_data/simulated_primates_5X.10.bootstrap.gene.tre
```
Here, the `-e` option is used to input a set of extra trees that ASTRAL uses to expand its search space. The file provided simply has 200 bootstrap replicates for each of the these 10 simulated genes.
A similar option `-f` can be used when input trees have species labels instead of gene labels (only consequential when for multi-individual datasets).



Miscellaneous 
---------------

### Memory
For big datasets (say more than 500 taxa) increasing the memory available to java might be necessary. Note that you should never give java more memory than what you have available on your machine. So, for example, if you have 4GB of free memory, you can invoke ASTRAL using the following command to make 3GB available to java:

```
java -Xmx3000M -jar __astral.jar__ -i in.tree
```

### Memory issues due to taxon names

Some times, you have different names for the same species in your input gene trees. ASTRAL has no way of knowing they are the same species, and will think they are different species. It will hence use a lot of memory and will fail. 

The obvious solution is to use the `-a` option described earlier to map several names in the gene trees to the same name in the species tree output. This way, ASTRAL at least knows about the mapping between gene tree names and species names. This works to some extent, but has a drawback. If you have thousands of names in the gene trees, ASTRAL still will take up a lot of memory. But we have made a workaround, which is given using `-R` option. 

Assume you have an input file `gt.tre` and you have created a mapping file called `mapping.txt` (same format as multi-individual mapping file) that maps the gene tree names to species tree names. Then, you can run:

~~~bash
java -jar __astral.jar__ -i gt.tre  -a mapping.txt   -o gt-renamed.tre -R > newmapping.txt
~~~

This will generate the new gene tree file called `gt-renamed.tre` and a new mapping file called `newmapping.txt ` (in case you do have multiple individuals *also*). You can now use these new files as input to ASTRAL. 

~~~bash 
java -jar __astral.jar__ -i gt-renamed.tre  -a newmapping.txt   -o astral.tre 
~~~

If you do not have multiple individuals, check `newmapping.txt`. It should show only one gene name per species name (second column should be always 1). In that case, you can ommit `-a`.

### Other options

* `-m [a number]`: removes genes with less that the specified number of leaves in them. Thus, this is useful for requiring a certain level of taxon occupancy. 
* `-k completed`: To build the set X (and *not* to score the species tree), ASTRAL internally completes the gene trees. To see these completed gene trees, run this option. This option is usable only when you also have `-o`.
* `-k bootstraps` and `-k bootstraps_norun`: these options output the bootstrap replicate inputs to ASTRAL. These are useful if you want to run ASTRAL separately on each bootstrap replicate on a cluster. 
* `-k searchspace_norun`: outputs the search space (constraint set X) and exits. 
* `--polylimit`: when ASTRAL adds new bipartitions to its search space, it partially does so based on a quadratic number of resolutions per each polytomy on a greedy consensus of gene trees. This could be slow, and therefore, we limit it to small polytomies. If you like to increase or decrease the search space, adjusting this option helps. 
* `--samplingrounds`: For multi-individual datasets, this option controls how many rounds of individual sampling is used in building the constraint set. Adjust to reduce/increase the search space for multi-individual datasets



### Acknowledgment

ASTRAL code uses bytecode and some reverse engineered code from PhyloNet package (with permission from the authors).


### Bug Reports
contact: ``astral-users@googlegroups.com``
