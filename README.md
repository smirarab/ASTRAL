DESCRIPTION:
-----------
ASTRAL is a tool for estimating an unrooted species tree given a set of unrooted gene trees.
ASTRAL is statistically consistent under the multi-species coalescent model (and thus is useful for handling incomplete lineage sorting, i.e., ILS).
ASTRAL finds the species tree that has the maximum number of shared induced quartet trees with the set of gene trees, subject to the constraint that the set of bipartitions in the species tree comes from a predefined set of bipartitions. This predefined set is empirically decided by ASTRAL (but see tutorial on how to expand it). 



The current code corresponds to **ASTRAL-III** (see below for the publication).
The algorithm was designed by Tandy Warnow and Siavash Mirarab originally. ASTRAL-III incorporates many ideas by Chao Zhang and Maryam Rabiee.
[Code developers](https://github.com/smirarab/ASTRAL/graphs/contributors) are mainly Siavash Mirarab, Chao Zhang, Maryam Rabiee, and Erfan Sayyari.

Email: `astral-users@googlegroups.com` for questions.



##### Publications:

- The original algorithm (ASTRAL-I) is described in:
    - Mirarab, Siavash, Rezwana Reaz, Md. Shamsuzzoha Bayzid, Theo Zimmermann, M Shel Swenson, and Tandy Warnow. “ASTRAL: Genome-Scale Coalescent-Based Species Tree.” Bioinformatics (ECCB special issue) 30 (17): i541–i548. 2014. [doi:10.1093/bioinformatics/btu462](doi.org/10.1093/bioinformatics/btu462).
- All the versions between 4.7.4  and 5.1.0 corresponds to ASTRAL-II, described in:
    * Mirarab, Siavash, and Tandy Warnow. “ASTRAL-II: Coalescent-Based Species Tree Estimation with Many Hundreds of Taxa and Thousands of Genes.”. Bioinformatics (ISMB special issue) 31 (12): i44–i52. 2015. [doi:10.1093/bioinformatics/btv234](http://bioinformatics.oxfordjournals.org/content/31/12/i44)
- Since version 5.1.1, the code corresponds to **ASTRAL-III**, described in:
    * Zhang, Chao, Maryam Rabiee, Erfan Sayyari, and Siavash Mirarab. 2018. “ASTRAL-III: Polynomial Time Species Tree Reconstruction from Partially Resolved Gene Trees.” BMC Bioinformatics 19 (S6): 153. [doi:10.1186/s12859-018-2129-y](https://doi.org/10.1186/s12859-018-2129-y).
- For **multi-individual** datasets, the relevant paper to cite is:
	* Rabiee, Maryam, Erfan Sayyari, and Siavash Mirarab. 2018. “Multi-Allele Species Reconstruction Using ASTRAL.” Molecular Phylogenetics and Evolution 130 (2018): 439489. [doi:10.1016/j.ympev.2018.10.033](https://www.sciencedirect.com/science/article/pii/S1055790317308424?dgcid=author).
- Since version 4.10.0, ASTRAL can also compute branch length (in coalescent units) and a measure of support called “local posterior probability”, described here:
    * Sayyari, Erfan, and Siavash Mirarab. “Fast Coalescent-Based Computation of Local Branch Support from Quartet Frequencies.” Molecular Biology and Evolution 33 (7): 1654–68. 2016. [doi:10.1093/molbev/msw079](http://mbe.oxfordjournals.org/content/early/2016/05/12/molbev.msw079.short?rss=1)
- ASTRAL can also perform a test of polytomies in species trees (`-t 10` option):
    * Sayyari, Erfan, and Siavash Mirarab. “Testing for Polytomies in Phylogenetic Species Trees Using Quartet Frequencies.” Genes 9 (3): 132. 2018. [doi:10.3390/genes9030132](http://www.mdpi.com/2073-4425/9/3/132).


Documentations
-----------

- The rest of this README file
- **Our [tutorial](astral-tutorial.md)**.
- The chapter of Siavash Mirarab's dissertation that describes ASTRAL in detail is provided [here](thesis-astral.pdf).
- Publications shown above have scientific details
- A [developer guide](developer-guide.md).

INSTALLATION:
-----------
There is no installation required to run ASTRAL.
You simply need to download the [zip file](https://github.com/smirarab/ASTRAL/raw/master/Astral.5.6.3.zip)
and extract the contents to a folder of your choice. Alternatively, you can clone the [github repository](https://github.com/smirarab/ASTRAL/). You can run `make.sh` to build the project or simply use the jar file that is included with the repository.

ASTRAL is a java-based application, and should run in any environment (Windows, Linux, Mac, etc.) as long as java is installed. Java 1.5 or later is required. We have tested ASTRAL only on Linux and MAC.

To test your installation, go to the place where you put the uncompressed ASTRAL, and run:

```
java -jar astral.5.6.3.jar -i test_data/song_primates.424.gene.tre
```

This should quickly finish. There are also other sample input files under `test_data/` that can be used.

ASTRAL can be run from any directory. You just need to run `java -jar /path/to/astral/astral.5.6.3.jar`.
Also, you can move `astral.5.6.3.jar` to any location you like and run it from there, but note that you need
to move the `lib` directory as well.

EXECUTION:
-----------
ASTRAL currently has no GUI. You need to run it through the command-line. In a terminal, go the location where you have downloaded the software, and issue the following command:

```
  java -jar astral.5.6.3.jar
```

This will give you a list of options available in ASTRAL.

To find the species tree given a set of gene trees in a file called `in.tree`, use:

```
java -jar astral.5.6.3.jar -i in.tree
```

The results will be outputted to the standard output. To save the results in a file use the `-o` option (**Strongly recommended**):

```
java -jar astral.5.6.3.jar -i in.tree -o out.tre
```
To save the logs (**also recommended**), run:

```
java -jar astral.5.6.3.jar -i in.tree -o out.tre 2>out.log
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
* branch supports measured as [local posterior probabilities](http://mbe.oxfordjournals.org/content/early/2016/05/12/molbev.msw079.short?rss=1). 
* It can also annotate branches with other quantities, such as quartet support, as described in the [tutorial](astral-tutorial.md).

The ASTRAL tree leaves the branch length of terminal branches empty. Some tools for visualization and tree editing do not like this (e.g., ape). In FigTree, if you open the tree several times, it eventually opens up (at least on our machines). In ape, if you ask it to ignore branch lengths all together, it works. In general, if you tool does not like the lack of terminal branches, you can add a dummy branch length, [as in this script](https://github.com/smirarab/global/blob/master/src/mirphyl/utils/add-bl.py). 

### Other features (local posterior, bootstrapping):

Please refer to the [tutorial](astral-tutorial.md) for all other features, including bootstrapping, branch annotation, and local posterior probability.

### Memory:
For big datasets (say more than 5000 taxa), increasing the memory available to Java can result in speedups. Note that you should give Java only as much free memory as you have available on your machine. So, for example, if you have 3GB of free memory, you can invoke ASTRAL using the following command to make all the 3GB available to Java:

```
java -Xmx3000M -jar astral.5.6.3.jar -i in.tree
```

Acknowledgment
-----------
ASTRAL code uses bytecode and some reverse engineered code from PhyloNet package (with permission from the authors).


Bug Reports:
-----------
contact ``astral-users@googlegroups.com``
