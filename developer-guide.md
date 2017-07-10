This document will be used to document the code. Algorithms are given in the papers and will not be repeated here. 

The project comes with an .project file that can be imported to Eclipse. 
Make sure you select at least Java 1.6 for this to compile fine in your Eclipse. 

## Versioning

Please update the version number in the Commandline.java file everytime the code changed. 
After updating the version, you need to commit to the repository the new file.
To do that, 

- run `make.sh`
- remove the old zip file
- git add the new zip file
- commit to git

## Design

The code is designed such that various phylogeny reconstruction methods
based on dynamic programming could be implemented. 
This *generality* creates some complications, because lots of functions
have to be abstracted out to base classes. 
Currently, the general algorithm is instantiated twice, 
once for ASTRAL, and once for DynaDup. 
While ASTRAL is under constant development, testing, and use, 
DynaDup has not been tested or used in a long time (it may or may not
be functioning at this point). 

By Java convention, we use an `I` prefix for interfaces, 
and an `Abstract` prefix for abstract classes that need to
be instantiated for each individual algorithm (e.g., for ASTARL). 
Instances for DynaDup are prefixed with `DL` and those relevant
to ASTRAL are prefixed by `WQ`. 

Another important point to consider is that the code was originally
based on phylonet and inherits some of phylonet's structure to date.
Many of the classes are used in binary form and as is from phylonet. 
A few of them are overwritten in ASTRAL. 

## Code components

The code is in several packages but [phylonet.coalescent](main/phylonet/coalescent/)  is where the majority of the code resides. When we don't mention a package, we mean this place. 


### General

In ASTRAL, each clade or partition is represented as a BitSet. 

*  [phylonet.util.BitSet](main/phylonet/util/BitSet.java): This simply modifies Java's Bitset class for improved efficiency and usability.
*  [phylonet.tree.model.sti.STITreeCluster](main/phylonet/tree/model/sti/STITreeCluster.java): This is a cluster, which is simply a subset of the set of leaves. Each Cluster should be associated with a taxon identifier, and inside includes a bitset.
	* [phylonet.tree.model.sti.STITreeCluster.Vertex](main/phylonet/tree/model/sti/STITreeCluster.java): This is an inner class inside STITreeCluster, which simply constitutes a node in the dynamic programming. It has various extra fields (in addition to the implicit reference to its outter class, which is the cluster in the dynamic programming) for backtracking of the dynamic programming. 

Then, there are interfaces and abstract classes:

* [IClusterCollection](main/phylonet/coalescent/IClusterCollection.java), [AbstractClusterCollection](main/phylonet/coalescent/AbstractClusterCollection.java):
Keep the set of clusters (used as the set `X` in ASTRAL)
	* Always has a *top* cluster; all other clusters are a subset of the top cluster
	* Has the important capacity to generate a subset of itself limited only to those clusters that are a subset of a given cluster; thus, it can generate a subset of itself with a different top element
	* Can find all resolutions of the top element among pairs of clusters in it. Useful in the dynamic programming. 
* [AbstractDataCollection](main/phylonet/coalescent/AbstractDataCollection.java):
Mostly just has methods for setting up the set `X`. Thus, this sets up and maintains an instance of IClusterCollection 
* [AbstractInference](main/phylonet/coalescent/AbstractInference.java):
The hub for all types of operations that ASTRAL can perform (dynamic programming inference, scoring, building set X, etc.).
Actual computations are delegated to other classes. This just "manages" stuff. Has access to input options, to data collections, and to weight
calculators. Sets everything up and starts computations. 
* [AbstractComputeMinCostTask](main/phylonet/coalescent/AbstractComputeMinCostTask.java): This is where the **dynamic programming** algorithm is actually implemented. Uses Vertex and IClusterCollection, in addition to the WeightCalculator to perform the dynamic programming. 
* [AbstractWeightCalculator](main/phylonet/coalescent/AbstractWeightCalculator.java): Subclasses of this are where the weights in the dynamic programming are computed. The abstract class does little useful work. Important methods are abstract. 



### Relevant to ASTRAL

Classes specific to ASTRAL are:

* [WQClusterCollection](main/phylonet/coalescent/WQClusterCollection.java): Doesn't do anything. Work done in the abstract class. 
* [WQComputeMinCostTask](main/phylonet/coalescent/WQComputeMinCostTask.java): Doesn't do anything. Work done in the abstract class.
* [WQDataCollection](main/phylonet/coalescent/WQDataCollection.java):
Has important functionality related to building the search space, including all the ASTRAL-II heuristics for augmenting the set `X`, the heuristics
for completing incomplete gene trees, heuristics for handling multi-individual datasets, etc. 
* [WQInference](main/phylonet/coalescent/WQInference.java):
has implementation of functions for setting up the data structures, 
and also:
	* Has functions for scoring a given species tree
	* Has functions (scoreBrances) for branch annotations using BipartitionWeightCalculator. 
* [WQWeightCalculator](main/phylonet/coalescent/WQWeightCalculator.java): This class knows how to find the score of a tripartition from the gene trees. Has two algorithms, implemented in two nested inner classes:
	* WQWeightCalculator.SetWeightCalculator: this is the old algorithm, used in ASTRAL-I and is based on set intersections
	* WQWeightCalculator.TraversalWeightCalculator: this is the new algorithm based on tree traversal, used in ASTRAL-II
* [BipartitionWeightCalculator](main/phylonet/coalescent/BipartitionWeightCalculator.java): This class computes a weight of Bipartition versus a set of input tripartitions. This is relevant to our local posterior probability paper and computation of branch lengths
* [Posterior](main/phylonet/coalescent/Posterior.java): This computes posterior probabilities given quartet frequencies
* [SimilarityMatrix](main/phylonet/coalescent/SimilarityMatrix.java): This implements simple functionalities for a distance matrix. The distance matrix is used in completion of gene trees and in adding to set X. 
* [Tripartition](main/phylonet/coalescent/Tripartition.java): This class  represents a tripartition


### Relevant to input/output

Classes relevant to the command-line GUI, including handling of input and output and taxon names are:

* [GlobalMaps](main/phylonet/coalescent/GlobalMaps.java): Static class, keeping instances of singleton classes that should be accessed only through here:
	* taxon identifier
	* random number generator
	* taxon name map 	
* [CommandLine](main/phylonet/coalescent/CommandLine.java): Reads command line input using the JSAP library and starts the correct inference based on the options provided by the user
	* Handles bootstrapping options (to be refactored out)
	* Handles creation of some singletons, such as taxon identifier, name maps, etc.  
* [Options](main/phylonet/coalescent/Options.java): Saves user input options. Inference class has an instance, which is created by the CommandLine. 

### Relevant to species names, individual names, etc.

* [TaxonIdentifier](main/phylonet/coalescent/TaxonIdentifier.java):
 Maps taxon names to taxon IDs.  
* [TaxonNameMap](main/phylonet/coalescent/TaxonNameMap.java):
Maps the names in the gene trees to the names in the species tree 
* [SpeciesMapper](main/phylonet/coalescent/SpeciesMapper.java): Maps IDs between gene trees and the species tree
* [SingleIndividualSample.java](main/phylonet/coalescent/SingleIndividualSample.java):  This class keeps track of a single-individual sample
 of a multi-individual dataset. For each species, 
we will include only one of its individuals in any instance of this class. 
* [Solution](main/phylonet/coalescent/Solution.java): Used for building and keeping the output (not important; to be removed maybe)
* [phylonet.tree.io.NewickWriter](main/phylonet/tree/io/NewickWriter.java): simple modifications to phylonet's newick writer




