This document will be used to document the code. Algorithms are given in the papers and will not be repeated here. 

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

The code is in several packages. [phylonet.coalescent](main/phylonet/coalescent/) this is where the majority of the code resides. When we don't mention a package, we mean this place. 


### General

In ASTRAL, each clade or partition is represented as a BitSet. 

*  [phylonet.util.BitSet](main/phylonet/util/BitSet.java): This simply modifies Java's Bitset class for improved efficiency and usability.
*  [phylonet.tree.model.sti.STITreeCluster](main/phylonet/tree/model/sti/STITreeCluster.java): This is a cluster, which is simply a subset of the set of leaves. Each Cluster should be associated with a taxon identifier, and inside includes a bitset.
	* [phylonet.tree.model.sti.STITreeCluster.Vertex](main/phylonet/tree/model/sti/STITreeCluster.java): This is an inner class inside Vertex, which simply constitutes a node in the dynamic programming. It has various extra fields (in addition to the implicit reference to its out class, which is the cluster in the dynamic programming) for backtracking of the dynamic programming. 

The, there are interfaces and abstract classes:

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
* [AbstractComputeMinCostTask](main/phylonet/coalescent/AbstractComputeMinCostTask.java): This is where the **dynamic programming** algorithm is actually implemented. Uses Vertex and IClusterCollection, in addition to the WeightCalculator to perform one run of the dynamic programming. 
* [AbstractWeightCalculator](main/phylonet/coalescent/AbstractWeightCalculator.java): This is where the weights in the dynamic programming are computed. 

### Relevant to input/output

Classes relevant to the command-line GUI, including handling of input and output, are:

* [CommandLine](main/phylonet/coalescent/CommandLine.java): Reads command line input using the JSAP library and starts the execution
* [Options](main/phylonet/coalescent/Options.java): Saves user input options 
* [GlobalMaps](main/phylonet/coalescent/GlobalMaps.java): Static class, keeping instances of singleton classes that should be accessed only through here:
	* taxon identifier
	* random number generator
	* taxon name map 	
* [TaxonIdentifier](main/phylonet/coalescent/TaxonIdentifier.java):
 Maps taxon names to taxon IDs.  
* [TaxonNameMap](main/phylonet/coalescent/TaxonNameMap.java):
Maps the names in the gene trees to the names in the species tree 
* [SpeciesMapper](main/phylonet/coalescent/SpeciesMapper.java): Maps IDs between gene trees and the species tree
* [Solution](main/phylonet/coalescent/Solution.java): Used for keeping the output
* [phylonet.tree.io.NewickWriter](main/phylonet/tree/io/NewickWriter.java): simple modifications to phylonet's newick writer


### Relevant to ASTRAL

Classes specific to ASTRAL are:

* [WQClusterCollection](main/phylonet/coalescent/WQClusterCollection.java): Doesn't do anything. Work done in the abstract class. 
* [WQComputeMinCostTask](main/phylonet/coalescent/WQComputeMinCostTask.java): Doesn't do anything. Work done in the abstract class.
* [WQDataCollection](main/phylonet/coalescent/WQDataCollection.java)
* [WQInference](main/phylonet/coalescent/WQInference.java)
* [WQWeightCalculator](main/phylonet/coalescent/WQWeightCalculator.java)
* [BipartitionWeightCalculator](main/phylonet/coalescent/BipartitionWeightCalculator.java)
* [Posterior](main/phylonet/coalescent/Posterior.java)
* [SimilarityMatrix](main/phylonet/coalescent/SimilarityMatrix.java)
* [Tripartition](main/phylonet/coalescent/Tripartition.java)







