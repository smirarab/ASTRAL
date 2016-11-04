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

The code is in several packages:

* `phylonet.coalescent`: this is where the majority of the code resides. 
* `phylonet.tree.io`, `phylonet.tree.model.sti`, `phylonet.util`: these each include a single class, overwriting parts of java or phylonet for our purposes. 

We focus on `phylonet.coalescent`. 

### General
First, some interfaces and abstract classes are:

* [IClusterCollection.java](), [AbstractClusterCollection.java](main/phylonet/coalescent/IClusterCollection.java](), [AbstractClusterCollection.java):
	* Keep the set of clusters (used as X in ASTRAL)
* [AbstractInference.java](main/phylonet/coalescent/AbstractInference.java):
	* The hub for all types of operations (dynamic programming inference, scoring, building set X, etc. ) 
* [AbstractDataCollection.java](main/phylonet/coalescent/AbstractDataCollection.java):
	* Methods for setting up the set X 
* [AbstractComputeMinCostTask.java](main/phylonet/coalescent/AbstractComputeMinCostTask.java):
* [AbstractWeightCalculator.java](main/phylonet/coalescent/AbstractWeightCalculator.java):

### Relevant to input/output

Classes relevant to the command-line GUI, including handling of input and output, are:

* [CommandLine.java](main/phylonet/coalescent/CommandLine.java): Reads command line input using the JSAP library and starts the execution
* [Options.java](main/phylonet/coalescent/Options.java): Saves user input options 
* [GlobalMaps.java](main/phylonet/coalescent/GlobalMaps.java): Static class, keeping instances of singleton classes that should be accessed only through here:
	* taxon identifier
	* random number generator
	* taxon name map 	
* [TaxonIdentifier.java](main/phylonet/coalescent/TaxonIdentifier.java):
 Maps taxon names to taxon IDs.  
* [TaxonNameMap.java](main/phylonet/coalescent/TaxonNameMap.java):
Maps the names in the gene trees to the names in the species tree 
* [SpeciesMapper.java](main/phylonet/coalescent/SpeciesMapper.java): Maps IDs between gene trees and the species tree
* [Solution.java](main/phylonet/coalescent/Solution.java): Used for keeping the output


### Relevant to ASTRAL

Classes specific to ASTRAL are:

* [WQClusterCollection.java](main/phylonet/coalescent/WQClusterCollection.java): Doesn't do anything. Work done in the abstract class. 
* [WQComputeMinCostTask.java](main/phylonet/coalescent/WQComputeMinCostTask.java): Doesn't do anything. Work done in the abstract class.
* [WQDataCollection.java](main/phylonet/coalescent/WQDataCollection.java)
* [WQInference.java](main/phylonet/coalescent/WQInference.java)
* [WQWeightCalculator.java](main/phylonet/coalescent/WQWeightCalculator.java)
* [BipartitionWeightCalculator.java](main/phylonet/coalescent/BipartitionWeightCalculator.java)
* [Posterior.java](main/phylonet/coalescent/Posterior.java)
* [SimilarityMatrix.java](main/phylonet/coalescent/SimilarityMatrix.java)
* [Tripartition.java](main/phylonet/coalescent/Tripartition.java)







