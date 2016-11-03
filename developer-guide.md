This document will be used to document the code. Algorithms are given in the papers and will not be repeated here. 

### Design

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

### Code components

The code is in several packages:

* `phylonet.coalescent`: this is very the majority of the code resides. 
* `phylonet.tree.io`, `phylonet.tree.model.sti`, `phylonet.util`: these each include a single class, overwriting parts of java or phylonet for our purposes. 

We focus on `phylonet.coalescent`. 

First, some interfaces and abstract classes are:

* [IClusterCollection.java](), [AbstractClusterCollection.java]():
* [AbstractInference.java]():
* [AbstractDataCollection.java]():
* [AbstractComputeMinCostTask.java]():
* [AbstractWeightCalculator.java]():


Classes relevant to the command-line GUI, including handling of input and output, are:

* [CommandLine.java]()
* [Options.java]()
* [GlobalMaps.java]()
* [TaxonIdentifier.java]()
* [TaxonNameMap.java]()
* [SpeciesMapper.java]()
* [Solution.java]()


Classes specific to ASTRAL are:

* [WQClusterCollection.java]()
* [WQComputeMinCostTask.java]()
* [WQDataCollection.java]()
* [WQInference.java]()
* [WQWeightCalculator.java]()
* [BipartitionWeightCalculator.java]()
* [Posterior.java]()
* [SimilarityMatrix.java]()
* [Tripartition.java]()







