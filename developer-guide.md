This document will be used to document the code. Algorithmic techniques are given in the paper and will not be repeated here. 

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

### Code components



