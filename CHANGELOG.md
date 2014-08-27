- Version 4.6.3:
  - Improved scalability of missing taxon completion algorithm (n^2 logn+n|X|^2)
  - 
  - A bit of code refactoring
  - Adding more stderr logs
  - Added option to prevent addition of extra taxa
  
- Version 4.6.2:
  - Merge in 4.4.3 and 4.4.4 to mutli-allele branch
  
- Version 4.6.1:
  - Merge in 4.4.2 changes to mutli-allele branch
  
- Version 4.6.0:
  - Incorporate bug fix from 4.4.1
  
- Version 4.6.0:
  - Handle unresolved gene trees
    
- Version 4.5.1:
    - This version automatically adds extra bipartitions using a calculations based on quartet distances between two taxa 	

- Version 4.5.0:
  - ASTRAL can now handle multi individual gene trees

- Version 4.4.4:
  - Added a script to fix support values on output files that were incorrect from version 4.3.1 to 4.4.1

- Version 4.4.3:
  - Print bootstrap support as internal branch labels instead of branch length values
  - ** Bug fix**: on some machines greedy used to choke (see commit 807edf)
    
- Version 4.4.2:
  - **Bug Fix (IMPORTANT):** Support values drawn on the main tree were incorrect in previous versions since 4.3.1 (related to rerooting of trees). 
  - Prompts changed slightly 
 
- Version 4.4.1:
  - Print a user-friendly error when extra trees have taxon not in main trees.

- Version 4.4.0:
  - **Bug fix:** when gene trees had extreme levels of missing taxa, and there existed two taxa that never appeared together in the same gene tree, an uncaught division by zero could leave to wrong placement of one or both taxa. 
    
- Version 4.3.1:
  - Output consensus bootstrap tree as well
  - Output support only for internal branches
  - Fix default seed number to 692
  - Root all output trees *arbitrarily* on one taxon

- Version 4.3.0 (**command-line options changed!** Notice some commands are not backward compatible!):
  - Performs bootstrapping (`-b`, `-r`, `-g`, `-s` added)
  - Uses a library for argument parsing. We had to change some of the option names (exact version: `-xt` --> `-x`, extra trees: `-ex` --> `-e`).
  - Does not output the number of quartets satisfied by each branch; that information was not interpretable. 

- Version 4.2.1:
  - Read input gene trees with internal node labels
  - Updated the prompts and messages (e.g. Score instead of Cost)
  - Software outputs the version number too

- Version 4.2.0:
  - **Bug fix:** overflow happened for large number of genes.  Use long instead of int
  - Improve the command-line help a bit

- Version 4.1.1:
  - Added a hidden option to force algorithm used for estimation (only for debugging) 

- Version 4.1:
  - Updated the command-line so that ASTRAL algorithm is the default and DynaDup specific options are removed (-wq omitted and is default now)

- Version 4.0:
  - Added support for handling missing data by completing bipartitions
