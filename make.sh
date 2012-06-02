#!/bin/sh
cd main
javac -classpath ../main.jar phylonet/util/BitSet.java phylonet/coalescent/DuplicationWeightCounter.java phylonet/coalescent/MGDInference_DP.java phylonet/coalescent/CannotResolveException.java phylonet/tree/model/sti/STITreeCluster.java
jar cvfm ../mgd.jar ../manifest.text phylonet/util/BitSet.class phylonet/coalescent/DuplicationWeightCounter*.class phylonet/coalescent/MGDInference_DP*.class phylonet/coalescent/CannotResolveException.class phylonet/tree/model/sti/STITreeCluster.class
