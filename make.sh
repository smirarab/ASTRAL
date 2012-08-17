#!/bin/sh
cd main
rm phylonet/util/BitSet.class phylonet/coalescent/*.class phylonet/util/BitSet.class phylonet/tree/model/sti/STITreeCluster.class 
javac -classpath ../main.jar:../lib/jsr166.jar phylonet/util/BitSet.java phylonet/coalescent/*ClusterCollection*java phylonet/coalescent/DuplicationWeightCounter.java phylonet/util/BitSet.java phylonet/coalescent/ComputeMinCostTask.java phylonet/coalescent/MGDInference_DP.java phylonet/coalescent/CannotResolveException.java phylonet/tree/model/sti/STITreeCluster.java phylonet/coalescent/DeepCoalescencesCounter.java
jar cvfm ../mgd.jar ../manifest.text phylonet/util/BitSet.class phylonet/coalescent/DuplicationWeightCounter*.class phylonet/coalescent/*ClusterCollection*class phylonet/coalescent/MGDInference_DP*.class phylonet/coalescent/CannotResolveException.class phylonet/tree/model/sti/STITreeCluster.class phylonet/coalescent/DeepCoalescencesCounter.class phylonet/coalescent/ComputeMinCostTask.class
chmod +x ../mgd.jar
