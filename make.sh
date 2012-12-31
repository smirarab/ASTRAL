#!/bin/sh
version=2.2.4
cd main
rm phylonet/coalescent/*.class phylonet/util/BitSet.class phylonet/tree/model/sti/STITreeCluster*.class 
javac -classpath ../main.jar:../lib/jsr166.jar phylonet/util/BitSet.java phylonet/coalescent/*ClusterCollection*java phylonet/coalescent/DuplicationWeightCounter.java phylonet/util/BitSet.java phylonet/coalescent/ComputeMinCostTask.java phylonet/coalescent/MGDInference_DP.java phylonet/coalescent/CannotResolveException.java phylonet/tree/model/sti/STITreeCluster.java phylonet/coalescent/DeepCoalescencesCounter.java phylonet/coalescent/STBipartition.java
jar cvfm ../mgd.$version.jar ../manifest.text phylonet/util/BitSet.* phylonet/coalescent/DuplicationWeightCounter*.* phylonet/coalescent/*ClusterCollection* phylonet/coalescent/MGDInference_DP*.* phylonet/coalescent/CannotResolveException.* phylonet/tree/model/sti/STITreeCluster*.class phylonet/tree/model/sti/STITreeCluster.java phylonet/coalescent/DeepCoalescencesCounter.* phylonet/coalescent/ComputeMinCostTask.* phylonet/coalescent/STBipartition.*
cd ..
chmod +x mgd.$version.jar
cp mgd.$version.jar mgd.jar
sed -e "s/__mgd.jar__/mgd.$version.jar/g" README.template > README
rm -r DynaDup/*
mkdir DynaDup
cd DynaDup
ln -s ../lib .
ln -s ../main.jar .
ln -s ../README .
ln -s ../mgd.$version.jar .
cd ..
zip -r DynaDup.$version.zip DynaDup 
