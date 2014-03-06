#!/bin/sh
version=3.1.1

cd main

rm phylonet/coalescent/*.class phylonet/util/BitSet.class phylonet/tree/model/sti/STITreeCluster*.class 

javac -classpath ../main.jar:../lib/jsr166.jar phylonet/util/BitSet.java phylonet/coalescent/*java phylonet/util/BitSet.java phylonet/tree/model/sti/STITreeCluster.java

jar cvfm ../astral.$version.jar ../manifest.text phylonet/util/BitSet.* phylonet/coalescent/*.* phylonet/tree/model/sti/STITreeCluster*.*

cd ..

chmod +x astral.$version.jar
cp astral.$version.jar astral.jar
sed -e "s/__astral.jar__/astral.$version.jar/g" README.template > README
rm -r Astral/*
mkdir Astral
cd Astral
ln -s ../lib .
ln -s ../main.jar .
ln -s ../README .
ln -s ../astral.$version.jar .
cd ..
zip -r Astral.$version.zip Astral 
