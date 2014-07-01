#!/bin/sh

version=`grep _versinon main/phylonet/coalescent/CommandLine.java|grep String|sed -e "s/.*= .//g" -e "s/.;//g"`
echo Version $version

cd main

rm phylonet/coalescent/*.class phylonet/util/BitSet.class phylonet/tree/model/sti/STITreeCluster*.class  

javac -source 1.5  -target 1.5 -classpath ../main.jar:../lib/jsr166.jar:../lib/JSAP-2.1.jar phylonet/util/BitSet.java phylonet/coalescent/*java phylonet/tree/model/sti/STITreeCluster.java

jar cvfm ../astral.$version.jar ../manifest.text phylonet/util/BitSet.* phylonet/coalescent/*.* phylonet/tree/model/sti/STITreeCluster*.*

cd ..

chmod +x astral.$version.jar
sed -e "s/__astral.jar__/astral.$version.jar/g" -e "s/__astral.zip__/Astral.$version.zip/g" README.template > README.md
rm -r Astral/*
mkdir Astral
cd Astral
ln -s ../lib .
ln -s ../main.jar .
ln -s ../README.md .
ln -s ../astral.$version.jar .
ln -s ../main/test_data .
ln -s ../astral-tutorial.pdf .
cd ..
zip -r Astral.$version.zip Astral 
