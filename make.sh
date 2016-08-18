#!/bin/bash

set -u
set -e
set -x
set -o pipefail

version=`grep _versinon main/phylonet/coalescent/CommandLine.java|grep String|sed -e "s/.*= .//g" -e "s/.;//g"`
echo Version $version

cd main

rm -f phylonet/coalescent/*.class phylonet/util/BitSet.class phylonet/tree/model/sti/STITreeCluster*.class phylonet/tree/io/NewickWriter.class

javac -source 1.7 -target 1.7 -classpath ../lib/main.jar:../lib/colt.jar:../lib/JSAP-2.1.jar:../lib/jocl-2.0.0.jar phylonet/util/BitSet.java phylonet/coalescent/*.java phylonet/tree/model/sti/STITreeCluster.java phylonet/tree/io/NewickWriter.java

jar cvfm ../astral.$version.jar ../manifest.text phylonet/util/BitSet.* phylonet/coalescent/*.* phylonet/tree/model/sti/STITreeCluster*.* phylonet/tree/io/NewickWriter.*

cd ..

chmod +x astral.$version.jar
sed -e "s/__astral.jar__/astral.$version.jar/g" -e "s/__astral.zip__/Astral.$version.zip/g" README.template.md > README.md
rm -fr Astral/*
mkdir -p  Astral
cd Astral
ln -s ../lib .
ln -s ../README.md .
ln -s ../astral.$version.jar .
ln -s ../main/test_data .
ln -s ../astral-tutorial.md .
ln -s ../thesis-astral.pdf .
cd ..
