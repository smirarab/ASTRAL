#!/bin/bash

set -u
set -e
set -x

version=`grep _versinon main/phylonet/coalescent/CommandLine.java|grep String|sed -e "s/.*= .//g" -e "s/.;//g"`
echo Version $version

cd main

rm -f phylonet/coalescent/*.class phylonet/util/BitSet.class phylonet/tree/model/sti/STITreeCluster*.class phylonet/tree/io/NewickWriter.class

javac -source 1.7  -target 1.7 -classpath ../lib/main.jar:../lib/colt.jar:../lib/JSAP-2.1.jar phylonet/util/BitSet.java phylonet/coalescent/*java phylonet/tree/model/sti/STITreeCluster.java phylonet/tree/io/NewickWriter.java

jar cvfm ../astral.$version.jar ../manifest.text phylonet/util/BitSet.* phylonet/coalescent/*.* phylonet/tree/model/sti/STITreeCluster*.* phylonet/tree/io/NewickWriter.*

cd ..

chmod +x astral.$version.jar
sed -e "s/__astral.jar__/astral.$version.jar/g" -e "s/__astral.zip__/Astral.$version.zip/g" README.template.md > README.md
sed -e "s/__astral.jar__/astral.$version.jar/g" -e "s/__astral.zip__/Astral.$version.zip/g" astral-tutorial-template.md > astral-tutorial.md
rm -fr Astral/*
mkdir -p  Astral
cd Astral
ln -s ../lib .
ln -s ../README.md .
ln -s ../astral.$version.jar .
ln -s ../main/test_data .
ln -s ../astral-tutorial.pdf .
ln -s ../thesis-astral.pdf .
cd ..
rm -f Astral.$version.zip
zip -r Astral.$version.zip Astral 

set +x
echo "
Build finished successfully. You can distribute Astral.$version.zip or simply run astral.$version.jar. 
  Note that if you are moving astral.$version.jar to some other location, you need to also move the lib directory."
