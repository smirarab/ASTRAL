#!/bin/sh

set -u
set -e
set -x
set -o pipefail

version=`grep _versinon main/phylonet/coalescent/CommandLine.java|grep String|sed -e "s/.*= .//g" -e "s/.;//g"`
echo Version $version

cd main

rm phylonet/coalescent/*.class phylonet/util/BitSet.class phylonet/tree/model/sti/STITreeCluster*.class  

javac -source 1.5  -target 1.5 -classpath ../lib/main.jar:../lib/JSAP-2.1.jar phylonet/util/BitSet.java phylonet/coalescent/*java phylonet/tree/model/sti/STITreeCluster.java

jar cvfm ../astral.$version.jar ../manifest.text phylonet/util/BitSet.* phylonet/coalescent/*.* phylonet/tree/model/sti/STITreeCluster*.*

cd ..

chmod +x astral.$version.jar
sed -e "s/__astral.jar__/astral.$version.jar/g" -e "s/__astral.zip__/Astral.$version.zip/g" README.template > README.md
rm -fr Astral/*
mkdir -p  Astral
cd Astral
ln -s ../lib .
ln -s ../README.md .
ln -s ../astral.$version.jar .
ln -s ../main/test_data .
ln -s ../astral-tutorial.pdf .
cd ..
rm -f Astral.$version.zip
zip -r Astral.$version.zip Astral 

set +x
echo "
Build finished successfully. You can distribute Astral.$version.zip or simply run astral.$version.jar. 
  Note that if you are moving astral.$version.jar to some other location, you need to also move the lib directory."
