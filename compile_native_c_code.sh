for word in `java -XshowSettings:properties -version 2>&1 | grep java.home`; do DIR=${word}/../include/; done && for word in `ls -d ${DIR}*/`; do SUBDIR=${word}; done || (( DIR=. && SUBDIR=. ))
if [ -z $JAVA_HOME ]; then DIR2=.; SUBDIR2=.; else DIR2=$JAVA_HOME; for word in `ls -d ${DIR2}/*/`; do SUBDIR2=${word}; done || (( DIR2=. && SUBDIR2=. )); fi
mv lib/libAstral.so lib/libAstral.so.old || echo "no default library found"
g++ -std=c++11 -I"$DIR" -I"$SUBDIR" -I"$DIR2" -I"$SUBDIR2" -march=native -Ofast -fPIC -o lib/libAstral.so -shared main/phylonet_coalescent_Polytree_PTNative.cpp || mv lib/libAstral.so.old lib/libAstral.so || echo "Native library not compiled and default version not found."
