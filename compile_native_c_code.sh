for word in `java -XshowSettings:properties -version 2>&1 | grep java.home`; do test -d ${word}/include/ && DIR=${word}/include/ || DIR=${word}/../include/; done && for word in `ls -d ${DIR}*/`; do SUBDIR=${word}; done || (( DIR=. && SUBDIR=. ))
test `uname` == "Darwin" && ext='dylib' || ext='so'
if [ -z $JAVA_HOME ]; then DIR2=.; SUBDIR2=.; else DIR2=$JAVA_HOME; for word in `ls -d ${DIR2}/*/`; do SUBDIR2=${word}; done || (( DIR2=. && SUBDIR2=. )); fi
mv lib/libAstral.${ext} lib/libAstral.${ext}.old || echo "no default library found"
g++ -std=c++11 -I"$DIR" -I"$SUBDIR" -I"$DIR2" -I"$SUBDIR2" -march=native -Ofast -fPIC -o lib/libAstral.${ext} -shared main/phylonet_coalescent_Polytree_PTNative.cpp || mv lib/libAstral.${ext}.old lib/libAstral.${ext} || echo "Native library not compiled and default version not found."
