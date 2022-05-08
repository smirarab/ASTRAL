word=`java -XshowSettings:properties -version 2>&1 | grep java.home|sed -e "s/.* =//g"`
DIR="${word}/../include/"; 
echo $DIR
for word in `ls -d "${DIR}*/"`; do SUBDIR=${word}; done || (( DIR=. && SUBDIR=. ))
if [ -z ${SUBDIR} ]; then DIR=.; SUBDIR=.; fi
if [ -z $JAVA_HOME ]; then DIR2=.; SUBDIR2=.; else DIR2=$JAVA_HOME/../include/; for word in `ls -d ${DIR2}*/`; do SUBDIR2=${word}; done || (( DIR2=. && SUBDIR2=. )); fi
if [ -z ${SUBDIR2} ]; then DIR2=.; SUBDIR2=.; fi
for word in `java -XshowSettings:properties -version 2>&1 | grep java.home`; do DIR3=${word}/include/; done && for word in `ls -d ${DIR3}*/`; do SUBDIR3=${word}; done || (( DIR3=. && SUBDIR3=. ))
if [ -z ${SUBDIR3} ]; then DIR3=.; SUBDIR3=.; fi
if [ -z $JAVA_HOME ]; then DIR4=.; SUBDIR4=.; else DIR4=$JAVA_HOME/include/; for word in `ls -d ${DIR4}*/`; do SUBDIR4=${word}; done || (( DIR4=. && SUBDIR4=. )); fi
if [ -z ${SUBDIR4} ]; then DIR4=.; SUBDIR4=.; fi
if [ `uname -s` == Darwin ]; then LIBNAME=libAstral.dylib; else LIBNAME=libAstral.so; fi
mv lib/${LIBNAME} lib/${LIBNAME}.old || echo "no default library found"
g++ -std=c++11 -I"$DIR" -I"$SUBDIR" -I"$DIR2" -I"$SUBDIR2" -I"$DIR3" -I"$SUBDIR3" -I"$DIR4" -I"$SUBDIR4" -march=native -Ofast -fPIC -o lib/${LIBNAME} -shared main/phylonet_coalescent_Polytree_PTNative.cpp || icc -std=c++11 -I"$DIR" -I"$SUBDIR" -I"$DIR2" -I"$SUBDIR2" -I"$DIR3" -I"$SUBDIR3" -I"$DIR4" -I"$SUBDIR4" -march=native -Ofast -fPIC -o lib/${LIBNAME} -shared main/phylonet_coalescent_Polytree_PTNative.cpp || mv lib/${LIBNAME}.old lib/${LIBNAME} || echo "Native library not compiled and default version not found." 
