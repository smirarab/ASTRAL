#!/bin/bash

cd main/gputest/

javac -Xlint:deprecation -cp ../../JOCL-0.2.0RC.jar GPUCall.java Main.java
java -Djava.library.path=../../ -cp ../../JOCL-0.2.0RC.jar:. Main

cd ../../
