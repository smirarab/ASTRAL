#!/bin/bash
i=1;
counter=0;
output=0;
sum=0
while [ $i -le 50 ]
do
rep=`printf "%02d\n" $i`
output=$(sed '3!d' ${rep}/astral-gpu-comparison.txt)
if ! [ -z $output ]
then
counter=$((counter+1))
sum=$(bc <<< "scale = 5; $sum + $output")
fi
i=$((i+1))
done
bc <<< "scale = 5; $sum/$counter" > stats.txt
