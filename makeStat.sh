#! /bin/bash

OUTPUT_DIR="stats"

for file in `ls /media/angelo/BigLinux/Documents/Dropbox/2017_SCHEDULER/inputTest/inputTest*`
do
	#echo $file
	
	Debug/SensorsClustering -f $file -k 7 -e 4 -l 16 | grep "Maximum MAX correlation"
done

#Debug/SensorsClustering -f inputScenario/inputTest.dat -g inputScenario/inputTest.dat -n 120 -o outputRis/outputTest.dot -k 9 -s 100 -i 50 -e 3 -r 1

#Debug/SensorsClustering -f inputScenario/inputTest.dat -o outputRis/outputTest.dot -k 7 -i 100 -e 3
