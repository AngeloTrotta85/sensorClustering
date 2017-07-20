#!/bin/bash

echo "make all"
echo "Building file: MyCoord.cpp"
echo "Invoking: Cross G++ Compiler"
g++ -std=c++0x -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"src/MyCoord.d" -MT"src/MyCoord.o" -o "src/MyCoord.o" "../src/MyCoord.cpp"
echo "Finished building: MyCoord.cpp"
 
echo "Building file: SensorsClustering.cpp"
echo "Invoking: Cross G++ Compiler"
g++ -std=c++0x -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"src/SensorsClustering.d" -MT"src/SensorsClustering.o" -o SensorsClustering.o SensorsClustering.cpp
echo "Finished building: SensorsClustering.cpp"
 
echo "Building target: SensorsClustering"
echo "Invoking: Cross G++ Linker"
g++  -o SensorsClustering  MyCoord.o SensorsClustering.o
