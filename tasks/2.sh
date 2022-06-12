#!/bin/bash
for i in {0..10}
do
    x=$(echo "1-0.5*$i" | bc)
    y=$(echo "3.5-0.5*$i" | bc)
    a="Thread count: 20; #number of computing threads the program should run on (this is one less than the total threads)\r\nChannels: 1024;\r\nStop when collisions reach: 5e6;\r\nReal time plotting with gnuplot: 0;\r\nDensity (g/cm3): 3.67;\r\nFWHM: 6;\r\nEnergy of the source: 661.7;\r\nUpdate real time every N collisions: 1e4;\r\ngnuplot executable: gnuplot;\r\nCross section file generated by xcom: ../crs_NaI;\r\nSave file:; #leave it blank if you dont want to save\r\nParticle source position: ($x,$y,2);\r\nCylinder (top,bottom,radius): (1.5,-1.5,2.5);"
    echo -e $a > tmp
    ./mtransport tmp
    rm tmp
done
