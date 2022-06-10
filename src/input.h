#ifndef input_H_
#define input_H_

#include <stdio.h>
#include <stdlib.h>
#include "tmath.h"


struct input_data{
    int channels;
    struct vector cylinder;
    struct vector sourcePos;
    double rho;
    int threadCount;
    unsigned long particleNum;
    unsigned long updateFreq;
    int realtime;
    double sourceEnergy;
    double sigma;
    char gnuplotExecutable[255];
    char xcomLocation[255];
    char savefile[255];
};

struct input_data getInput(int argc,char** argv);


#endif
