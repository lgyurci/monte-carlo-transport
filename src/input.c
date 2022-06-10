#include <stdio.h>
#include <stdlib.h>
#include "tmath.h"
#include "input.h"
#include <string.h>

void copySFromInto(char* a, char* b){
    int i = 0;
    for (; i < strlen(a); i++){
        b[i] = a[i];
    }
    b[i] = '\0';
}

void setDefaultInputs(struct input_data * id){
    id->channels = 1024;
    id->cylinder.x = 2;
    id->cylinder.y = -2;
    id->cylinder.z = 2;
    copySFromInto("gnuplot\0",id->gnuplotExecutable);
    id->particleNum = 0;
    id->realtime = 1;
    id->rho = 3.67;
    copySFromInto("\0",id->savefile);
    id->sigma = 15;
    id->sourceEnergy = 4000;
    id->threadCount = 20;
    id->updateFreq = 1e4;
    id->sourcePos.x = 3;
    id->sourcePos.y = 3;
    id->sourcePos.z = 3;
    copySFromInto("crs_NaI",id->xcomLocation);
}

void configure(struct input_data * id,char *option,char *params){

}


struct input_data getInput(int argc,char** argv){
    struct input_data id;
    setDefaultInputs(&id);
    char config[255] = "config\0";

    if (argc > 1) {
        copySFromInto(argv[1],config);
    }

    FILE *f;
    f = fopen(config,"r");

    if (f == NULL) {
        printf("Opening config file named '%s' failed.\n",config);
        exit(1);
    } else {
        char line[200];
        int linecnt = 0;

        while (fgets(line,200,f) != NULL){
            linecnt++;
            char option[50];
            int ol = 0;
            for (; line[ol] != ':' && line[ol] != '\0';ol++){
                option[ol] = line[ol];
            }
            option[ol] = '\0';
            ol++;
            configure(&id,option,&(line[ol]));
        }
    }

    return id;
}
