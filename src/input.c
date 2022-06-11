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

double parseDouble(char *line){
    int ol = 0;
    char tc[20] = {0};
    int i;
    for(i = 0; line[ol] != ';' && line[ol] != '\0'; ol++){
        tc[i++] = line[ol];
    }
    tc[i] = '\0';
    return atof(tc);
}

int parseInt(char *line){
    int ol = 0;
    char tc[20] = {0};
    int i;
    for(i = 0; line[ol] != ';' && line[ol] != '\0'; ol++){
        tc[i++] = line[ol];
    }
    tc[i] = '\0';
    return atoi(tc);
}

long parseLong(char *line){
    int ol = 0;
    char tc[20] = {0};
    int i;
    for(i = 0; line[ol] != ';' && line[ol] != '\0'; ol++){
        tc[i++] = line[ol];
    }
    tc[i] = '\0';
    return strtol(tc);
}



void parseString(char *line,char *tc){
    int ol = 0;
    int i;
    for(i = 0; line[ol] != ';' && line[ol] != '\0'; ol++){
        tc[i++] = line[ol];
    }
    tc[i] = '\0';
}

char *parseVector(char *tc, int ol){
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
            char tc[20] = {0};

            if (strcmp(option,"Thread count") == 0) {
                id->threadCount = parseInt(line[ol++]);
                if (id->threadCount <= 0) {
                    exit(0);
                }
            }

            if (strcmp(option,"Channels") == 0) {
                id->channels = parseInt(line[ol++]);
                if (id->channels <= 0) {
                    exit(0);
                }
            }

            if (strcmp(option,"Stop when collisions reach") == 0) {
                id->particleNum = parseLong(line[ol++]);
                if (id->particleNum < 0) {
                    exit(0);
                }
            }

            if (strcmp(option,"Real time plotting with gnuplot") == 0) {
                int b = parseInt(line[ol++]);
                if (b == 0 || b == 1) id->realtime = b;
            }

            if (strcmp(option,"Density (g/cm3)") == 0) {
                id->rho = parseDouble(line[ol++]);
                if (id->rho <= 0) exit(0);
            }

            if (strcmp(option,"FWHM") == 0) {
                id->sigma = parseDouble(line[ol++]);
                if (id->sigma < 0) exit(0);
            }

            if (strcmp(option,"Energy of the source") == 0) {
                id->sourceEnergy = parseDouble(line[ol++]);
                if (id->sourceEnergy <= 0) exit(0);
            }

            if (strcmp(option,"Update real time every N collisions") == 0) {
                id->updateFreq = parseLong(line[ol++]);
                if (id->updateFreq < 0) exit(0);
            }

            if (strcmp(option,"gnuplot executable") == 0) {
                char tc [255] = {0};
                parseString(line,tc);
                copySFromInto(tc,id->gnuplotExecutable);
            }
        }
    }

    return id;
}
