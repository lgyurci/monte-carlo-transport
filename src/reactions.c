#include <stdio.h>
#include <stdlib.h>
#include "transport.h"

double **xcom;
int x;
int y;

void initReactions(){
    FILE *f;

    f = fopen("../data.pl","r");


    if (f){
        int rows = 0;
        int columns = 0;
        int rfound = 0;
        char c = 0;
        char c0 = 0;
        int cc = 0;
        while ((c = getc(f)) != EOF){
            if (cc == 1){
                c0 = c;
                cc = 0;
            }
            if (rfound == 0 && c == '\t'){
                columns++;
            }
            if (c == '\n') {
                if (c0 >= 48 && c0 <= 57){
                    rows++;
                }
                rfound = 1;
                cc = 1;
            }

        }

        xcom = malloc(sizeof(double)*columns);
        for (int i = 0; i < columns; i++){
            xcom[i] = malloc(sizeof(double)*rows);
        }

        x = columns;
        y = rows;

        rewind(f);

        char buff[255];
        int bp = 0;
        int xcp = 0;
        int xrp = 0;

        cc = 0;
        c0 = 0;

        while ((c = getc(f)) != EOF){
            if (cc == 1){
                c0 = c;
                cc = 0;
            }

            if (c == '\n') {
                cc = 1;

            }

            if (c0 >= 48 && c0 <= 57){
                if (c == '\t'){
                    if (xcp >= x){
                        xcp = 0;
                        xrp++;
                    }
                    buff[bp++] = '\0';
                    bp = 0;
                    xcom[xcp++][xrp] = atof(buff);
                } else {
                    buff[bp++] = c;
                }
            }
        }

        fclose(f);
    }
}

double getCrossection(double energy,int effect){ //effect: 1 - compton scattering, 2 - photoelectric effect, 3 - pair production in nucleus, 4 - pair production in electron bands
    int i;
    energy = energy/1000;
    for (i = 0; i < y && energy > xcom[0][i]; i++);
    double energy_low = xcom[0][i-1];
    double energy_high = xcom[0][i];
    double crossection_low = xcom[effect][i-1];
    double crossection_high = xcom[effect][i];
    double crossection = ((energy - energy_low)/(energy_high - energy_low))*(crossection_high - crossection_low) + crossection_low;
    return crossection;
}

double shuffle_freeway_length(double energy){
    double sumcross = 0;
    for (int i = 1; i < x; i++){
        sumcross += getCrossection(energy,i);
    }
    double lambda = -1*log(drand())/sumcross;
    return lambda;
}

int shuffle_reaction(double energy){
    double min = -1*log(drand())/getCrossection(energy,1);
    int mini = 1;
    for (int i = 2; i < x;i++){
        double c = -1*log(drand())/getCrossection(energy,i);
        if (c < min){
            min = c;
            mini = i;
        }
    }
    return mini;
}

void freeReactions(){
    for (int i = 0; i < x; i++){
        free(xcom[i]);
    }
    free(xcom);
}
/*

int main(){
    initReactions();
    freeReactions();
}*/
