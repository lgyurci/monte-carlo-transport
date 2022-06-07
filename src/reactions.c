#include <stdio.h>
#include <stdlib.h>
#include "transport.h"
#include "tmath.h"
#include "wyhash.h"

double **xcom;
int x;
int y;
double ro;
double sourceEnergy;
double *crossectionsAtSourceEnergy;
double sumCrossectionsAtSourceEnergy;

double getCrossection(double energy,int effect){ //effect: 0 - total, 1 - compton scattering, 2 - photoelectric effect, 3 - pair production in nucleus, 4 - pair production in electron bands
    if (energy == sourceEnergy) { //cache a forrás energiájára, mert ez gyakran kell
        if (effect == 0){
            return sumCrossectionsAtSourceEnergy;
        } else {
            return crossectionsAtSourceEnergy[effect];
        }
    }
    if (energy <= xcom[0][0]) return 0;
    if (effect == 0){
        double sum = 0;
        for (int i = 1; i < x; i++){
            sum += getCrossection(energy,i);
        }
        return sum;
    }
    int i;
    energy = energy/1000;
    for (i = 0; i < y && energy > xcom[0][i]; i++);
    double energy_low = xcom[0][i-1];
    double energy_high = xcom[0][i];
    double crossection_low = xcom[effect][i-1];
    double crossection_high = xcom[effect][i];
    double crossection = ((energy - energy_low)/(energy_high - energy_low))*(crossection_high - crossection_low) + crossection_low;
    return crossection*ro;
}

void getCrs(double energy, double *crs){

    energy = energy/1000;
    if (energy <= xcom[0][0] || energy > xcom[0][y-1]) {
        for (int j = 0; j < x; j++){
            crs[j] = 0;
        }
        return;
    }

    crs[0] = 0;

    int j;
    for (j = 0; j < y && energy > xcom[0][j]; j++);

    double energy_low = xcom[0][j-1];
    double energy_high = xcom[0][j];

 /*   int eil = 0;
    int eih = y-1;

    while (eih-eil != 1){
        int a = eil+(eih-eil)/2;
        if (xcom[0][a] < energy){
            eil = a;
        } else {
            eih = a;
        }

    }

    double energy_low = xcom[0][eil];
    double energy_high = xcom[0][eih];*/

    for (int i = 1; i < x; i++){
        double crossection_low = xcom[i][j-1];
        double crossection_high = xcom[i][j];
        double crossection = ((energy - energy_low)/(energy_high - energy_low))*(crossection_high - crossection_low) + crossection_low;
        crs[i] = crossection*ro;
        crs[0] += crs[i];
    }
}

int getX(){
    return x;
}

void initReactions(double roo,double source){
    sourceEnergy = -1;
    ro = roo;
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

        crossectionsAtSourceEnergy = malloc(sizeof(double)*x);
        sumCrossectionsAtSourceEnergy = 0;
        for (int i = 1; i < x; i++){
            crossectionsAtSourceEnergy[i] = getCrossection(source,0);
            sumCrossectionsAtSourceEnergy += crossectionsAtSourceEnergy[i];
        }
        sourceEnergy = source;
    }
}

/*double calcTotalCrossection(double energy,double *crossections){
    double total = 0;
    for (int i = 1; i < x; i++){
        crs[i] = getCrossection(energy,i);
        total += crs[i];
    }
    return total;
}*/

double shuffle_freeway_length(uint64_t *random, double *crs){
    //double lambda = -1*log(drandt(random))/getCrossection(energy,0);
    double lambda = -1*log(drandt(random))/crs[0];
    //double lambda = -1*log(drandt(random))/getCrossection(energy,0);
    return lambda;
}

int shuffle_reaction(uint64_t *random, double *crs){
    double r = drandt(random)*crs[0];
    double s = crs[1];
    for (int i = 1; i <= x; i++){
        if (r <= s){
            return i;
        }
        if (i >= 4) return -1;
        s += crs[i+1];
    }
    //printf("asdff somethingiswrong\n");
    return 0;
}

/*int shuffle_reaction(double energy,uint64_t *random,double totalCrossection){
    double r = drandt(random)*getCrossection(energy,0);
    double s = getCrossection(energy,1);
    for (int i = 1; i <= x; i++){
        if (r <= s){
            return i;
        }
        if (i >= 4) return -1;
        s += getCrossection(energy,i+1);
    }
    return 0;
}*/

void freeReactions(){
    for (int i = 0; i < x; i++){
        free(xcom[i]);
    }
    free(xcom);
    free(crossectionsAtSourceEnergy);
}

struct vector comptonScatter(struct vector direction,double *energy,uint64_t *random){

    double lambda = 511/(*energy);

    double R = 0;
    double X = 0;
    int foundx = 0;

    while (foundx == 0){
        double r1 = drandt(random);
        double r2 = drandt(random);
        double r3 = drandt(random);
        if (r1 <= (1+2/lambda)/(9+2/lambda)){
            R = 1+2*r2/lambda;
            if (r3 <= 4*(1/R-1/(R*R))){
                X = R;
                foundx = 1;
            }
        } else {
            R = (1+2/lambda)/(1+2*r2/lambda);
            double a = lambda-R*lambda+1;
            if (r3<=0.5*(a*a+1/R)){
                X = R;
                foundx = 1;
            }
        }
    }

    double E2 = 511/(X*lambda);

    double costheta = 1+lambda-lambda*X;

    struct vector ret = isotropicScatter(direction,costheta,random);
    *energy = E2;
    return ret;
}
/*

int main(){
    initReactions();
    freeReactions();
}*/
