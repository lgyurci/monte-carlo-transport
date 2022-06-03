#include "transport.h"
#include "tmath.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "reactions.h"

#define inf 1e9
#define pi 3.14159
#define channels 1024

double pztop = 2;
double pzbottom = -2;
double radius = 2;

int ch[channels] = { 0 };

double maxenergy = 600;

void traceParticle(struct vector pos, struct vector direction, double particle_energy){
    int particle_is_alive = 1;
    double sumen = 0;
    while (particle_is_alive == 1){
        double lastPoint = intersect_cylinder_in(pos,direction,pztop,pzbottom,radius);
        double freeway = shuffle_freeway_length(particle_energy);
        if (freeway < lastPoint){
            pos = transloc(pos,vMult(direction,freeway));
            int reac = shuffle_reaction(particle_energy);
            double eenergy = 0;
            if (reac == 1){ //Compton
                eenergy = particle_energy;
                comptonScatter(direction,&particle_energy);
                eenergy = eenergy - particle_energy;
                sumen += eenergy;
            }
            if (reac == 2){ //Photoelectric effect
                particle_is_alive = 0;
                eenergy = particle_energy;
                sumen += eenergy;
            }
        }
        else {
            particle_is_alive = 0;
        }
    }

    if (sumen > 0) ch[((int)(sumen*channels/maxenergy))]++;
}

int main(){
    initRandrand();
    initReactions();
    int particleNum = 1e6;

    struct vector sourcePos;
    sourcePos.x = 3;
    sourcePos.y = 3;
    sourcePos.z = 3;

    double sourceEnergy = 500;
    int colled = 0;

    for (int i = 0; i < particleNum && colled == 0; i++){
        struct vector direction = isotropicDirection(); //TODO kúpban iránysorsolás
        double firstPoint = intersect_cylinder_out(sourcePos,direction,pztop,pzbottom,radius);
        if (firstPoint < inf){
            traceParticle(transloc(sourcePos,vMult(direction,firstPoint)),direction,sourceEnergy);
        }
    }

    for (int i = 0; i < channels; i++){
        printf("%d\n",ch[i]);
    }
    freeReactions();
    return 0;
}

/*void testrand(){
    struct timeval ct; //idő lekérésere szolgáló struct
    gettimeofday(&ct,NULL);
    long peTi = ct.tv_sec * (int)1e6 + ct.tv_usec;
    for (int i=0;i<1e9;i++){
        drand();
    }
    gettimeofday(&ct,NULL);
    printf ("stdlib rand: %f\n",(double)(ct.tv_sec * (int)1e6 + ct.tv_usec - peTi)/1e6);

    MTRand r = seedRand(1337);

    gettimeofday(&ct,NULL);
    peTi = ct.tv_sec * (int)1e6 + ct.tv_usec;
    for(int i=0; i<1e9; i++) {
        genRand(&r);
    }
    gettimeofday(&ct,NULL);
    printf ("Mersenne-twister: %f\n",(double)(ct.tv_sec * (int)1e6 + ct.tv_usec - peTi)/1e6);
}*/
