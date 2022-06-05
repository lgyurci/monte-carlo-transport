#include "transport.h"
#include "tmath.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "reactions.h"
#include "mtwister.h"

#define inf 1e9
#define pi 3.14159
#define channels 1024

double pztop = 2;
double pzbottom = -2;
double radius = 2;

int ch[channels] = { 0 };

double maxenergy = 4100;

//double ro = 3.67;

double traceParticle(struct vector pos, struct vector direction, double particle_energy, MTRand *random){ //returns: absorbed energy by the crystal
    int particle_is_alive = 1;
    double sumen = 0;
    while (particle_is_alive == 1){
        double lastPoint = intersect_cylinder_in(pos,direction,pztop,pzbottom,radius);
        double freeway = shuffle_freeway_length(particle_energy,random);
        if (freeway < lastPoint){
            pos = transloc(pos,vMult(direction,freeway));
            int reac = shuffle_reaction(particle_energy,random);
            double eenergy = 0;
            if (reac == 1){ //Compton scattering
                eenergy = particle_energy;
                direction = comptonScatter(direction,&particle_energy,random);
                eenergy = eenergy - particle_energy;
            }
            if (reac == 2){ //Photoelectric effect
                particle_is_alive = 0;
                eenergy = particle_energy;
            }
            if (reac == 3 || reac == 4){ //Pair production in nucleus/electron band
                particle_is_alive = 0;
                eenergy = particle_energy - 1022;
                struct vector annihil1 = isotropicDirection(random);
                struct vector annihil2 = vMult(annihil1,-1);
                eenergy += traceParticle(pos,annihil1,511,random);
                eenergy += traceParticle(pos,annihil2,511,random);
              //  printf("%f\n",sumen);
            }

            sumen += eenergy;
        }
        else {
            particle_is_alive = 0;
        }
    }

    return sumen;
}

int main(){
  //  initRandrand();
    double roh = 3.67;
    initReactions(roh);

    int particleNum = 5e5;
    double sourceEnergy = 4000;
    double fwhm = 30;
    double sigma = fwhm/2.355;

    MTRand mainrand = initRandrandt();

    struct vector sourcePos;
    sourcePos.x = 3;
    sourcePos.y = 3;
    sourcePos.z = 3;

    double sinalpha = sqrt(radius*radius+(pztop-pzbottom)*(pztop-pzbottom))/vAbs(sourcePos);
    double cosalpha = sqrt(1-sinalpha*sinalpha);

    //struct timeval ct; //idő lekérésere szolgáló struct
    //gettimeofday(&ct,NULL);
    //long peTi = ct.tv_sec * (int)1e6 + ct.tv_usec;
    int colledParticles = 0;

    for (int i = 0; colledParticles < particleNum; i++){
        //struct vector direction = isotropicDirection();
        struct vector direction = coneDirection(cosalpha,vMult(sourcePos,-1),&mainrand);
        double firstPoint = intersect_cylinder_out(sourcePos,direction,pztop,pzbottom,radius);
        if (firstPoint < inf){
            double energyAbsorbed = traceParticle(transloc(sourcePos,vMult(direction,firstPoint)),direction,sourceEnergy,&mainrand);
            if (energyAbsorbed > 0) {
                colledParticles++;
                energyAbsorbed = energyAbsorbed + sigma*drandnt(&mainrand);
                if (maxenergy < energyAbsorbed) energyAbsorbed = maxenergy;
                ch[((int)(energyAbsorbed*(channels-1)/maxenergy))]++;
            }
        }
    }

    //gettimeofday(&ct,NULL);
    //printf ("Average speed: %f Mp/s\n",colledParticles/1e6/((double)(ct.tv_sec * (int)1e6 + ct.tv_usec - peTi)/1e6));

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
