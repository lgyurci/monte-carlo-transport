#include "transport.h"
#include "tmath.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "reactions.h"
#include "mtwister.h"
#include <pthread.h>

#define inf 1e9
#define pi 3.14159
#define channels 1024

double pztop = 2;
double pzbottom = -2;
double radius = 2;

//double maxenergy = 4100;

struct tracingThreadArgs{
    MTRand *random;
    int *chan;
    int stop_at_colls;
    double sourceEnergy;
    double cosalpha;
    double sigma;
    struct vector sourcePos;
    int particlesTraced;
    int detectorParticles;
    double sumEnergy;
    double maxEnergy;
};

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

void plot(FILE *gp_pipe,int *chan, int chansize){
    fprintf(gp_pipe,"plot '-'\n");
    for (int i = 0; i < chansize; i++){
        fprintf(gp_pipe,"%d %d\n",i,chan[i]);
        //printf("%d %f\n",i,ch[i]);
    }
    fprintf(gp_pipe,"e\n");
}

void *tracingThread(void* arg){
    struct tracingThreadArgs *targs = (struct tracingThreadArgs*) arg;

    int colledParticles = 0;
    int i = 0;

    for (; colledParticles < targs->stop_at_colls; i++){
        //struct vector direction = isotropicDirection();
        struct vector direction = coneDirection(targs->cosalpha,vMult(targs->sourcePos,-1),targs->random);
        double firstPoint = intersect_cylinder_out(targs->sourcePos,direction,pztop,pzbottom,radius);
        if (firstPoint < inf){
            targs->detectorParticles++;
            double energyAbsorbed = traceParticle(transloc(targs->sourcePos,vMult(direction,firstPoint)),direction,targs->sourceEnergy,targs->random);
            targs->sumEnergy += energyAbsorbed;
            if (energyAbsorbed > 0) {
                colledParticles++;
                energyAbsorbed = energyAbsorbed + targs->sigma*drandnt(targs->random);
                if (targs->maxEnergy < energyAbsorbed) energyAbsorbed = targs->maxEnergy;
                (targs->chan)[((int)(energyAbsorbed*(channels-1)/(targs->maxEnergy)))]++;
            }
        }
    }
    targs->particlesTraced += i;
}

int main(){
  //  initRandrand();
    double roh = 3.67;
    initReactions(roh);

    int particleNum = 5e8;
    int update = 1e4;
    double sourceEnergy = 4000;
    double fwhm = 30;
    double sigma = fwhm/2.355;

    double maxenergy = sourceEnergy+sigma*7;
    double energyPerChannel = maxenergy/channels;

    MTRand mainrand = initRandrandt();

    FILE* gp_pipe = popen ("gnuplot -persistent", "w");
 //   fprintf(gp_pipe,"set title 'Spectrum'\n");
    fprintf(gp_pipe,"set tmargin 4\n");
    fprintf(gp_pipe,"set logscale y\n");

    struct vector sourcePos;
    sourcePos.x = 3;
    sourcePos.y = 3;
    sourcePos.z = 3;

    double sinalpha = sqrt(radius*radius+(pztop-pzbottom)*(pztop-pzbottom))/vAbs(sourcePos);
    double cosalpha = sqrt(1-sinalpha*sinalpha);
    double ang = 2*pi*(1-cosalpha);
    double particleMultiplier = 4*pi/ang;

    int *sumTChannels = malloc(sizeof(int)*channels);
  /*  for (int i = 0; i < channels; i++){
        sumTChannels[i] = 0;
    }*/

    int threadCount = 4;

    int **threadChannels = malloc(sizeof(int*)*threadCount);
    struct tracingThreadArgs *targs = malloc(sizeof(struct tracingThreadArgs)*threadCount);
    MTRand *threadRands = malloc(sizeof(MTRand)*threadCount);
    pthread_t *threads = malloc(sizeof(pthread_t)*threadCount);

    struct timeval ct;
    gettimeofday(&ct,NULL);
    long peTi = ct.tv_sec * (int)1e6 + ct.tv_usec;

    for (int i = 0; i < threadCount; i++){
        threadChannels[i] = malloc(sizeof(int)*channels);
        for (int j = 0; j < channels; j++){
            threadChannels[i][j] = 0;
        }
        threadRands[i] = initRandt((long)(drandt(&mainrand)*1e6));
        targs[i].chan = threadChannels[i];
        targs[i].random =  &(threadRands[i]);
        targs[i].stop_at_colls = update;
        targs[i].sourceEnergy = sourceEnergy;
        targs[i].cosalpha = cosalpha;
        targs[i].sigma = sigma;
        targs[i].sourcePos = sourcePos;
        targs[i].particlesTraced = 0;
        targs[i].detectorParticles = 0;
        targs[i].sumEnergy = 0;
        targs[i].maxEnergy = maxenergy;
        pthread_create(&(threads[i]),NULL,tracingThread,&(targs[i]));
    }

    int sumcoll = 0;
    int sumpart = 0;
    int sumdet = 0;
    double sumEnergy = 0;

    for (int i = 0; i < particleNum/update; i++){

        for (int j = 0; j < threadCount; j++){
            pthread_join(threads[j],NULL);
        }

        for (int j = 0; j < threadCount; j++){
            sumpart += targs[j].particlesTraced;
            sumdet += targs[j].detectorParticles;
            sumEnergy += targs[j].sumEnergy;
        }

        for (int j = 0; j < channels; j++){
            sumTChannels[j] = 0;
            for (int k = 0; k < threadCount; k++){
                sumTChannels[j] += threadChannels[k][j];
                sumcoll += threadChannels[k][j];
            }
        }
        for (int j = 0; j < threadCount; j++){
            pthread_create(&(threads[j]),NULL,tracingThread,&(targs[j]));
        }

        //Minden más kód ide jön!!!!

        gettimeofday(&ct,NULL);

  //      printf ("Average speed: %f kp/s\n",sumcoll/1e3/((double)(ct.tv_sec * (int)1e6 + ct.tv_usec - peTi)/1e6));
  //      printf ("Total collisions: %d\n",sumcoll);
        fprintf(gp_pipe,"set label 1 'Average detected particle speed: %.1f kp/s' at screen 0.5,screen 0.99 center\n",sumcoll/1e3/((double)(ct.tv_sec * (int)1e6 + ct.tv_usec - peTi)/1e6));
        fprintf(gp_pipe,"set label 2 'Total number of hits: %.3f Mp' at screen 0.05,screen 0.99\n",sumcoll/1e6);
        fprintf(gp_pipe,"set label 3 'Total particles traced: %.3f Mp' at screen 0.05,screen 0.96 \n",(double) sumpart/1e6);
        fprintf(gp_pipe,"set label 4 'Max efficiency: %.5f' at screen 0.95, screen 0.99 right\n",sumEnergy/(sumdet*sourceEnergy));
        fprintf(gp_pipe,"set label 5 'Total efficiency: %.5f' at screen 0.95, screen 0.96 right\n",sumEnergy/(sumpart*particleMultiplier*sourceEnergy));
        fprintf(gp_pipe,"set label 6 'Average particle speed: %.2f Mp/s' at screen 0.5, screen 0.96 center\n",sumpart/1e6/((double)(ct.tv_sec * (int)1e6 + ct.tv_usec - peTi)/1e6));
        fprintf(gp_pipe,"set label 7 'Average source activity: %.2f MBq' at screen 0.05, screen 0.93 left\n",(sumpart*particleMultiplier)/1e6/((double)(ct.tv_sec * (int)1e6 + ct.tv_usec - peTi)/1e6));
       // fprintf(gp_pipe,"set label 4 'Total number of hits: %.3f M' at screen 0.95,screen 0.98 right\n",sumcoll/1e6);
//        fprintf(gp_pipe,"set logscale y\n");

        plot(gp_pipe,sumTChannels,channels);

        sumcoll = 0;
        sumpart = 0;
        sumdet = 0;
        sumEnergy = 0;

    }


    //gettimeofday(&ct,NULL);
    //printf ("Average speed: %f Mp/s\n",colledParticles/1e6/((double)(ct.tv_sec * (int)1e6 + ct.tv_usec - peTi)/1e6));

 /*   for (int i = 0; i < channels; i++){
        printf("%d\n",ch[i]);
    }*/
   // plot(gp_pipe);
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
