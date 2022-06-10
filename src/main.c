#include "transport.h"
#include "tmath.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "reactions.h"
#include "wyhash.h"
#include <pthread.h>

#define inf 1e9
#define pi 3.14159
#define channels 1024
#define devw 25

double pztop = 2;
double pzbottom = -2;
double radius = 2;

struct tracingThreadArgs{
    uint64_t random[4];
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


double traceParticle(struct vector pos, struct vector direction, double particle_energy, uint64_t *random, double *crs){ //returns: absorbed energy by the crystal
    int particle_is_alive = 1;
    double sumen = 0;
    while (particle_is_alive == 1){
        double lastPoint = intersect_cylinder_in(pos,direction,pztop,pzbottom,radius);
        getCrs(particle_energy,crs);
        double freeway = shuffle_freeway_length(random,crs);
        if (freeway < lastPoint){
            pos = transloc(pos,vMult(direction,freeway));
            int reac = shuffle_reaction(random,crs);
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
                eenergy += traceParticle(pos,annihil1,511,random,crs);
                eenergy += traceParticle(pos,annihil2,511,random,crs);

            }
            sumen += eenergy;
        }
        else {
            particle_is_alive = 0;
        }
    }

    return sumen;
}

void plot(FILE *gp_pipe,int *chan, int chansize,double energyPerChannel){
    fprintf(gp_pipe,"plot '-'\n");
    for (int i = 0; i < chansize; i++){
        fprintf(gp_pipe,"%f %d\n",i*energyPerChannel,chan[i]);
    }
    fprintf(gp_pipe,"e\n");
}

void *tracingThread(void* arg){
    struct tracingThreadArgs *targs = (struct tracingThreadArgs*) arg;

    double *crs = malloc(getX()*sizeof(double));

    int colledParticles = 0;
    int i = 0;

    for (; colledParticles < targs->stop_at_colls; i++){
        //struct vector direction = isotropicDirection();
        struct vector direction = coneDirection(targs->cosalpha,vMult(targs->sourcePos,-1),targs->random);
        double firstPoint = intersect_cylinder_out(targs->sourcePos,direction,pztop,pzbottom,radius);
        if (firstPoint < inf){
            targs->detectorParticles++;
            double energyAbsorbed = traceParticle(transloc(targs->sourcePos,vMult(direction,firstPoint)),direction,targs->sourceEnergy,targs->random,crs);
            targs->sumEnergy += energyAbsorbed;
            if (energyAbsorbed > 0) {
                colledParticles++;
                energyAbsorbed = energyAbsorbed + targs->sigma*drandnt(targs->random);
                if (energyAbsorbed < 0) energyAbsorbed = 0;
                if (targs->maxEnergy < energyAbsorbed) energyAbsorbed = targs->maxEnergy;
                (targs->chan)[((int)(energyAbsorbed*(channels-1)/(targs->maxEnergy)))]++;
            }
        }
    }
    targs->particlesTraced += i;
    free(crs);

    return 0;
}

int main(){

    double roh = 3.67;
    initReactions(roh);

    int threadCount = 20;
    unsigned long particleNum = 1e7;//5e8;
    int update = 1e4;
    double sourceEnergy = 4000;
    double fwhm = 30;
    double sigma = fwhm/2.355;

    double maxenergy = sourceEnergy+sigma*7;
    double energyPerChannel = maxenergy/channels;

    //uint64_t mainrand = initRandrandt();
    uint64_t mainrand[4];
    make_secret(time(NULL),mainrand);

    FILE* gp_pipe = popen ("gnuplot -persistent", "w");

    fprintf(gp_pipe,"set tmargin 4\n");
    fprintf(gp_pipe,"set logscale y\n");
    fprintf(gp_pipe,"set ylabel 'Beütésszám'\n");
    fprintf(gp_pipe,"set xlabel 'Energia (keV)'\n");

    struct vector sourcePos;
    sourcePos.x = 3;
    sourcePos.y = 3;
    sourcePos.z = 3;

    double sinalpha = sqrt(radius*radius+(pztop-pzbottom)*(pztop-pzbottom))/vAbs(sourcePos);
    double cosalpha = sqrt(1-sinalpha*sinalpha);
    double ang = 2*pi*(1-cosalpha);
    double particleMultiplier = 4*pi/ang;

    int *sumTChannels = malloc(sizeof(int)*channels);

    int **threadChannels = malloc(sizeof(int*)*threadCount);
    struct tracingThreadArgs *targs = malloc(sizeof(struct tracingThreadArgs)*threadCount);
    //uint64_t **threadRands = malloc(sizeof(uint64_t)*threadCount);
    pthread_t *threads = malloc(sizeof(pthread_t)*threadCount);

    struct timeval ct;
    gettimeofday(&ct,NULL);
    long peTi = ct.tv_sec * (int)1e6 + ct.tv_usec;

    for (int i = 0; i < threadCount; i++){
        threadChannels[i] = malloc(sizeof(int)*channels);
        for (int j = 0; j < channels; j++){
            threadChannels[i][j] = 0;
        }
        //threadRands[i] = initRandt((long)(drandt(&mainrand)*1e6));
        make_secret(wyrand(mainrand),targs->random);
        targs[i].chan = threadChannels[i];
      //  targs[i].random =  &(threadRands[i]);
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

    double teffs[devw] = {0};
    double meffs[devw] = {0};
    int effsp = 0;
    int effsf = 0;

    for (int i = 0; i < particleNum/(update*threadCount); i++){

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

        //Everything else comes here!!!! The threads are already running, so a little heavier code won't have a (big) impact on the simulation

        gettimeofday(&ct,NULL);

        fprintf(gp_pipe,"set label 1 'Average detected particle speed: %.1f kp/s' at screen 0.5,screen 0.99 center\n",sumcoll/1e3/((double)(ct.tv_sec * (int)1e6 + ct.tv_usec - peTi)/1e6));
        fprintf(gp_pipe,"set label 2 'Total number of hits: %.3f Mp' at screen 0.05,screen 0.99\n",sumcoll/1e6);
        fprintf(gp_pipe,"set label 3 'Total particles traced: %.3f Mp' at screen 0.05,screen 0.96 \n",(double) sumpart/1e6);
        double meff = sumEnergy/(sumdet*sourceEnergy);
        double teff = sumEnergy/(sumpart*particleMultiplier*sourceEnergy);
        meffs[effsp] = meff;
        teffs[effsp] = teff;
        effsp++;
        if (effsp > devw-1) {
            effsp = 0;
            effsf = 1;
        }
        int ec = 0;
        if (effsf == 0){
            ec = effsp;
        } else {
            ec = devw;
        }
        double avgm = 0;
        double avgt = 0;

        for (int k = 0; k < ec; k++){
            avgm += meffs[k];
            avgt += teffs[k];
        }
        avgm = avgm/ec;
        avgt = avgt/ec;

        double devm = 0;
        double devt = 0;

        for (int k = 0; k < ec; k++){
            devm += (meffs[k] - avgm)*(meffs[k] - avgm);
            devt += (teffs[k] - avgt)*(teffs[k] - avgt);
        }

        devm = sqrt(devm/ec);
        devt = sqrt(devt/ec);

        fprintf(gp_pipe,"set label 4 'Max efficiency: %.5f pm %.7f' at screen 0.95, screen 0.99 right\n",meff,devm);
        fprintf(gp_pipe,"set label 5 'Total efficiency: %.5f pm %.7f' at screen 0.95, screen 0.96 right\n",teff,devt);
        fprintf(gp_pipe,"set label 6 'Average particle speed: %.2f Mp/s' at screen 0.5, screen 0.96 center\n",sumpart/1e6/((double)(ct.tv_sec * (int)1e6 + ct.tv_usec - peTi)/1e6));
        fprintf(gp_pipe,"set label 7 'Average source activity: %.2f MBq' at screen 0.05, screen 0.93 left\n",(sumpart*particleMultiplier)/1e6/((double)(ct.tv_sec * (int)1e6 + ct.tv_usec - peTi)/1e6));

        plot(gp_pipe,sumTChannels,channels,energyPerChannel);

        sumcoll = 0;
        sumpart = 0;
        sumdet = 0;
        sumEnergy = 0;

    }

   // freeReactions();
    return 0;
}
