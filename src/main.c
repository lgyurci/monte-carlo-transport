#include "transport.h"
#include "tmath.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#define inf 1e9
#define channels 1024

int main(){
    initRandrand();
    int ch[channels] = { 0 };
    int particleNum = 1e5;

    struct vector sourcePos;
    sourcePos.x = 3;
    sourcePos.y = 3;
    sourcePos.z = 3;

    double pztop = 2;
    double pzbottom = -2;
    double radius = 2;

    double sourceEnergy = 500;

    for (int i = 0; i < particleNum; i++){
        struct vector direction = isotropicDirection();
        double firstPoint = intersect_cylinder_out(sourcePos,direction,pztop,pzbottom,radius);
        if (firstPoint < inf){
            struct vector pos = transloc(sourcePos,vMult(direction,firstPoint));
            double lastPoint = intersect_cylinder_in(pos,direction,pztop,pzbottom,radius);
            pos = transloc(pos,vMult(direction,lastPoint));


        }
    }
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
