#include <math.h>
#include <stdio.h>
#include "tmath.h"

struct vector crossProduct(struct vector v1, struct vector v2){ //vektoriális szorzat (egyéb indoklás nélkül)
    struct vector v;
    v.x = v1.y*v2.z-v1.z*v2.y;
    v.y = v1.z*v2.x-v1.x*v2.z;
    v.z = v1.x*v2.y-v1.y*v2.x;
    return v;
}
double vAbs(struct vector v){ //egy vektor hossza
    return sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
}
double vAngle(struct vector v1, struct vector v2){ //két vektor által bezárt szög
    double invangle = (v1.x*v2.x+v1.y*v2.y+v1.z*v2.z)/(vAbs(v1)*vAbs(v2)); //skaláris szorzatból számolva
    if (invangle > 1) invangle = 1; //néha a double típus határai miatt itt jelent meg 1-nél nagyobb vagy kisebb érték (matematikailag nem lehetne)
    if (invangle < -1) invangle = -1;
    return acos(invangle);
}


struct vector vAdd(struct vector v1, struct vector v2){ //két vektor összege
    struct vector v;
    v.x = v1.x+v2.x;
    v.y = v1.y+v2.y;
    v.z = v1.z+v2.z;
    return v;
}
struct vector vMult(struct vector v1, double lambda){ //vektor szorzása skalárral
    struct vector v;
    v.x = v1.x*lambda;
    v.y = v1.y*lambda;
    v.z = v1.z*lambda;
    return v;
}
struct vector reduceTo(struct vector v, double targetLength){ //egy vektor hosszának beállítása az adott értékre, irányának megtartásával
    double lambda = sqrt(targetLength*targetLength/(v.x*v.x+v.y*v.y+v.z*v.z));
    v.x = v.x*lambda;
    v.y = v.y*lambda;
    v.z = v.z*lambda;
    if (targetLength < 0) v = vMult(v,-1);
    return v;
}

struct vector transloc(struct vector p, struct vector v){ //pont eltolása vektorral
    struct vector mp;
    mp.x = p.x+v.x;
    mp.y = p.y+v.y;
    mp.z = p.z+v.z;
    return mp;
}
struct vector pToV(struct vector p1,struct vector p2){ //két pont által meghatározott vektor
    //p2-be mutat
    struct vector v;
    v.x = p2.x-p1.x;
    v.y = p2.y-p1.y;
    v.z = p2.z-p1.z;
    return v;
}
double pDist(struct vector p1, struct vector p2){ //két pont távolsága
    return sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)+(p1.z-p2.z)*(p1.z-p2.z));
}

double dotProduct(struct vector v1, struct vector v2){ //skaláris szorzat
    return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
}

void printVector(struct vector v){
    printf("%f\t%f\t%f\n",v.x,v.y,v.z);
}
