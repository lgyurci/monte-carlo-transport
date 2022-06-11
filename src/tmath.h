#ifndef tmath_H_
#define tmath_H_
#include <math.h>

struct vector{
    double x,y,z;
};

double dotProduct(struct vector v1, struct vector v2); //skaláris szorzat
double pDist(struct vector p1, struct vector p2); //két pont távolsága
double vAbs(struct vector v); //vektor abszolút értéke
struct vector crossProduct(struct vector v1, struct vector v2); //vektoriális szorzat
double vAngle(struct vector v1, struct vector v2); //két vektor által bezárt szög
struct vector vAdd(struct vector v1, struct vector v2); //két vektor összege
struct vector reduceTo(struct vector v, double targetLength); //vektor hosszát egy adott értékre állítja, az írányát nem változtatja
struct vector vMult(struct vector v1, double lambda); //vektor szorzása skalárral
struct vector transloc(struct vector p, struct vector v); //pont eltolása egy adott vektorral
struct vector pToV(struct vector p1,struct vector p2); ///két pontból csinál vektort
void printVector(struct vector v); //egy vektor kiírása a standard outra (csak debug funkció)

#endif
