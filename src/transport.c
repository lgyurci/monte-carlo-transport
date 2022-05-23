#include <stdio.h>
#include <stdlib.h>
#include "tmath.h"
#include <math.h>
#include <time.h>
#define pi 3.14159
#define inf 1e9

double drand(){
    return ((double) rand())/RAND_MAX;
}

struct vector isotropicDirection(){
    struct vector direction;
    direction.z = drand()*2-1;
    double fi = 2*pi*drand();
    double r = sqrt(1-direction.z*direction.z);
    direction.x = r*cos(fi);
    direction.y = r*sin(fi);
    return direction;
}

double intersect_plane(struct vector position, struct vector direction, double planeZ){
    double dist = 0;
    if (vAbs(direction) < 1e-3) {
        dist = inf;
    } else {
        dist = (planeZ - position.z)/direction.z;
    }
    return dist;
}

struct vector intersect_cylinder(struct vector position, struct vector direction, double radius){
    double a = direction.x*direction.x + direction.y*direction.y;
    double b = 2*(position.x*direction.x+position.y*direction.y);
    double c = position.x*position.x + position.y*position.y-radius*radius;
    double d = b*b-4*a*c;
    struct vector dist;
    dist.x = 0;
    dist.y = 0;
    if (d < 0){
        dist.x = inf;
        dist.y = inf;
    }
    else {
        dist.x = (-1*b+sqrt(d))/(2*a);
        dist.y = (-1*b-sqrt(d))/(2*a);
    }
//    printf("ic %f %f\n",dist.x,dist.y);
    return dist;
}

double intersect_cylinder_in(struct vector position,struct vector direction,double pztop,double pzbottom,double radius){
    struct vector d = intersect_cylinder(position,direction,radius);
    double dp = intersect_plane(position,direction,pztop);
    double dn = intersect_plane(position,direction,pzbottom);
    double fidesz = 0;
    double dk = 0;
    if (d.x > d.y){
        fidesz = d.x;
    } else {
        fidesz = d.y;
    }

    if (dp > dn){
        dk = dp;
    } else {
        dk = dn;
    }

    if (fidesz < dk){
        return fidesz;
    }

    return dk;
}

double intersect_cylinder_out(struct vector position,struct vector direction,double pztop,double pzbottom,double radius){
    struct vector de = intersect_cylinder(position,direction,radius);
    double d1 = de.x;
    double d2 = de.y;
    double dmin = 0;
    double dmax = 0;
    double dist = 0;
    if (d1 >= d2) {
        dmax = d1;
        dmin = d2;
    } else {
        dmax = d2;
        dmin = d1;
    }

    if (dmin >= inf || dmax < 0 || (position.z > pztop && direction.z > 0) || (position.z < pzbottom && direction.z < 0)){
       // printf("asd %f %f\n",de.x,de.y);
        dist = inf;
    } else {
        double dcyl = 0;

        if (d1*d2 > 0){
            dcyl = dmin;
        } else {
            dcyl = dmax;
        }

        double zcyl = position.z + dcyl*direction.z;

        if (zcyl < pzbottom || zcyl > pztop){
            dcyl = inf;
        }

        double dp = intersect_plane(position,direction,pztop);
        double xp = position.x + dp*direction.x;
        double yp = position.y + dp*direction.y;

        if (xp*xp + yp*yp > radius*radius){
            dp = inf;
        }

        double dn = intersect_plane(position,direction,pzbottom);
        double xn = position.x + dn*direction.x;
        double yn = position.y + dn*direction.y;

        if (xn*xn+yn*yn > radius*radius){
            dn = inf;
        }

        if (dn < dp && dn < dcyl) return dn;
        if (dp < dn && dp < dcyl) return dp;
        return dcyl;
    }
    return dist;
}

void initRandrand(){
    time_t t;
    srand((unsigned) time(&t));
}

/*int main(){
    initRandrand();
    for (int i = 0; i < 10000; i++){
        struct vector asd = isotropicDirection();
        struct vector pos;
        pos.x = 0;
        pos.y = 0;
        pos.z = 0;
        double asd2 = intersect_cylinder_in(pos,asd,2,-2,1);
        if (asd2 < inf) {
            asd = vMult(asd,asd2);
            asd = transloc(pos,asd);
            printf("%f\t%f\t%f\n",asd.x,asd.y,asd.z);
        }

      //  printf("%f\n",asd2);
    }
    return 0;
}*/
