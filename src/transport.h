#ifndef transport_H_
#define transport_H_
#include <stdlib.h>
#include "tmath.h"
#include <math.h>
#include <time.h>
#include "random.h"

double drandt(uint64_t *random);
struct vector isotropicDirection(uint64_t *random);
double intersect_plane(struct vector position, struct vector direction, double planeZ);
struct vector intersect_cylinder(struct vector position, struct vector direction, double radius);
double intersect_cylinder_in(struct vector position,struct vector direction,double pztop,double pzbottom,double radius);
double intersect_cylinder_out(struct vector position,struct vector direction,double pztop,double pzbottom,double radius);
//MTRand initRandrandt();
//MTRand initRandt(long seed);
struct vector isotropicScatter(struct vector direction,double angle,uint64_t *random);
struct vector coneDirection(double cosalpha,struct vector direction,uint64_t *random);
double drandnt(uint64_t *random);

#endif
