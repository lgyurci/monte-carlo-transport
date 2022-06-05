#ifndef transport_H_
#define transport_H_
#include <stdlib.h>
#include "tmath.h"
#include "mtwister.h"
#include <math.h>
#include <time.h>

double drandt(MTRand *random);
struct vector isotropicDirection(MTRand *random);
double intersect_plane(struct vector position, struct vector direction, double planeZ);
struct vector intersect_cylinder(struct vector position, struct vector direction, double radius);
double intersect_cylinder_in(struct vector position,struct vector direction,double pztop,double pzbottom,double radius);
double intersect_cylinder_out(struct vector position,struct vector direction,double pztop,double pzbottom,double radius);
MTRand initRandrandt();
MTRand initRandt(long seed);
struct vector isotropicScatter(struct vector direction,double angle,MTRand *random);
struct vector coneDirection(double cosalpha,struct vector direction,MTRand *random);
double drandnt(MTRand *random);

#endif
