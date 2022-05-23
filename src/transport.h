#ifndef transport_H_
#define transport_H_
#include <stdlib.h>
#include "tmath.h"
#include <math.h>
#include <time.h>

double drand();
struct vector isotropicDirection();
double intersect_plane(struct vector position, struct vector direction, double planeZ);
struct vector intersect_cylinder(struct vector position, struct vector direction, double radius);
double intersect_cylinder_in(struct vector position,struct vector direction,double pztop,double pzbottom,double radius);
double intersect_cylinder_out(struct vector position,struct vector direction,double pztop,double pzbottom,double radius);
void initRandrand();

#endif
