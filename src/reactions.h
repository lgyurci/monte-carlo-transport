#ifndef reactions_H_
#define reactions_H_
#include <stdio.h>
#include <stdlib.h>
#include "transport.h"
#include "wyhash.h"

void initReactions(double roo,double source,char *xcomfile);
double getCrossection(double energy,int effect); //effect: 1 - compton scattering, 2 - photoelectric effect, 3 - pair production in nucleus, 4 - pair production in electron bands
double shuffle_freeway_length(uint64_t *random,double *crs);
int shuffle_reaction(uint64_t *random,double *crs);
void getCrs(double energy, double *crs);
void freeReactions();
double *allocCrossections();
int getX();
struct vector comptonScatter(struct vector direction,double *energy,uint64_t *random);


#endif
