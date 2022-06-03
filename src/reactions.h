#ifndef reactions_H_
#define reactions_H_
#include <stdio.h>
#include <stdlib.h>
#include "transport.h"

void initReactions();
double getCrossection(double energy,int effect); //effect: 1 - compton scattering, 2 - photoelectric effect, 3 - pair production in nucleus, 4 - pair production in electron bands
double shuffle_freeway_length(double energy);
int shuffle_reaction(double energy);
void freeReactions();
struct vector comptonScatter(struct vector direction,double *energy);


#endif
