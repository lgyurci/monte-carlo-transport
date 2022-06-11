#ifndef random_H_
#define random_H_
#include "wyhash.h"


static inline double drandt(uint64_t *random){
    return wy2u01(wyrand(random));
}

static inline double drandnt(uint64_t *random){
    return wy2gau(wyrand(random));
}

#endif
