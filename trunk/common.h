#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

#include <frei0r.hpp>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <cmath>

#define TRI_IND(y,x) (x+(y+y*y)/2)
#define TRI_IND1(y,x) (x-1+(y*y-y)/2)
#define CLIP0(x) (x < 0 ? 0 : x )
#define CLIP(x,y) (x > y ? y : x )

#define POW_COMP(x,t,s) ( x <= t ? 0.5*pow(x/t,s) : 1.0-0.5*pow((1-x)/(1-t),s*(1-t)/t) )

// BGRA in memory (= big-endian), equals ARGB in bit-order on little-endian (=cairo on x86)
typedef struct col32bit_s {
    uint8_t b;
    uint8_t g;
    uint8_t r;
    uint8_t a;
} col32bit;

typedef struct col128bit_s {
    uint32_t b;
    uint32_t g;
    uint32_t r;
    uint32_t a;
} col128bit;

typedef struct {
  int w;
  int h;
  int stride;
  int bpp;
  int size;
} ScreenGeometry;

int static inline fCLIP(int x,int y){
    return (x > y ? y : x );
}

int static inline fCLIP0(int x){
    return (x < 0 ? 0 : x );
}

uint8_t static inline fCLIP8(uint32_t x){
    return ( x>255 ) ? 255 : (uint8_t)(x);
}

#endif // COMMON_H_INCLUDED
