# Introduction #

Due to the non-linear nature of the RGB to HSV (HSI) color space conversion, it's a well known fact that most conversion routines are utilizing floating-point arithmetic and store HSV values in floating point numbers, respectively. Furthermore, it's impossible to store full color infomation of RGB color space with a given numeric precicion in HSV space with identical precision. The non-linear transformation requires higher precision in HSV space to get a lossless, neutral RGB->HSV->RGB conversion chain.

For this reason, a purely integer based implementation was developed, which allows to add a variable number of bits on the HSV side compared to the 3x8-bit resolution of the RGB component values.

# Details #

The developed conversion routines, which are part of the color\_compress0r _frei0r_ effect plugin, utilize a quadruple of 4x8-bit integer (struct of 4 times uint8) to store RGBA pixel information, and a quadruple of 4x16-bit integer (struct if 4 time uint16) to store HSVA values. The conversion itself relies on 32-bit integer arithmetic. Effective numeric precision for HSV values can be adjusted at complie time.

Sweeping the HSV numeric precision showed that 4 additional bits are sufficient to store 8-bit RGB colorspace information without any loss. Going RGB(8-bit)->HSV(12-bit)->RGB(8-bit) builds a neutral transformation. The following table lists the remaining bit error versus  HSV precision and clearly shows that 8-bit HSV precision is way too low.

| **HSV resolution** | **0-bit error** | **1-bit error** | **2-bit error** | **>=3-bit error**| **max. error = sum()**|
|:-------------------|:----------------|:----------------|:----------------|:-----------------|:|
| 12-bit | 100% |  |  |  | 0 |
| 11-bit | 98.6% | 1.4% |  |  | 1 |
| 10-bit | 90% | 10% |  |  | 2 |
| 9-bit  | 80% | 16% | 4% |  | 3 |
| 8-bit  | 60% | 12% | 8% | 20% | 7 |

# C/C++ Code #

```
#define ABITS 4
#define HSCALE 256

typedef struct col32bit {
    union {
        uint8_t r;
        uint8_t h;
    };
    union {
        uint8_t g;
        uint8_t s;
    };
    union {
        uint8_t b;
        uint8_t v;
    };
    uint8_t a;
} col32bit_t;

typedef struct col64bit {
    union {
        uint16_t r;
        uint16_t h;
    };
    union {
        uint16_t g;
        uint16_t s;
    };
    union {
        uint16_t b;
        uint16_t v;
    };
    uint16_t a;
} col64bit_t;

inline col64bit rgba2hsva(const col32bit c){
  col64bit r;
  uint32_t iMin,iMax,chroma;
  const uint32_t k1=255 << ABITS;
  const uint32_t k2=HSCALE << ABITS;
  
  // there's no need to touch alpha, it remains 8-bit
  r.a=c.a;

  if (c.r > c.g) {
      iMax = max (c.r, c.b);
      iMin = min (c.g, c.b);
  } else {
      iMax = max (c.g, c.b);
      iMin = min (c.r, c.b);
  }

  chroma = iMax - iMin;
  // set value
  r.v = iMax << ABITS;

  // set saturation
  if (r.v == 0)
    r.s = 0;
  else
    r.s = (k1*chroma)/iMax;

  // set hue 
  if (r.s == 0)
      r.h = 0;
  else {
      if ( c.r == iMax ) {
        r.h  =  (k2*(6*chroma+c.g - c.b))/(6*chroma);
        if (r.h >= k2) r.h -= k2;
      } else if (c.g  == iMax)
        r.h  =  (k2*(2*chroma+c.b - c.r )) / (6*chroma);
      else // (c.b == iMax )
        r.h  =  (k2*(4*chroma+c.r - c.g )) / (6*chroma);
  }
  return r;
}

inline col32bit hsva2rgba(const col64bit c){
  col32bit r;
  uint32_t m;
  int32_t H,X,ih,is,iv;
  const uint32_t k1=255 << ABITS;
  const int32_t k2=HSCALE << ABITS;
  const int32_t k3=1<<(ABITS-1);
  r.a=c.a;

  // set chroma and min component value m
  //chroma = ( c.v * c.s )/k1;
  //m = c.v - chroma;
  m = ((uint32_t)c.v*(k1 - (uint32_t )c.s ))/k1;

  // chroma  == 0 <-> c.s == 0 --> m=c.v
  if (c.s == 0) {
      r.b = ( r.g = ( r.r = c.v >> ABITS ));
  } else {
    ih=(int32_t)c.h;
    is=(int32_t)c.s;
    iv=(int32_t)c.v;

    H = (6*ih)/k2;
    X = ((iv*is)/k2)*(k2-abs(6*ih- 2*(H>>1)*k2 - k2)) ;

    // removing additional bits --> unit8
    X=((X+iv*(k1 - is))/k1 + k3) >>ABITS;
    m=m >> ABITS;

    // ( chroma + m ) --> c.v ;
    switch (H) {
        case 0:
          r.r = c.v >> ABITS ;
          r.g = X;
          r.b = m ;
          break;
        case 1:
          r.r = X;
          r.g = c.v >> ABITS;
          r.b = m ;
          break;
        case 2:
          r.r = m ;
          r.g = c.v >> ABITS;
          r.b = X;
          break;
        case 3:
          r.r = m ;
          r.g = X;
          r.b = c.v >> ABITS;
          break;
        case 4:
          r.r = X;
          r.g = m ;
          r.b = c.v >> ABITS;
          break;
        case 5:
          r.r = c.v >> ABITS;
          r.g = m ;
          r.b = X;
          break;
    }
  }
  return r;
}




```