#ifndef FREI0R_COMMON_H_INCLUDED
#define FREI0R_COMMON_H_INCLUDED

#include <algorithm>
#include <cmath>
#include <inttypes.h>
#include <cstdlib>

#define TRI_IND(y,x) (x+(y+y*y)/2)
#define TRI_IND1(y,x) (x-1+(y*y-y)/2)
#define CLIP0(x) (x < 0 ? 0 : x )
#define CLIP(x,y) (x > y ? y : x )

#define CDIST_RGB 0
#define CDIST_L1 0
#define CDIST_L2 1
#define CDIST_HUE 2
#define CDIST_VWEIGHT 4
#define CDIST_SWEIGHT 8

#define ABITS 4
#define VSCALE 255
#define HSCALE 256
#define AVSCALE (VSCALE << ABITS)
#define AHSCALE (HSCALE << ABITS)

#define POW_COMP(x,tx,ty,s) ( x <= tx ? ty*pow(x/tx,s*tx/ty) : 1.0-(1-ty)*pow((1-x)/(1-tx),s*(1-tx)/(1-ty)))

typedef bool f0r_param_boolcxx;

// BGRA in memory (= big-endian), equals ARGB in bit-order on little-endian (=cairo on x86)
typedef struct col32bitBGRA {
    union {
        uint8_t b;
        uint8_t v;
    };
    union {
        uint8_t g;
        uint8_t s;
    };
    union {
        uint8_t r;
        uint8_t h;
    };
    uint8_t a;
} col32bitBGRA_t;

// BGRA in memory (= big-endian), equals ARGB in bit-order on little-endian (=cairo on x86)
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

typedef struct col128bit {
    union {
        uint32_t r;
        uint32_t h;
    };
    union {
        uint32_t g;
        uint32_t s;
    };
    union {
        uint32_t b;
        uint32_t v;
    };
    uint32_t a;
} col128bit_t;

typedef struct colfloat {
    union {
        float r;
        float h;
    };
    union {
        float g;
        float s;
    };
    union {
        float b;
        float v;
    };
    uint8_t a;
} colfloat_t;

typedef struct ScreenGeometry {
  int w;
  int h;
  int stride;
  int bpp;
  int size;
} ScreenGeometry_t;

inline int fCLIP(int x,int y){
    return (x > y ? y : x );
}

inline int fCLIP0(int x){
    return (x < 0 ? 0 : x );
}

inline uint8_t fCLIP8(uint32_t x){
    return ( x>255 ) ? 255 : (uint8_t)(x);
}

inline col32bit inline_rgba2hsva(const col32bit c){

  col32bit r;
  int iMin,delta,delta360;
  r.a=c.a;

  if (c.r > c.g) {
      r.b = std::max (c.r, c.b);
      iMin = std::min (c.g, c.b);
  } else {
      r.b = std::max (c.g, c.b);
      iMin = std::min (c.r, c.b);
  }

  delta = r.b - iMin;

  if (r.b == 0)
    r.g = 0;
  else
    r.g = (uint8_t)((255*delta)/r.b);

  if (r.g == 0)
      r.r = 0;
  else {
      int h;
      delta360=delta*360;
      if (c.r == r.b) {
        //h = (60 * ((int)c.g - (int)c.b)) / delta;
        //h = (60 * (c.g - c.b)) / delta;
        h = (360*delta +  60 * (c.g - c.b));
        //if (h < 0) h += 360;
        if (h >= delta360) h -= delta360;
      } else if (c.g == r.b)
        h = (120*delta + 60 * (c.b - c.r));
      else
        h = (240*delta + 60 * (c.r - c.g));

      r.r = (uint8_t) ((HSCALE*h)/delta360);
  }

  return r;
}

inline col32bit inline_hsva2rgba(const col32bit c){
  col32bit r;
  int i,f,s;
  uint8_t m,iq,it;
// R G B
// H S V

  r.a=c.a;

  if (c.g == 0) {
      r.r=c.b;
      r.g=c.b;
      r.b=c.b;
  } else {
    m = (uint8_t) (((255 - c.g) * c.b )/255  ) ; // =m

    i=(6*c.r)/HSCALE;
    f = ( 6 * c.r - HSCALE*i);

    s=HSCALE*255;
    iq = (uint8_t) (( (s - c.g *f        ) * c.b )/s) ;
    it = (uint8_t) (( (s - c.g *(255-f)  ) * c.b )/s) ;

    switch (i) {
        case 0:
          r.r = c.b ;
          r.g = it;
          r.b = m;
          break;

        case 1:
          r.r = iq;
          r.g = c.b;
          r.b = m;
          break;

        case 2:
          r.r = m;
          r.g = c.b ;
          r.b = it;
          break;

        case 3:
          r.r = m;
          r.g = iq;
          r.b = c.b;
          break;

        case 4:
          r.r = it;
          r.g = m;
          r.b = c.b;
          break;

        case 5:
          r.r = c.b;
          r.g = m;
          r.b = iq;
          break;
    }
  }
  return r;
}

inline col64bit inline_rgba2hsva64(const col32bit c){

  col64bit r;
  uint32_t iMin,iMax,chroma,k1,k2;
  r.a=c.a;

  k1=AVSCALE;
  k2=AHSCALE;

  if (c.r > c.g) {
      iMax = std::max (c.r, c.b);
      //r.b = std::max (c.r, c.b) << ABITS ;
      iMin = std::min (c.g, c.b);
  } else {
      iMax = std::max (c.g, c.b);
      //r.b = std::max (c.g, c.b) << ABITS;
      iMin = std::min (c.r, c.b);
  }

  //chroma = r.b - iMin;
  chroma = iMax - iMin;
  // set value V --> B
  r.b = iMax << ABITS;

  // set saturation S --> G
  if (r.b == 0)
    r.g = 0;
  else
    r.g = (k1*chroma)/iMax;

  // set Hue H --> R
  if (r.g == 0)
      r.r = 0;
  else {
      if ( c.r == iMax ) {
        // h= 60 * ((c.g-c.b)/chroma mod 6 )
        //r.r= (((k1*(c.g-c.b))/chroma) % (k1*6) )/6;
        // r.r /k1 =   ( 60*(c.g - c.b) / chroma) % 360 ) / 360 ;
        r.r  =  (k2*(6*chroma+c.g - c.b))/(6*chroma)    ;
        if (r.r >= k2) r.r -= k2;
      } else if (c.g  == iMax)
        // h= 60 * ((c.b-c.r)/chroma + 2 )
        // r.r /k1 =   ( 60*(c.b - c.r) / chroma) + 120 ) / 360 ;
        // r.r /k1 =   ( (c.b - c.r) / chroma) + 2 ) / 6 ;
        r.r  =  (k2*(2*chroma+c.b - c.r )) / (6*chroma);
      else // (c.b == iMax )
        // h= 60 * ((c.r-c.g)/chroma + 4  )
        r.r  =  (k2*(4*chroma+c.r - c.g )) / (6*chroma);
  }
  return r;
}

inline col32bit inline_hsva2rgba64(const col64bit c){
  col32bit r;

  //uint32_t chroma;
  uint32_t m,k1;
  int H,X,k2,k3,ir,ig,ib;
  // R G B
  // H S V
  r.a=c.a;

  k1=AVSCALE;
  k2=AHSCALE;
  k3=1<<(ABITS-1);

  // set chroma
  //chroma = ( c.b * c.g )/k1;
  //m = c.b - chroma;
  m = ((uint32_t)c.b*(k1 - (uint32_t )c.g ))/k1;

  // chroma  == 0 <-> c.g == 0 --> m=c.b
  if (c.g == 0) {
      r.b = ( r.g = ( r.r = c.b >> ABITS ));
  } else {
    ir=(int)c.r;
    ig=(int)c.g;
    ib=(int)c.b;

    H = (6*ir)/k2;
    //X = chroma*(1-abs(X' -1.0  )) ;
    //X' = (6*c.r/k2)- 2*H
    //X' = (6*c.r- 2*H*k2)/k2;
    //X = chroma*(1-abs(   (6*c.r- 2*H*k2)/k2 -k2/k2  )) ;
    //X = chroma*(1-abs((6*c.r- 2*H*k2- k2)/k2 )) ;
    //X = k1*chroma*(k2-abs(6*((int)c.r)- 2*H*k2 - k2))/(k1*k2) ;

    X = ((ib*ig)/k2)*(k2-abs(6*ir- 2*(H>>1)*k2 - k2)) ;

    // removing additional bits --> unit8
    X=( (X+ib*(k1 - ig ))/k1 + k3 ) >>ABITS;
    m=m >> ABITS;

    // ( chroma + m ) --> c.b ;
    switch (H) {
        case 0:
          r.r = c.b >> ABITS ;
          r.g = X;
          r.b = m ;
          break;

        case 1:
          r.r = X;
          r.g = c.b >> ABITS;
          r.b = m ;
          break;

        case 2:
          r.r = m ;
          r.g = c.b >> ABITS;
          r.b = X;
          break;

        case 3:
          r.r = m ;
          r.g = X;
          r.b = c.b >> ABITS;
          break;

        case 4:
          r.r = X;
          r.g = m ;
          r.b = c.b >> ABITS;
          break;

        case 5:
          r.r = c.b >> ABITS;
          r.g = m ;
          r.b = X;
          break;
    }
  }

  return r;
}

inline int hsdist(uint16_t h,uint16_t s){
    int r;
    // f(sf,hf)=0.5*(1-sf)+sf*hf
    // = 0.5*(1-s/AVSCALE)+s*h/(AHSCALE*AVSCALE)
    // AVSCALE*f= 0.5*(AVSCALE-s)+s*h/AHSCALE
    const int inv_s_weight=2;
    s=AVSCALE-s;
    r=(int)s*h;
    // = (AVSCALE-s)/inv_s_weight +r/AHSCALE
    return r/AHSCALE+(AVSCALE-s)/inv_s_weight;
}

inline uint8_t inline_rgba2cdist(const col32bit c1,const col32bit c2, int mode){
    int db,dg,dr;
    int dh;
    col64bit c1t,c2t;

    const int inv_v_weight=2;
    const int k1= (1+inv_v_weight);

    if (mode & CDIST_HUE) {
        c1t=inline_rgba2hsva64(c1);
        c2t=inline_rgba2hsva64(c2);
        // h360 = 360*h/HSCALE
        // h1 > h2 -> /\ :: d=nvh1-h2 > 180 ?
        dh=std::abs((int)c1t.h-c2t.h);
        if( 2*dh > AHSCALE )
            dh=AHSCALE-dh;
        //are going to pay respect to saturation ?
        if(mode & CDIST_SWEIGHT)
            dh=std::abs(hsdist(dh,c1t.s)-hsdist(dh,c2t.s));
        // are going to weight the value, too ?
        if(mode & CDIST_VWEIGHT)
            dh= (dh*inv_v_weight+std::abs((int)c1t.v - c2t.v))/k1;
        return (uint8_t) (dh >> ABITS);

    } else { // RGB
        db=c1.b-c2.b;
        dg=c1.g-c2.g;
        dr=c1.r-c2.r;

        return (uint8_t) (mode & CDIST_L2 ? \
            std::sqrt((db*db+dg*dg+dr*dr)/3) : \
            (std::abs(db)+std::abs(dg)+std::abs(dr))/3 );
    }
}

inline col32bit alpha_composite(const col32bit top,const col32bit bot){
    col32bit r;
    uint32_t t1,t2;

    if (bot.a == 0 ) return top;
    if (top.a == 0 ) return bot;

    t1=bot.a*(255-top.a);
    t2=(t1+255*top.a);
    r.a = (uint8_t)(t2/255);
    r.b = (uint8_t) ((255*top.b*top.a + t1* bot.b) / t2);
    r.g = (uint8_t) ((255*top.g*top.a + t1* bot.g) / t2);
    r.r = (uint8_t) ((255*top.r*top.a + t1* bot.r) / t2);
    return r;
}

inline col32bit alpha_composite(const col32bitBGRA top,const col32bit bot){
    col32bit r;
    uint32_t t1,t2;

    if (bot.a == 0 ) return {{top.r},{top.g},{top.b},top.a};
    if (top.a == 0 ) return bot;

    t1=bot.a*(255-top.a);
    t2=(t1+255*top.a);
    r.a = (uint8_t)(t2/255);
    r.b = (uint8_t) ((255*top.b*top.a + t1* bot.b) / t2);
    r.g = (uint8_t) ((255*top.g*top.a + t1* bot.g) / t2);
    r.r = (uint8_t) ((255*top.r*top.a + t1* bot.r) / t2);
    return r;
}

uint8_t rgba2cdist(const col32bit c1,const col32bit c2, int mode);
col32bit rgba2hsva(const col32bit c);
col32bit hsva2rgba(const col32bit c);
col64bit rgba2hsva64(const col32bit c);
col32bit hsva2rgba64(const col64bit c);

#endif // FREI0R_COMMON_H_INCLUDED
