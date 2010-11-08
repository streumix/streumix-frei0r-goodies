/*
 *
 * color_compress0r Plugin
 * Copyright (C) 2010 Toby Mangold
 *
 * This source code  is free software; you can  redistribute it and/or
 * modify it under the terms of the GNU Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This source code is distributed in the hope that it will be useful,
 * but  WITHOUT ANY  WARRANTY; without  even the  implied  warranty of
 * MERCHANTABILITY or FITNESS FOR  A PARTICULAR PURPOSE.  Please refer
 * to the GNU Public License for more details.
 *
 * You should  have received  a copy of  the GNU Public  License along
 * with this source code; if  not, write to: Free Software Foundation,
 * Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <cstdlib>
#include <cassert>
#include <cstring>
#include <cmath>
#include <arpa/inet.h>

#include <vector>
#include <iostream>

#define CDIST 0
#define ABITS 4
#define HSCALE 256

/*
 * ABITS STAT:
 * 4 : 100%/+-0
 * 3 : 98.6%/+-0 1.4%/+-1 (max total=1))
 * 2 : 90%/+-0 10%/+-1 (max total=2)
 * 1 : 80%/+-0 16%/+-1 4%/+-2 (max total=3)
 * 0 : 60%/+-0 12%/+-1 8%/+-2 20%/ >3 (max total=7)
 */

using namespace std;

// BGRA in memory (= big-endian), equals ARGB in bit-order on little-endian (=cairo on x86)
typedef struct col32bit_s {
    uint8_t b;
    uint8_t g;
    uint8_t r;
    uint8_t a;
} col32bit;

typedef struct col64bit_s {
    uint16_t b;
    uint16_t g;
    uint16_t r;
    uint16_t a;
} col64bit;

typedef struct col128bit_s {
    uint32_t b;
    uint32_t g;
    uint32_t r;
    uint32_t a;
} col128bit;

class col32bite {
public:
    uint8_t b;
    uint8_t g;
    uint8_t r;
    uint8_t a;

    float db,dg,dr;

    int h;
//void operator<<(){
//   cout << "["<< b << ", "<< g << ", "<< r "," << a << "]";
//    }

void print(){
    cout << "[" << b << "," << g << "," << r << "," << a << "]";
    }


int dist(col32bit y){
    return abs(b-y.b)+abs(g-y.g)+abs(r-y.r)+abs(a-y.a);
    }

int dist(col128bit y){
    return abs(b-uint8_t(y.b>>ABITS))+abs(g-uint8_t(y.g>>ABITS))+abs(r-uint8_t(y.r>>ABITS))+abs(a-uint8_t(y.a));
    }

uint8_t max(){
    if (b > g) return (r>b) ? r : b;
    else return (r>g) ? r : g;
    }

uint8_t min(){
    if (b < g) return (r<b) ? r : b;
    else return (r<g) ? r : g;
    }

}; // class

col32bite rgba2hsvaf(const col32bit c){

  col32bite r;
  float iMin,delta;
  r.a=c.a;

  if (c.r > c.g) {
      r.b = std::max (c.r, c.b);
      iMin = std::min (c.g, c.b);
  } else {
      r.b = std::max (c.g, c.b);
      iMin = std::min (c.r, c.b);
  }
  delta = r.b - iMin;

  r.db=r.b/255.0;
  if (r.b == 0){
    r.g = 0;
    r.dg =0.0;
  } else
    r.g = (uint8_t)((255.0*delta)/r.b);
    r.dg = delta/r.b;

  if (r.g == 0) {
      r.r = 0;
      r.dr = 0.0;
  } else {
      float h;
      if (c.r == r.b)
        h = (60.0 * (c.g - c.b)) / delta;
      else if (c.g == r.b)
        h = 120.0 + ( 60.0 * (c.b - c.r)) / delta;
      else
        h = 240.0 + ( 60.0 * (c.r - c.g)) / delta;

      if (h < 0.0) h += 360.0;
      else if (h >= 360.0) h -= 360.0;
      r.r = (uint8_t) floor((HSCALE*h)/360.0);
      r.dr = h/360.0;
  }
  return r;
}

col32bite hsva2rgbaf(const col32bite c){
  col32bite r;
  float chroma,X,m;
  int H;
  // R G B
  // H S V
  r.a=c.a;

  // set chroma
  chroma = ( c.db * c.dg );
  m = c.db - chroma;

  // chroma  == 0 <-> c.g == 0 --> m=c.b
  if (c.dg == 0) {
      r.dr=c.db;
      r.dg=c.db;
      r.db=c.db;
  } else {

    H = floor(3*c.dr);
    X = chroma*(1-abs( (6*c.dr)- 2*H - 1.0  )) + m ;

    H = floor(6*c.dr);

    // ( chroma + m ) --> c.b ;
    switch (H) {
        case 0:
          r.dr = c.db  ;
          r.dg = X;
          r.db = m ;
          break;

        case 1:
          r.dr = X;
          r.dg = c.db ;
          r.db = m ;
          break;

        case 2:
          r.dr = m ;
          r.dg = c.db ;
          r.db = X;
          break;

        case 3:
          r.dr = m ;
          r.dg = X;
          r.db = c.db ;
          break;

        case 4:
          r.dr = X;
          r.dg = m ;
          r.db = c.db ;
          break;

        case 5:
          r.dr = c.db ;
          r.dg = m ;
          r.db = X;
          break;
    }
  }
  r.b=floor(r.db*255.0);
  r.g=floor(r.dg*255.0);
  r.r=floor(r.dr*255.0);

  return r;
}

col32bit rgba2hsva(const col32bit c){

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

col32bit hsva2rgba(const col32bit c){
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

col64bit rgba2hsva64(const col32bit c){

  col64bit r;
  uint32_t iMin,iMax,chroma,k1,k2;
  r.a=c.a;

  k1=255 << ABITS;
  k2=HSCALE << ABITS;

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

col32bit hsva2rgba64(const col64bit c){
  col32bit r;

  //uint32_t chroma;
  uint32_t m,k1;
  int H,X,k2,k3,ir,ig,ib;
  // R G B
  // H S V
  r.a=c.a;

  k1=255<<ABITS;
  k2=HSCALE<<ABITS;
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

col128bit rgba2hsva128(const col32bit c){

  col128bit r;
  uint32_t iMin,iMax,chroma,k1,k2;
  r.a=c.a;

  k1=255 << ABITS;
  k2=HSCALE << ABITS;

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

col32bit hsva2rgba128(const col128bit c){
  col32bit r;

  //uint32_t chroma;
  uint32_t m,k1,k3;
  int H,X,k2;
  // R G B
  // H S V
  r.a=c.a;

  k1=255<<ABITS;
  k2=HSCALE<<ABITS;
  k3=1<<(ABITS-1);

  // set chroma
  //chroma = ( c.b * c.g )/k1;
  //m = c.b - chroma;
  m = (c.b*(k1 - c.g ))/k1;

  // chroma  == 0 <-> c.g == 0 --> m=c.b
  if (c.g == 0) {
      r.b = ( r.g = ( r.r = c.b >> ABITS ));
  } else {

    H = (int)((3*c.r)/k2);
    //X = chroma*(1-abs(X' -1.0  )) ;
    //X' = (6*c.r/k2)- 2*H
    //X' = (6*c.r- 2*H*k2)/k2;
    //X = chroma*(1-abs(   (6*c.r- 2*H*k2)/k2 -k2/k2  )) ;
    //X = chroma*(1-abs((6*c.r- 2*H*k2- k2)/k2 )) ;
    //X = k1*chroma*(k2-abs(6*((int)c.r)- 2*H*k2 - k2))/(k1*k2) ;

    X = ((c.b*c.g)/k2)*(k2-abs(6*((int)c.r)- 2*H*k2 - k2)) ;

    // removing additional bits --> unit8
    X=( (X+c.b*(k1 - c.g ))/k1 + k3 ) >>ABITS;
    m=m >> ABITS;

    H = (int)(6*c.r)/k2;
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


int main()
{
    int ir,ig,ib,d1,d2,d3;
    int d,dmax,count=0;
    float df,dmaxf,kf1,kf2;
    vector<int> dr_stat(5,0);
    vector<int> dg_stat(5,0);
    vector<int> db_stat(5,0);
    col32bit ipix,c8pix,opix,o8pix;
    //col32bite cfpix,ofpix;
    col64bit cpix;

    const int N=(256*256*256)/1000;
    dmax=0;
    dmaxf=0.0;
    kf1=1.0*(255<<ABITS);
    kf2=1.0*(HSCALE<<ABITS);

    for (ir=0; ir<256;ir++){
        for (ig=0; ig<256;ig++){
            for (ib=0; ib<256;ib++){
              ipix={(uint8_t)ib,(uint8_t)ig,(uint8_t)ir,0};
              cpix=rgba2hsva64(ipix);
              //c8pix=rgba2hsva(ipix);
              //cfpix=rgba2hsvaf(ipix);
              //c8pix={(uint8_t)(cpix.b>>ABITS),(uint8_t)(cpix.g>>ABITS),(uint8_t)(cpix.r>>ABITS),0};

              opix=hsva2rgba64(cpix);
              //ofpix=hsva2rgbaf(cfpix);

              //opix={(uint8_t)(cpix.b>>ABITS),(uint8_t)(cpix.g>>ABITS),(uint8_t)(cpix.r>>ABITS),0};

              //if((df=abs(255.0*o8pix.db - ipix.b)) > dmaxf ) dmaxf=df;
              //if((df=abs(255.0*o8pix.dg - ipix.g)) > dmaxf ) dmaxf=df;
              //if((df=abs(255.0*o8pix.dr - ipix.r)) > dmaxf ) dmaxf=df;

              //d1=abs((int)c8pix.r-(int)cpix.r);
              //d2=abs((int)c8pix.g-(int)cpix.g);
              //d3=abs((int)c8pix.b-(int)cpix.b);

              d1=(int)opix.r-(int)ipix.r;
              d2=(int)opix.g-(int)ipix.g;
              d3=(int)opix.b-(int)ipix.b;
              d=abs(d1)+abs(d2)+abs(d3);
              if (d<=2){
                  dr_stat[d1+2]++;
                  dg_stat[d2+2]++;
                  db_stat[d3+2]++;
              } else count ++;
              if (d > dmax ){
                        dmax=d;
                        //df=ofpix.dist(opix);
              }
            }
        }
    }
    //cout << "found " << count << " conversion errors (max="<< dmax << "|"<< dmaxf <<")" << endl;
    cout << "found " << count << " large conversion errors (max="<< dmax <<")" << endl;
    cout << "R_stat " << dr_stat[0]/N << ", "<< dr_stat[1]/N << ", "<< dr_stat[2]/N \
         << ", "<< dr_stat[3]/N << ", "<< dr_stat[4]/N << endl;
    cout << "G_stat " << dg_stat[0]/N << ", "<< dg_stat[1]/N << ", "<< dg_stat[2]/N \
         << ", "<< dg_stat[3]/N << ", "<< dg_stat[4]/N << endl;
    cout << "B_stat " << db_stat[0]/N << ", "<< db_stat[1]/N << ", "<< db_stat[2]/N \
         << ", "<< db_stat[3]/N << ", "<< db_stat[4]/N << endl;
}
