/*
 *
 * frei0r_common.cpp
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

#include "frei0r_common.hpp"

uint8_t rgba2cdist(const col32bit c1,const col32bit c2, int mode){
    int db,dg,dr;
    int dh;
    col64bit c1t,c2t;

    const int inv_v_weight=2;
    const int k1= (1+inv_v_weight);

    if (mode & CDIST_HUE) {
        c1t=rgba2hsva64(c1);
        c2t=rgba2hsva64(c2);
        // h360 = 360*h/HSCALE
        // h1 > h2 -> /\ :: d=nvh1-h2 > 180 ?
        dh=std::abs((int)c1t.h-(int)c2t.h);

        if( 2*dh >= AHSCALE )
            dh=AHSCALE-dh-1;
        dh*=2;
// TODO (tmangold#1#): There's some strange behaviour for l2 + hue mode. No idea yet ...        if (mode & CDIST_L2) dh=(dh*dh)/AHSCALE;
        // are we going to pay respect to saturation ?
        if(mode & CDIST_SWEIGHT)
            dh=hsdist(dh,std::abs((int)c1t.s-(int)c2t.s));
        //  are we going to weight the value, too ?
        if(mode & CDIST_VWEIGHT)
            dh= (dh*inv_v_weight+std::abs((int)c1t.v - (int)c2t.v))/k1;
        // double-check
        // (dh*iw + |v1-v2|)/(1+iw)
        // (dh + w*|v1-v2|)/(w+1)

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

col32bit hsva2rgba64(const col64bit c){
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
