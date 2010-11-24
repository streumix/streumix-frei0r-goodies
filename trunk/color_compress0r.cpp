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
#include <cairo.h>
#include <arpa/inet.h>

#include <frei0r.hpp>
#include "frei0r_common.hpp"
#include "hist_plot.hpp"

//typedef bool f0r_param_boolcxx;

class color_compress0r: public frei0r::filter , hist_plot {
  uint8_t map[256];
  f0r_param_double map_xth;
  f0r_param_double map_yth;
  f0r_param_double map_slope;
  ScreenGeometry fscreen;

  uint16_t map64(const uint16_t v);

  public:
  color_compress0r(int wdt, int hgt);
// ~color_compress0r();

  virtual void update();
  f0r_param_double param_xth;
  f0r_param_double param_yth;
  f0r_param_double param_slope;
  f0r_param_boolcxx param_hist;
  f0r_param_boolcxx param_hsv;
};

uint16_t color_compress0r::map64(const uint16_t v){
    uint16_t a16,sm16,m16,n16;
    uint8_t m8;
    const uint16_t k= 1 << ABITS;
    const uint16_t amask= k-1;

    m8 = v >> ABITS;
    m16 = ((uint16_t)map[m8]);
    sm16=m16 << ABITS;
    if (m8 == 255) return sm16;
    else { // linear interpolation based on ABITS
        a16=v & amask;
        if (a16 == 0) return sm16;
        n16=((uint16_t)map[m8+1]);
        return (m16*(k-a16)+n16*a16);
    }
}

color_compress0r::color_compress0r(int wdt, int hgt) : hist_plot(4) {
    // no alpha mapping
    hdata[3].active=false;
    // color defaults are set to RGBA
    hdata[0].c={0,0,1};
    hdata[1].c={0,1,0};
    hdata[2].c={1,0,0};
    hdata[3].c={0.5,0.5,0.5};

    // set defaults
    param_xth=0.5;
    param_yth=0.5;
    param_slope=0.8;
    param_hist=false;
    param_hsv=false;

    register_param(param_xth, "x-threshold", "x value of I/O mapping point with given slope");
    register_param(param_yth, "y-threshold", "y value of I/O mapping point with given slope");
    register_param(param_slope, "slope", "slope at (x/y) threshold point");
    register_param(param_hsv, "hsv", "HSVA color model");
    register_param(hdata[0].active, "channel_0", "activate mapping for channel 1 (B/V)");
    register_param(hdata[1].active, "channel_1", "activate mapping for channel 2 (G/S)");
    register_param(hdata[2].active, "channel_2", "activate mapping for channel 3 (R/H)");
    register_param(hdata[3].active, "channel_3", "activate mapping for channel 4 (A)");
    register_param(param_hist, "histogram", "overlay histogram info");
    register_param(rhsize, "hsize", "size of histogram relative to frame size");
    register_param(rhx, "hposition_x", "relative position of histogram in x direction");
    register_param(rhy, "hposition_y", "relative position of histogram in y direction");

    map_xth=-1.0;
    map_yth=-1.0;
    map_slope=-1.0;

    fscreen.w=wdt;
    fscreen.h=hgt;
    fscreen.bpp=32;
    fscreen.size = fscreen.w * fscreen.h;
    fscreen.stride = fscreen.w;
    hist_init(fscreen);
 }

void color_compress0r::update() {
  int i,j,m,n;
  unsigned int l;
  union {
      col32bitBGRA c;
      uint32_t u;
  } cc;

  //param_hsv=true;
  //param_hist=false;

  //uint32_t * surface_buf;
  uint32_t a8b,t1,t2;

  col64bit chsv64;
  col32bit chsv;
  col32bit *colin = (col32bit*)in;
  col32bit *colout = (col32bit*)out;

  //re-adjust mapping table
  if ((map_xth != param_xth) || (map_yth != param_yth) || ( map_slope != param_slope)) {
      map_xth = param_xth;
      map_yth = param_yth;
      map_slope = param_slope;
      double a;
      if (param_slope < 0.01)
        a=tan(M_PI_2*0.01);
      else if (param_slope > 0.99)
        a=tan(M_PI_2*0.99);
      else
        a=tan(M_PI_2*param_slope);
      for (i=0; i<256;i++){
        map[i]=int(0.5+255.0*POW_COMP(i/255.0,map_xth,map_yth,a));
      }
  }

  // overlay histogram
  if (param_hist){  // there's a histogram to overlay
    // zero data
    for(l=0;l<hist_data_n;l++){
        memset( hdata[l].v, '\0',4*256 );
    }
    if ( param_hsv ) {
      for(i=0; i<fscreen.size; i++){
        chsv=rgba2hsva(colin[i]);
        (hdata[0].v[ chsv.b ])++;
        (hdata[1].v[ chsv.g ])++;
        (hdata[2].v[ chsv.r ])++;
        (hdata[3].v[ chsv.a ])++;
      }
    } else {
      for(i=0; i<fscreen.size; i++){
        (hdata[0].v[ colin[i].b ])++;
        (hdata[1].v[ colin[i].g ])++;
        (hdata[2].v[ colin[i].r ])++;
        (hdata[3].v[ colin[i].a ])++;
      }
    }
    hist_set_range();
    hist_draw(map);

    // map & composite overlay
    if(param_hsv){
        for(n=0; n<fscreen.h; n++){
          for(m=0; m<fscreen.w; m++){
            j=m+n*fscreen.stride;

            chsv64=rgba2hsva64(colin[j]);
            if( hdata[0].active ) chsv64.b=map64(chsv64.b);
            if( hdata[1].active ) chsv64.g=map64(chsv64.g);
            if( hdata[2].active ) chsv64.r=map64(chsv64.r);
            if( hdata[3].active ) chsv64.a=map[chsv64.a];
            chsv=hsva2rgba64(chsv64);

            cc.u=surface_buf[m+n*hist_screen.stride];

            a8b=hdata[3].active ? map[chsv.a] : chsv.a;
            t1=a8b*(255-cc.c.a);
            t2=(t1+255*cc.c.a);
            chsv.a= (uint8_t)(t2/255);
            if ( chsv.a > 0 ){
                colout[j].b = (uint8_t) ((255*cc.c.b*cc.c.a + t1*chsv.b) / t2);
                colout[j].g = (uint8_t) ((255*cc.c.g*cc.c.a + t1*chsv.g) / t2);
                colout[j].r = (uint8_t) ((255*cc.c.r*cc.c.a + t1*chsv.r) / t2);
                colout[j].a=chsv.a;
            } else {
                colout[j]=chsv;
            }
          }
        }
    } else { // !param_hsv
        for(n=0; n<fscreen.h; n++){
          for(m=0; m<fscreen.w; m++){
            j=m+n*fscreen.stride;

            cc.u=surface_buf[m+n*hist_screen.stride];

            a8b=hdata[3].active ? map[colin[j].a] : colin[j].a;
            t1=a8b*(255-cc.c.a);
            t2=(t1+255*cc.c.a);
            colout[j].a= (uint8_t)(t2/255);
            if (colout[j].a > 0 ){
                colout[j].b = (uint8_t) ((255*cc.c.b*cc.c.a + t1*(hdata[0].active ? map[colin[j].b] : colin[j].b)) / t2);
                colout[j].g = (uint8_t) ((255*cc.c.g*cc.c.a + t1*(hdata[1].active ? map[colin[j].g] : colin[j].g)) / t2);
                colout[j].r = (uint8_t) ((255*cc.c.r*cc.c.a + t1*(hdata[2].active ? map[colin[j].r] : colin[j].r)) / t2);
            } else {
                // no alpha -> no need to overlay, just map colors as required
                colout[j].b= hdata[0].active ? map[colin[j].b] : colin[j].b;
                colout[j].g= hdata[1].active ? map[colin[j].g] : colin[j].g;
                colout[j].r= hdata[2].active ? map[colin[j].r] : colin[j].r;
            }
          }
        }
    }
  } else { // no histogram -> mapping only
    // pixels are F0R_COLOR _MODEL_BGRA8888 format
    if( param_hsv ){
        for(i=0; i<fscreen.size; i++){
            chsv64=rgba2hsva64(colin[i]);
            if( hdata[0].active ) chsv64.b=map64(chsv64.b);
            if( hdata[1].active ) chsv64.g=map64(chsv64.g);
            if( hdata[2].active ) chsv64.r=map64(chsv64.r);
            if( hdata[3].active ) chsv64.a=map[chsv64.a];
            colout[i]=hsva2rgba64(chsv64);
        }
    } else { // !param_hsv
        for(i=0; i<fscreen.size; i++){
            colout[i].b= hdata[0].active ? map[colin[i].b] : colin[i].b;
            colout[i].g= hdata[1].active ? map[colin[i].g] : colin[i].g;
            colout[i].r= hdata[2].active ? map[colin[i].r] : colin[i].r;
            colout[i].a= hdata[3].active ? map[colin[i].a] : colin[i].a;
        }
    }

  }
}

frei0r::construct<color_compress0r> plugin("color_compress0r",
				  "color space / levels compressor based on pow() function",
				  "Mangold, Toby",
				  1,0,F0R_COLOR_MODEL_RGBA8888
				  );
