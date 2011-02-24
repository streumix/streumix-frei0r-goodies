/*
 *
 * color_distanc0r Plugin
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

class color_distanc0r: public frei0r::filter,hist_plot {
  uint8_t map[256];
  f0r_param_double map_xth;
  f0r_param_double map_yth;
  f0r_param_double map_slope;
  ScreenGeometry fscreen;

  col32bit refcol;

  public:
  color_distanc0r(int wdt, int hgt);
// ~color_distanc0r();

  virtual void update();
  f0r_param_color param_color;
  f0r_param_double param_xth;
  f0r_param_double param_yth;
  f0r_param_double param_slope;
  f0r_param_boolcxx param_hist;
  f0r_param_boolcxx param_hue;
  f0r_param_boolcxx param_vmult;
  f0r_param_boolcxx param_smult;
  f0r_param_boolcxx param_l2;
  f0r_param_boolcxx param_alpha;
};

color_distanc0r::color_distanc0r(int wdt, int hgt): hist_plot(1) {
    // set defaults
    param_color={0,0,0};
    param_xth=0.5;
    param_yth=0.5;
    param_slope=0.8;
    param_hist=false;
    param_hue=false;
    param_vmult=false;
    param_smult=false;
    param_l2=false;
    param_alpha=false;

    hdata[0].c={0.8,0.3,0.3};

    register_param(param_color, "color", "ref color for distance");
    register_param(param_xth, "x-threshold", "x value of I/O mapping point with given slope");
    register_param(param_yth, "y-threshold", "y value of I/O mapping point with given slope");
    register_param(param_slope, "slope", "slope at (x/y) threshold point");
    register_param(param_hue, "hue", "HUE based color distance");
    register_param(param_vmult, "valueweight", "hue distance is weighted with value");
    register_param(param_smult, "satweight", "hue distance is weighted with saturation");
    register_param(param_l2, "l2norm", "utilize a l2-distance norm instead of l1");
    register_param(param_alpha, "alpha_mode", "store distance into alpha channel");
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

void color_distanc0r::update() {
  int i,j,m,n;
  union {
      col32bitBGRA c;
      uint32_t u;
  } cc;

  int cmode=0;
  uint8_t tdist;
  col32bit *colin = (col32bit*)in;
  col32bit *colout = (col32bit*)out;


  // there's atwist !!!
  refcol.b= (uint8_t)floor(param_color.b*255);
  refcol.g= (uint8_t)floor(param_color.g*255);
  refcol.r= (uint8_t)floor(param_color.r*255);
  refcol.a= 255;

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

  // construct mode
  cmode  = (param_l2 ? CDIST_L2 : 0) \
        | (param_hue ? CDIST_HUE : 0) \
        | (param_vmult ? CDIST_VWEIGHT : 0) \
        | (param_smult ? CDIST_SWEIGHT : 0);

  // calc cdists
  if (param_hist) {
    // there's a histogram to overlay -> zero data
    memset( hdata[0].v, '\0',256*sizeof(unsigned int) );

    for(n=0; n<fscreen.size; n++){
        tdist=rgba2cdist(colin[n],refcol, cmode);
        hdata[0].v[tdist]++;
        if (param_alpha){
            colout[n]=colin[n];
            colout[n].a=255-map[tdist];
        } else {
            colout[n].b=colout[n].g=colout[n].r=map[tdist];
            colout[n].a=colin[n].a;
        }
    }
  } else {
    if(param_alpha){
        for(n=0; n<fscreen.size; n++){
            colout[n]=colin[n];
            colout[n].a=255-map[rgba2cdist(colin[n],refcol,cmode)];
        }
    } else {
        for(n=0; n<fscreen.size; n++){
            colout[n].b=colout[n].g=colout[n].r=map[rgba2cdist(colin[n],refcol,cmode)];
            //colout[n].b=colout[n].g=colout[n].r=rgba2cdist(colin[n],refcol,cmode);
            colout[n].a=colin[n].a;
        }
    }
  }

  // overlay histogram
  if (param_hist){  // there's a histogram to overlay
    hist_set_range();
    hist_draw(map);

    //composite overlay
    j=0;
    for(n=0; n<fscreen.h; n++){
        for(m=0; m<fscreen.w; m++){
            cc.u=surface_buf[m+n*hist_screen.stride];
            colout[j]=alpha_composite(cc.c,colout[j]);
            j++;
        }
    }
  }
}

frei0r::construct<color_distanc0r> plugin("color_distanc0r",
				  "color distance incl. non-linear distance compression",
				  "Mangold, Toby",
				  1,0,
				  F0R_COLOR_MODEL_RGBA8888);
