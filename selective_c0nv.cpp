/*
 *
 * selective_c0nv Plugin
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

#include "hist_plot.hpp"
//#include "common.h"

// number of frames in ringbuffer (has to be 2^x ), defines max. kernel time duration
#define PLANES_MAX 8
#define KERNELS_REDUCED_PREC 2
#define KERNELS_PER_FRAME (256>>KERNELS_REDUCED_PREC)
// max. radius of kernel
#define R_MAX 10

#define RB_INCR(x) ((++x) & (PLANES_MAX-1))
#define RB_OFFSET(x,y) ((x+y) & (PLANES_MAX-1))

class convkernel {
  int kernels_per_frame;
  int kernel_size;

public:
  uint8_t* data;
  int radius;
  int frames;
  int frame_size;

  convkernel(){
    radius=0;
    frames=0;
    data=NULL;
  }

  ~convkernel(){
    if (data) delete data;
  }

  int get_offset(int n){
    if (data)
        return (n-1)*kernel_size;
    else
        return 0;
  }

  void update(double kr){
    int i,j,k,l,m,n,p0;
    float r,w_sum;
    float w[TRI_IND(R_MAX,R_MAX)];

    // clipping requested radius/kernel size to match preserved memory
    radius = kr > R_MAX ? R_MAX : floor(kr) ;
    // no time invariant convolution at the moment
    frames=1;
    kernels_per_frame = KERNELS_PER_FRAME ;
    kernel_size=TRI_IND(radius,radius)+1;
    frame_size = kernel_size*kernels_per_frame;

    if (data) delete data ;
    data = new uint8_t[frames*frame_size] ;
    assert (data);

    for (n=1;n<=kernels_per_frame;n++){
      w_sum=0.0;
      // in a first step, we calculate the unscaled coeffs
      for (k=0;k<=radius;k++){
        for (l=0;l<=k;l++){
            r=sqrt(k*k+l*l);
            i=TRI_IND(k,l);
            if (r <= radius*n/kernels_per_frame) {
                r= r*kernels_per_frame/n/ radius;
                r=r*r;
                r=r*r-2.0*r+1;
                w[i]=r;
                // Depending on its local coordinates, a coeff will be reused x times
                // when doing a full convolution.
                // Remember: We're storing just half of a quadrant with "overlapping" boundaries
                if (k==l || l==0){
                  if (k==0) w_sum+=r;
                  else w_sum+=4*r;
                } else w_sum+=8*r;
            } else w[i]=0.0;
        }
      }
      // now we can normalize to an unit integrale sum
      p0=get_offset(n);
      i=0;
      for (k=0;k<=radius;k++){
        for (l=0;l<=k;l++){
            j=TRI_IND(k,l);
            m=p0+j;
            data[m] = (uint8_t) (0.5+255.0*w[j]/w_sum);
            // intergrating the final integer vlues to collect rounding errors
            if (k==l || l==0){
              if (k==0) i+=data[m];
              else i+=4*data[m];
            } else i+=8*data[m];
        }
      }
      // Since the convolution kernel must convserve "energy",
      // a final adjustment fixes rounding errors
      data[p0] += 255-i;
//      for (l=0;l<kernel_size;l++){
//          data[p0+l] = (uint8_t) (255.0*w[l]/w_sum);
//      }
    }
  }
  // end of class convkernel
};

class selective_c0nv: public frei0r::filter , hist_plot {
  uint8_t map[256];
  f0r_param_double map_xth;
  f0r_param_double map_slope;
  ScreenGeometry fscreen;
  int max_distance;
  convkernel kernel;

  col128bit *planebuf;
  col128bit *planetable[PLANES_MAX];
  int plane;

  col32bit rgba2hsva(const col32bit);
  col32bit hsva2rgba(const col32bit);

  public:
  selective_c0nv(int wdt, int hgt);
  ~selective_c0nv();

  virtual void update();
  f0r_param_double param_xth;
  f0r_param_double param_slope;
  f0r_param_double param_r;
  f0r_param_color param_color;
  f0r_param_boolcxx param_plot;
  f0r_param_boolcxx param_inv;
};

selective_c0nv::selective_c0nv(int wdt, int hgt) {
    // set defaults
    int i;

    param_xth=0.5;
    param_slope=0.8;
    param_plot=false;
    param_inv=false;
    param_r=5;
    param_color.r=1.0;
    param_color.g=1.0;
    param_color.b=1.0;

    register_param(param_xth, "threshold", "input value of mapping, where output level is 50% at the specified slope");
    register_param(param_slope, "slope", "slope at threshold value, where output level is 50%");
    register_param(param_color, "color", "ref color for distance");
    register_param(param_inv, "inverse", "inverse the regular distance->radius mapping");
    register_param(param_r, "radius", "Maximum radius of dynamic convolution kernel = spatial extension");

    register_param(param_plot, "map", "activate plot of mapping function");
    register_param(rhsize, "hsize", "size of histogram relative to frame size");
    register_param(rhx, "hposition_x", "relative position of histogram in x direction");
    register_param(rhy, "hposition_y", "relative position of histogram in y direction");

    map_xth=-1.0;
    map_slope=-1.0;

    fscreen.w=wdt;
    fscreen.h=hgt;
    fscreen.bpp=32;
    fscreen.size = fscreen.w * fscreen.h;
    fscreen.stride = fscreen.w;

    planebuf =  new col128bit[PLANES_MAX*fscreen.size];
    //memset(planetable[RB_OFFSET(plane,kernel.frames-1)],'\0',sizeof(col128bit)*fscreen.size);
    for(i=0;i<PLANES_MAX;i++)
      planetable[i] = &planebuf[fscreen.size *i];
    plane = 0;

    block_histplot=true;
    hist_init(fscreen);
 }

selective_c0nv::~selective_c0nv() {
  delete planebuf;
}

void selective_c0nv::update() {
  int i,j,j0,m,n;
  int w0,h0,w,h,wl,hl;
  int w1,w2,h1,h2;
  int lr;
  int p0;
  int d;

  union {
      col32bit c;
      uint32_t u;
  } cc;
  uint32_t t1,t2;

  col32bit *colin = (col32bit*)in;
  col32bit *colout = (col32bit*)out;
  col128bit *oplane;

  // updating the scaling base
  int ir,ig,ib;

  ir=255*param_color.r;
  ig=255*param_color.g;
  ib=255*param_color.b;
  max_distance= std::max(ir,255-ir) + std::max(ig,255-ig) + std::max(ib,255-ib);

  //re-adjust convolution kernel
  if (kernel.radius != param_r) kernel.update(param_r);

  oplane=planetable[plane];
  // clear output buffer and re-use, most advanced frame in ring buffer
  memset((char *)colout,'\0',sizeof(col32bit)*fscreen.size);
  memset(planetable[RB_OFFSET(plane,kernel.frames-1)],'\0',sizeof(col128bit)*fscreen.size);

  // re-set mapping table
  if ((map_xth != param_xth) || ( map_slope != param_slope)) {
      map_xth = param_xth;
      map_slope = param_slope;
      double a;
      if (param_slope < 0.01)
        a=tan(M_PI_2*0.01);
      else if (param_slope > 0.99)
        a=tan(M_PI_2*0.99);
      else
        a=tan(M_PI_2*param_slope);
      for (i=0; i<256;i++){
        map[i]=int(0.5+255.0*POW_COMP(i/255.0,map_xth,a));
      }
  }

  // iterating over all source pixels from input
  for(h0=0; h0<fscreen.h; h0++) {
   for(w0=0; w0<fscreen.w; w0++) {
    // deriving local scaling factor v=[0;255] from color c of current pixel
    j0=h0*fscreen.stride+w0;
    //Calculating color distance, which is a fast l1-norm for now to avoid
    //the slow sqrt() function call. l2 might be added ... optinally ... later
    d=map[255*(std::abs(ir-colin[j0].r) + std::abs(ig-colin[j0].g) + std::abs(ib-colin[j0].b))/max_distance];
    // What are we going to blur ? The far-distant or close colors ?
    if(param_inv) d =  1 + (d >> KERNELS_REDUCED_PREC);
    else d = KERNELS_PER_FRAME - (d >> KERNELS_REDUCED_PREC);

    t2=0;
    p0=kernel.get_offset(d);
    // according to profiler, fCLIP would get re-eval every loop without the following pre-calc ....
    h1=CLIP0(h0-kernel.radius);
    h2=CLIP(h0+kernel.radius,fscreen.h-1);
    w1=CLIP0(w0-kernel.radius);
    w2=CLIP(w0+kernel.radius,fscreen.w-1);
    // building convolution of source pixel, color c @ w0,h0 with kernel
    for(h=h1; h<=h2; h++) {
      for(w=w1; w<=w2; w++) {
        j=w+h*fscreen.stride;
        wl = std::abs(w-w0);
        hl = std::abs(h-h0);
        // calculate relative radius to grab scaling factor = elem of convolution kernel
        lr = (hl <= wl) ? TRI_IND(wl,hl) : TRI_IND(hl,wl);
        t1 = kernel.data[p0+lr];
        oplane[j0].b += t1*colin[j].b;
        oplane[j0].g += t1*colin[j].g;
        oplane[j0].r += t1*colin[j].r;
        // Alpha channel is used to sum up "energie" contributions from blurred pixels.
        // This value is required for normalization afterwards, otherwise "energy-less"
        // black cannot be blurred ...
        oplane[j0].a += t1;
        //t2+=t1;
      }
    }
    t2=oplane[j0].a;
    t2=255-t2;
    // done with source pixel @ w,h
    // ensuring pixel energy conservation
    oplane[j0].b += t2*colin[j].b;
    oplane[j0].g += t2*colin[j].g;
    oplane[j0].r += t2*colin[j].r;
    oplane[j0].a += t2;
   }
  }

  // transfering plane "0" to output
  for(h0=0; h0<fscreen.h; h0++) {
    for(w0=0; w0<fscreen.w; w0++) {
        j0=h0*fscreen.stride+w0;
//        colout[j0].b = fCLIP8(oplane[j0].b/255);
//        colout[j0].g = fCLIP8(oplane[j0].g/255);
//        colout[j0].r = fCLIP8(oplane[j0].r/255);
//        colout[j0].a = fCLIP8(oplane[j0].a/255);
        t1=oplane[j0].a;
        colout[j0].b = (uint8_t)(oplane[j0].b/t1);
        colout[j0].g = (uint8_t)(oplane[j0].g/t1);
        colout[j0].r = (uint8_t)(oplane[j0].r/t1);
        colout[j0].a = colin[j0].a;
    }
  }

  //param_plot=false;
  // overlay histogram
  if (param_plot){
    hist_draw(map);
    // composite overlay
    for(n=0; n<fscreen.h; n++){
      for(m=0; m<fscreen.w; m++){
        j=m+n*fscreen.stride;

        // the cairo pixel to overlay
        cc.u=surface_buf[m+n*hist_screen.stride];

        t1= ((uint32_t)colout[j].a)*(255-cc.c.a);
        t2= (t1+255*cc.c.a);
        colout[j].a= (uint8_t)(t2/255);
        if (colout[j].a > 0 ){
            colout[j].b = (uint8_t) ((255*cc.c.b*cc.c.a + t1*colout[j].b) / t2);
            colout[j].g = (uint8_t) ((255*cc.c.g*cc.c.a + t1*colout[j].g) / t2);
            colout[j].r = (uint8_t) ((255*cc.c.r*cc.c.a + t1*colout[j].r) / t2);
        }
      } // for m
    } // for n
  } // if param_plot

}

frei0r::construct<selective_c0nv> plugin("selective_c0nv",
				  "dynamic convolution based blur filter",
				  "Mangold, Toby",
				  1,0);
