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

#define ABITS 4
#define HSCALE 256

#define HIST_DATA_N 4
#define POW_COMP(x,tx,ty,s) ( x <= tx ? ty*pow(x/tx,s*tx/ty) : 1.0-(1-ty)*pow((1-x)/(1-tx),s*(1-tx)/(1-ty)))

typedef bool f0r_param_boolcxx;

typedef struct {
  int w;
  int h;
  int stride;
  int bpp;
  int size;
} ScreenGeometry;

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

typedef struct {
    unsigned int min;
    unsigned int max;
} urange;

typedef struct {
    f0r_param_boolcxx active;
    f0r_param_color c;
    urange range;
    unsigned int v[256];
} hist_data_channel;

typedef struct {
    double width_grid;
    double gray_grid;
    double width_hist;
    double alpha_hist;
    double width_map;
    double gray_map;
    double width_shadow_incr;
    double offset_shadow;
    double alpha_shadow;
} hist_plot_config;

class hist_plot
{
    urange gminmax;
    cairo_surface_t *surface;
    cairo_t *cr;

    int draw_hdata(unsigned int *h,urange mm,double r,double g,double b,double a){
        int i;
        double s=100.0/(log2((double) (mm.max-mm.min+1) ));
        double x;
        cairo_path_t * t;

        if(!cr) return 0;
        cairo_set_source_rgba(cr, r, g,b,a);
        x= h[0] > mm.min ? log2((double) (h[0]-mm.min+1) ) : 0;
        cairo_move_to (cr, 0, (100.0-s*x));
        for (i=1; i<256;i++){
            x= h[i] > mm.min ?  log2((double)(h[i]-mm.min+1)) : 0;
            cairo_line_to (cr, i, (100.0-s*x));
        }
        t=cairo_copy_path(cr);
        cairo_line_to (cr, 255.0, 100.0);
        cairo_line_to (cr, 0, 100.0);
        cairo_close_path(cr);
        cairo_fill(cr);
        cairo_set_source_rgb (cr, r,g,b);
        cairo_append_path(cr,t);
        cairo_stroke(cr);
        cairo_path_destroy(t);
        return 1;
    }

protected:
    ScreenGeometry hist_screen;
    uint32_t * surface_buf;
    bool block_histplot;

    void dump_png(const char *s){
        if(surface) cairo_surface_write_to_png (surface,s);
    }

public:
    f0r_param_double rhsize;
    f0r_param_double rhx;
    f0r_param_double rhy;
    hist_plot_config hconf;
    hist_data_channel hdata[HIST_DATA_N];

    hist_plot(int wdt, int hgt){
        hist_plot();
    }

    hist_plot(){
        int i;
        block_histplot=false;
        cr=NULL;
        surface=NULL;
        // init of max data
        gminmax= { 2*hist_screen.size,0};

        for(i=0;i<HIST_DATA_N;i++){
            hdata[i].range=gminmax;
            hdata[i].active=true;
        }
        // no alpha mapping
        hdata[3].active=false;
        // color defaults are set to RGBA
        hdata[0].c={0,0,1};
        hdata[1].c={0,1,0};
        hdata[2].c={1,0,0};
        hdata[3].c={0.5,0.5,0.5};
        // relative size
        rhsize=0.5;
        // position on screen
        rhx=0.5;
        rhy=0.95;

        hist_screen.w=0;
        // drawing defaults
        hconf.width_grid=1.75;
        hconf.gray_grid=0.8;
        hconf.width_hist=1.75;
        hconf.alpha_hist=0.2;
        hconf.width_map=2.0;
        hconf.gray_map=0.9;
        hconf.width_shadow_incr=1;
        hconf.offset_shadow=1.0;
        hconf.alpha_shadow=0.5;
    }

    hist_plot(ScreenGeometry sg){
        hist_plot();
        hist_init(sg);
    }

    ~hist_plot(){
        // DEBUG only !!!!!
        // cairo_surface_write_to_png (surface,"test_cairo_draw.png");
        // DEBUG
        if (cr) cairo_destroy (cr);
        if (surface) cairo_surface_destroy (surface);
    }

    void hist_init(ScreenGeometry sg){
        // cairo init
        surface = NULL;
        cr = NULL;
        surface_buf= NULL;
        hist_screen=sg;
        hist_screen.stride=0;
    }

    uint32_t * hist_get_data(){
        if (surface) return (uint32_t *) cairo_image_surface_get_data(surface);
        else return NULL;
    }

    int hist_draw(uint8_t *map){
        unsigned int i;
        double dx,dy;
        cairo_path_t * tmp_path;

        if (hist_screen.w < 1 ) return 0;
        else {
            // reusing a surface/context seems to cause serious trouble
            // therefore, let's destroy the existing stuff and create a brand new context / surface
            if (cr) cairo_destroy (cr);
            if (surface) cairo_surface_destroy (surface);
            surface = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, hist_screen.w,hist_screen.h);
            cr = cairo_create (surface);
            hist_screen.stride=cairo_image_surface_get_stride(surface) / sizeof(uint32_t);
            surface_buf=(uint32_t *) hist_get_data();
        }

        // clearing surface is not required due to the recreation of the surface
        //cairo_save (cr);
        //cairo_set_operator (cr, CAIRO_OPERATOR_CLEAR);
        //cairo_set_source_rgba (cr, 0.0, 0.0, 0.0,0.0);
        //cairo_set_operator (cr, CAIRO_OPERATOR_SOURCE);
        //cairo_paint (cr);
        //cairo_restore (cr);

        cairo_set_source_rgba (cr, 0.0, 0.0, 0.0,hconf.alpha_shadow);
        cairo_set_line_width (cr,hconf.width_grid+hconf.width_shadow_incr);

        dx= hist_screen.w*rhsize/255.0;
        cairo_translate (cr, hist_screen.w*(1-rhsize)*rhx, (hist_screen.h-100*dx)*rhy);
        cairo_scale (cr, dx, dx);

        cairo_translate (cr, hconf.offset_shadow,hconf.offset_shadow);
        for (i=0; i<=10;i++){
            cairo_move_to (cr, i*255.0/10, 0);
            cairo_rel_line_to (cr, 0, 100.0);
        }
        for (i=0; i<=4;i++){
            cairo_move_to (cr, 0, i*25.0);
            cairo_rel_line_to (cr, 255.0, 0);
        }
        tmp_path=cairo_copy_path(cr);
        cairo_stroke(cr);
        cairo_translate (cr, -hconf.offset_shadow,-hconf.offset_shadow);

        dx=hconf.width_grid;
        dy=hconf.width_hist;
        cairo_device_to_user_distance(cr,&dx,&dy);
        cairo_set_line_width (cr,dx);
        cairo_set_source_rgb (cr, hconf.gray_grid, hconf.gray_grid, hconf.gray_grid);
        cairo_append_path(cr,tmp_path);
        cairo_stroke(cr);
        cairo_path_destroy(tmp_path);

        // drawing of hdata
        if(!block_histplot){
            cairo_set_line_width (cr,dy);
            for(i=0;i<HIST_DATA_N;i++){
                if (hdata[i].active){
                    draw_hdata(hdata[i].v,gminmax,
                        hdata[i].c.r,hdata[i].c.g,hdata[i].c.b,hconf.alpha_hist);
                }
            }
        }
        // now let's overlay the mapping
        if (map) {
            cairo_set_source_rgba (cr, 0.0, 0.0, 0.0,hconf.alpha_shadow);
            dx=hconf.width_map+hconf.width_shadow_incr;
            dy=hconf.width_map;
            cairo_device_to_user_distance(cr,&dx,&dy);
            cairo_set_line_width (cr,dx);
            dx=100.0/255.0;
            cairo_translate (cr, hconf.offset_shadow,hconf.offset_shadow);
            cairo_move_to (cr, 0, 100.0);
            for (i=0; i<256;i++){
                cairo_line_to (cr, i, 100.0-dx*map[i]);
            }
            tmp_path=cairo_copy_path(cr);
            cairo_stroke(cr);
            cairo_translate (cr,-hconf.offset_shadow,-hconf.offset_shadow);
            cairo_set_source_rgb (cr, hconf.gray_map, hconf.gray_map, hconf.gray_map);
            cairo_set_line_width (cr,dy);
            cairo_append_path(cr,tmp_path);
            cairo_stroke(cr);
            cairo_path_destroy(tmp_path);
        }
        return 1;
    }

    urange hist_set_range(int n){
        hdata[n].range={2*hist_screen.size,0};
        return hist_get_range(n);
    }

    urange hist_set_range(){
        int i;
        gminmax= {2*hist_screen.size,0};
        for(i=0;i<HIST_DATA_N;i++){
            hdata[i].range=gminmax;
        }
        return hist_get_range();
    }

    urange hist_get_range(int n){
        unsigned int i;
        if ( hdata[n].range.max == 0) {
            for (i=0; i<256;i++){
                 if (hdata[n].v[i] < hdata[n].range.min) hdata[n].range.min = hdata[n].v[i];
                 if (hdata[n].v[i] > hdata[n].range.max) hdata[n].range.max = hdata[n].v[i];
            }
        }
        return hdata[n].range;
    }

    urange hist_get_range(){
        int i;
        urange r;
        if (!gminmax.max) {
            for (i=0; i<HIST_DATA_N;i++){
                r=hist_get_range(i);
                if (hdata[i].active) gminmax.min=std::min(gminmax.min,r.min);
                if (hdata[i].active) gminmax.max=std::max(gminmax.max,r.max);
            }
        }
        return gminmax;
    }
};

class color_compress0r: public frei0r::filter , hist_plot {
  uint8_t map[256];
  f0r_param_double map_xth;
  f0r_param_double map_yth;
  f0r_param_double map_slope;
  ScreenGeometry fscreen;

  col32bit rgba2hsva(const col32bit);
  col64bit rgba2hsva64(const col32bit);
  col32bit hsva2rgba64(const col64bit);
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


col32bit color_compress0r::rgba2hsva(const col32bit c){

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

col64bit color_compress0r::rgba2hsva64(const col32bit c){

  col64bit r;
  uint32_t iMin,iMax,chroma;
  const uint32_t k1=255 << ABITS;
  const uint32_t k2=HSCALE << ABITS;

  r.a=c.a;

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

col32bit color_compress0r::hsva2rgba64(const col64bit c){
  col32bit r;

  //uint32_t chroma;
  uint32_t m;
  int H,X,ir,ig,ib;
  const uint32_t k1=255<<ABITS;
  const int k2=HSCALE<<ABITS;
  const int k3=1<<(ABITS-1);

  // there's no need to touch alpha, it remains 8-bit
  r.a=c.a;

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

color_compress0r::color_compress0r(int wdt, int hgt) {
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

//color_compress0r::~color_compress0r() {
//  delete planebuf ;
//}

void color_compress0r::update() {
  int i,j,m,n;
  union {
      col32bit c;
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
    for(i=0;i<HIST_DATA_N;i++){
        memset( hdata[i].v, '\0',4*256 );
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
				  1,0);
