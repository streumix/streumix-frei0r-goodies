#include <cairo.h>
#include <arpa/inet.h>

#include "frei0r_common.hpp"

#define HIST_DATA_N 4

typedef struct {
    unsigned int min;
    unsigned int max;
} urange;

typedef struct {
    bool active;
    struct { double b,g,r,a; } c;
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

class hist_plot {
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
    unsigned int hist_data_n;

    void dump_png(const char *s){
        if(surface) cairo_surface_write_to_png (surface,s);
    }

public:
    f0r_param_double rhsize;
    f0r_param_double rhx;
    f0r_param_double rhy;
    hist_plot_config hconf;
    hist_data_channel *hdata;

    hist_plot(){
        hist_plot(HIST_DATA_N);
    }

    hist_plot(int data_channels){
        unsigned int i;
        block_histplot=false;
        cr=NULL;
        surface=NULL;
        hist_data_n=data_channels;
        // init of max data
        gminmax= { 2*hist_screen.size,0};

        hdata = new hist_data_channel[hist_data_n];

        for(i=0;i<hist_data_n;i++){
            hdata[i].range=gminmax;
            hdata[i].active=true;
            hdata[i].c={1,1,1};
        }
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
        //if (cr) cairo_destroy (cr);
        //if (surface) cairo_surface_destroy (surface);
        //if (hdata) delete hdata;
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
            for(i=0;i<hist_data_n;i++){
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
        unsigned int i;
        gminmax= {2*hist_screen.size,0};
        for(i=0;i<hist_data_n;i++){
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
        unsigned int i;
        urange r;
        if (!gminmax.max) {
            for (i=0; i<hist_data_n;i++){
                r=hist_get_range(i);
                if (hdata[i].active) gminmax.min=std::min(gminmax.min,r.min);
                if (hdata[i].active) gminmax.max=std::max(gminmax.max,r.max);
            }
        }
        return gminmax;
    }
};

