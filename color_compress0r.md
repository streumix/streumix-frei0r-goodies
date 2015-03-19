# Introduction #

A frei0r video effect plugin for compressing RGBA or HSVA colorspace with optional histogram and transfer curve overlay display based on Cairo (scaleable/antialiased).


# Details #

The plugin, called _color\_compress0r_, can do the following:
  * work in either RGBA or HSVA color space
  * colorspace conversion is based on 32-bit integer arithmetic (No floating point!)
  * 16-bit (12-bit effective) resolution is utilized to store temporal HSVA pixel to allow for a loss-less RGBY->HSVA-RGBA conversion chain.
  * selectively map the color channels based on a common compression curve, which is based on two pow() functions and controlled by just 3 parameters: The slope(1) at a specified intersection point x,y (2,3). ( -> smooth, strictly monotonic mapping / flexible and easy to control)
  * the current mapping curve and histrogram data of the color channels (in log scale) can optionally be display as an overlay graphic, which is freely adjustable in position and size. (Current implementaion does not allow real-time video at high fps. Optimizations are still possible, but maybe Cairo is too slow in general. We'll see ...)
  * a XML effect file for nice integration into Kdenlive does exist, too.