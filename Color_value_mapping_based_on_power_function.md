# Introduction #
Adjusting color channels values in a common task in image processing. Often a local (de-)compression of the limited value range is required, which is typically achived by so called _color curve_ effects. These kind of effects map the color component value range via an adjustable mapping curve. When is comes to the specification of the transfer functions, different approaches do exist. One of the most flexible ways to define the mapping is based on piecewise defined spline functions. But this flexibility comes with the cost of a typically large, and worst case undefined number of parameters dependent on the number of spline sections.
To overcome this issue, a different approach based on two power functions is presented. It covers a wide range of common application scenarios, whereas it requires just three scalar parameters to define the mapping curve.

# Details #
Since we focus on a local compression or decompression of the normalized component value range [0,1], the start and end point of the range are not going to be shifted. Therefore, we get two fixed points _(x,y)_for the required mapping function _x_ -> _y_ : _(0,0)_ and _(1,1)_. A third point of mapping, which is called the threshold point _(x<sub>th</sub>,y<sub>th</sub>)_, is specified as a parameter. The third parameter of the mapping function represents the slope _s<sub>th</sub>_of the function at the threshold point.
Based on these three input parameters the following mapping curve can be defined.

<wiki:gadget url="http://streumix-frei0r-goodies.googlecode.com/svn/wiki/local.mathml-gadget.xml" border="0" width="100%" height="100ex" up\_content="f(x) = {(0,x <= 0),(y\_(t h)\*(x/x\_(t h))^(s\_(t h)\*y\_(t h)/x\_(t h)),0 <= x <= x\_(t h)),(1-(1-y\_(t h))\*((1-x)/(1-x\_(t h)))^(s\_(t h)\*(1-x\_(t h))/(1-y\_(t h))),x\_(t h) < x <= 1),(1,1 < x ):}"/>


![http://streumix-frei0r-goodies.googlecode.com/svn/wiki/pow-map-image.png](http://streumix-frei0r-goodies.googlecode.com/svn/wiki/pow-map-image.png)