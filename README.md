# Kepler's Goat Herd

Code for solving Kepler's equation using contour integration, following Philcox et al. (2021, [arXiv](https://arxiv.org/abs/2103.15829)). This uses a method originally proposed by Ullisch (2020) to solve the "geometric goat problem".

The C++ code contains implementations of a variety of solution methods:
- Newton-Raphson: The quadratic Newton-Raphson root finder.
- Danby: The quartic root finder described in Danby (1988).
- Series: An elliptical series method, as described in Murray & Dermott.
- Contour: A new method based on contour integration.

Given an array of mean anomalies, an eccentricity and a desired precision, the code will estimate the eccentricity using each method. The accuracy of each approach is increased until the desired precision is reached, and timing is performed using the C++ `chrono` package.

To compile the code using g++, simply run ```g++ -o kepler keplers_goat_herd.cpp -std=c++17 -ffast-math -Wall -O3```.  The code can be run using ```./kepler```. The individual functions, e.g. ```compute_contour``` can also be used outside of this script, given an input array of mean anomalies and an eccentricity.

We also provide a pure numpy version of the contour integration function in ```keplers_goat_herd.py```. This is around 9 times slower than the C++ code. Python bindings for the C++ code can be added if these would be of use.

**Authors**:
- Oliver Philcox (Princeton, [ohep2@cantab.ac.uk](mailto:ohep2@cantab.ac.uk))


# Java Version
The java version uses the standard math library.  A FastMath implementation would probably speed things up.  Add a Bessel implementation to test series method.
To compile ```javac keplers_goat_herd.java```.  To run type ```java -classpath . keplers_goat_herd```

# Utility Functions
The util directory contains cpp, python, and java versions of mToE.  This function is very easy to include in your own code in order to use Philcox's method for solving Kepler's equation.
