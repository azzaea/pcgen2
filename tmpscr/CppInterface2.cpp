#include "ConcreteShapes.h"

// This is just to show we can split up the
// interface into multiple .cpp files.  Interface
// functions that create instances of the Square and
// Circle classes are written here.  Those that call
// the nonmember functions in the Standard C++
// "code base" are located in the CppInterface.cpp file.


// Class Member Function Interfaces:

// Interface to Square member
// function area(.):
// [[Rcpp::export]]
double squareArea(double side)
{
  Square sq(side);
  return sq.area();
}

// Interface to Circle member
// function area(.):
// [[Rcpp::export]]
double circleArea(double radius)
{
  Circle circ(radius);
  return circ.area();
}
