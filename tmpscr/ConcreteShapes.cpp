#include "ConcreteShapes.h"

// Circle implementation:
Circle::Circle(double radius) :radius_(radius) {}
double Circle::area() const
{
  double pi = std::acos(-1.0);
  return pi * radius_ * radius_;
}

// Square implementation:
Square::Square(double side) :side_(side) {}
double Square::area() const
{
  return side_ * side_;
}

