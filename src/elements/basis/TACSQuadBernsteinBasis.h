/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef TACS_QUAD_BERNSTEIN_BASIS_H
#define TACS_QUAD_BERNSTEIN_BASIS_H

#include "TACSElementBasis.h"

/**
   Basis class for a quadratic quad element
*/
class TACSQuadraticQuadBernsteinBasis : public TACSElementBasis {
 public:
  ElementLayout getLayoutType();
  void getVisPoint( int n, double pt[] );
  int getNumNodes();
  int getNumParameters();
  int getNumQuadraturePoints();
  double getQuadratureWeight( int n );
  double getQuadraturePoint( int n, double pt[] );
  int getNumElementFaces();
  int getNumFaceQuadraturePoints( int face );
  double getFaceQuadraturePoint( int face, int n, double pt[], double t[] );
  void computeBasis( const double pt[], double N[] );
  void computeBasisGradient( const double pt[], double N[], double Nxi[] );
};

/**
   Basis class for a cubic quad element
*/
class TACSCubicQuadBernsteinBasis : public TACSElementBasis {
 public:
  ElementLayout getLayoutType();
  void getVisPoint( int n, double pt[] );
  int getNumNodes();
  int getNumParameters();
  int getNumQuadraturePoints();
  double getQuadratureWeight( int n );
  double getQuadraturePoint( int n, double pt[] );
  int getNumElementFaces();
  int getNumFaceQuadraturePoints( int face );
  double getFaceQuadraturePoint( int face, int n, double pt[], double t[] );
  void computeBasis( const double pt[], double N[] );
  void computeBasisGradient( const double pt[], double N[], double Nxi[] );
};

/**
   Basis class for a quintic quad element
*/
class TACSQuarticQuadBernsteinBasis : public TACSElementBasis {
 public:
  ElementLayout getLayoutType();
  void getVisPoint( int n, double pt[] );
  int getNumNodes();
  int getNumParameters();
  int getNumQuadraturePoints();
  double getQuadratureWeight( int n );
  double getQuadraturePoint( int n, double pt[] );
  int getNumElementFaces();
  int getNumFaceQuadraturePoints( int face );
  double getFaceQuadraturePoint( int face, int n, double pt[], double t[] );
  void computeBasis( const double pt[], double N[] );
  void computeBasisGradient( const double pt[], double N[], double Nxi[] );
};

/**
   Basis class for a quintic quad element
*/
class TACSQuinticQuadBernsteinBasis : public TACSElementBasis {
 public:
  ElementLayout getLayoutType();
  void getVisPoint( int n, double pt[] );
  int getNumNodes();
  int getNumParameters();
  int getNumQuadraturePoints();
  double getQuadratureWeight( int n );
  double getQuadraturePoint( int n, double pt[] );
  int getNumElementFaces();
  int getNumFaceQuadraturePoints( int face );
  double getFaceQuadraturePoint( int face, int n, double pt[], double t[] );
  void computeBasis( const double pt[], double N[] );
  void computeBasisGradient( const double pt[], double N[], double Nxi[] );
};

#endif // TACS_QUAD_BERNSTEIN_BASIS_H
