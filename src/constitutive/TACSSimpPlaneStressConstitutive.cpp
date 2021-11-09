/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#include "TACSSimpPlaneStressConstitutive.h"

const char* TACSSimpPlaneStressConstitutive::psName = "TACSSimpPlaneStressConstitutive";

const char* TACSSimpPlaneStressConstitutive::getObjectName(){
  return psName;
}
/*
  PlaneStressStiffness member function definitions
*/
TACSSimpPlaneStressConstitutive::TACSSimpPlaneStressConstitutive( TACSMaterialProperties *props,
								  TacsScalar _t,
								  int _tNum,
								  TacsScalar _q,
								  TacsScalar _tlb,
								  TacsScalar _tub ):
  TACSPlaneStressConstitutive(NULL){
  properties = props;
  if (properties){
    properties->incref();
  }
  t = _t;
  tNum = _tNum;
  q = _q;
  tlb = _tlb;
  tub = _tub;
}

TACSSimpPlaneStressConstitutive::~TACSSimpPlaneStressConstitutive(){
  if (properties){
    properties->decref();
  }
}


// Retrieve the global design variable numbers
int TACSSimpPlaneStressConstitutive::getDesignVarNums( int elemIndex,
						       int dvLen, int dvNums[] ){
  if (tNum >= 0){
    if (dvNums && dvLen >= 1){
      dvNums[0] = tNum;
    }
    return 1;
  }
  return 0;
}

// Set the element design variable from the design vector
int TACSSimpPlaneStressConstitutive::setDesignVars( int elemIndex,
						    int dvLen,
						    const TacsScalar dvs[] ){
  if (tNum >= 0 && dvLen >= 1){
    t = dvs[0];
    return 1;
  }
  return 0;
}

// Get the element design variables values
int TACSSimpPlaneStressConstitutive::getDesignVars( int elemIndex,
						    int dvLen,
						    TacsScalar dvs[] ){
  if (tNum >= 0 && dvLen >= 1){
    dvs[0] = t;
    return 1;
  }
  return 0;
}

// Get the lower and upper bounds for the design variable values
int TACSSimpPlaneStressConstitutive::getDesignVarRange( int elemIndex,
							int dvLen,
							TacsScalar lb[],
							TacsScalar ub[] ){
  if (tNum >= 0 && dvLen >= 1){
    if (lb){ lb[0] = tlb; }
    if (ub){ ub[0] = tub; }
    return 1;
  }
  return 0;
}

// Evaluate the material density
TacsScalar TACSSimpPlaneStressConstitutive::evalDensity( int elemIndex,
							 const double pt[],
							 const TacsScalar X[] ){
  if (properties){
    //printf("elem: %d %e %e\n", elemIndex, t, properties->getDensity());
    return t*properties->getDensity();
  }
  return 0.0;
}

// Add the derivative of the density
void TACSSimpPlaneStressConstitutive::addDensityDVSens( int elemIndex,
							TacsScalar scale,
							const double pt[],
							const TacsScalar X[],
							int dvLen,
							TacsScalar dfdx[] ){
  if (properties && tNum >= 0){
    dfdx[0] += scale*properties->getDensity();
  }
}

// Evaluate the stresss
void TACSSimpPlaneStressConstitutive::evalStress( int elemIndex,
						  const double pt[],
						  const TacsScalar X[],
						  const TacsScalar e[],
						  TacsScalar s[] ){
  TacsScalar C[6];
  if (properties){
    evalTangentStiffness(elemIndex, pt, X, C);
    s[0] = (C[0]*e[0] + C[1]*e[1] + C[2]*e[2]);
    s[1] = (C[1]*e[0] + C[3]*e[1] + C[4]*e[2]);
    s[2] = (C[2]*e[0] + C[4]*e[1] + C[5]*e[2]);
  }
  else {
    s[0] = s[1] = s[2] = 0.0;
  }
}

// Evaluate the tangent stiffness
void TACSSimpPlaneStressConstitutive::evalTangentStiffness( int elemIndex,
							    const double pt[],
							    const TacsScalar X[],
							    TacsScalar C[] ){
  if (properties){
    properties->evalTangentStiffness2D(C);
    // Get the penalized stiffness value
    auto penalty = pow(t,q);
    C[0] *= penalty;  C[1] *= penalty;  C[2] *= penalty;
    C[3] *= penalty;  C[4] *= penalty;  C[5] *= penalty;
  }
  else {
    C[0] = C[1] = C[2] = C[3] = C[4] = C[5] = 0.0;
  }
}

// Add the derivative of the stress w.r.t. design variables
void TACSSimpPlaneStressConstitutive::addStressDVSens( int elemIndex,
						       TacsScalar scale,
						       const double pt[],
						       const TacsScalar X[],
						       const TacsScalar e[],
						       const TacsScalar psi[],
						       int dvLen,
						       TacsScalar dfdx[] ){
  if (properties && tNum >= 0){
    TacsScalar C[6];
    properties->evalTangentStiffness2D(C);
    // Get the penalized density
    auto dpenalty = q*pow(t, q-1.0);
    TacsScalar s[3];
    s[0] = (C[0]*e[0] + C[1]*e[1] + C[2]*e[2]);
    s[1] = (C[1]*e[0] + C[3]*e[1] + C[4]*e[2]);
    s[2] = (C[2]*e[0] + C[4]*e[1] + C[5]*e[2]);

    // Compute the derivative w.r.t. the design vector
    dfdx[0] += dpenalty*scale*(s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2]);
  }
}
