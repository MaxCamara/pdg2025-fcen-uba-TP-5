//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-07 20:22:23 taubin>
//------------------------------------------------------------------------
//
// Optimization.cpp
//
// Software developed for the course
// Digital Geometry Processing
// Copyright (c) 2025, Gabriel Taubin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//     * Redistributions of source code must retain the above
//       copyright notice, this list of conditions and the following
//       disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials
//       provided with the distribution.
//     * Neither the name of the Brown University nor the names of its
//       contributors may be used to endorse or promote products
//       derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GABRIEL
// TAUBIN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
// OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
// USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

// ASSIGNMENT 5 - TODO
//
// get your implementations various functions from Assignment 4
// - search for TODO in this file
//
// complete implementation of
// void Optimization::clusterVerticesApply();

#include <iostream>
#include <iomanip>
#include <cmath>
#include "Optimization.hpp"
#include "Geometry.hpp"
#include "Partition.hpp"
#include "HexGridPartition.hpp"
#include "Heap.hpp"
#include "wrl/IndexedFaceSetVariables.hpp"

//////////////////////////////////////////////////////////////////////
Optimization::Optimization():
  _ifsInput(nullptr),
  _ifsOptimized(nullptr),
  /* _stop(false), */
  _steps(2),
  _lambda(0.5f),
  _mu(0.525f),
  _jacobiWeightSmoothing(0.5f),
  _edgeLengths(),
  _minEdgeLength(0.0f),
  _maxEdgeLength(0.0f),
  _targetEdgeLength(0.0f),
  _clusterMin(),
  _clusterMax(),
  _clusterResolution(0),
  _nClusters(0),
  _vSelIndex(0),
  _eSelIndex(0),
  _fSelIndex(0) {
}

//////////////////////////////////////////////////////////////////////
void Optimization::clear() {
  _ifsInput     = (IndexedFaceSet*)0;
  _ifsOptimized = (IndexedFaceSet*)0;
  _steps = 2;
  _lambda = 0.5f;
  _mu = 0.525f;
  _jacobiWeightSmoothing = 0.5f;
  _edgeLengths.clear();
  _minEdgeLength = 0.0f;
  _maxEdgeLength = 0.0f;
  _targetEdgeLength = 0.0f;
  _nClusters = 0;
  _vSelIndex = 0;
  _eSelIndex = 0;
  _fSelIndex = 0;
}

//////////////////////////////////////////////////////////////////////
IndexedFaceSet* Optimization::getInput() {
  return _ifsInput;
}

//////////////////////////////////////////////////////////////////////
IndexedFaceSet* Optimization::getOptimized() {
  return _ifsOptimized;
}

//////////////////////////////////////////////////////////////////////
void Optimization::setInput(IndexedFaceSet* ifsInput) {
  _ifsInput = ifsInput;
}

//////////////////////////////////////////////////////////////////////
void Optimization::saveOptimized() {
  if(_ifsInput!=(IndexedFaceSet*)0 && _ifsOptimized!=(IndexedFaceSet*)0) {
    (*_ifsInput) = (*_ifsOptimized);
  }
}

void Optimization::setSelectedVertexIndex(const int value) {
  _vSelIndex = (value<0)?0:(value%64);
}

void Optimization::setSelectedEdgeIndex(const int value) {
  _eSelIndex = (value<0)?0:(value%64);
}

void Optimization::setSelectedFaceIndex(const int value) {
  _fSelIndex = (value<0)?0:(value%64);
}

//////////////////////////////////////////////////////////////////////
void Optimization::setOptimized(IndexedFaceSet* ifsOptimized, const bool reset) {
  _ifsOptimized = ifsOptimized;
  if(_ifsOptimized!=(IndexedFaceSet*)0 && _ifsInput!=(IndexedFaceSet*)0) {

    vector<int>&   coordIndex      = _ifsOptimized->getCoordIndex();
    vector<float>& coord           = _ifsOptimized->getCoord();

    if(reset || coord.size()==0) {
      _ifsOptimized->clear();

      vector<int>&   coordIndexInput = _ifsInput->getCoordIndex();
      vector<float>& coordInput      = _ifsInput->getCoord();
        
      coordIndex.insert
        (coordIndex.end(), coordIndexInput.begin(),coordIndexInput.end());
      coord.insert
        (coord.end(), coordInput.begin(),coordInput.end());
    }

    if(Geometry::isTriangulated(coordIndex)==false) {
      Geometry::triangulate(coord,coordIndex);
      if(_ifsOptimized->hasColorPerFace() || _ifsOptimized->hasColorPerCorner())
        _ifsOptimized->clearColor();
      if(_ifsOptimized->hasNormalPerFace() || _ifsOptimized->hasNormalPerCorner())
        _ifsOptimized->clearNormal();
      if(_ifsOptimized->hasTexCoordPerCorner())
        _ifsOptimized->clearTexCoord();
    }

    if(_ifsOptimized->hasNormalPerFace()==false) {
      _ifsOptimized->clearNormal();
      _ifsOptimized->setNormalPerVertex(false); // default value is TRUE
      Geometry::computeNormalsPerFace
        (coord,coordIndex,_ifsOptimized->getNormal());
    }

    IndexedFaceSetVariables ifsv(*_ifsOptimized);
    PolygonMesh* pmesh = ifsv.getPolygonMesh(true);

    Geometry::computeEdgeLengths(coord,*pmesh,_edgeLengths);

    int nE = pmesh->getNumberOfEdges();

    _minEdgeLength = 0.0f;
    _maxEdgeLength = 0.0f;
    for(int iE=0;iE<nE;iE++) {
      float eLength = _edgeLengths[iE];
      if(iE==0 || eLength<_minEdgeLength) _minEdgeLength = eLength;
      if(iE==0 || eLength>_maxEdgeLength) _maxEdgeLength = eLength;
    }
    _targetEdgeLength = (_minEdgeLength+_maxEdgeLength)/2.0f;
  }
}

//////////////////////////////////////////////////////////////////////
int Optimization::getSteps() {
  return _steps;
}

void Optimization::setSteps(const int value) {
  _steps = (value>0)?value:0;
}

float Optimization::getLambda() {
  return _lambda;
}

void  Optimization::setLambda(const float value) {
  _lambda = value;
}

float Optimization::getMu() {
  return _mu;
}

void  Optimization::setMu(const float value) {
  _mu = value;
}

// kappa : pass-band frequency
float Optimization::getKappa() {
  float kappa = (_lambda+_mu)/(_lambda*_mu);
  return kappa;
}

// kappa : pass-band frequency
void Optimization::setKappa(const float value) {
  float kappa = (0.0f<value)?value:0.001;
  if(_lambda*kappa>=1.0f) _lambda = (1.0f/kappa)-0.001f;
  _mu = _lambda/(_lambda*kappa-1.0f);
}

//////////////////////////////////////////////////////////////////////
float Optimization::getJacobiWeightData() {
  return 1.0f-_jacobiWeightSmoothing;
}

//////////////////////////////////////////////////////////////////////
void  Optimization::setJacobiWeightData(const float value) {
  _jacobiWeightSmoothing = (value<0.0f)?1.0f:(value>1.0f)?0.0f:(1.0f-value);
}

//////////////////////////////////////////////////////////////////////
float Optimization::getJacobiWeightSmoothing() {
  return _jacobiWeightSmoothing;
}

//////////////////////////////////////////////////////////////////////
void  Optimization::setJacobiWeightSmoothing(const float value) {
  _jacobiWeightSmoothing = (value<0.0f)?0.0f:(value>1.0f)?1.0f:value;
}

//////////////////////////////////////////////////////////////////////
void Optimization::laplacianSmoothingVertexCoordinatesRun() {

    // Jacobi iteration for energy function
    // E(x) = \sum_{0<=iE<nE} {1} \| x_{iV0}-x_{iV1}^0\|^2

    std::cout << "void Optimization::laplacianSmoothingVertexCoordinatesRun() {\n";
    std::cout << "  lambda = " << _lambda << "\n";
    std::cout << "  mu     = " << _mu     << "\n";
    std::cout << "  steps  = " << _steps  << "\n";

    IndexedFaceSetVariables ifsv(*_ifsOptimized);
    PolygonMesh* pmesh = ifsv.getPolygonMesh(true);
    int nV = pmesh->getNumberOfVertices();
    int nE = pmesh->getNumberOfEdges();

    // this is the signal that we want to smooth over the primal graph
    vector<float>& x = _ifsOptimized->getCoord();

    // allocate local arrays to accumulate displacements and weights
    vector<float> dx(3*nV,0.0f);
    vector<float> wx(nV,0.0f);

    for(int step=0;step<_steps;step++) {

        // zero accumulators dx[] and wx[]
        for (int iV=0; iV<nV; iV++) {
            wx[iV] = 0;
            dx[iV*3] = 0;
            dx[iV*3 + 1] = 0;
            dx[iV*3 + 2] = 0;
        }

        // accumulate displacement vectors and weights
        //
        // for each edge iE (iV0,iV1) {
        //   let dx_iV0_iV1  be the displacement vector from x_iV0 to x_iV1
        //   let  w_iE>0 be the weight assigned to the edge (such as 1.0)
        //   add  dx_iV0_iV1 multiplied by w_iE to the accumulator dx_iV0
        //   add  w_iE to the accumulator wx_iV0
        //   add -dx_iV0_iV1 multiplied by w_iE to the accumulator dx_iV1
        //   add  w_iE to the accumulator wx_iV1
        // }
        for (int iE=0; iE<nE; iE++) {
            int iV0 = pmesh->getVertex0(iE);
            int iV1 = pmesh->getVertex1(iE);
            vector<float> dx_iV0_iV1 = {x[iV1*3]-x[iV0*3], x[iV1*3+1]-x[iV0*3+1], x[iV1*3+2]-x[iV0*3+2]};
            float w_iE = 1.0;

            dx[iV0*3] += w_iE * dx_iV0_iV1[0];
            dx[iV0*3+1] += w_iE * dx_iV0_iV1[1];
            dx[iV0*3+2] += w_iE * dx_iV0_iV1[2];
            wx[iV0] += w_iE;

            dx[iV1*3] += w_iE * -dx_iV0_iV1[0];
            dx[iV1*3+1] += w_iE * -dx_iV0_iV1[1];
            dx[iV1*3+2] += w_iE * -dx_iV0_iV1[2];
            wx[iV1] += w_iE;
        }

        // normalize the displacement vectors
        // (make sure that you do not divide by zero!)
        //
        // for each vertex iV {
        //   dx_iV /= wx_iV
        // }
        for (int iV=0; iV<nV; iV++) {
            if (wx[iV]!=0) {
                dx[iV*3] /= wx[iV];
                dx[iV*3+1] /= wx[iV];
                dx[iV*3+2] /= wx[iV];
            }
        }

        // alternate between _lambda and _mu for even and odd step values
        float lambda = (step%2==0)?_lambda:_mu;

        // apply displacement
        //
        // for each vertex iV {
        //   x_iV += lambda * dx_iV
        // }
        for (int iV=0; iV<nV; iV++) {
            x[iV*3] += lambda * dx[iV*3];
            x[iV*3+1] += lambda * dx[iV*3+1];
            x[iV*3+2] += lambda * dx[iV*3+2];
        }

    }

    // since the vertex coordinates have changed, recompute face normals
    _ifsOptimized->clearNormal();
    _ifsOptimized->setNormalPerVertex(false);
    vector<float>& coord      = _ifsOptimized->getCoord(); // == x
    vector<float>& normal     = _ifsOptimized->getNormal();
    vector<int>&   coordIndex = _ifsOptimized->getCoordIndex();
    Geometry::computeNormalsPerFace(coord,coordIndex,normal);

    // Remember that none of the smoothing methods modify the
    // connectivity of the mesh
    // - as a result, if the mesh had any other properties, they are
    //   still valid after smoothing
    // - same thing with respect to selections

    // it would be useful to preserve the face colors set by some of the
    // remeshing methods
    // _ifsOptimized->clearColor();

    // it is very unlikely that any of the meshes that we use will have
    // texture coordinates
    // _ifsOptimized->clearTexCoord();

    // it would be useful to preserve the selections set by some of the
    // remeshing methods
    // ifsv.clearAllSelection();

    //Las longitudes de las aristas cambian
    Geometry::computeEdgeLengths(coord,*pmesh,_edgeLengths);

    std::cout << "}\n";
}

//////////////////////////////////////////////////////////////////////
void Optimization::laplacianSmoothingFaceNormalsRun(const bool normalize) {

    // Jacobi iteration for energy function
    // E(n) = \sum_{0<=iE<nED} {1} \| n_{iF0}-n_{iF1}^0\|^2

    std::cout << "void Optimization::laplacianSmoothingFaceNormalsRun() {\n";
    std::cout << "  normalize         = " << ((normalize)?"true":"false") << "\n";
    std::cout << "  lambda            = " << _lambda << "\n";
    std::cout << "  mu                = " << _mu     << "\n";
    std::cout << "  steps             = " << _steps  << "\n";

    IndexedFaceSetVariables ifsv(*_ifsOptimized);
    PolygonMesh* pmesh = ifsv.getPolygonMesh(true);
    int nF = pmesh->getNumberOfFaces();
    int nE = pmesh->getNumberOfEdges();

    // we need normals per face not indexed
    if(_ifsOptimized->getNormalBinding()!=IndexedFaceSet::Binding::PB_PER_FACE) {
        // recompute face normals
        _ifsOptimized->clearNormal();
        _ifsOptimized->setNormalPerVertex(false);
        vector<float>& coord      = _ifsOptimized->getCoord();
        vector<float>& normal     = _ifsOptimized->getNormal();
        vector<int>&   coordIndex = _ifsOptimized->getCoordIndex();
        Geometry::computeNormalsPerFace(coord,coordIndex,normal);
    }

    // this is the signal that we want to smooth over the dual graph
    vector<float>& n = _ifsOptimized->getNormal();

    // allocate local arrays to accumulate the normal displacements and weights
    vector<float> dn(3*nF,0.0f);
    vector<float> wn(nF,0.0f);

    for(int step=0;step<_steps;step++) {

        // zero accumulators dn[] and wn
        for (int iF=0; iF<nF; iF++) {
            wn[iF] = 0;
            dn[iF*3] = 0;
            dn[iF*3 + 1] = 0;
            dn[iF*3 + 2] = 0;
        }

        // accumulate displacement vectors
        //
        // for each regular edge iE {
        //   find the two faces iF0 and iF1 incident to the edge
        //   let dn_iF0_iF1  be the displacement vector from n_iF0 to n_iF1
        //   let w_iE>0 be the weight assigned to the edge (such as 1.0)
        //   add  dn_iF0_iF1 multiplied by w_iE to the accumulator dn_iF0
        //   add  w_iE to the accumulator wn_iF0
        //   add -dx_iF0_iF1 multiplied by w_iE to the accumulator dn_iF1
        //   add  w_iE to the accumulator wn_iF1
        // }
        for (int iE=0; iE<nE; iE++) {
            if (!pmesh->isRegularEdge(iE)) continue;
            int iF0 = pmesh->getEdgeFace(iE, 0);
            int iF1 = pmesh->getEdgeFace(iE, 1);
            vector<float> dn_iF0_iF1 = {n[iF1*3]-n[iF0*3], n[iF1*3+1]-n[iF0*3+1], n[iF1*3+2]-n[iF0*3+2]};
            float w_iE = 1.0;

            dn[iF0*3] += w_iE * dn_iF0_iF1[0];
            dn[iF0*3+1] += w_iE * dn_iF0_iF1[1];
            dn[iF0*3+2] += w_iE * dn_iF0_iF1[2];
            wn[iF0] += w_iE;

            dn[iF1*3] += w_iE * -dn_iF0_iF1[0];
            dn[iF1*3+1] += w_iE * -dn_iF0_iF1[1];
            dn[iF1*3+2] += w_iE * -dn_iF0_iF1[2];
            wn[iF1] += w_iE;
        }

        // normalize the displacement vectors
        // (make sure that you do not divide by zero!)
        //
        // for each face iF {
        //   dn_iF /= wn_iF
        // }
        for (int iF=0; iF<nF; iF++) {
            if (wn[iF]!=0) {
                dn[iF*3] /= wn[iF];
                dn[iF*3+1] /= wn[iF];
                dn[iF*3+2] /= wn[iF];
            }
        }

        // alternate between _lambda and _mu for even and odd step values
        float lambda = (step%2==0)?_lambda:_mu;

        // apply displacement
        //
        // for each face iF {
        //   n_iF += lambda * dn_iF
        // }
        for (int iF=0; iF<nF; iF++) {
            n[iF*3] += lambda * dn[iF*3];
            n[iF*3+1] += lambda * dn[iF*3+1];
            n[iF*3+2] += lambda * dn[iF*3+2];
        }

    }

    if(normalize) {
        // normalize face normals to unit length
        for (int iF=0; iF<nF; iF++) {
            float normalLength = sqrt(pow(n[iF*3], 2) + pow(n[iF*3+1], 2) + pow(n[iF*3+2], 2));
            n[iF*3] /= normalLength;
            n[iF*3+1] /= normalLength;
            n[iF*3+2] /= normalLength;
        }
    }

    // clear colors (or maybe not?)
    _ifsOptimized->clearColor();

    // clear all selection buffers (or maybe not?)
    ifsv.clearAllSelection();

    std::cout << "}\n";
}

//////////////////////////////////////////////////////////////////////
void  Optimization::jacobiRun() {

    // This method is very similar to
    // Optimization::laplacianSmoothingVertexCoordinatesRun()

    // Jacobi iteration for energy function
    // E(x) = (1-sigma)*\sum_{0<=iV<nV} \| x_{iV}-x0_{iV}\|^2 +
    //        (  sigma)*\sum_{0<=iE<nE} \| x_{iV0}-x_{iV1}\|^2

    std::cout << "void Optimization::jacobiSmoothingRun() {\n";
    std::cout << "  lambda = " << _lambda << "\n";
    std::cout << "  mu     = " << _mu     << "\n";
    std::cout << "  steps  = " << _steps  << "\n";
    std::cout << "  sigma  = " << _jacobiWeightSmoothing  << "\n";

    IndexedFaceSetVariables ifsv(*_ifsOptimized);
    PolygonMesh* pmesh = ifsv.getPolygonMesh(true);
    int nV = pmesh->getNumberOfVertices();
    int nE = pmesh->getNumberOfEdges();

    // these are the target coordinates
    vector<float>& x0 = _ifsInput->getCoord();

    // this is the signal that we want to smooth over the primal graph
    vector<float>& x = _ifsOptimized->getCoord();

    // allocate local arrays to accumulate displacements and weights
    vector<float> dx(3*nV,0.0f);
    vector<float> wx(nV,0.0f);

    float sigma = _jacobiWeightSmoothing;
    // setJacobiWeightSmoothing(const float value);
    // constraints this value to 0<=sigma<=1
    // note that if sigma==1.0f then this method produces the same result as
    // Optimization::laplacianSmoothingVertexCoordinatesRun()
    // and if sigma==0 the regularization term is ignored

    for(int step=0;step<_steps;step++) {

        // zero accumulators dx[] and wx[]
        for (int iV=0; iV<nV; iV++) {
            wx[iV] = 0;
            dx[iV*3] = 0;
            dx[iV*3 + 1] = 0;
            dx[iV*3 + 2] = 0;
        }

        // accumulate displacement vectors and weights
        //
        // 1) accumulate the contribution of the data term
        //
        // w = 1.0f-sigma;
        // for each vertex iV {
        //   let dx_iV be the displacement vector from x_iV to x0_iV
        //   add dx_iV multiplied by w to the accumulator dx_iV
        //   add the weight w to the accumulator wx_iV
        //
        // }
        float w = 1.0 - sigma;
        for (int iV=0; iV<nV; iV++) {
            vector<float> dx_iV = {x0[iV*3]-x[iV*3], x0[iV*3+1]-x[iV*3+1], x0[iV*3+2]-x[iV*3+2]};
            dx[iV*3] += w * dx_iV[0];
            dx[iV*3+1] += w * dx_iV[1];
            dx[iV*3+2] += w * dx_iV[2];
            wx[iV] += w;
        }

        // 2) accumulate the contribution of the regularization term
        //    - note that this part can be copied from
        //      Optimization::laplacianSmoothingVertexCoordinatesRun()
        //
        // w = sigma;
        // for each edge iE (iV0,iV1) {
        //   let dx_iV0_iV1  be the displacement vector from x_iV0 to x_iV1
        //   add  dx_iV0_iV1 multiplied by w to the accumulator dx_iV0
        //   add  w to the accumulator wx_iV0
        //   add -dx_iV0_iV1 multiplied by w to the accumulator dx_iV1
        //   add  w to the accumulator wx_iV1
        // }
        // normalize the displacement vectors
        // (make sure that you do not divide by zero!)
        //
        // for each vertex iV {
        //   dx_iV /= wx_iV
        // }
        w = sigma;
        for (int iE=0; iE<nE; iE++) {
            int iV0 = pmesh->getVertex0(iE);
            int iV1 = pmesh->getVertex1(iE);
            vector<float> dx_iV0_iV1 = {x[iV1*3]-x[iV0*3], x[iV1*3+1]-x[iV0*3+1], x[iV1*3+2]-x[iV0*3+2]};

            dx[iV0*3] += w * dx_iV0_iV1[0];
            dx[iV0*3+1] += w * dx_iV0_iV1[1];
            dx[iV0*3+2] += w * dx_iV0_iV1[2];
            wx[iV0] += w;

            dx[iV1*3] += w * -dx_iV0_iV1[0];
            dx[iV1*3+1] += w * -dx_iV0_iV1[1];
            dx[iV1*3+2] += w * -dx_iV0_iV1[2];
            wx[iV1] += w;
        }
        for (int iV=0; iV<nV; iV++) {
            if (wx[iV]!=0) {
                dx[iV*3] /= wx[iV];
                dx[iV*3+1] /= wx[iV];
                dx[iV*3+2] /= wx[iV];
            }
        }

        // alternate between _lambda and _mu for even and odd step values
        float lambda = (step%2==0)?_lambda:_mu;

        // apply displacement
        //
        // for each vertex iV {
        //   x_iV += lambda * dx_iV
        // }
        for (int iV=0; iV<nV; iV++) {
            x[iV*3] += lambda * dx[iV*3];
            x[iV*3+1] += lambda * dx[iV*3+1];
            x[iV*3+2] += lambda * dx[iV*3+2];
        }
    }

    // since the vertex coordinates have changed, recompute face normals
    _ifsOptimized->clearNormal();
    _ifsOptimized->setNormalPerVertex(false);
    vector<float>& coord      = _ifsOptimized->getCoord();
    vector<float>& normal     = _ifsOptimized->getNormal();
    vector<int>&   coordIndex = _ifsOptimized->getCoordIndex();
    Geometry::computeNormalsPerFace(coord,coordIndex,normal);

    // clear colors ???
    _ifsOptimized->clearColor();

    // clear all selection buffers ???
    ifsv.clearAllSelection();

    //Las longitudes de las aristas cambian
    Geometry::computeEdgeLengths(coord,*pmesh,_edgeLengths);

    std::cout << "}\n";
}

vector<float>& Optimization::getEdgeLengths() {
  return _edgeLengths;
}

float Optimization::getMinEdgeLength() {
  return _minEdgeLength;
}

float Optimization::getMaxEdgeLength() {
  return _maxEdgeLength;
}

float Optimization::getTargetEdgeLength() {
  return _targetEdgeLength;
}

void Optimization::setMinEdgeLength(const float value) {
  _minEdgeLength = (value<0.0f)?0.0f:value;
}

void Optimization::setMaxEdgeLength(const float value) {
  _maxEdgeLength = (value<_minEdgeLength)?_minEdgeLength:value;
}

void Optimization::setTargetEdgeLength(const float value) {
  _targetEdgeLength = 
    (value<_minEdgeLength)?_minEdgeLength:
    (value>_maxEdgeLength)?_maxEdgeLength:
    value;
}

int Optimization::getQuantizationResolution() {
  return _clusterResolution;
}

void Optimization::setQuantizationResolution
(const int resolution) {
  _clusterResolution = resolution;
}

void Optimization::setQuantizationBox
(const Vec3f& min, const Vec3f& max) {
  _clusterMin        = min;
  _clusterMax        = max;
}

int Optimization::getNumberOfClusters() {
  return _nClusters;
}

void Optimization::clusterVerticesApply() {

  if(_clusterResolution<=0) return;

  IndexedFaceSetVariables ifsv(*_ifsOptimized);
  PolygonMesh* pmesh = ifsv.getPolygonMesh(true);

  // int nV = pmesh->getNumberOfVertices();
  // int nE = pmesh->getNumberOfEdges();
  int nF = pmesh->getNumberOfFaces();
  int nC = pmesh->getNumberOfCorners();

  vector<int>&   coordIndex  = _ifsOptimized->getCoordIndex();
  vector<float>& coord       = _ifsOptimized->getCoord();

  HexGridPartition hgp(_clusterMin,_clusterMax,_clusterResolution);
  hgp.insertPoints(coord);

  _nClusters = hgp.getNumberOfCells();

  // sample the point cloud
  vector<float> newCoord;
  vector<int>   vertexMap;
  hgp.sample(newCoord,&vertexMap);
  // int nPoints  = static_cast<int>(coord.size()/3)
  // int nSamples = static_cast<int>(newCoord.size()/3)
  // vertexMap maps [0:nPoints) onto [0:nSamples)
  
  int iF,iC0,iC1;

  // determine which faces have to be deleted, if any, because edges
  // have collapsed
  //
  // int nFacesDeleted = 0;
  vector<bool> deleteFace(nF,false);

  // TODO
  for(iF=iC0=iC1=0;iC1<nC;iC1++) {
    if(coordIndex[iC1]>=0) continue;
    int iV0, iV1;
    for(int iC=iC0; iC<iC1; iC++){
      iV0 = pmesh->getSrc(iC);
      if(iC==iC0){
        iV1 = pmesh->getSrc(iC1);
      } else {
        iV1 = pmesh->getSrc(iC-1);
      }
      //Si la arista fue colapsada (es decir, sus dos vértices se mapean al mismo vértice nuevo), la cara se borra
      if(vertexMap[iV0]==vertexMap[iV1]){
        deleteFace[iF] = true;
        break;
      }
    }

    iC0 = iC1+1; iF++;
  }

  // generate output coordIndex array
  vector<int> newCoordIndex;
  for(iF=iC0=iC1=0;iC1<nC;iC1++) {
    if(coordIndex[iC1]>=0) continue;
    // mesh is triangulated
    // assert(iC1-iC0==3);
    if(deleteFace[iF]==false) {
      
      // TODO

      // - for each corner of the face in coordIndex, get the input vertex index
      // - map the input vertex input onto a sample index
      // - push back the sample index onto the newCoordIndex array
      // - don't forget to push back a -1
      for(int iC=iC0; iC<iC1; iC++){
        int iV = pmesh->getSrc(iC);
        int iS = vertexMap[iV];
        newCoordIndex.push_back(iS);
      }
      newCoordIndex.push_back(-1);
    }
    iC0 = iC1+1; iF++;
  }

  // fix _ifsRendering in place
  coord.clear();
  coord.insert(coord.end(),newCoord.begin(),newCoord.end());
  coordIndex.clear();
  coordIndex.insert(coordIndex.end(),newCoordIndex.begin(),newCoordIndex.end());
  // rebuild PolygonMesh
  ifsv.deletePolygonMesh();
  // pmesh = ifsv.getPolygonMesh(true);

  // clear colors and normals
  _ifsOptimized->clearColor();
  _ifsOptimized->clearNormal();

  // recompute face normals
  _ifsOptimized->setNormalPerVertex(false);
  Geometry::computeNormalsPerFace(coord,coordIndex,_ifsOptimized->getNormal());

  // clear all selection buffers
  ifsv.clearAllSelection();
}

void Optimization::clusterVerticesApply
(const Vec3f& min, const Vec3f& max, const int resolution) {
  setQuantizationBox(min,max);
  setQuantizationResolution(resolution);
  clusterVerticesApply();
}

void Optimization::_collapseEdgesSelect
(const EdgeCollapseIndependentSet indepSet,
 vector<int>&                     edgeSelection,
 bool                             colorIncidentFaces) {

    int   i,iV,iF,iE,iV0,iV1,iV2,iV3,nEF;
    int   iC,iC0,iC1,iC0n,iC0p,iC1n,iC1p;

    IndexedFaceSetVariables ifsv(*_ifsOptimized);
    PolygonMesh* pmesh = ifsv.getPolygonMesh(true);

    int nV = pmesh->getNumberOfVertices();
    int nE = pmesh->getNumberOfEdges();
    int nF = pmesh->getNumberOfFaces();

    // clear output edge selection array
    edgeSelection.clear();
    edgeSelection.resize(nE,-1); // no edge is initially selected

    // color faces for visualization purposes
    vector<int>*   colorIndex = (vector<int>*)0;
    vector<float>* color      = (vector<float>*)0;
    int noChangesColorIndex     = 0;
    int collapseEdgesColorIndex = 1;
    if(colorIncidentFaces) {

        // set _ifsOptimized
        // as color per face indexed with two colors
        colorIndex = &(_ifsOptimized->getColorIndex());
        colorIndex->clear();
        color      = &(_ifsOptimized->getColor());
        color->clear();
        _ifsOptimized->setColorPerVertex(false);

        // color 0
        Color noChangesColor(0.8f,0.8f,0.8f);
        noChangesColorIndex =
            static_cast<int>(color->size()/3); // ==0
        color->push_back(noChangesColor.r);
        color->push_back(noChangesColor.g);
        color->push_back(noChangesColor.b);

        // color 1
        Color collapseEdgesColor(0.9f,0.5f,0.7f);
        collapseEdgesColorIndex =
            static_cast<int>(color->size()/3); // ==1
        color->push_back(collapseEdgesColor.r);
        color->push_back(collapseEdgesColor.g);
        color->push_back(collapseEdgesColor.b);

        // initially paint all the faces with color noChangesColorIndex
        for(int iF=0;iF<nF;iF++)
            colorIndex->push_back(noChangesColorIndex);
    }

    // compute error metric for all edges and insert them into a
    Heap heap;
    for(iE=0;iE<nE;iE++) {
        heap.add(_edgeLengths[iE],iE);
    }

    // creating an independent set ef edges to split ...

    // initially mark all the vertices as not used
    vector<bool> usedVertex(nV,false);

    while(heap.length()>0) {
        /* iH = */  heap.delMin();
        iE      =  heap.getLastIKey();

        // you may want to set up a threshold so that only edges shorter
        // than the threshold are collapsed
        //
        // float eLength = heap.getLastFKey();
        // if(eLength>=_lowEdgeLength) break;
        float eLength = heap.getLastFKey();
        if(eLength>=_targetEdgeLength) break;

        // if edge is not regular, continue
        if (!pmesh->isRegularEdge(iE)) continue;

        // if either one of the two vertices iV0 and iV1 of the edge iE
        // are marked as used, continue
        int iV0 = pmesh->getVertex0(iE);
        int iV1 = pmesh->getVertex1(iE);
        if (usedVertex[iV0] || usedVertex[iV1]) continue;

        // mark the edge as selected
        edgeSelection[iE] = _eSelIndex;

        // mark the two vertices, and depending on the value of indepSet,
        // some neighboring vertices as used
        //
        // do you need to take relative orientation into account?
        //
        //        iV3                      iV3
        //     /       \                /       \
        //    / <-iC1-- \              / --iC1-> \
        // iV0 --------- iV1   or   iV0 --------- iV1
        //    \ --iC0-> /              \ --iC0-> /
        //     \       /                \       /
        //        iV2                      iV2
        //
        // mark the x vertices as used
        switch(indepSet) {
        case EdgeCollapseIndependentSet::VERTICES_8:
            //     x - x - x
            //      \ / \ /
            // iV0-> x - x <- iV1
            //      / \ / \
            //     x - x - x
            usedVertex[iV0] = true;
            usedVertex[iV1] = true;

            iC0 = pmesh->getEdgeHalfEdge(iE, 0);
            iC0p = pmesh->getPrev(iC0);
            iC0n = pmesh->getNext(iC0);
            iC1 = pmesh->getEdgeHalfEdge(iE, 1);
            iC1p = pmesh->getPrev(iC1);
            iC1n = pmesh->getNext(iC1);
            iV2 = pmesh->getDst(iC0n);
            iV3 = pmesh->getDst(iC1n);
            usedVertex[iV2] = true;
            usedVertex[iV3] = true;

            //Si las aristas son regulares, los twins de iC0p, iC0n, iC1p e iC1n se encuentran en las caras que contienen los otros vértices
            iC = pmesh->getTwin(iC0p);
            if (iC>=0) {
                //Como iC es el half-edge twin, el half-edge que lo precede en la cara es el que se origina en el vértice que buscamos
                iC = pmesh->getPrev(iC);
                iV = pmesh->getSrc(iC);
                usedVertex[iV] = true;
            }
            iC = pmesh->getTwin(iC0n);
            if (iC>=0) {
                iC = pmesh->getPrev(iC);
                iV = pmesh->getSrc(iC);
                usedVertex[iV] = true;
            }
            iC = pmesh->getTwin(iC1p);
            if (iC>=0) {
                iC = pmesh->getPrev(iC);
                iV = pmesh->getSrc(iC);
                usedVertex[iV] = true;
            }
            iC = pmesh->getTwin(iC1n);
            if (iC>=0) {
                iC = pmesh->getPrev(iC);
                iV = pmesh->getSrc(iC);
                usedVertex[iV] = true;
            }
            break;
        case EdgeCollapseIndependentSet::VERTICES_4:
            //          x
            //         / \
            // iV0 -> x - x <- iV1
            //         \ /
            //          x
            usedVertex[iV0] = true;
            usedVertex[iV1] = true;

            iC0 = pmesh->getEdgeHalfEdge(iE, 0);
            iC1 = pmesh->getEdgeHalfEdge(iE, 1);
            iC0n = pmesh->getNext(iC0);
            iC1n = pmesh->getNext(iC1);
            iV2 = pmesh->getDst(iC0n);
            iV3 = pmesh->getDst(iC1n);
            usedVertex[iV2] = true;
            usedVertex[iV3] = true;
            break;
        case EdgeCollapseIndependentSet::VERTICES_2:
            //          o
            //         / \
            // iV0 -> x - x <- iV1
            //         \ /
            //          o
            usedVertex[iV0] = true;
            usedVertex[iV1] = true;
            break;
        }

        // for visualization purposes
        // set the colorIndex of the two faces incident to the edge to
        // collapseEdgesColorIndex;
        nEF = pmesh->getNumberOfEdgeFaces(iE);
        for(i=0;i<nEF;i++) {
            iF = pmesh->getEdgeFace(iE,i);
            if(colorIncidentFaces) {
                (*colorIndex)[iF] = collapseEdgesColorIndex;
            }
        }

    }
}

void Optimization::collapseEdgesShow
(const EdgeCollapseIndependentSet indepSet) {
    if(_ifsOptimized==nullptr) return;
  IndexedFaceSetVariables ifsv(*_ifsOptimized);
  vector<int>& edgeSelection = ifsv.getEdgeSelection();
  _collapseEdgesSelect(indepSet,edgeSelection,true);
}

void Optimization::collapseEdgesApply
(const EdgeCollapseIndependentSet indepSet) {
    if(_ifsOptimized==nullptr) return;

    // select an independent set of edges
    //
    // use the _ifsOptimized edgeSelection array instead ?
    // vector<int>& edgeSelection = _ifsOptimized->getEdgeSelection();
    vector<int> edgeSelection;
    _collapseEdgesSelect(indepSet,edgeSelection,false);

    IndexedFaceSetVariables ifsv(*_ifsOptimized);
    PolygonMesh* pmesh = ifsv.getPolygonMesh(true);
    vector<float>& coord      = _ifsOptimized->getCoord();
    vector<int>&   coordIndex = _ifsOptimized->getCoordIndex();

    int nV = pmesh->getNumberOfVertices();
    int nE = pmesh->getNumberOfEdges();
    int nF = pmesh->getNumberOfFaces();
    int nC = pmesh->getNumberOfCorners();

    // map old vertices to new vertices and create new coord array
    vector<int>    vertexMap(nV,-1);
    vector<float>  newCoord;

    // iVnew = 0
    //
    // assign a new vertex index to each collapsed edge
    //
    // for each selected edge iE {
    //    assign new index iVnew to iV0 and iV1
    //    compute newCoord_iVnew as the midpoint of coord_iV0 and coord_iV1
    //    if(errorMetric==EdgeCollapseErrorMetric::GARLAND_HECKBERT) {
    //      save old values of qVMatrix4x4
    //    }
    //    increment iVnew
    // }
    int iVnew = 0;
    vector<bool> selectedEdgeVertex(nV, false);
    vector<bool> selectedEdgeFace(nF, false);

    for (int iE=0; iE<nE; iE++) {
        if (edgeSelection[iE] != _eSelIndex) continue;

        int iV0 = pmesh->getVertex0(iE);
        int iV1 = pmesh->getVertex1(iE);

        selectedEdgeVertex[iV0] = true;
        selectedEdgeVertex[iV1] = true;
        int incidentFaces = pmesh->getNumberOfEdgeFaces(iE);
        for (int i=0; i<incidentFaces; i++) {
            int iF = pmesh->getEdgeFace(iE, i);
            selectedEdgeFace[iF] = true;
        }

        vertexMap[iV0] = iVnew;
        vertexMap[iV1] = iVnew;

        vector<float> coord_iV0 = {coord[iV0*3], coord[iV0*3+1], coord[iV0*3+2]};
        vector<float> coord_iV1 = {coord[iV1*3], coord[iV1*3+1], coord[iV1*3+2]};
        vector<float> coord_iVnew = {
            (coord_iV0[0] + coord_iV1[0]) / 2,
            (coord_iV0[1] + coord_iV1[1]) / 2,
            (coord_iV0[2] + coord_iV1[2]) / 2
        };
        newCoord.push_back(coord_iVnew[0]);
        newCoord.push_back(coord_iVnew[1]);
        newCoord.push_back(coord_iVnew[2]);

        iVnew++;
    }

    //
    // assign new vertex indices to vertices which are not ends of
    // collapsed edges
    //
    // for each vertex iV {
    //   if vertex iV is not end of a selected edge {
    //      assign new index iVnew to iV
    //      copy coord_iV onto newCoord_iVnew
    //      if(errorMetric==EdgeCollapseErrorMetric::GARLAND_HECKBERT) {
    //        use values stored in newQMatrix4x4 and _errorGarlandHeckbert()
    //        to generate new matrix for the new vertex and to compute
    //        coordinates of the new vertex
    //      } else {
    //        set coordinates of the new vertex as edge midpoint
    //      }
    //      increment iVnew
    //   }
    // }
    for (int iV=0; iV<nV; iV++) {
        if (!selectedEdgeVertex[iV]) {
            vertexMap[iV] = iVnew;
            newCoord.push_back(coord[iV*3]);
            newCoord.push_back(coord[iV*3+1]);
            newCoord.push_back(coord[iV*3+2]);
            iVnew++;
        }
    }

    // since faces incident to selected edges must be deleted, you may
    // want to create a data structure to identify them before the next
    // loop

    // create the new coordIndex array
    vector<int> newCoordIndex;
    // traverse the coordIndex array
    int iF,iC0,iC1,iC;
    for(iF=iC0=iC1=0;iC1<nC;iC1++) {
        if(coordIndex[iC1]>=0) continue;
        //
        // if face iF is incident to a selected edge, it has to be
        // deleted, continue
        if (selectedEdgeFace[iF]) {
            iC0 = iC1+1; iF++;
            continue;
        }

        // for each corner iC of face iF {
        //   using the vertexMap array determine the new vertex index iVnew
        //   push_back that value onto the new coordIndex array
        // }
        // don't forget to add a face separator
        for(iC=iC0; iC<iC1; iC++){
            iVnew = vertexMap[pmesh->getSrc(iC)];
            newCoordIndex.push_back(iVnew);
        }
        newCoordIndex.push_back(-1);

        // advance to next face
        iC0 = iC1+1; iF++;
    }

    // fix _ifsRendering
    coord.swap(newCoord);
    coordIndex.swap(newCoordIndex);

    // the PolygonMesh is no longer valid
    ifsv.deletePolygonMesh();

    // since the number of vertices and faces has changed colors and
    // are no longer valid
    _ifsOptimized->clearColor();
    // if the mesh had colors and you want to preserve them, a better
    // solution would be to map old colors onto new colors

    // since a lot of vertices and faces may have changed, recompute
    // face normals
    _ifsOptimized->clearNormal();
    _ifsOptimized->setNormalPerVertex(false);
    Geometry::computeNormalsPerFace(coord,coordIndex,_ifsOptimized->getNormal());

    // clear all selection buffers
    ifsv.clearAllSelection();

    //Actualizo la variable _edgeLengths que se usa en este algoritmo para poder usarlo varias veces sobre un mismo IndexedFaceSet sin que haya errores
    pmesh = ifsv.getPolygonMesh(true);
    Geometry::computeEdgeLengths(coord,*pmesh,_edgeLengths);

    _ifsOptimized->printInfo("  ");
}

void Optimization::_adaptiveSubdivisionSelect
(vector<int>& vertexSelection,
 vector<int>& edgeSelection,
 const SplitEdgesMode mode,
 bool colorIncidentFaces) {

    if(_ifsOptimized==nullptr) return;
    vector<int>& coordIndex = _ifsOptimized->getCoordIndex();
    if(Geometry::isTriangulated(coordIndex)==false) return;

    IndexedFaceSetVariables ifsv(*_ifsOptimized);
    PolygonMesh* pmesh = ifsv.getPolygonMesh(true);

    int nV = pmesh->getNumberOfVertices();
    int nE = pmesh->getNumberOfEdges();
    int nF = pmesh->getNumberOfFaces();
    int nC = pmesh->getNumberOfCorners();

    // for visualization purposes
    // set color per face indexed with three colors
    vector<int>*   colorIndex             = (vector<int>*)0;
    vector<float>* color                  = (vector<float>*)0;
    int            noSplitColorIndex      = 0;
    int            split2ColorIndex       = 1;
    int            split4ColorIndex       = 2;
    if(colorIncidentFaces) {

        colorIndex = &(_ifsOptimized->getColorIndex());
        colorIndex->clear();
        colorIndex->resize(nF,noSplitColorIndex);

        color = &(_ifsOptimized->getColor());
        color->clear();

        _ifsOptimized->setColorPerVertex(false);

        Color noSplitColor(0.8f,0.8f,0.8f); // index 0
        color->push_back(noSplitColor.r);
        color->push_back(noSplitColor.g);
        color->push_back(noSplitColor.b);

        Color split2EdgesColor(0.6f,0.9f,0.3f); // index 1
        color->push_back(split2EdgesColor.r);
        color->push_back(split2EdgesColor.g);
        color->push_back(split2EdgesColor.b);

        Color split4EdgesColor(0.3f, 0.6f,0.9f); // index 2
        color->push_back(split4EdgesColor.r);
        color->push_back(split4EdgesColor.g);
        color->push_back(split4EdgesColor.b);
    }

    // select vertices of edges which have to be split

    // int nVselected = 0;
    vertexSelection.clear();
    vertexSelection.resize(nV,-1);
    // use
    // vertexSelection[iV] = _vSelIndex;
    // to select vertex iV

    if(mode==SplitEdgesMode::ALL) {
        // select al the vertices
        for (int iV=0; iV<nV; iV++) vertexSelection[iV] = _vSelIndex;

    } else if(mode==SplitEdgesMode::SELECTED) {

        vector<int>& eSel = ifsv.getEdgeSelection();
        // select the vertices iV0 and iV1 of each selected edge (eSel[iE]!=-1)
        for (int iE=0; iE<nE; iE++) {
            if (eSel[iE] != _eSelIndex) continue;
            int iV0 = pmesh->getVertex0(iE);
            int iV1 = pmesh->getVertex1(iE);
            vertexSelection[iV0] = _vSelIndex;
            vertexSelection[iV1] = _vSelIndex;
        }

    } else if(mode==SplitEdgesMode::LONG) {

        // using the array _edgeLengths of pre-computed edge lengths
        // select all vertices of edges longer than highEdgeLength
        float highEdgeLength = _targetEdgeLength;
        for (int iE=0; iE<nE; iE++) {
            if (_edgeLengths[iE] <= highEdgeLength) continue;
            int iV0 = pmesh->getVertex0(iE);
            int iV1 = pmesh->getVertex1(iE);
            vertexSelection[iV0] = _vSelIndex;
            vertexSelection[iV1] = _vSelIndex;
        }
    }

    // std::cout << "nVselected = " << nVselected << std::endl;

    // select all the edges with the two ends selected
    int iE,iF,iV;

    edgeSelection.clear();
    edgeSelection.resize(nE,-1);
    for(iE=0;iE<nE;iE++) {
        // if(_edgeLengths[iE]<=_highEdgeLength) continue;
        iV = pmesh->getVertex0(iE);
        if(vertexSelection[iV]==-1) continue;
        iV = pmesh->getVertex1(iE);
        if(vertexSelection[iV]==-1) continue;
        edgeSelection[iE] = _eSelIndex;
    }

    // color faces for visualization purposes
    if(colorIncidentFaces) {
        int iF,iC0,iC1,nVsel;

        for(iF=iC0=iC1=0;iC1<nC;iC1++) {
            if(coordIndex[iC1]>=0) continue;
            // remember that the faces are all triangles
            // assert(iC1-iC0==3);
            assert(iC1-iC0==3);
            //
            // count number of selected vertices of the triangle
            nVsel = 0;
            // ...
            for(int iC=iC0; iC<iC1; iC++){
                iV = pmesh->getSrc(iC);
                if (vertexSelection[iV]==_vSelIndex) nVsel++;
            }
            // painT the face accordingly
            (*colorIndex)[iF] =
                (nVsel<=1)?noSplitColorIndex:
                    (nVsel==2)?split2ColorIndex:
                    split4ColorIndex;

            // advance to next face
            iC0=iC1+1; iF++;
        }
    }
}

void Optimization::adaptiveSubdivisionShow(const SplitEdgesMode mode) {
  IndexedFaceSetVariables ifsv(*_ifsOptimized);
  vector<int>& vertexSelection = ifsv.getVertexSelection();
  vector<int>& edgeSelection   = ifsv.getEdgeSelection();
  _adaptiveSubdivisionSelect(vertexSelection,edgeSelection,mode,true);
}

void Optimization::adaptiveSubdivisionApply
(const SplitEdgesMode mode, const bool colorFaces) {

    vector<int> vertexSelection;
    vector<int> edgeSelection;
    _adaptiveSubdivisionSelect(vertexSelection,edgeSelection,mode,true);

    vector<float>& coord      = _ifsOptimized->getCoord();
    vector<int>&   coordIndex = _ifsOptimized->getCoordIndex();

    // _ifsOptimized should have color per face indexed if colorFaces==true
    if(colorFaces==false) _ifsOptimized->clearColor();
    // vector<float>& color = _ifsOptimized->getColor();
    vector<int>&   colorIndex = _ifsOptimized->getColorIndex();

    // remove normal vectors & texCoords
    // normal vectors will be recomputed at the end
    _ifsOptimized->clearNormal();
    _ifsOptimized->clearTexCoord();

    IndexedFaceSetVariables ifsv(*_ifsOptimized);
    PolygonMesh* pmesh = ifsv.getPolygonMesh(true);
    int nV = pmesh->getNumberOfVertices();
    int nE = pmesh->getNumberOfEdges();
    int nF = pmesh->getNumberOfFaces();
    int nC = pmesh->getNumberOfCorners();

    // a new vertex must be created for each selected edge

    // assign new vertex indices to edges and compute coordinates of new
    // vertices
    //
    // since all the original vertices are preserved, new vertex indices
    // start with iVnew=nV
    //
    // when you create new coordinates for a new vertices, you can just
    // push them back onto the original coord array.

    vector<int> edgeToNewVertex(nE,-1);
    // for each edge iE {
    //   if either iV0 or iV1 is not selected, continue;
    //   create new vertex at edge midpoint and save them
    //   save the new vertex index into the edgeToNewVertex array
    // }
    int iVnew = nV;
    for (int iE=0; iE<nE; iE++) {
        int iV0 = pmesh->getVertex0(iE);
        int iV1 = pmesh->getVertex1(iE);
        if (vertexSelection[iV0] != _vSelIndex || vertexSelection[iV1] != _vSelIndex) continue;

        vector<float> coord_iV0 = {coord[iV0*3], coord[iV0*3+1], coord[iV0*3+2]};
        vector<float> coord_iV1 = {coord[iV1*3], coord[iV1*3+1], coord[iV1*3+2]};
        vector<float> coord_iVnew = {
            (coord_iV0[0] + coord_iV1[0]) / 2,
            (coord_iV0[1] + coord_iV1[1]) / 2,
            (coord_iV0[2] + coord_iV1[2]) / 2
        };
        coord.push_back(coord_iVnew[0]);
        coord.push_back(coord_iVnew[1]);
        coord.push_back(coord_iVnew[2]);

        edgeToNewVertex[iE] = iVnew;
        iVnew++;
    }


    // now split the triangles
    // for each triangle iF connecting vertices iV0, iV1, and iV2 {
    //   count the number of selected vertices nVsel=0,1,2,or 3
    //   if nVsel==0 or nVsel==1 {
    //     keep the triangle unchanged
    //   } else if nVsel==2 {
    //     for example, if the selected vertices are iV1 and iV2
    //     find the edge iE defined by iV1 and iV2
    //     find the new vertex index iV12 associated with the edge
    //
    //     split the input triangle into two output triangles
    //
    //           iV0
    //        /   |    \
    //       /    |     \
    //     iV1 - iV12 - iV2
    //
    //   } else if nVsel==3 {
    //
    //     split into 3 triangles
    //
    //           iV2
    //         /      \
    //       iV02 -- iV12
    //       /  \    /   \
    //     iV0 - iV01 - iV1
    //   }
    // }
    int iF,iC0,iC1,iV0,iV1,iV2,nVsel;

    for(iF=iC0=iC1=0;iC1<nC;iC1++) {
        if(coordIndex[iC1]>=0) continue;

        nVsel = 0;
        iV0 = pmesh->getSrc(iC0);
        iV1 = pmesh->getSrc(iC0+1);
        iV2 = pmesh->getSrc(iC0+2);
        if (vertexSelection[iV0] == _vSelIndex) nVsel++;
        if (vertexSelection[iV1] == _vSelIndex) nVsel++;
        if (vertexSelection[iV2] == _vSelIndex) nVsel++;

        if (nVsel==0 || nVsel==1) {
            iC0=iC1+1; iF++;
            continue;
        }

        if (nVsel==2) {
            //Creo una variable para almacenar los vértices en las puntas de la arista que se va a dividir en dos, y otra para almacenar el vértice restante
            vector<int> splitEdge;
            int remainingVertex;

            if (vertexSelection[iV0] == _vSelIndex) {
                splitEdge.push_back(iV0);
            } else {
                remainingVertex = iV0;
            }

            if (vertexSelection[iV1] == _vSelIndex) {
                splitEdge.push_back(iV1);
            } else {
                remainingVertex = iV1;
            }

            if (vertexSelection[iV2] == _vSelIndex) {
                splitEdge.push_back(iV2);
            } else {
                remainingVertex = iV2;
            }

            int iE = pmesh->getEdge(min(splitEdge[0], splitEdge[1]), max(splitEdge[0], splitEdge[1]));
            iVnew = edgeToNewVertex[iE];

            //Genero dos caras nuevas: la primera la ubico al final de coordIndex, y la otra reemplaza la cara que fue dividida en dos
            //La orientación de las caras generadas debe mantener la orientación de la original respecto de las caras vecinas
            //Para saber cuál orientación usar basta verificar si la arista partida es (iV0,iV2), ya que esta es la única arista que está orientada del segundo componente al primero

            if (splitEdge[0]==iV0 && splitEdge[1]==iV2) {
                //Si la arista partida es (iV0, iV2), doy vuelta el vector splitEdge para que splitEdge[0] siempre sea el que apunta al nuevo vértice
                int tmp = splitEdge[0];
                splitEdge[0] = splitEdge[1];
                splitEdge[1] = tmp;
            }

            coordIndex.push_back(splitEdge[0]);
            coordIndex.push_back(iVnew);
            coordIndex.push_back(remainingVertex);
            coordIndex.push_back(-1);

            coordIndex[iC0] = remainingVertex;
            coordIndex[iC0+1] = iVnew;
            coordIndex[iC0+2] = splitEdge[1];

            //Las nuevas caras tienen el mismo color que la cara original (la cara nueva que reemplaza a la original hereda su iF y por lo tanto ya cumple esto)
            if (colorFaces) colorIndex.push_back(colorIndex[iF]);
        }

        if (nVsel==3) {
            int iE01 = pmesh->getEdge(min(iV0, iV1), max(iV0, iV1));
            int iE12 = pmesh->getEdge(min(iV1, iV2), max(iV1, iV2));
            int iE02 = pmesh->getEdge(min(iV0, iV2), max(iV0, iV2));

            int iV01 = edgeToNewVertex[iE01];
            int iV12 = edgeToNewVertex[iE12];
            int iV02 = edgeToNewVertex[iE02];

            //Ubico tres de las caras nuevas al final de coordIndex, y la cuarta reemplaza la cara que fue dividida
            //Sabiendo que por definición la orientación es iV0 -> iV1 -> iV2 -> iV0, las caras nuevas se definen para mantener esta orientación respecto de las caras vecinas

            coordIndex.push_back(iV0);
            coordIndex.push_back(iV01);
            coordIndex.push_back(iV02);
            coordIndex.push_back(-1);

            coordIndex.push_back(iV1);
            coordIndex.push_back(iV12);
            coordIndex.push_back(iV01);
            coordIndex.push_back(-1);

            coordIndex.push_back(iV2);
            coordIndex.push_back(iV02);
            coordIndex.push_back(iV12);
            coordIndex.push_back(-1);

            coordIndex[iC0] = iV01;
            coordIndex[iC0+1] = iV12;
            coordIndex[iC0+2] = iV02;

            if (colorFaces) {
                //A las tres caras nuevas añadidas al final de coordIndex les asigno el color de la cara original
                colorIndex.push_back(colorIndex[iF]);
                colorIndex.push_back(colorIndex[iF]);
                colorIndex.push_back(colorIndex[iF]);
            }
        }

        // advance to next face
        iC0=iC1+1; iF++;
    }

    // clear all selection buffers
    ifsv.clearAllSelection();

    // delete pmesh
    ifsv.deletePolygonMesh();
    pmesh = nullptr;

    // recompute face normals
    _ifsOptimized->setNormalPerVertex(false);
    Geometry::computeNormalsPerFace(coord,coordIndex,_ifsOptimized->getNormal());

    //Actualizo las longitudes en _edgeLengths ahora que cambiaron las aristas
    pmesh = ifsv.getPolygonMesh(true);
    Geometry::computeEdgeLengths(coord,*pmesh,_edgeLengths);

}

void Optimization::_equalizeValencesSelect
(vector<int>& edgeSelection, bool colorIncidentFaces) {

  // for each edge e do
  // let a, b, c, d be the vertices of the two triangles adjacent to e
  // deviation_pre =
  //   abs(valence(a)-target_val(a)) +
  //   abs(valence(b)-target_val(b)) +
  //   abs(valence(c)-target_val(c)) +
  //   abs(valence(d)-target_val(d))
  // flip(e)
  // deviation_post =
  //   abs(valence(a)-target_val(a)) +
  //   abs(valence(b)-target_val(b)) +
  //   abs(valence(c)-target_val(c)) +
  //   abs(valence(d)-target_val(d))
  // if deviation_pre <= deviation_post do
  // flip(e)

  // vector<float>& coord      = _ifsOptimized->getCoord();
  // vector<int>&   coordIndex = _ifsOptimized->getCoordIndex();
  IndexedFaceSetVariables ifsv(*_ifsOptimized);
  PolygonMesh* pmesh = ifsv.getPolygonMesh(true);

  int nV,nE,nF,nEF,iE,iV0,iV1,iV2,iV3,iC0,iC1,iC0n,iC1n,iF,i;
  int targetValence0,targetValence1,targetValence2,targetValence3;
  int errBefore,errAfter,errDiff;

  nV = pmesh->getNumberOfVertices();
  nE = pmesh->getNumberOfEdges();
  nF = pmesh->getNumberOfFaces();

  vector<int> valence;
  Geometry::computeValences(*pmesh,valence);

  // clear edge selection
  edgeSelection.clear();
  edgeSelection.resize(nE,-1);

  // set color per face
  vector<int>*   colorIndex               = (vector<int>*)0;
  vector<float>* color                    = (vector<float>*)0;
  int            noChangesColorIndex      = 0;
  int            flipEdgesColorIndex = 1;
  if(colorIncidentFaces) {
    colorIndex = &(_ifsOptimized->getColorIndex());
    colorIndex->clear();
    color = &(_ifsOptimized->getColor());
    color->clear();
    _ifsOptimized->setColorPerVertex(false);
    Color noChangesColor(0.8f,0.8f,0.8f); // index 0
    noChangesColorIndex = static_cast<int>(color->size()/3);
    color->push_back(noChangesColor.r);
    color->push_back(noChangesColor.g);
    color->push_back(noChangesColor.b);
    Color flipEdgesColor(0.7f,0.5f,0.9f); // index 1
    flipEdgesColorIndex = static_cast<int>(color->size()/3);
    color->push_back(flipEdgesColor.r);
    color->push_back(flipEdgesColor.g);
    color->push_back(flipEdgesColor.b);
    for(iF=0;iF<nF;iF++)
      colorIndex->push_back(noChangesColorIndex);
  }

  Heap heap;
  for(iE=0;iE<nE;iE++) {
    if(pmesh->isRegularEdge(iE)) {
      // pmesh->getNumberOfEdgeFaces(iE)==2
      // pmesh->getNumberOfEdgeHalfEdges(iE)==2

      // iV0 = pmesh->getVertex0(iE);
      // iV1 = pmesh->getVertex1(iE);

      iC0  = pmesh->getEdgeHalfEdge(iE,0);
      iC0n = pmesh->getNext(iC0);
      iC1  = pmesh->getEdgeHalfEdge(iE,1);
      iC1n = pmesh->getNext(iC1);

      //      iV1
      //    /     \
      // iV0 ----- iV2
      //    \     /
      //      iV3

      iV0  = pmesh->getSrc(iC0n);
      iV1  = pmesh->getDst(iC0n);
      iV2  = pmesh->getSrc(iC1n);
      iV3  = pmesh->getDst(iC1n);

      targetValence0 = (pmesh->isBoundaryVertex(iV0))?4:6;
      targetValence1 = (pmesh->isBoundaryVertex(iV1))?4:6;
      targetValence2 = (pmesh->isBoundaryVertex(iV2))?4:6;
      targetValence3 = (pmesh->isBoundaryVertex(iV3))?4:6;

      // current
      // (iV0,iV1,iV2) & (iV2,iV3,iV0)
      errBefore =
        abs(valence[iV0]-targetValence0) +
        abs(valence[iV1]-targetValence1) +
        abs(valence[iV2]-targetValence2) +
        abs(valence[iV3]-targetValence3);
      
      // swap
      // (iV0,iV1,iV3) & (iV1,iV2,iV3)
      errAfter =
        abs(valence[iV0]-1-targetValence0) +
        abs(valence[iV1]+1-targetValence1) +
        abs(valence[iV2]-1-targetValence2) +
        abs(valence[iV3]+1-targetValence3);

      if((errDiff=errAfter-errBefore)<0) {
        heap.add((float)errDiff,iE);
      }

    }
  }

  // create an independent set ef edges to flip
  vector<bool> usedVertex(nV,false);
  while(heap.length()>0) {
    /* iH = */ heap.delMin();
    iE      =  heap.getLastIKey();
    // errDiff = -heap.getLastFKey();

    iC0  = pmesh->getEdgeHalfEdge(iE,0);  // iV2->iV0
    iC0n = pmesh->getNext(iC0);           // iV0->iV1
    iC1  = pmesh->getEdgeHalfEdge(iE,1);  // iV0->iV2
    iC1n = pmesh->getNext(iC1);           // iV2->iV3

    iV0  = pmesh->getSrc(iC0n);
    iV1  = pmesh->getDst(iC0n);
    iV2  = pmesh->getSrc(iC1n);
    iV3  = pmesh->getDst(iC1n);

    if(usedVertex[iV0] || usedVertex[iV1]
       /* || usedVertex[iV2] || usedVertex[iV3] */ ) {
      // iE REJECTED
      continue;
    }

    if(pmesh->getEdge(iV1,iV3)>=0) {
      // iE REJECTED (DIAGONAL)
      continue;
    }

    usedVertex[iV0]   = true;
    usedVertex[iV1]   = true;
    usedVertex[iV2]   = true;
    usedVertex[iV3]   = true;
    edgeSelection[iE] = _eSelIndex;

    if(colorIncidentFaces) {
      nEF = pmesh->getNumberOfEdgeFaces(iE);
      for(i=0;i<nEF;i++) {
        iF = pmesh->getEdgeFace(iE,i);
        (*colorIndex)[iF] = flipEdgesColorIndex;
      }
    }
  }
}

void Optimization::equalizeValencesShow() {
  IndexedFaceSetVariables ifsv(*_ifsOptimized);
  vector<int>& edgeSelection = ifsv.getEdgeSelection();
  _equalizeValencesSelect(edgeSelection,true);
}

void Optimization::equalizeValencesApply() {

  vector<int> edgeSelection;
  _equalizeValencesSelect(edgeSelection,false);

  vector<float>& coord      = _ifsOptimized->getCoord();
  vector<int>&   coordIndex = _ifsOptimized->getCoordIndex();

  IndexedFaceSetVariables ifsv(*_ifsOptimized);
  PolygonMesh* pmesh = ifsv.getPolygonMesh(true);

  int nE,iV0,iV1,iV2,iV3,iE,iC0,iC1,iC0n,iC1n,iF0,iF1;

  // nV = pmesh->getNumberOfVertices();
  nE = pmesh->getNumberOfEdges();
  // nF = pmesh->getNumberOfFaces();

  // ERROR | it should not happen
  if(static_cast<int>(edgeSelection.size())!=nE) return;

  for(iE=0;iE<nE;iE++) {
    if(edgeSelection[iE]<0) continue; // edge not selected for flipping

    iC0  = pmesh->getEdgeHalfEdge(iE,0);  // iV2->iV0
    iC0n = pmesh->getNext(iC0);           // iV0->iV1
    iC1  = pmesh->getEdgeHalfEdge(iE,1);  // iV0->iV2
    iC1n = pmesh->getNext(iC1);           // iV2->iV3

    iV0  = pmesh->getSrc(iC0n); // == pmesh->getDst(iC0) == pmesh->getSrc(iC1);
    iV1  = pmesh->getDst(iC0n);
    iV2  = pmesh->getSrc(iC1n); // == pmesh->getDst(iC1) == pmesh->getSrc(iC0);
    iV3  = pmesh->getDst(iC1n);

    // we know that the edge is regular
    iF0  = pmesh->getFace(iC0); // (iV0,iV1,iV2)
    iF1  = pmesh->getFace(iC1); // (iV2,iV3,iV0)

    //      iV1
    //    / iF0 \
    // iV0 ----- iV2
    //    \ iF1 /
    //      iV3

    // iF0 : (iV0,iV1,iV2) -> (iV0,iV1,iV3)
    coordIndex[4*iF0+0] = iV0;
    coordIndex[4*iF0+1] = iV1;
    coordIndex[4*iF0+2] = iV3;

    // iF1 : (iV2,iV3,iV0) -> (iV2,iV3,iV1)
    coordIndex[4*iF1+0] = iV2;
    coordIndex[4*iF1+1] = iV3;
    coordIndex[4*iF1+2] = iV1;
  }

  // fix _ifsRendering

  // rebuild PolygonMesh
  ifsv.deletePolygonMesh();
  // pmesh = ifsv.getPolygonMesh(true);

  // clear colors ???
  // _ifsOptimized->clearColor();

  // recompute face normals
  _ifsOptimized->clearNormal();
  _ifsOptimized->setNormalPerVertex(false);
  Geometry::computeNormalsPerFace(coord,coordIndex,_ifsOptimized->getNormal());

  // clear all selection buffers ???
  ifsv.clearAllSelection();

  _ifsOptimized->printInfo("  ");
}
