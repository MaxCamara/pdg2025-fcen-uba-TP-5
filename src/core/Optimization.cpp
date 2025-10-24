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

  // TODO - get from Assignment 4
  
}

//////////////////////////////////////////////////////////////////////
void Optimization::laplacianSmoothingFaceNormalsRun(const bool normalize) {

  // TODO - get from Assignment 4

}

//////////////////////////////////////////////////////////////////////
void  Optimization::jacobiRun() {

  // TODO - get from Assignment 4

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

  // TODO - get from Assignment 4

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

  // TODO - get from Assignment 4

}

void Optimization::_adaptiveSubdivisionSelect
(vector<int>& vertexSelection,
 vector<int>& edgeSelection,
 const SplitEdgesMode mode,
 bool colorIncidentFaces) {

  // TODO - get from Assignment 4

}

void Optimization::adaptiveSubdivisionShow(const SplitEdgesMode mode) {
  IndexedFaceSetVariables ifsv(*_ifsOptimized);
  vector<int>& vertexSelection = ifsv.getVertexSelection();
  vector<int>& edgeSelection   = ifsv.getEdgeSelection();
  _adaptiveSubdivisionSelect(vertexSelection,edgeSelection,mode,true);
}

void Optimization::adaptiveSubdivisionApply
(const SplitEdgesMode mode, const bool colorFaces) {

  // TODO - get from Assignment 4

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
