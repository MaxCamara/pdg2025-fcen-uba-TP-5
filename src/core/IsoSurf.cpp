//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-07 20:22:20 taubin>
//------------------------------------------------------------------------
//
// IsoSurf.cpp
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

#include "IsoSurf.hpp"
#include "Geometry.hpp"
#include "core/SimpleGraphMap.hpp"

const int IsoSurf::_faceTable[12][12] = {
  // i=0                       // i=1
  //   1-(01)-5                //   3-(03)-7
  // (08)    (10)              // (05)    (07)
  //   0-(00)-4                //   1-(01)-5
  // (04)    (06)              // (08)    (10)
  //   2-(02)-6                //   0-(00)-4
  { 0,8,1,1,5,10,4,6,6,2,2,4}, { 1,5,3,3,7,7,5,10,4,0,0,8},
  // i=2                       // i=3
  //   0-(00)-4                //   2-(02)-6
  // (04)    (06)              // (09)    (11)
  //   2-(02)-6                //   3-(03)-7
  // (09)    (11)              // (05)    (07)
  //   3-(03)-7                //   1-(01)-5
  { 2,4,0,0,4,6,6,11,7,3,3,9}, { 3,9,2,2,6,11,7,7,5,1,1,5},
  // i=4                       // i=5
  //   4-(06)-6                //   0-(04)-2
  // (00)    (02)              // (08)    (09)
  //   0-(04)-2                //   1-(05)-3
  // (08)    (09)              // (01)    (03)
  //   1-(05)-3                //   5-(07)-7
  { 0,0,4,6,6,2,2,9,3,5,1,8},  { 1,8,0,4,2,9,3,3,7,7,5,1},
  // i=6                       // i=7
  //   5-(07)-7                //   1-(05)-3
  // (10)    (11)              // (01)    (03)
  //   4-(06)-6                //   5-(07)-7
  // (00)    (02)              // (10)    (11)
  //   0-(04)-2                //   4-(06)-6
  { 4,10,5,7,7,11,6,2,2,4,0,0},{ 5,1,1,5,3,3,7,11,6,6,4,10},
  // i=8                       // i=9
  //   2-(09)-3                //   6-(11)-7
  // (04)    (05)              // (02)    (03)
  //   0-(08)-1                //   2-(09)-3
  // (00)    (01)              // (04)    (05)
  //   4-(10)-5                //   0-(08)-1
  { 0,4,2,9,3,5,1,1,5,10,4,0}, { 2,2,6,11,7,3,3,5,1,8,0,4},
  // i=10                      // i=11
  //   0-(08)-1                //   4-(10)-5
  // (00)    (01)              // (06)    (07)
  //   4-(10)-5                //   6-(11)-7
  // (06)    (07)              // (02)    (03)
  //   6-(11)-7                //   2-(09)-3
  { 4,0,0,8,1,1,5,7,7,11,6,6}, { 6,6,4,10,5,7,7,3,3,9,2,2}
};

//    vertices      //    edges                 //    faces
//      6-----7     //        [6]---11---[7]    //        1
//     /|    /|     //        /|         /|     //        | 3
//    4-----5 |     //       6 2        7 3     //        |/
//    | 2---|-3     //      /  |       /  |     //    4---+---5
//    |/    |/      //    [4]---10---[5]  |     //       /|
//    0-----1       //     |   |      |   |     //      2 |
//                  //     |  [2]--9--|--[3]    //        0
//     i            //     0  /       1  /      //
//     | j          //     | 4        | 5       //
//     |/           //     |/         |/        //
//     +---k        //    [0]---8----[1]        //

const int IsoSurf::_edgeTable[12][2] = {
  /* _edgeTable[ 0] = */ {0,4}, // +4 | [i,i+1],j  ,k
  /* _edgeTable[ 1] = */ {1,5}, // +4 | [i,i+1],j  ,k+1
  /* _edgeTable[ 2] = */ {2,6}, // +4 | [i,i+1],j+1,k
  /* _edgeTable[ 3] = */ {3,7}, // +4 | [i,i+1],j+1,k+1
  /* _edgeTable[ 4] = */ {0,2}, // +2 | i  ,[j,j+1]  ,k
  /* _edgeTable[ 5] = */ {1,3}, // +2 | i  ,[j,j+1]  ,k+1
  /* _edgeTable[ 6] = */ {4,6}, // +2 | i+1,[j,j+1]  ,k
  /* _edgeTable[ 7] = */ {5,7}, // +2 | i+1,[j,j+1]  ,k+1
  /* _edgeTable[ 8] = */ {0,1}, // +1 | i  ,j  ,[k,k+1]
  /* _edgeTable[ 9] = */ {2,3}, // +1 | i  ,j+1,[k,k+1]
  /* _edgeTable[10] = */ {4,5}, // +1 | i+1,j  ,[k,k+1]
  /* _edgeTable[11] = */ {6,7}  // +1 | i+1,j+1,[k,k+1]
};


const int (*IsoSurf::getEdgeTable())[2] {
  return _edgeTable;
}

//////////////////////////////////////////////////////////////////////
// static
int IsoSurf::makeCellFaces
(bool b[/*8*/], int iE[/*12*/], vector<int>& coordIndex) {

  // if(b==null || b.length<8) return;
  // if(iE==null || iE.length<12) return;
  // if(coordIndex==null) return;

  // link vertices
  int i,j,k,j0,j1;
  int next[12];
  for(i=0;i<12;i++)
    next[i] = -1;
  for(i=j=0;j<3;j++) { // j=0,1,2
    for(j0=0;j0<2;j0++) {
      for(j1=0;j1<2;j1++,i++) {
        if(iE[i]>=0) {
          k =
            (b[_faceTable[i][0]])?
            ((!b[_faceTable[i][2]])?1:(!b[_faceTable[i][4]])?3:5):
            ((!b[_faceTable[i][8]])?7:(!b[_faceTable[i][10]])?9:11);
          next[i] = _faceTable[i][k];
        }
      }
    }
  }

  // traverse linked list and output faces
  int j_prev,j_next,nCorners;
  int nFaces = 0;
  for(i=0;i<12;i++) {
    if(next[i]>=0) {
      nCorners = 0;
      j_prev = -1;
      j = i;
      do {
        j_next = next[j];
        // skip repeated vertex indices
        if(j_prev<0 || iE[j]!=iE[j_prev] || (j_next==i && iE[j]!=iE[i])) {
          coordIndex.push_back(iE[j]); nCorners++;
        }
        j_prev = j; j=j_next; next[j_prev] = -1;
      } while(j!=i);
      if(nCorners>=3) {
        // faces with 3 or more corners are OK
        coordIndex.push_back(-1);
        nFaces++;
      } else {
        // remove faces with less than 3 corners
        for(;nCorners>0;nCorners--)
          coordIndex.pop_back();
      }
    }
  }
  return nFaces;
}


//////////////////////////////////////////////////////////////////////
//
// compute isosurface for dense regular grid scalar field
//
// static
void IsoSurf::computeIsosurface
(const Vec3f& center, const Vec3f& size,
 const int depth, const float scale, const bool isCube,
 vector<float>& fGrid, IndexedFaceSet& isoSurface) {

  float dx = size.x/2.0f;
  float dy = size.y/2.0f;
  float dz = size.z/2.0f;
  float dMax = dx; if(dy>dMax) dMax=dy; if(dz>dMax) dMax=dz;
  if(isCube) {
    dx = dMax; dy = dMax; dz = dMax;
  }
  if(scale>0.0f) {
    dx *= scale; dy *= scale; dz *= scale;
  }
  float x,y,z;
  float x0 = center.x-dx; float y0 = center.y-dy; float z0 = center.z-dz;
  float x1 = center.x+dx; float y1 = center.y+dy; float z1 = center.z+dz;
  Vec3f min(x0,y0,z0);
  Vec3f max(x1,y1,z1);

  int iCell,ix,iy,iz,jx,jy,jz;

  // create output surface

  IndexedFaceSet*   surface       = &isoSurface; // 
  surface->clear();
  vector<float>&    coordIfs      = surface->getCoord();
  vector<int>&      coordIndexIfs = surface->getCoordIndex();

  int   iVB[8]; // bbox corner vertex indices
  Vec3f v[8];   // bbox corner coordinates
  float F[8];   // function values at bbox corners
  bool  b[8];   // function is positive or negative ?
  float tj,tk;
  int   iE[12],i,j,k,iV,N = 1<<depth;
  
  SimpleGraphMap graph; // stores map from grid edges to isovertices

  Vec4f f;
  Vec3f minCell,maxCell,nCell;
  vector<float> coordCell;
  vector<float> normalCell;
  // int nFacesCell    = 0;

  for(iCell=iz=0;iz<N;iz++) {
    minCell.z = (((float)(N-iz  ))*z0+((float)(iz  ))*z1)/((float)N);
    maxCell.z = (((float)(N-iz-1))*z0+((float)(iz+1))*z1)/((float)N);
    for(iy=0;iy<N;iy++) {
      minCell.y = (((float)(N-iy  ))*y0+((float)(iy  ))*y1)/((float)N);
      maxCell.y = (((float)(N-iy-1))*y0+((float)(iy+1))*y1)/((float)N);
      for(ix=0;ix<N;ix++,iCell++) {
        minCell.x = (((float)(N-ix  ))*x0+((float)(ix  ))*x1)/((float)N);
        maxCell.x = (((float)(N-ix-1))*x0+((float)(ix+1))*x1)/((float)N);
        
        for(i=0;i<8;i++) {
          v[i].z = (((i>>0)&0x1)==0)?minCell.z:maxCell.z;
          v[i].y = (((i>>1)&0x1)==0)?minCell.y:maxCell.y;
          v[i].x = (((i>>2)&0x1)==0)?minCell.x:maxCell.x;
          jz = iz+((i>>0)&0x1);
          jy = iy+((i>>1)&0x1);
          jx = ix+((i>>2)&0x1);
          iVB[i] = iV = jx+(N+1)*(jy+(N+1)*jz);
          b[i] = ((F[i]=fGrid[static_cast<size_t>(iV)])<0.0f);
        }
        
        for(i=0;i<12;i++) {
          iV   = -1;
          j    = _edgeTable[i][0];
          k    = _edgeTable[i][1];
          if(b[j]!=b[k]) {
            // look for the index of the isovertex associated with this grid edge
            iV = graph.get(iVB[j],iVB[k]);
            if(iV<0) { // need to create a new vertex
              iV = static_cast<int>((coordIfs.size()/3));
              // isovertex coordinates
              tk = F[j]/(F[j]-F[k]);
              tj = F[k]/(F[k]-F[j]);
              x  = tj*v[j].x+tk*v[k].x;
              y  = tj*v[j].y+tk*v[k].y;
              z  = tj*v[j].z+tk*v[k].z;
              coordIfs.push_back(x);
              coordIfs.push_back(y);
              coordIfs.push_back(z);
              // save the association
              graph.insert(iVB[j],iVB[k],iV);
            }
            
          }
          iE[i] = iV;
        }
        
        // make the faces of the isosurface within this cell
        /* nFacesCell = */ IsoSurf::makeCellFaces(b,iE,coordIndexIfs);
      }
    }
  }

  // add normal per face
  surface->clearNormal();
  vector<float>& normalIfs = surface->getNormal();
  surface->setNormalPerVertex(false);
  Geometry::computeNormalsPerFace(coordIfs,coordIndexIfs,normalIfs);
  
}

//////////////////////////////////////////////////////////////////////
//
// compute isosurface for sparse regular grid scalar field
//
// static
void IsoSurf::computeIsosurface
(HexGridPartition& hgp,
 vector<float>& fGrid, IndexedFaceSet& isoSurface) {
  isoSurface.clear();
  int nGridVertices = hgp.getNumberOfVertices();
  if(static_cast<int>(fGrid.size())!=nGridVertices) return;

  IndexedFaceSet*   surface       = &isoSurface; // 
  vector<float>&    coordIfs      = surface->getCoord();
  vector<int>&      coordIndexIfs = surface->getCoordIndex();

  Vec3f& gridMin = hgp.getMin();
  float x0=gridMin.x, y0=gridMin.y, z0=gridMin.z;
  Vec3f& gridMax = hgp.getMax();
  float x1=gridMax.x, y1=gridMax.y, z1=gridMax.z;

  int   iVB[8]; // cell corner vertex indices
  Vec3f v[8];   // cell corner coordinates
  float F[8];   // function values at bbox corners
  bool  b[8];   // function is positive or negative ?
  float tj,tk,x,y,z;
  int   iE[12],h,h0,h1,h2,i,j,k,iV,iVmap,iC,ix,iy,iz,jx,jy,jz,iCell;
  int   N = hgp.getResolution();

  Vec4f f;
  Vec3f minCell,maxCell,nCell;
  vector<float> coordCell;
  vector<float> normalCell;

  // used to map grid edges, represented as grid vertex pairs, onto isovertices
  SimpleGraphMap graphMap;

  // maps regular grid scan order vertex indices onto sparce grid vertex indices
  map<int,int> vMap;
  
  map<int,int> first =  hgp.getFirstMap();
  map<int,int>::iterator iMap;

  int nVmap = 0;
  for(iMap=first.begin();iMap!=first.end();iMap++) {
    iCell = iMap->first;
    iC = iCell; ix = iC%N; iC/=N; iy = iC%N; iC/=N; iz = iC;

    // cell min x,y,z coords
    minCell.x = (((float)(N-ix  ))*x0+((float)(ix  ))*x1)/((float)N);
    minCell.y = (((float)(N-iy  ))*y0+((float)(iy  ))*y1)/((float)N);
    minCell.z = (((float)(N-iz  ))*z0+((float)(iz  ))*z1)/((float)N);

    // cell max x,y,z coords
    maxCell.x = (((float)(N-ix-1))*x0+((float)(ix+1))*x1)/((float)N);
    maxCell.y = (((float)(N-iy-1))*y0+((float)(iy+1))*y1)/((float)N);
    maxCell.z = (((float)(N-iz-1))*z0+((float)(iz+1))*z1)/((float)N);

    // fill the arrays iVB, F, and b, with the 8 cell vertex values
    for(h=0;h<8;h++) {
      h0 = (h  )%2; x = v[h].x = (h0==0)?minCell.x:maxCell.x; jx = ix+h0;
      h1 = (h/2)%2; y = v[h].y = (h1==0)?minCell.y:maxCell.y; jy = iy+h1;
      h2 = (h/4)%2; z = v[h].z = (h2==0)?minCell.z:maxCell.z; jz = iz+h2;
      
      iV = jx+(N+1)*(jy+(N+1)*jz);
      if(vMap.count(iV)==0) {
        vMap[iV] = nVmap++;
      }
      iVB[h] = iVmap = vMap[iV];
      F[h] = fGrid[iVmap];
      b[h] = (F[h]<0.0f);
    }

    // create the isovertices 
    for(i=0;i<12;i++) {
      iV   = -1;
      j    = _edgeTable[i][0];
      k    = _edgeTable[i][1];
      if(b[j]!=b[k]) {
        // look for the isovertex index associated with this grid edge
        iV = graphMap.get(iVB[j],iVB[k]);
        if(iV<0) {
          // create a new isovertex on first visit to this regular grid edge
          iV = static_cast<int>((coordIfs.size()/3));
          // compute the isovertex coordinates
          tk = F[j]/(F[j]-F[k]);
          tj = F[k]/(F[k]-F[j]);
          x  = tj*v[j].x+tk*v[k].x;
          y  = tj*v[j].y+tk*v[k].y;
          z  = tj*v[j].z+tk*v[k].z;
          coordIfs.push_back(x);
          coordIfs.push_back(y);
          coordIfs.push_back(z);
          // save the association in the graphMap
          graphMap.insert(iVB[j],iVB[k],iV);
        }
      }
      iE[i] = iV;
    }

    // make the faces of the isosurface within this cell
    /* nFacesCell = */ IsoSurf::makeCellFaces(b,iE,coordIndexIfs);
  }

  // add normal per face
  vector<float>& normalIfs = surface->getNormal();
  surface->setNormalPerVertex(false);
  Geometry::computeNormalsPerFace(coordIfs,coordIndexIfs,normalIfs);
  
}
