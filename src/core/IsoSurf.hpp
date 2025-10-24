//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-07 20:22:21 taubin>
//------------------------------------------------------------------------
//
// IsoSurf.hpp
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

#ifndef _ISO_SURF_HPP_
#define _ISO_SURF_HPP_

#include<vector>
#include "wrl/Types.hpp"
#include "wrl/IndexedFaceSet.hpp"
#include "core/HexGridPartition.hpp"

using namespace std;

class IsoSurf {

public:

  static const int (*getEdgeTable())[2];

  static int makeCellFaces
    (bool b[/*8*/], int iE[/*12*/], vector<int>& coordIndex);

  static void computeIsosurface
  (const Vec3f& center, const Vec3f& size,
   const int depth, const float scale, const bool cube,
   vector<float>& fGrid, IndexedFaceSet& isoSurface);

  static void computeIsosurface
  (HexGridPartition& hgp,
   vector<float>& fGrid, IndexedFaceSet& isoSurface);

private:

  static const int _faceTable[12][12];
  static const int _edgeTable[12][2];

};

#endif // _ISO_SURF_HPP_
