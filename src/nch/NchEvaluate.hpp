//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-07 20:24:34 taubin>
//------------------------------------------------------------------------
//
// NchEvaluate.hpp
//
// Software developed for the course
// Digital Geometry Processing
// Copyright (c) 2025, Gabriel Taubin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the Brown University nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL GABRIEL TAUBIN BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef NCH_EVALUATE_HPP
#define NCH_EVALUATE_HPP

#include<vector>
// #include <QObject>
#include <QThread>
#include <QProgressBar>
// #include "wrl/Types.hpp"
#include "wrl/IndexedFaceSet.hpp"
#include "core/HexGridPartition.hpp"

using namespace std;

class NchEvaluate : public QThread {
  Q_OBJECT
  
public:

  static void setProgressBar(QProgressBar* progressBar);
  
  NchEvaluate
  (const vector<float>& coord,
   const vector<float>& normal,
   vector<float>& rho,
   const Vec3f& center, const Vec3f& size,
   const int depth, const float scale, const bool cube,
   vector<float>& fGrid,
   IndexedFaceSet* surface,
   const bool multiThreaded = false,
   QObject *parent = nullptr);
  
  NchEvaluate
  (const vector<float>& coord,
   const vector<float>& normal,
   vector<float>& rho,
   HexGridPartition& hgp,
   vector<float>& fGrid,
   IndexedFaceSet* surface,
   const bool multiThreaded = false,
   QObject *parent = nullptr);
  
  // virtual ~NchEvaluate();
  
  virtual void run();

private:
    
  void _runRegular();
  void _runAdaptive();

signals:
  
  void progressReset();
  void progressSetRange(int min, int max);
  void progressSetValue(int value);
  
private:
  static QProgressBar* _progressBar;

  const vector<float>& _coord;
  const vector<float>& _normal;
  vector<float>&       _rho;
  // for regular grid
  const Vec3f          _center;
  const Vec3f          _size;
  const int            _depth;
  const float          _scale;
  const bool           _cube;
  // for adaptive grid
  HexGridPartition*    _hgpPtr;
  //
  vector<float>&       _fGrid;
  IndexedFaceSet*      _surface;
  const bool           _multiThreaded;
};

#endif // NCH_EVALUATE_HPP
