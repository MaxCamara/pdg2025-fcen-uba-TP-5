//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin / 3D Shape Tech LLC
//  Time-stamp: <2025-08-07 20:24:28 taubin>
//------------------------------------------------------------------------
//
// TokenizerFile.cpp
//
// Software developed for the course
// Digital Geometry Processing
// Copyright (c) 2025, Gabriel Taubin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
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

#include <stdio.h>
#include "TokenizerFile.hpp"

TokenizerFile::TokenizerFile(FILE* fp):
  Tokenizer(),
  _fp(fp) {
}

char TokenizerFile::getc() {
  return static_cast<char>(std::getc(_fp)); // c library function
}

// #define LINE_BUFFER_LENGTH 1024
//
// bool TokenizerFile::getline() {
//   static char buffer[LINE_BUFFER_LENGTH];
//   bool success = false;
//   clear();
//   if(fgets(buffer,LINE_BUFFER_LENGTH,_fp)!=nullptr) {
//     for(size_t i=0;i<LINE_BUFFER_LENGTH-1;i++) {
//       if(buffer[i]=='\0' || buffer[i]=='\n') break;
//       push_back(buffer[i]);
//     }
//     success = true;
//   }
//   return success;
// }
