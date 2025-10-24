//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-07 20:23:47 taubin>
//------------------------------------------------------------------------
//
// GuiPanelRendering.hpp
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

#ifndef _GUI_PANEL_RENDERING_HPP_
#define _GUI_PANEL_RENDERING_HPP_

#include "ui_GuiPanelRendering.h"

#include <QResizeEvent>
#include <QLabel>
#include <QLineEdit>
#include <QSpinBox>
#include <QPushButton>
#include <QCheckBox>
#include "GuiPanel.hpp"

class GuiPanelRendering : public GuiPanel, public Ui::GuiPanelRendering {

  Q_OBJECT;

  static bool  _registered;
  static bool  registerPanel();

public:

               GuiPanelRendering(QWidget *parent = 0);
  virtual     ~GuiPanelRendering();

protected:

  virtual void mousePressEvent(QMouseEvent * event) Q_DECL_OVERRIDE;

private:

  void _setCheckBoxRow
    (QString title, QFont& font, QCheckBox*& cbName);
  void _setCheckBoxRowColor
    (QString title, QFont& font, QCheckBox*& cbName ,QLabel*& lHex, QLabel*& lColor);
  void _setCheckBoxRowColor
    (QString title, QFont& font, QLabel*& lName ,QLabel*& lHex, QLabel*& lColor);
  void _updateColor
    (QRgb rgb, QLabel* labelHex, QLabel* labelColor);
  QRgb _selectColor
    (QString title, QRgb rgb, QLabel* labelHex, QLabel* labelColor);

public slots:

  virtual void updateState() override; // c++11

private slots:

  void on_cbAllVertices_stateChanged(int state);
  void on_cbSelectedVertices_stateChanged(int state);
  void on_cbAllEdges_stateChanged(int state);
  void on_cbSelectedEdges_stateChanged(int state);
  void on_cbAllFaces_stateChanged(int state);
  void on_cbSelectedFaces_stateChanged(int state);
  void on_cbNormals_stateChanged(int state);
  void on_cbAllPolylines_stateChanged(int state);
  void on_cbSelectedPolylines_stateChanged(int state);
  void on_spinPointSize_valueChanged(int value);
  void on_spinNormalFactor_valueChanged(int value);

private:

  static const char* bg255;
  static const char* bg100;
  static const char* bg200;

};

#endif // _GUI_PANEL_RENDERING_HPP_
