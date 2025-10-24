//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-10-22 15:16:28 taubin>
//------------------------------------------------------------------------
//
// GuiPanelImplicit.hpp
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

#ifndef _GUI_PANEL_IMPLICIT_HPP_
#define _GUI_PANEL_IMPLICIT_HPP_

#include "ui_GuiPanelImplicit.h"

#include <QResizeEvent>
#include <QLabel>
#include <QLineEdit>
#include <QSpinBox>
#include <QPushButton>
#include <QCheckBox>
#include "GuiPanel.hpp"
#include "core/HexGridPartition.hpp"

class GuiMainWindow;

class GuiPanelImplicit : public GuiPanel, public Ui::GuiPanelImplicit {

  Q_OBJECT;

  static bool _registered;
  static bool registerPanel();

public:

  GuiPanelImplicit(QWidget *parent = 0);
  virtual ~GuiPanelImplicit();

public slots:

  virtual void updateState() override; // overrides virtual GuiPanel::updateState()

protected:

  virtual void mousePressEvent(QMouseEvent * event) Q_DECL_OVERRIDE;


private slots:

  // Hexahedral Grid
  void hexGridSetDepth(const int depth);
  void hexGridDepthUp();
  void hexGridDepthDown();
  void on_spinBoxHexGridDepth_valueChanged(int depth);
  void on_pushButtonHexGridAdd_clicked();
  void on_pushButtonHexGridRemove_clicked();
  void on_editHexGridScale_returnPressed();
  void on_checkBoxHexGridCube_stateChanged(int state);
  void on_checkBoxHexGridOccupiedCells_stateChanged(int state);

  // Points
  // bool pointsHasNormals();
  // bool pointsHasColors();
  void on_pushButtonPointsRemove_clicked();
  void on_pushButtonPointsShow_clicked();
  void on_pushButtonPointsHide_clicked();
  void on_checkBoxPointsHasNormals_stateChanged();
  void on_checkBoxPointsHasColors_stateChanged();

  // Samples
  // void on_comboSampleSurface_currentIndexChanged(int index);
  void on_pushButtonSampleSurface_clicked();
  void on_pushButtonSamplePoints_clicked();
  void on_pushButtonSamplesSwap_clicked();
  void on_pushButtonSamplesRemove_clicked();
  void on_pushButtonSamplesShow_clicked();
  void on_pushButtonSamplesHide_clicked();

  // Surface
  void on_pushButtonSurfaceRemove_clicked();
  void on_pushButtonSurfaceShow_clicked();
  void on_pushButtonSurfaceHide_clicked();

  // Plane
  void on_pushButtonPlaneRemove_clicked();
  void on_pushButtonPlaneShow_clicked();
  void on_pushButtonPlaneHide_clicked();
  void on_pushButtonPlaneFit_clicked();

  // Surface Reconstruction
  void on_pushButtonSurfRecMultiplePlanes_clicked();
  void on_pushButtonSurfRecContinuous_clicked();
  void on_pushButtonSurfRecWatertight_clicked();
  void on_pushButtonSurfRecJacobi_clicked();

  // Non Convex Hull
  void on_comboNchSource_currentIndexChanged(int index);
  void on_pushButtonNchFit_clicked();
  void on_pushButtonNchReconstruct_clicked();
  void on_pushButtonNchFitMultiThreaded_clicked();
  void on_pushButtonNchReconstructMultiThreaded_clicked();

private:
  HexGridPartition *_hexGridPartition;
};

#endif // _GUI_PANEL_IMPLICIT_HPP_
