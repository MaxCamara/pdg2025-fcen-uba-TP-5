//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-08-07 20:23:43 taubin>
//------------------------------------------------------------------------
//
// GuiPanelOptimization.cpp
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

#include <iostream>

#include "GuiApplication.hpp"
#include "GuiMainWindow.hpp"
#include "GuiPanelOptimization.hpp"
#include "wrl/SceneGraphProcessor.hpp"
#include "core/Geometry.hpp"

// static
Color GuiPanelOptimization::_diffuseColor(0.5f,0.7f,0.9f);

//////////////////////////////////////////////////////////////////////

// static variable evaluated at start-up time
bool GuiPanelOptimization::_registered = GuiPanelOptimization::registerPanel();

// static method run once at start-up time
bool GuiPanelOptimization::registerPanel() {
  return
    GuiMainWindow::registerPanelFactory
    ("Optimization",
     [](QWidget *parent)->GuiPanel*{ return new GuiPanelOptimization(parent);});
}

//////////////////////////////////////////////////////////////////////
GuiPanelOptimization::GuiPanelOptimization(QWidget* parent):
  GuiPanel(parent),
  _optimization()
{
  setupUi(this);

  GuiViewerData& data = getApp()->getMainWindow()->getData();
  _optimization.setSelectedVertexIndex(data.getSelectedVertexIndex());
  _optimization.setSelectedEdgeIndex(data.getSelectedEdgeIndex());
  _optimization.setSelectedFaceIndex(data.getSelectedFaceIndex());
  
  comboLaplacianSmoothing->blockSignals(true);
  comboLaplacianSmoothing->addItem("VERTEX COORDINATES");
  comboLaplacianSmoothing->addItem("FACE NORMALS");
  comboLaplacianSmoothing->blockSignals(false);
  
  comboCollapseEdgesIndependentSet->blockSignals(true);
  comboCollapseEdgesIndependentSet->addItem("2 VERTICES");  // index==0
  comboCollapseEdgesIndependentSet->addItem("4 VERTICES");  // index==1
  comboCollapseEdgesIndependentSet->addItem("8 VERTICES");  // index==2
  comboCollapseEdgesIndependentSet->blockSignals(false);

  comboSplitEdgesMode->blockSignals(true);
  comboSplitEdgesMode->addItem("LONG");
  comboSplitEdgesMode->addItem("ALL");
  comboSplitEdgesMode->addItem("SELECTED");
  comboSplitEdgesMode->blockSignals(false);

  IndexedFaceSet* ifsInput =
    _getNamedIndexedFaceSet("SURFACE",false,_diffuseColor);
  _optimization.setInput(ifsInput);

  IndexedFaceSet* ifsOptimized =
    _getNamedIndexedFaceSet("OPTIMIZED",false,_diffuseColor);
  if(ifsOptimized==(IndexedFaceSet*)0)
    ifsOptimized = _getNamedIndexedFaceSet("OPTIMIZED",true,_diffuseColor);
  _optimization.setOptimized(ifsOptimized);

  checkBoxShowInput->blockSignals(true);
  checkBoxShowInput->setChecked(false);
  checkBoxShowInput->blockSignals(false);
  _showNamed("SURFACE",false);

  checkBoxShowOutput->blockSignals(true);
  checkBoxShowOutput->setChecked(true);
  checkBoxShowOutput->blockSignals(false);
  _showNamed("OPTIMIZED",true);
}

//////////////////////////////////////////////////////////////////////
GuiPanelOptimization::~GuiPanelOptimization() {
}

//////////////////////////////////////////////////////////////////////
// overrides pure virtual GuiPanel::updateState()

void GuiPanelOptimization::updateState() {
  std::cout << "GuiPanelOptimization::updateState() {\n";

  // if OPTIMIZED Shape has been removed from the SceneGraph
  // somewhere else, reset _optimization
  IndexedFaceSet* ifsOptimized =
    _getNamedIndexedFaceSet("OPTIMIZED",false,_diffuseColor);
  if(ifsOptimized==(IndexedFaceSet*)0) {
    _optimization.clear();
    checkBoxShowInput->blockSignals(true);
    checkBoxShowInput->setChecked(true);
    checkBoxShowInput->blockSignals(false);
    checkBoxShowOutput->blockSignals(true);
    checkBoxShowOutput->setChecked(false);
    checkBoxShowOutput->blockSignals(false);
  }

  bool hasInput    = (_optimization.getInput()!=(IndexedFaceSet*)0);
  bool hasOptimized = (_optimization.getOptimized()!=(IndexedFaceSet*)0);
  if(hasOptimized) {
    pushButtonReset->setText("RESET");
  } else {
    pushButtonReset->setText("ADD");
  }
  pushButtonSave->setEnabled(hasInput && hasOptimized);
  pushButtonRemove->setEnabled(hasOptimized);
  pushButtonShowToggle->setEnabled(hasInput && hasOptimized);

  pushButtonLaplacianSmoothingRun->setEnabled(hasOptimized);
  pushButtonJacobiRun->setEnabled(hasOptimized);

  spinBoxSmoothingSteps->blockSignals(true);
  spinBoxSmoothingSteps->setValue(_optimization.getSteps());
  spinBoxSmoothingSteps->blockSignals(false);

  editLambda->setText
    ("  "+QString::number(_optimization.getLambda(),'f',6));
  editMu->setText
    ("  "+QString::number(_optimization.getMu(),'f',6));
  editKappa->setText
    ("  "+QString::number(_optimization.getKappa(),'f',6));
  editJacobiWeightData->setText
    ("  "+QString::number(_optimization.getJacobiWeightData(),'f',6));
  editJacobiWeightSmoothing->setText
    ("  "+QString::number(_optimization.getJacobiWeightSmoothing(),'f',6));

  pushButtonSplitEdgesShow->setEnabled(hasOptimized);
  pushButtonSplitEdgesApply->setEnabled(hasOptimized);
  pushButtonCollapseEdgesShow->setEnabled(hasOptimized);
  pushButtonCollapseEdgesApply->setEnabled(hasOptimized);
  pushButtonEqualizeValencesShow->setEnabled(hasOptimized);
  pushButtonEqualizeValencesApply->setEnabled(hasOptimized);
  pushButtonClusterVerticesApply->setEnabled(hasOptimized);
    
  if(hasInput) {

    float minEdgeLength = _optimization.getMinEdgeLength();
    
    editMinEdgeLength->setText
      ("  "+QString::number(minEdgeLength,'f',6));
    editMaxEdgeLength->setText
      ("  "+QString::number(_optimization.getMaxEdgeLength(),'f',6));
    editTargetEdgeLength->setText
      ("  "+QString::number(_optimization.getTargetEdgeLength(),'f',6));

    // vertex clustering

    GuiViewerData& data       = getApp()->getMainWindow()->getData();
    float          bboxScale  = data.getBBoxScale();
    SceneGraph*    pWrl       = data.getSceneGraph();
    Vec3f&         bboxCenter = pWrl->getBBoxCenter();
    Vec3f&         bboxSide   = pWrl->getBBoxSize();

    float side = bboxSide.x;
    if(bboxSide.y>side) side = bboxSide.y;
    if(bboxSide.z>side) side = bboxSide.z;
    side *= bboxScale;
    
    float clusterSize = 0.0f;
    float clusterResolution = 0.0f;

    if(_optimization.getQuantizationResolution()<=0) {
      clusterSize = minEdgeLength/3.0f;
      clusterResolution = std::ceil(side/clusterSize);
    } else {
      int value = _optimization.getQuantizationResolution();
      clusterResolution = static_cast<float>(value);
      clusterSize = side/clusterResolution;
    }

    Vec3f clusterMin(bboxCenter);
    clusterMin.x -= side/2;
    clusterMin.y -= side/2;
    clusterMin.z -= side/2;

    Vec3f clusterMax(bboxCenter);
    clusterMax.x += side/2;
    clusterMax.y += side/2;
    clusterMax.z += side/2;

    _optimization.setQuantizationBox
        (clusterMin,clusterMax);
    _optimization.setQuantizationResolution
        (clusterResolution);

    editClusterSize->setText
      ("  "+QString::number(clusterSize,'f',6));

    int value = static_cast<int>(clusterResolution);
    spinBoxClusterResolution->blockSignals(true);
    spinBoxClusterResolution->setValue(value);
    spinBoxClusterResolution->blockSignals(false);

    int nClusters = _optimization.getNumberOfClusters();
    editClusterNumber->setText("  "+QString("%1").arg(nClusters));

  } else {
    editMinEdgeLength->setText("  "+QString::number(0.0f,'f',6));
    editMaxEdgeLength->setText("  "+QString::number(0.0f,'f',6));
    editTargetEdgeLength->setText("  "+QString::number(0.0f,'f',6));

    editClusterSize->setText("  "+QString::number(0.0f,'f',6));
    spinBoxClusterResolution->blockSignals(true);
    spinBoxClusterResolution->setValue(1);
    spinBoxClusterResolution->blockSignals(true);

    editClusterNumber->setText("  0");
  }


  std::cout << "}\n";
}

//////////////////////////////////////////////////////////////////////
// panelShow

void GuiPanelOptimization::on_pushButtonReset_clicked() {
  on_pushButtonRemove_clicked();
  IndexedFaceSet* ifsInput =
    _getNamedIndexedFaceSet("SURFACE",false,_diffuseColor);
  _optimization.setInput(ifsInput);
  IndexedFaceSet* ifsOptimized =
    _getNamedIndexedFaceSet("OPTIMIZED",false,_diffuseColor);
  if(ifsOptimized == (IndexedFaceSet*)0)
    ifsOptimized = _getNamedIndexedFaceSet("OPTIMIZED",true,_diffuseColor);
  _optimization.setOptimized(ifsOptimized);
  _refresh3DView();
  checkBoxShowInput->setChecked(false);
  checkBoxShowOutput->setChecked(true);
}

void GuiPanelOptimization::on_pushButtonSave_clicked() {
  _optimization.saveOptimized();
  _refresh3DView();
}

void GuiPanelOptimization::on_pushButtonRemove_clicked() {
  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl!=(SceneGraph*)0) {
    _optimization.setOptimized(nullptr);
    SceneGraphProcessor processor(*pWrl);
    processor.removeSceneGraphChild("OPTIMIZED");
    mainWindow->setSceneGraph(pWrl,false);
    mainWindow->refresh();
    checkBoxShowInput->setChecked(true);
  }
}

void GuiPanelOptimization::on_checkBoxShowInput_stateChanged(int state) {
  bool surfaceChecked = (state!=0);
  checkBoxShowOutput->blockSignals(true);
  checkBoxShowOutput->setChecked(!surfaceChecked);
  checkBoxShowOutput->blockSignals(false);
  _showNamed("SURFACE",surfaceChecked);
  _showNamed("OPTIMIZED",!surfaceChecked);
}

void GuiPanelOptimization::on_checkBoxShowOutput_stateChanged(int state) {
  bool optimizedChecked = (state!=0);
  checkBoxShowInput->blockSignals(true);
  checkBoxShowInput->setChecked(!optimizedChecked);
  checkBoxShowInput->blockSignals(false);
  _showNamed("OPTIMIZED",optimizedChecked);
  _showNamed("SURFACE",!optimizedChecked);
}

void GuiPanelOptimization::on_pushButtonShowToggle_clicked() {
  bool value = checkBoxShowInput->isChecked();
  checkBoxShowInput->blockSignals(true);
  checkBoxShowInput->setChecked(!value);
  checkBoxShowInput->blockSignals(false);
  checkBoxShowOutput->blockSignals(true);
  checkBoxShowOutput->setChecked(value);
  checkBoxShowOutput->blockSignals(false);
  _showNamed("SURFACE",!value);
  _showNamed("OPTIMIZED",value);
}

void GuiPanelOptimization::on_checkBoxShowNormals_stateChanged(int state) {
  GuiMainWindow* mainWin = getApp()->getMainWindow();
  GuiViewerData& data = mainWin->getData();
  data.setPaintNormals(state!=0);
  mainWin->resetSceneGraph();
}

void GuiPanelOptimization::on_pushButtonNormalsRecompute_clicked() {
  IndexedFaceSet* ifsOptimized = _optimization.getOptimized();
  if(ifsOptimized==nullptr) return;
  ifsOptimized->clearNormal();
  ifsOptimized->setNormalPerVertex(false);
  Geometry::computeNormalsPerFace
    (ifsOptimized->getCoord(),
     ifsOptimized->getCoordIndex(),
     ifsOptimized->getNormal());
  _refresh3DView();
}

//////////////////////////////////////////////////////////////////////
// panelSmoothingParameters

void GuiPanelOptimization::on_spinBoxSmoothingSteps_valueChanged(int steps) {
  _optimization.setSteps(steps);
  // updateState();
}

void GuiPanelOptimization::on_editLambda_returnPressed() {
  QString str = editLambda->text();
  float value = str.toFloat();
  _optimization.setLambda(value);
  updateState();
}

void GuiPanelOptimization::on_editMu_returnPressed() {
  QString str = editMu->text();
  float value = str.toFloat();
  _optimization.setMu(value);
  updateState();
}

void GuiPanelOptimization::on_editKappa_returnPressed() {
  QString str = editKappa->text();
  float value = str.toFloat();
  _optimization.setKappa(value);
  updateState();
}

void GuiPanelOptimization::on_pushButtonMuToggleSign_clicked() {
  QString str = editMu->text();
  float value = str.toFloat();
  _optimization.setMu(-value);
  updateState();
}

//////////////////////////////////////////////////////////////////////
// panelLaplacianSmoothing

void GuiPanelOptimization::on_pushButtonLaplacianSmoothingRun_clicked() {
  QString type = comboLaplacianSmoothing->currentText();
  if(type=="VERTEX COORDINATES") {
    _optimization.laplacianSmoothingVertexCoordinatesRun();
  } else if(type=="FACE NORMALS") {
    const bool normalize = true;
    _optimization.laplacianSmoothingFaceNormalsRun(normalize);
  }
  _refresh3DView();
}

//////////////////////////////////////////////////////////////////////
// panelJacobi

void GuiPanelOptimization::on_editJacobiWeightData_returnPressed() {
  QString str = editJacobiWeightData->text();
  float value = str.toFloat();
  _optimization.setJacobiWeightData(value);
  updateState();
}

void GuiPanelOptimization::on_editJacobiWeightSmoothing_returnPressed() {
  QString str = editJacobiWeightSmoothing->text();
  float value = str.toFloat();
  _optimization.setJacobiWeightSmoothing(value);
  updateState();
}

void GuiPanelOptimization::on_pushButtonJacobiRun_clicked() {
  _optimization.jacobiRun();
  _refresh3DView();
}

//////////////////////////////////////////////////////////////////////
// panelEdgeLengths

void GuiPanelOptimization::on_editTargetEdgeLength_returnPressed() {
  QString str = editTargetEdgeLength->text();
  float value = str.toFloat();
  _optimization.setTargetEdgeLength(value);
  updateState();
}

//////////////////////////////////////////////////////////////////////
// panelClusterVertices

void GuiPanelOptimization::on_spinBoxClusterResolution_valueChanged() {
  int value = spinBoxClusterResolution->value();
  _optimization.setQuantizationResolution(value);
  float clusterResolution = static_cast<float>(value);

  GuiViewerData& data       = getApp()->getMainWindow()->getData();
  float          bboxScale  = data.getBBoxScale();
  SceneGraph*    pWrl       = data.getSceneGraph();
  Vec3f&         bboxSide   = pWrl->getBBoxSize();
  float          cubeSide   = bboxSide.x;
  if(bboxSide.y>cubeSide) cubeSide = bboxSide.y;
  if(bboxSide.z>cubeSide) cubeSide = bboxSide.z;
  cubeSide *= bboxScale;
  
  float clusterSize = cubeSide/clusterResolution;
  editClusterSize->setText
    ("  "+QString::number(clusterSize,'f',6));
}

void GuiPanelOptimization::on_pushButtonClusterVerticesApply_clicked() {
  _optimization.clusterVerticesApply();
  _refresh3DView();
}

//////////////////////////////////////////////////////////////////////
// panelCollapseEdges

Optimization::EdgeCollapseIndependentSet
GuiPanelOptimization::_getIndependentSet() {
  Optimization::EdgeCollapseIndependentSet indepSet =
    Optimization::EdgeCollapseIndependentSet::VERTICES_2;
  QString txt = comboCollapseEdgesIndependentSet->currentText();
  if(txt=="4 VERTICES")
    indepSet = Optimization::EdgeCollapseIndependentSet::VERTICES_4;
  else if(txt=="8 VERTICES")
    indepSet = Optimization::EdgeCollapseIndependentSet::VERTICES_8;
  return indepSet;
}

void GuiPanelOptimization::on_comboCollapseEdgesIndependentSet_currentIndexChanged
(int index) {
  (void)index;
  on_pushButtonCollapseEdgesShow_clicked();
}

void GuiPanelOptimization::on_pushButtonCollapseEdgesShow_clicked() {
  GuiViewerData& data = getApp()->getMainWindow()->getData();
  _optimization.setSelectedVertexIndex(data.getSelectedVertexIndex());
  _optimization.setSelectedEdgeIndex(data.getSelectedEdgeIndex());
  _optimization.setSelectedFaceIndex(data.getSelectedFaceIndex());
  _optimization.collapseEdgesShow(_getIndependentSet());
  _refresh3DView();
}

void GuiPanelOptimization::on_pushButtonCollapseEdgesApply_clicked() {
  GuiViewerData& data = getApp()->getMainWindow()->getData();
  _optimization.setSelectedVertexIndex(data.getSelectedVertexIndex());
  _optimization.setSelectedEdgeIndex(data.getSelectedEdgeIndex());
  _optimization.setSelectedFaceIndex(data.getSelectedFaceIndex());
  _optimization.collapseEdgesApply(_getIndependentSet());
  _refresh3DView();
}

//////////////////////////////////////////////////////////////////////
// panelSplitEdges

void GuiPanelOptimization::on_comboSplitEdgesMode_currentIndexChanged(int index) {
  (void) index;
  on_pushButtonSplitEdgesShow_clicked();
}

Optimization::SplitEdgesMode GuiPanelOptimization::_getSplitEdgesMode() {
  Optimization::SplitEdgesMode mode = Optimization::SplitEdgesMode::ALL;
  QString txt = comboSplitEdgesMode->currentText();
  if(txt=="SELECTED")
    mode = Optimization::SplitEdgesMode::SELECTED;
  else if(txt=="LONG")
    mode = Optimization::SplitEdgesMode::LONG;
  return mode;
}

void GuiPanelOptimization::on_pushButtonSplitEdgesShow_clicked() {
  Optimization::SplitEdgesMode mode = _getSplitEdgesMode();
  GuiViewerData& data = getApp()->getMainWindow()->getData();
  _optimization.setSelectedVertexIndex(data.getSelectedVertexIndex());
  _optimization.setSelectedEdgeIndex(data.getSelectedEdgeIndex());
  _optimization.setSelectedFaceIndex(data.getSelectedFaceIndex());
  _optimization.adaptiveSubdivisionShow(mode);
  _refresh3DView();
}

void GuiPanelOptimization::on_pushButtonSplitEdgesApply_clicked() {
  Optimization::SplitEdgesMode mode = _getSplitEdgesMode();
  GuiViewerData& data = getApp()->getMainWindow()->getData();
  _optimization.setSelectedVertexIndex(data.getSelectedVertexIndex());
  _optimization.setSelectedEdgeIndex(data.getSelectedEdgeIndex());
  _optimization.setSelectedFaceIndex(data.getSelectedFaceIndex());
  _optimization.adaptiveSubdivisionApply(mode,true /* paintFaces */);
  _refresh3DView();
}

//////////////////////////////////////////////////////////////////////
// panelEqualizeValences

void GuiPanelOptimization::on_pushButtonEqualizeValencesShow_clicked() {
  _optimization.equalizeValencesShow();
  _refresh3DView();
}

void GuiPanelOptimization::on_pushButtonEqualizeValencesApply_clicked() {
  _optimization.equalizeValencesApply();
  _refresh3DView();
}

//////////////////////////////////////////////////////////////////////
void GuiPanelOptimization::mousePressEvent(QMouseEvent * event) {

  int x = event->position().x();
  int y = event->position().y();

  bool clickedOnLabel = false;

  if(labelShow->geometry().contains(x,y,true)) {
    panelShow->setVisible(panelShow->isHidden());
    clickedOnLabel = true;
  } else if(labelSmoothingParameters->geometry().contains(x,y,true)) {
    panelSmoothingParameters->setVisible(panelSmoothingParameters->isHidden());
    clickedOnLabel = true;
  } else if(labelLaplacianSmoothing->geometry().contains(x,y,true)) {
    panelLaplacianSmoothing->setVisible(panelLaplacianSmoothing->isHidden());
    clickedOnLabel = true;
  } else if(labelJacobi->geometry().contains(x,y,true)) {
    panelJacobi->setVisible(panelJacobi->isHidden());
    clickedOnLabel = true;
  } else if(labelEdgeLengths->geometry().contains(x,y,true)) {
    panelEdgeLengths->setVisible(panelEdgeLengths->isHidden());
    clickedOnLabel = true;
  } else if(labelClusterVertices->geometry().contains(x,y,true)) {
    panelClusterVertices->setVisible(panelClusterVertices->isHidden());
    clickedOnLabel = true;
  } else if(labelCollapseEdges->geometry().contains(x,y,true)) {
    panelCollapseEdges->setVisible(panelCollapseEdges->isHidden());
    clickedOnLabel = true;
  } else if(labelSplitEdges->geometry().contains(x,y,true)) {
    panelSplitEdges->setVisible(panelSplitEdges->isHidden());
    clickedOnLabel = true;
  } else if(labelEqualizeValences->geometry().contains(x,y,true)) {
    panelEqualizeValences->setVisible(panelEqualizeValences->isHidden());
    clickedOnLabel = true;
  }
  
  if(clickedOnLabel) updateState();

}

