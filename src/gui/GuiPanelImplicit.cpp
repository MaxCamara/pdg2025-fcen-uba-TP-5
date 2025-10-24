//------------------------------------------------------------------------
//  Copyright (C) Gabriel Taubin
//  Time-stamp: <2025-10-22 15:17:53 taubin>
//------------------------------------------------------------------------
//
// GuiPanelImplicit.cpp
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

#include <iostream>

#include "GuiApplication.hpp"
#include "GuiMainWindow.hpp"
#include "GuiPanelImplicit.hpp"
#include "nch/NchProcessor.hpp"
#include "wrl/SceneGraphProcessor.hpp"
#include "core/Geometry.hpp"

//////////////////////////////////////////////////////////////////////

// static variable evaluated at start-up time
bool GuiPanelImplicit::_registered = GuiPanelImplicit::registerPanel();

// static method run once at start-up time
bool GuiPanelImplicit::registerPanel() {
  return
    GuiMainWindow::registerPanelFactory
    ("Implicit",
     [](QWidget *parent)->GuiPanel*{ return new GuiPanelImplicit(parent);});
}

//////////////////////////////////////////////////////////////////////
GuiPanelImplicit::GuiPanelImplicit(QWidget* parent):
  GuiPanel(parent),
  _hexGridPartition(nullptr) {
  setupUi(this);

  comboSampleSurface->blockSignals(true);
  comboSampleSurface->addItem("FACES");
  comboSampleSurface->addItem("VERTICES");
  comboSampleSurface->blockSignals(false);

  comboNchSource->blockSignals(true);
  comboNchSource->addItem("POINTS");
  comboNchSource->addItem("SAMPLES");
  comboNchSource->blockSignals(false);
}

//////////////////////////////////////////////////////////////////////
GuiPanelImplicit::~GuiPanelImplicit() {
  if(_hexGridPartition!=(HexGridPartition*)0) {
    delete _hexGridPartition;
  }
}

//////////////////////////////////////////////////////////////////////
void GuiPanelImplicit::updateState() {

  // std::cout << "GuiPanelImplicit::updateState() {\n";

  auto mainWindow = getApp()->getMainWindow();
  if(mainWindow) {

    GuiViewerData& data  = mainWindow->getData();

    int   bboxDepth = data.getBBoxDepth();
    bool  bboxCube  = data.getBBoxCube();
    float bboxScale = data.getBBoxScale();

    SceneGraph* pWrl = data.getSceneGraph();
    if(pWrl!=(SceneGraph*)0) {
      Vec3f& bboxCenter = pWrl->getBBoxCenter();
      editBBoxCenterX->setText("  "+QString::number(bboxCenter.x,'f',4));
      editBBoxCenterY->setText("  "+QString::number(bboxCenter.y,'f',4));
      editBBoxCenterZ->setText("  "+QString::number(bboxCenter.z,'f',4));
      Vec3f& bboxSide   = pWrl->getBBoxSize();
      editBBoxSideX->setText("  "+QString::number(bboxSide.x,'f',4));
      editBBoxSideY->setText("  "+QString::number(bboxSide.y,'f',4));
      editBBoxSideZ->setText("  "+QString::number(bboxSide.z,'f',4));
    } else {
      editBBoxCenterX->setText("  "+QString::number(0.0f,'f',2));
      editBBoxCenterY->setText("  "+QString::number(0.0f,'f',2));
      editBBoxCenterZ->setText("  "+QString::number(0.0f,'f',2));
      editBBoxSideX->setText("  "+QString::number(0.0f,'f',2));
      editBBoxSideY->setText("  "+QString::number(0.0f,'f',2));
      editBBoxSideZ->setText("  "+QString::number(0.0f,'f',2));
    }

    spinBoxHexGridDepth->setValue(bboxDepth);
    checkBoxHexGridCube->setChecked(bboxCube);
    editHexGridScale->setText("  "+QString::number(bboxScale,'f',2));

    Vec4f& plane = data.getPlane();
    editPlaneCoefficient0->setText("  "+QString::number(plane.x));
    editPlaneCoefficient1->setText("  "+QString::number(plane.y));
    editPlaneCoefficient2->setText("  "+QString::number(plane.z));
    editPlaneCoefficient3->setText("  "+QString::number(plane.w));

    int nHexGridCells = 0;
    int nHexGridVertices = 0;
    if(checkBoxHexGridOccupiedCells->isChecked()) {
      if(_hexGridPartition!=(HexGridPartition*)0) {
        nHexGridCells = _hexGridPartition->getNumberOfCells();
        nHexGridVertices = _hexGridPartition->getNumberOfVertices();
      }
    } else {
      int N = 1<<bboxDepth;
      nHexGridCells = N*N*N;
      nHexGridVertices = (N+1)*(N+1)*(N+1);
    }
    editHexGridCells->setText("  "+QString::number(nHexGridCells));
    editHexGridVertices->setText("  "+QString::number(nHexGridVertices));

    SceneGraph* wrl = data.getSceneGraph();
    if(wrl==(SceneGraph*)0) {

      pushButtonHexGridAdd->setEnabled(false);
      pushButtonHexGridRemove->setEnabled(false);
      
      editPointsStatus->setText("SCENE GRAPH IS EMPTY");
      pushButtonPointsRemove->setEnabled(false);
      pushButtonPointsShow->setEnabled(false);
      pushButtonPointsHide->setEnabled(false);
      editPointsStatus->setText("SCENE GRAPH IS EMPTY");
      pushButtonSurfaceRemove->setEnabled(false);
      pushButtonSurfaceShow->setEnabled(false);
      pushButtonSurfaceHide->setEnabled(false);

      pushButtonPlaneRemove->setEnabled(false);
      pushButtonPlaneShow->setEnabled(false);
      pushButtonPlaneHide->setEnabled(false);

      pushButtonSampleSurface->setEnabled(false);
      pushButtonSamplesSwap->setEnabled(false);
      pushButtonSamplesRemove->setEnabled(false);
      pushButtonSamplesShow->setEnabled(false);
      pushButtonSamplesHide->setEnabled(false);
      
      pushButtonNchFit->setEnabled(false);
      pushButtonNchReconstruct->setEnabled(false);
      pushButtonNchFitMultiThreaded->setEnabled(false);
      pushButtonNchReconstructMultiThreaded->setEnabled(false);

      pushButtonPlaneFit->setEnabled(false);
      pushButtonSurfRecMultiplePlanes->setEnabled(false);
      pushButtonSurfRecContinuous->setEnabled(false);
      pushButtonSurfRecWatertight->setEnabled(false);
      pushButtonSurfRecJacobi->setEnabled(false);
      
    } else /* if(wrl!=(SceneGraph*)0) */ {

      SceneGraphProcessor processor(*wrl);

      bool hasGrid   = processor.hasGrid();
      pushButtonHexGridAdd->setEnabled(!hasGrid);
      pushButtonHexGridRemove->setEnabled(hasGrid);

      Node* plane = wrl->find("PLANE"); // should be a Shape node
      bool  hasPlane = (plane!=(Node*)0 && plane->isShape());
      if(hasPlane) {
        pushButtonPlaneRemove->setEnabled(true);
        bool show = plane->getShow();
        pushButtonPlaneShow->setEnabled(!show);
        pushButtonPlaneHide->setEnabled(show);
      } else {
        pushButtonPlaneRemove->setEnabled(false);
        pushButtonPlaneShow->setEnabled(false);
        pushButtonPlaneHide->setEnabled(false);
      }

      Node* points = wrl->find("POINTS"); // should be a Shape node
      bool  hasPoints = (points!=(Node*)0 && points->isShape());
      bool  hasNchRhosPoints = false;

      Node* samples = wrl->find("SAMPLES"); // should be a Shape node
      bool  hasSamples = (samples!=(Node*)0 && samples->isShape());
      bool  hasNchRhosSamples = false;

      if(hasPoints) {
        editPointsStatus->setText("SCENE GRAPH HAS A \"POINTS\" NODE");

        pushButtonPointsRemove->setEnabled(true);
        bool show = points->getShow();
        pushButtonPointsShow->setEnabled(!show);
        pushButtonPointsHide->setEnabled(show);

        bool hasPointNormals = false;
        bool hasPointColors  = false;

        Shape* shape = (Shape*)points;
        Node*  node = shape->getGeometry();
        if(node!=(Node*)0 && node->isIndexedFaceSet()) {
          IndexedFaceSet* ifs = (IndexedFaceSet*)node;
          hasPointNormals = (ifs->getNormalBinding()==IndexedFaceSet::PB_PER_VERTEX);
          hasPointColors  = (ifs->getColorBinding() ==IndexedFaceSet::PB_PER_VERTEX);
          hasNchRhosPoints     = (ifs->getVariable("nchRhos")!=nullptr);
        }

        checkBoxPointsHasNormals->blockSignals(true);
        checkBoxPointsHasNormals->setChecked(hasPointNormals);
        checkBoxPointsHasNormals->blockSignals(false);

        checkBoxPointsHasColors->blockSignals(true);
        checkBoxPointsHasColors->setChecked(hasPointColors);
        checkBoxPointsHasColors->blockSignals(false);

      } else {
        editPointsStatus->setText("NO \"POINTS\" NODE FOUND IN SCENE GRAPH");
        pushButtonPointsRemove->setEnabled(false);
        pushButtonPointsShow->setEnabled(false);
        pushButtonPointsHide->setEnabled(false);
      }

      if(hasSamples) {
        editSamplesStatus->setText("SCENE GRAPH HAS A \"SAMPLES\" NODE");
        pushButtonSamplesRemove->setEnabled(true);
        bool show = samples->getShow();
        pushButtonSamplesShow->setEnabled(!show);
        pushButtonSamplesHide->setEnabled(show);
        Shape* shape = (Shape*)samples;
        Node*  node = shape->getGeometry();
        if(node!=(Node*)0 && node->isIndexedFaceSet()) {
          IndexedFaceSet* ifs = (IndexedFaceSet*)node;
          hasNchRhosSamples     = (ifs->getVariable("nchRhos")!=nullptr);
        }
      } else {
        editSamplesStatus->setText("NO \"SAMPLES\" NODE FOUND IN SCENE GRAPH");
        pushButtonSamplesRemove->setEnabled(false);
        pushButtonSamplesShow->setEnabled(false);
        pushButtonSamplesHide->setEnabled(false);
      }
      
      pushButtonSamplesSwap->setEnabled(hasPoints && hasSamples);

      QString source = comboNchSource->currentText();
      if(source=="POINTS") {
        pushButtonNchFit->setEnabled(hasPoints);
        pushButtonNchReconstruct->setEnabled(hasNchRhosPoints);
        pushButtonNchFitMultiThreaded->setEnabled(hasPoints);
        pushButtonNchReconstructMultiThreaded->setEnabled(hasNchRhosPoints);
      } else /* if(source=="SAMPLES") */ {
        pushButtonNchFit->setEnabled(hasSamples);
        pushButtonNchReconstruct->setEnabled(hasNchRhosSamples);
        pushButtonNchFitMultiThreaded->setEnabled(hasSamples);
        pushButtonNchReconstructMultiThreaded->setEnabled(hasNchRhosSamples);
      }

      pushButtonPlaneFit->setEnabled(hasPoints);
      pushButtonSurfRecMultiplePlanes->setEnabled(hasPoints);
      pushButtonSurfRecContinuous->setEnabled(hasPoints);
      pushButtonSurfRecWatertight->setEnabled(hasPoints);
      pushButtonSurfRecJacobi->setEnabled(hasPoints);

      Node* surface    = wrl->find("SURFACE"); // should be a Shape node
      bool  hasSurface = (surface!=(Node*)0 && surface->isShape());
      if(hasSurface) {
        editSurfaceStatus->setText("SCENE GRAPH HAS A \"SURFACE\" NODE");
        pushButtonSurfaceRemove->setEnabled(true);
        bool show = surface->getShow();
        pushButtonSurfaceShow->setEnabled(!show);;
        pushButtonSurfaceHide->setEnabled(show);;
        pushButtonSampleSurface->setEnabled(true);
      } else {
        editSurfaceStatus->setText("NO \"SURFACE\" NODE FOUND IN SCENE GRAPH");
        pushButtonSurfaceRemove->setEnabled(false);
        pushButtonSurfaceShow->setEnabled(false);;
        pushButtonSurfaceHide->setEnabled(false);;
        pushButtonSampleSurface->setEnabled(false);
      }

      // get number of points
      int nPts = 0;
      if(hasPoints) {
        Shape* shape = (Shape*)points;
        Node* node = shape->getGeometry();
        if(node!=0 && node->isIndexedFaceSet()) {
          IndexedFaceSet* ifs = (IndexedFaceSet*)node;
          nPts = ifs->getNumberOfCoord();
        }
      }
      editPointsNumber->setText("  "+QString::number(nPts));

      // get number of samples
      int nSamples = 0;
      if(hasSamples) {
        Shape* shape = (Shape*)samples;
        Node* node = shape->getGeometry();
        if(node!=0 && node->isIndexedFaceSet()) {
          IndexedFaceSet* ifs = (IndexedFaceSet*)node;
          nSamples = ifs->getNumberOfCoord();
        }
      }
      editSamplesNumber->setText("  "+QString::number(nSamples));

      // get number of surface vertices and faces
      int nVertices = 0;
      int nFaces = 0;
      if(hasSurface) {
        Shape* shape = (Shape*)surface;
        Node* node = shape->getGeometry();
        if(node!=0 && node->isIndexedFaceSet()) {
          IndexedFaceSet* ifs = (IndexedFaceSet*)node;
          nVertices = ifs->getNumberOfCoord();
          nFaces = ifs->getNumberOfFaces();
        }
      }
      editSurfaceVertices->setText("  "+QString::number(nVertices));
      editSurfaceFaces->setText("  "+QString::number(nFaces));

    }
  }
  
  // std::cout << "}\n";
}

//////////////////////////////////////////////////////////////////////
void GuiPanelImplicit::mousePressEvent(QMouseEvent * event) {

  int x = event->position().x();
  int y = event->position().y();

  bool clickedOnLabel = false;

  if(labelBBox->geometry().contains(x,y,true)) {
    panelBBox->setVisible(panelBBox->isHidden());
    clickedOnLabel = true;
  } else if(labelHexGrid->geometry().contains(x,y,true)) {
    panelHexGrid->setVisible(panelHexGrid->isHidden());
    clickedOnLabel = true;
  } else if(labelPoints->geometry().contains(x,y,true)) {
    panelPoints->setVisible(panelPoints->isHidden());
    clickedOnLabel = true;
  } else if(labelSamples->geometry().contains(x,y,true)) {
    panelSamples->setVisible(panelSamples->isHidden());
    clickedOnLabel = true;
  } else if(labelSurface->geometry().contains(x,y,true)) {
    panelSurface->setVisible(panelSurface->isHidden());
    clickedOnLabel = true;
  } else if(labelPlane->geometry().contains(x,y,true)) {
    panelPlane->setVisible(panelPlane->isHidden());
    clickedOnLabel = true;
  } else if(labelNonConvexHull->geometry().contains(x,y,true)) {
    panelNonConvexHull->setVisible(panelNonConvexHull->isHidden());
    clickedOnLabel = true;
  } else if(labelSurfaceReconstruction->geometry().contains(x,y,true)) {
    panelSurfaceReconstruction->setVisible
      (panelSurfaceReconstruction->isHidden());
    clickedOnLabel = true;
  }

  if(clickedOnLabel) updateState();
}

void GuiPanelImplicit::on_spinBoxHexGridDepth_valueChanged(int depth) {
  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  int newDepth = (depth<0)?0:(depth>10)?10:depth;
  int oldDepth = data.getBBoxDepth();
  if(newDepth!=oldDepth) {
    data.setBBoxDepth(depth);
    // is scene graph has a GRID node => rebuild it
    SceneGraph* pWrl = data.getSceneGraph();
    SceneGraphProcessor processor(*pWrl);
    if(processor.hasGrid() && _hexGridPartition==nullptr) {
      float scale = data.getBBoxScale();
      bool  cube  = data.getBBoxCube();
      processor.gridAdd(newDepth,scale,cube);
      mainWindow->setSceneGraph(pWrl,false);
      mainWindow->refresh();
    }
  }
  updateState();
}

void GuiPanelImplicit::hexGridSetDepth(const int depth) {
  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  int newDepth = (depth<0)?0:(depth>10)?10:depth;
  int oldDepth = data.getBBoxDepth();
  if(newDepth!=oldDepth) {
    data.setBBoxDepth(newDepth);
    // is scene graph has a GRID node => rebuild it
    SceneGraph* pWrl = data.getSceneGraph();
    SceneGraphProcessor processor(*pWrl);
    if(processor.hasGrid() && _hexGridPartition==nullptr) {
      float scale = data.getBBoxScale();
      bool  cube  = data.getBBoxCube();
      processor.gridAdd(newDepth,scale,cube);
      mainWindow->setSceneGraph(pWrl,false);
      mainWindow->refresh();
    }
    updateState();
  }
}

void GuiPanelImplicit::hexGridDepthUp() {
  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  on_spinBoxHexGridDepth_valueChanged(data.getBBoxDepth()+1);
}

void GuiPanelImplicit::hexGridDepthDown() {
  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  on_spinBoxHexGridDepth_valueChanged(data.getBBoxDepth()-1);
}

void GuiPanelImplicit::on_pushButtonHexGridAdd_clicked() {

  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  int depth = data.getBBoxDepth();
  SceneGraph* pWrl = data.getSceneGraph();
  SceneGraphProcessor processor(*pWrl);
  data.setBBoxDepth(depth);
  float scale = data.getBBoxScale();
  bool  cube  = data.getBBoxCube();

  bool gridOccupiedCells = checkBoxHexGridOccupiedCells->isChecked();
  if(gridOccupiedCells && _hexGridPartition!=nullptr)
    processor.gridAdd(*_hexGridPartition); // occupied cells of _hexGridPartition
  else
    processor.gridAdd(depth,scale,cube); // full regular grid
  
  mainWindow->setSceneGraph(pWrl,false);
  mainWindow->refresh();
  updateState();
}

void GuiPanelImplicit::on_pushButtonHexGridRemove_clicked() {
  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  SceneGraph* pWrl = data.getSceneGraph();
  if(pWrl!=(SceneGraph*)0) {
    SceneGraphProcessor processor(*pWrl);

    processor.gridRemove();

    mainWindow->setSceneGraph(pWrl,false);
    mainWindow->refresh();
    updateState();
  }
}

void GuiPanelImplicit::on_editHexGridScale_returnPressed() {
  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  float scale = data.getBBoxScale();
  QString str = "1.05"; // _editHexGridScale->text();
  float value = str.toFloat();
  if(value!=scale) {
    data.setBBoxScale(value);
    SceneGraphProcessor processor(*(data.getSceneGraph()));
    if(processor.hasGrid()) {
      int   depth = data.getBBoxDepth();
      float scale = data.getBBoxScale();
      bool  cube  = data.getBBoxCube();
      processor.gridAdd(depth,scale,cube);
      mainWindow->setSceneGraph(data.getSceneGraph(),false);
      mainWindow->refresh();
      updateState();
    }
  }
}

void GuiPanelImplicit::on_checkBoxHexGridCube_stateChanged(int state) {
  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  data.setBBoxCube((state!=0));
  SceneGraphProcessor processor(*(data.getSceneGraph()));
  if(processor.hasGrid()) {
    int   depth = data.getBBoxDepth();
    float scale = data.getBBoxScale();
    bool  cube  = data.getBBoxCube();
    processor.gridAdd(depth,scale,cube);
    mainWindow->setSceneGraph(data.getSceneGraph(),false);
    mainWindow->refresh();
    updateState();
  }
}

void GuiPanelImplicit::on_checkBoxHexGridOccupiedCells_stateChanged(int state) {
  (void) state;
  updateState();
}


void GuiPanelImplicit::on_pushButtonPointsRemove_clicked() {
  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl!=(SceneGraph*)0) {
    SceneGraphProcessor processor(*pWrl);
    processor.pointsRemove();
    mainWindow->setSceneGraph(pWrl,false);
    mainWindow->refresh();
    updateState();
  }
}

void GuiPanelImplicit::on_pushButtonPointsShow_clicked() {
  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl==(SceneGraph*)0) return;
  Node* node = pWrl->find("POINTS");
  if(node==(Node*)0) return;
  node->setShow(true);
  mainWindow->setSceneGraph(pWrl,false);
  mainWindow->refresh();
  updateState();
}

void GuiPanelImplicit::on_pushButtonPointsHide_clicked() {
  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl==(SceneGraph*)0) return;
  Node* node = pWrl->find("POINTS");
  if(node==(Node*)0) return;
  node->setShow(false);
  mainWindow->setSceneGraph(pWrl,false);
  mainWindow->refresh();
  updateState();
}

void GuiPanelImplicit::on_checkBoxPointsHasNormals_stateChanged() {
  // bool state = checkBoxPointsHasNormals->isChecked();
  // checkBoxPointsHasNormals->blockSignals(true);
  // checkBoxPointsHasNormals->setChecked(!state);
  // checkBoxPointsHasNormals->blockSignals(false);
  updateState();
}

void GuiPanelImplicit::on_checkBoxPointsHasColors_stateChanged() {
  // bool state = checkBoxPointsHasColors->isChecked();
  // checkBoxPointsHasColors->blockSignals(true);
  // checkBoxPointsHasColors->setChecked(!state);
  // checkBoxPointsHasColors->blockSignals(false);
  updateState();
}


// void GuiPanelImplicit::on_comboSampleSurface_currentIndexChanged(int index) {
// }

void GuiPanelImplicit::on_pushButtonSampleSurface_clicked() {
  QString sampleElement = comboSampleSurface->currentText();

  // TODO Thu Mar 23 12:28:50 2023
  // - move to a core/Sampling class

  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl==(SceneGraph*)0) return;
  Shape* shapeSurface = dynamic_cast<Shape*>(pWrl->find("SURFACE"));
  if(shapeSurface==nullptr) return;
  Shape* shapePoints  = dynamic_cast<Shape*>(pWrl->find("POINTS"));
  if(shapePoints==nullptr) {
    Group* parent = dynamic_cast<Group*>(shapeSurface->getParent());
    if(parent==nullptr) return;
    shapePoints = new Shape();
    shapePoints->setName("POINTS");
    parent->addChild(shapePoints);
    Appearance* a = new Appearance();
    shapePoints->setAppearance(a);
    Material* m = new Material();
    a->setMaterial(m);
    QRgb qVC = data.getVertexColor();
    Color d(static_cast<uint>(qVC));
    m->setDiffuseColor(d);
    shapePoints->setGeometry(new IndexedFaceSet());
  }

  IndexedFaceSet* ifsP =
    dynamic_cast<IndexedFaceSet*>(shapePoints->getGeometry());
  ifsP->clear();
  ifsP->setNormalPerVertex(true);
  vector<float>& coordP  = ifsP->getCoord();
  vector<float>& normalP = ifsP->getNormal();

  IndexedFaceSet* ifsS =
    dynamic_cast<IndexedFaceSet*>(shapeSurface->getGeometry());
  if(sampleElement=="VERTICES") {
    vector<float>& coordS = ifsS->getCoord();
    coordP.insert(coordP.end(),coordS.begin(),coordS.end());
    if(ifsS->hasNormalPerVertex()) {
      vector<float>& normalS = ifsS->getNormal();
      normalP.insert(normalP.end(),normalS.begin(),normalS.end());
    } else {
      vector<int>& coordIndexS = ifsS->getCoordIndex();
      Geometry::computeNormalsPerVertex(coordS,coordIndexS,normalP);
    }
  } else /* if(sampleElement=="FACES") */ {
    vector<float>& coordS      = ifsS->getCoord();
    vector<int>&   coordIndexS = ifsS->getCoordIndex();
    Geometry::computeFaceCentroids(coordS,coordIndexS,coordP);
    if(ifsS->hasNormalPerFace()) {
      vector<float>& normalS      = ifsS->getNormal();
      vector<int>&   normalIndexS = ifsS->getNormalIndex();
      if(normalIndexS.size()==0) {
        normalP.insert(normalP.end(),normalS.begin(),normalS.end());        
      } else {
        int nF = ifsS->getNumberOfFaces();
        for(int iF=0;iF<nF;iF++) {
          int iN = normalIndexS[iF];
          normalP.push_back(normalS[3*iN  ]);
          normalP.push_back(normalS[3*iN+1]);
          normalP.push_back(normalS[3*iN+2]);
        }
      }
    } else {
      Geometry::computeNormalsPerFace(coordS,coordIndexS,normalP);
    }
  }

  shapeSurface->setShow(false);
  shapePoints->setShow(true);

  mainWindow->setSceneGraph(pWrl,false);
  mainWindow->refresh();
  updateState();

}

void GuiPanelImplicit::on_pushButtonSamplePoints_clicked() {
  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl!=(SceneGraph*)0) {
    // SceneGraphProcessor processor(*pWrl);

    Shape* shapePoints  = dynamic_cast<Shape*>(pWrl->find("POINTS"));
    if(shapePoints==nullptr) return;
    IndexedFaceSet* ifsPoints =
      dynamic_cast<IndexedFaceSet*>(shapePoints->getGeometry());
    if(ifsPoints==nullptr) return;
    if(ifsPoints->hasNormalPerVertex()==false) return;
    vector<float>& coordSrc  = ifsPoints->getCoord();
    vector<float>& normalSrc = ifsPoints->getNormal();

    Shape* shapeSample  = dynamic_cast<Shape*>(pWrl->find("SAMPLES"));
    if(shapeSample==nullptr) {
      Group* parent = dynamic_cast<Group*>(shapePoints->getParent());
      if(parent==nullptr) return;
      shapeSample = new Shape();
      shapeSample->setName("SAMPLES");
      parent->addChild(shapeSample);
      Appearance* a = new Appearance();
      shapeSample->setAppearance(a);
      Material* m = new Material();
      a->setMaterial(m);
      QRgb qVC = data.getVertexColor();
      Color d(static_cast<uint>(qVC));
      m->setDiffuseColor(d);
      shapeSample->setGeometry(new IndexedFaceSet());
    }
    IndexedFaceSet* ifsSample =
      dynamic_cast<IndexedFaceSet*>(shapeSample->getGeometry());
    ifsSample->clear();
    ifsSample->setNormalPerVertex(true);
    vector<float>& coordSample  = ifsSample->getCoord();
    vector<float>& normalSample = ifsSample->getNormal();

    Vec3f& center = pWrl->getBBoxCenter();
    Vec3f& size   = pWrl->getBBoxSize();
    float  scale  = data.getBBoxScale();
    bool   cube   = data.getBBoxCube();
    int    depth  = data.getBBoxDepth();

    if(_hexGridPartition!=(HexGridPartition*)0) {
      delete _hexGridPartition;
      _hexGridPartition = nullptr;
    }

    _hexGridPartition = new HexGridPartition(center,size,1<<depth,scale,cube);
    _hexGridPartition->insertPoints(coordSrc,normalSrc);
    _hexGridPartition->sample(coordSample,normalSample);

    shapePoints->setShow(false);
    shapeSample->setShow(true);

    mainWindow->setSceneGraph(pWrl,false);
    mainWindow->refresh();
    updateState();
  }
}

void GuiPanelImplicit::on_pushButtonSamplesSwap_clicked() {
  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl==(SceneGraph*)0) return;
  Node* nodeSample = pWrl->find("SAMPLES");
  Node* nodePoints = pWrl->find("POINTS");
  if(nodeSample==(Node*)0 || nodePoints==(Node*)0) return;
  bool showPoints = nodeSample->getShow();
  bool showSample = !showPoints;
  nodeSample->setShow(showSample);
  nodePoints->setShow(showPoints);
  mainWindow->setSceneGraph(pWrl,false);
  mainWindow->refresh();
  updateState();
}

void GuiPanelImplicit::on_pushButtonSamplesRemove_clicked() {
  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl!=(SceneGraph*)0) {
    SceneGraphProcessor processor(*pWrl);
    processor.samplesRemove();
    mainWindow->setSceneGraph(pWrl,false);
    mainWindow->refresh();
    updateState();
  }
}

void GuiPanelImplicit::on_pushButtonSamplesShow_clicked() {
  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl==(SceneGraph*)0) return;
  Node* node = pWrl->find("SAMPLES");
  if(node==(Node*)0) return;
  node->setShow(true);
  mainWindow->setSceneGraph(pWrl,false);
  mainWindow->refresh();
  updateState();
}

void GuiPanelImplicit::on_pushButtonSamplesHide_clicked() {
  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl==(SceneGraph*)0) return;
  Node* node = pWrl->find("SAMPLES");
  if(node==(Node*)0) return;
  node->setShow(false);
  mainWindow->setSceneGraph(pWrl,false);
  mainWindow->refresh();
  updateState();
}

void GuiPanelImplicit::on_pushButtonSurfaceRemove_clicked() {
  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl!=(SceneGraph*)0) {
    SceneGraphProcessor processor(*pWrl);
    processor.surfaceRemove();
    mainWindow->setSceneGraph(pWrl,false);
    mainWindow->refresh();
    updateState();
  }
}

void GuiPanelImplicit::on_pushButtonSurfaceShow_clicked() {
  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl==(SceneGraph*)0) return;
  Node* node = pWrl->find("SURFACE");
  if(node==(Node*)0) return;
  node->setShow(true);
  mainWindow->setSceneGraph(pWrl,false);
  mainWindow->refresh();
  updateState();
}

void GuiPanelImplicit::on_pushButtonSurfaceHide_clicked() {
  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl==(SceneGraph*)0) return;
  Node* node = pWrl->find("SURFACE");
  if(node==(Node*)0) return;
  node->setShow(false);
  mainWindow->setSceneGraph(pWrl,false);
  mainWindow->refresh();
  updateState();
}

// Plane
void GuiPanelImplicit::on_pushButtonPlaneRemove_clicked() {
  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl!=(SceneGraph*)0) {
    SceneGraphProcessor processor(*pWrl);
    processor.planeRemove();
    mainWindow->setSceneGraph(pWrl,false);
    mainWindow->refresh();
    updateState();
  }
}

void GuiPanelImplicit::on_pushButtonPlaneShow_clicked() {
  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl==(SceneGraph*)0) return;
  Node* node = pWrl->find("PLANE");
  if(node==(Node*)0) return;
  node->setShow(true);
  mainWindow->setSceneGraph(pWrl,false);
  mainWindow->refresh();
  updateState();
}

void GuiPanelImplicit::on_pushButtonPlaneHide_clicked() {
  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl==(SceneGraph*)0) return;
  Node* node = pWrl->find("PLANE");
  if(node==(Node*)0) return;
  node->setShow(false);
  mainWindow->setSceneGraph(pWrl,false);
  mainWindow->refresh();
  updateState();
}  

void GuiPanelImplicit::on_pushButtonPlaneFit_clicked() {
  auto mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl!=(SceneGraph*)0) {
    SceneGraphProcessor processor(*pWrl);
    Vec3f& center = pWrl->getBBoxCenter();
    Vec3f& size   = pWrl->getBBoxSize();
    float  scale  = data.getBBoxScale();
    bool   cube   = data.getBBoxCube();
    // int    depth = data.getBBoxDepth();
    Vec4f  f;
    processor.fitSinglePlane(center,size,scale,cube,f);
    data.setPlane(f);
    mainWindow->setSceneGraph(pWrl,false);
    mainWindow->refresh();
    updateState();
  }
}

void GuiPanelImplicit::on_pushButtonSurfRecMultiplePlanes_clicked() {
  GuiMainWindow* mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl!=(SceneGraph*)0) {

    SceneGraphProcessor processor(*pWrl);
    Vec3f& center = pWrl->getBBoxCenter();
    Vec3f& size   = pWrl->getBBoxSize();
    float  scale  = data.getBBoxScale();
    bool   cube   = data.getBBoxCube();
    int    depth = data.getBBoxDepth();
    vector<float> fVec;
    processor.fitMultiplePlanes(center,size,depth,scale,cube,fVec);

    mainWindow->setSceneGraph(pWrl,false);
    mainWindow->refresh();
    updateState();
  }
}

void GuiPanelImplicit::on_pushButtonSurfRecContinuous_clicked() {
  GuiMainWindow* mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl!=(SceneGraph*)0) {

    SceneGraphProcessor processor(*pWrl);
    Vec3f& center = pWrl->getBBoxCenter();
    Vec3f& size   = pWrl->getBBoxSize();
    float  scale  = data.getBBoxScale();
    bool   cube   = data.getBBoxCube();
    int    depth = data.getBBoxDepth();
    vector<float>& fGrid = data.getFunctionVertices();
    fGrid.clear();
    processor.fitContinuous(center,size,depth,scale,cube,fGrid);

    mainWindow->setSceneGraph(pWrl,false);
    mainWindow->refresh();
    updateState();
  }
}

void GuiPanelImplicit::on_pushButtonSurfRecWatertight_clicked() {
  GuiMainWindow* mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl!=(SceneGraph*)0) {

    SceneGraphProcessor processor(*pWrl);
    Vec3f& center = pWrl->getBBoxCenter();
    Vec3f& size   = pWrl->getBBoxSize();
    float  scale  = data.getBBoxScale();
    bool   cube   = data.getBBoxCube();
    int    depth = data.getBBoxDepth();
    vector<float>& fGrid = data.getFunctionVertices();
    fGrid.clear();
    processor.fitWatertight(center,size,depth,scale,cube,fGrid);

    mainWindow->setSceneGraph(pWrl,false);
    mainWindow->refresh();
    updateState();
  }
}

void GuiPanelImplicit::on_pushButtonSurfRecJacobi_clicked() {
  GuiMainWindow* mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl!=(SceneGraph*)0) {

    SceneGraphProcessor processor(*pWrl);
    Vec3f& center = pWrl->getBBoxCenter();
    Vec3f& size   = pWrl->getBBoxSize();
    float  scale  = data.getBBoxScale();
    bool   cube   = data.getBBoxCube();
    int    depth = data.getBBoxDepth();
    vector<float>& fGrid = data.getFunctionVertices();
    processor.fitOptimalJacobi(center,size,depth,scale,cube,fGrid);

    mainWindow->setSceneGraph(pWrl,false);
    mainWindow->refresh();
    updateState();
  }
  std::cout << "}\n";
}

// Non Convex Hull functions
void GuiPanelImplicit::on_comboNchSource_currentIndexChanged(int index) {
    (void)index;
  updateState();
}

void GuiPanelImplicit::on_pushButtonNchFit_clicked() {
  QString source = comboNchSource->currentText();

  GuiMainWindow* mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  QProgressBar*  progress = mainWindow->getProgressBar();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl!=(SceneGraph*)0) {
    NchProcessor::setProgressBar(progress);
    NchProcessor processor(*pWrl);
    processor.setSource(source);

    bool multiThreaded = false;
    processor.nchEstimateRhos(multiThreaded);

    mainWindow->setSceneGraph(pWrl,false);
    mainWindow->refresh();
    updateState();
  }
}

void GuiPanelImplicit::on_pushButtonNchFitMultiThreaded_clicked() {
  QString source = comboNchSource->currentText();

  GuiMainWindow* mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  QProgressBar*  progress = mainWindow->getProgressBar();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl!=(SceneGraph*)0) {
    NchProcessor::setProgressBar(progress);
    NchProcessor processor(*pWrl);
    processor.setSource(source);

    bool multiThreaded = true;
    processor.nchEstimateRhos(multiThreaded);

    mainWindow->setSceneGraph(pWrl,false);
    mainWindow->refresh();
    updateState();
  }
}

void GuiPanelImplicit::on_pushButtonNchReconstruct_clicked() {
  QString source = comboNchSource->currentText();

  GuiMainWindow* mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  QProgressBar*  progress = mainWindow->getProgressBar();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl!=(SceneGraph*)0) {

    NchProcessor::setProgressBar(progress);
    NchProcessor processor(*pWrl);
    processor.setSource(source);

    Vec3f& center = pWrl->getBBoxCenter();
    Vec3f& size   = pWrl->getBBoxSize();
    float  scale  = data.getBBoxScale();
    bool   cube   = data.getBBoxCube();
    int    depth  = data.getBBoxDepth();

    vector<float>& fGrid = data.getFunctionVertices();
    fGrid.clear();

    bool multiThreaded = false;
    if(checkBoxNchReconstructOccupiedCells->isChecked() && _hexGridPartition!=nullptr) {
      processor.nchReconstruct(multiThreaded,*_hexGridPartition,fGrid);
    } else {
      processor.nchReconstruct(multiThreaded,center,size,depth,scale,cube,fGrid);
    }

    mainWindow->setSceneGraph(pWrl,false);
    mainWindow->refresh();
    updateState();
  }
}

void GuiPanelImplicit::on_pushButtonNchReconstructMultiThreaded_clicked() {
  QString source = comboNchSource->currentText();

  GuiMainWindow* mainWindow = getApp()->getMainWindow();
  GuiViewerData& data = mainWindow->getData();
  QProgressBar*  progress = mainWindow->getProgressBar();
  SceneGraph*    pWrl = data.getSceneGraph();
  if(pWrl!=(SceneGraph*)0) {

    NchProcessor::setProgressBar(progress);
    NchProcessor processor(*pWrl);
    processor.setSource(source);
    
    Vec3f& center = pWrl->getBBoxCenter();
    Vec3f& size   = pWrl->getBBoxSize();
    float  scale  = data.getBBoxScale();
    bool   cube   = data.getBBoxCube();
    int    depth  = data.getBBoxDepth();
    
    vector<float>& fGrid = data.getFunctionVertices();
    fGrid.clear();

    bool multiThreaded = true;
    if(checkBoxNchReconstructOccupiedCells->isChecked() && _hexGridPartition!=nullptr) {
      processor.nchReconstruct(multiThreaded,*_hexGridPartition,fGrid);
    } else {
      processor.nchReconstruct(multiThreaded,center,size,depth,scale,cube,fGrid);
    }

    mainWindow->setSceneGraph(pWrl,false);
    mainWindow->refresh();
    updateState();
  }

}

