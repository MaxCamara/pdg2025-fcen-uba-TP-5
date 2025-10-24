** DGP2025 Digital Geometry Processing **

Assignment 5

The files are organized as in Assignment 1, and the compilation
process is the same.

Before compiling the application for the first time, edit the file
src/gui/GuiStrings.hpp or copy it from previous assignments

You have to implement several methods in new classes

The command line test program directory tst has been removed.

Additional instruction are included as comments in the source files

A number of files for the classes developed in previous assignments
are missing or empty in this repository. You should copy your version
of those files from assignment 4 in the core directory, because
otherwise the application will not compile.

A number of changes have been made to other classes, and new classes
have been added to the repository. Do not merge any other older files,
because the application may not compile.

This is the list of tasks in this assignment. It may look like a lot,
but the methods are closely related, and a lot of hints are
provided. I suggest that you make the minimal necessary changes to the
files so that the application compiles and runs, and then go work on
the tasks one at a time. Commit your chages as you make progress, at a
minimum after completing each task.

// class core/HexGridPartition
//
// This class is used to partition a set of 3D points into buckets
// corresponding to the cells of a regular hexahedral grid.
// The coordinates of the points are quantized to the desired
// resolution with respect to a bounding box, and inserted into
// the corresponding buckets. The buckets are represented as lists
// using an internal vector<int> next array. The first element of
// the bucket lists are accessed throup a map<int,int>. where the
// first index of the bucket (ix,iy,iz) is the scan order index
// ix+N*(iy+N*iz), where N is the number of cells each axis of
// the bounding box is subdivided into. The second index is the
// point index of the first point in the bucket list. Subsequent
// points are found using the next array. The last point in each
// bucket list has a -1 as its next value.
//
// We will use it to implement a simple version of the vertex clustering
// simplification method, which was pending since we worked on loading
// STL files. The Optimization panel in the GUI has a new section named
// "CLUSTER VERTICES", and the class core/Optimization has some new
// methods to support this function, which you also have to
// implement. Details below.
//
// We will also use the core/HexGridPartition to subsample a point cloud
// by selecting one sample point from each occupied cell. Since choosing
// the original point closest to the mean of the points contained in each
// cell has been shown to produce good results in practice, you will have
// to implement this sample() method of the class. There is a new
// Implicit panel in the GUI where this sampling method is used to
// generate a new node named SAMPLE in the scene graph. The Implicit
// panel allows you to chose to run NCH Surface reconstruction method
// either on the full point cloud, or on the subsampled point cloud to
// compare speeds. Again, more details below.
//
// Finally, we will also use the occupied cells of the
// core/HexGridPartition class to implicitly generate sparse regular
// grids in an attempt to save time in the computation of isosurfaces for
// surface reconstruction. More details below.

// class core/Optimization
//
// The following public methods are added to the class
//
//  void clusterVerticesApply
//  (const Vec3f& min, const Vec3f& max, const int resolution);
//  void clusterVerticesApply();
//  int  getQuantizationResolution();
//  void setQuantizationResolution
//  (const int resolution);
//  void setQuantizationBox
//  (const Vec3f& min, const Vec3f& max);
//  int getNumberOfClusters();
//
// The following private variables should also be added to the class
//
//  Vec3f           _clusterMin;
//  Vec3f           _clusterMax;
//  int             _clusterResolution;
//  int             _nClusters;
//
// More details about what needs to be done in this class in the
// source files

// class core/IsoSurf
//
// This class is fully implemented due to lack of time.
// However, I encourage you to take a deep look at the source code.
// This is my implementation of marching cubes, which produces
// polygon meshes with face ranging from 3 to 7 vertices.
// 
// The method
//   static void IsoSurf::computeIsosurface
//   (const Vec3f& center, const Vec3f& size,
//    const int depth, const float scale, const bool cube,
//    vector<float>& fGrid, IndexedFaceSet& isoSurface);
//
// is used to generate an isosurface from scalar data provided for all
// the vertices of a regular hexahedral grid, represented in the fGrid
// array.  The center,size,scale, and cube variables describe the
// bounding box of the regular grid, which is then split int N=1<<depth
// cells along each of the three axes. This regular grid has
// nGridCells=N*N*N, and nGridVertices=(N+1)*(N+1)*(N+1). The fGrid array
// should have nGridVertices values arranged in scan order
// ix+N*(iy+N*iz), where 0<=ix,iy,iz<=N

// The method
//   static void IsoSurf::computeIsosurface
//   (HexGridPartition& hgp,
//    vector<float>& fGrid, IndexedFaceSet& isoSurface);
//
// is used to generate an isosurface from scalar data provided for only
// the vertices of the sparse regular hexahedral grid corresponding to
// the occupied cells of the hgp variable. The scalar values
// corresponding to the vertices of these cells is represented in the
// fGrid array. The diffulty here is to figure exacly which vertices are
// involved here, and to be able to traverse them in scan order.  The
// information required to generate the underlaying regular grid is
// stored in the hgp variable.

// The Implicit panel has a last section named "SURFACE RECONSTRUCTION"
// The following methods have been added to the class
// wrl/SceneGraphProcessor, which can be run within this section
//
// void fitSinglePlane
// (const Vec3f& center, const Vec3f& size, const float scale, const bool cube,
//  Vec4f& f /* temporary output */ );
//
// The method fitSinglePlane() is fully implemented. You have to
// implement the following TWO methods. Detailed instructions are
// included as comments in the wrl/SceneGraphProcessor.cpp file, and
// will be discussed in class
//
// void fitMultiplePlanes
// (const Vec3f& center, const Vec3f& size,
//  const int depth, const float scale, const bool cube,
//  vector<float>& fVec );
//
// void fitContinuous
// (const Vec3f& center, const Vec3f& size,
//  const int depth, const float scale, const bool cube);

// The implementation of the following two methods is optional.
// If you implement them you wwill get extra credit
//
// void fitWatertight
// (const Vec3f& center, const Vec3f& size,
//  const int depth, const float scale, const bool cube,
//  vector<float>& fGrid);
//
// void fitOptimalJacobi
// (const Vec3f& center, const Vec3f& size,
//  const int depth, const float scale, const bool cube,
//  vector<float>& fGrid /* input & output */);
  
// The Implicit panel has a last section named "NON-CONVEX HULL", and
// there is a new src/nch directory containing the following classes
// 
// nch/NchProcessor
// nch/NchThread
// nch/NchEstimate
// nch/NchEvaluate
//
// These classes implement the non-convex hull surface reconstruction
// method that we discussed in class. The estimation of the rho
// parameters, which define the nch approximate signed distance function,
// is separate from the evaluation of the funcion on the vertices of the
// grid, and subsequent extraction of isosurfaces. The reason for this
// decision is that, since the estimation of the rho parameters is only
// function of the points, and independent of the grid, you can
// experiment with different grid sizes, and also with dense vs adaptive
// grids, to develop the intuition towards additional steps.
//
// Since the naive implementation of the non-convex hull method, using
// all the points in the input data set, and evaluating the nch function
// on all the vertices of a regular grid is time consuming, there are
// heuristic options to experiment with different ways to speed up the
// computation.
//
// On one hand there is the option of estimating the rho parameters from
// a downsampled point cloud. This option results in a nch function which
// requires much less time to estimate its rho parameters, and to
// evaluate it, but it is different from the nch function obtained using
// all the points. In particular, the nch function obtained from the
// subsample interpolates the points contained in the subsample, but it
// does not necessarily interpolates any of the other original points,
// and it is not clear before hand how much error is introduces in those
// points. An subsequent option here, which is not implemented, is to set
// an error threshold, and add to the subsample all the original points
// for which the nch function obtained from the subsample result in a
// value exceeding the threshold.
// 
// Since the non-convex hull methods is in principle massively
// paralelizable, both the estimation of the rho parameters and the
// evaluation of the nch function, another option that is implemented is
// to multi-threding. Here we the implementation will be based on the
// QThread class. If you haven't worked with multi-threading before, here
// you have a way to learn something. All the classes in the nch library
// are partially implemented. More details about what needs to be done in
// the source files

Submission:

  As usual,
  a) locate and remove any temporary files created by your editor
  b) delete the contents of the top subdirectories bin and build
  c) delete the file qt/DGP2023-A4.pro.user created by QtCreator
  d) delete the top build directory created by QtCreator
  e) commit all your changes at regular intervals and at milestones
  f) push your commits to your GitHub repository
  g) follow the assignment submission instructions in GitHub Classroom
     and submit the assignment by the deadline
