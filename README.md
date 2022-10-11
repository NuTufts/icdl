# icdl

Code dedicated to the reconstruction of particle interactions from images produced by Liquid Argon Time Projection Chambers.

Algorithms target the data coming from detectors in the Short Baseline Experiment: SBND, MicroBooNE, and ICARUS.

The reco algorithms are built on top of libraries that define the data format along with interfaces to the data.

Key components and external dependencies:

* ROOT: a library for storing large datasets, commonly used in High Energy Physics [external dependency]
* larlite: defines the data format for data extracted from SBN experiments and the Larsoft framework. Also provides interfaces to key geometry and detector properties. [dependencies ROOT]
* larcv: defines image and voxel formats along with metadata used to make connections to the detector geometry [dependencies: ROOT, larlite]
* laropencv: defines algorithms for LArTPCs that make use of the OpenCV library [dep: larlite, OpenCV (external)]
* ublarcvapp: utilities commonly needed by algorithms and analyses [dep: larlite, larcv, laropencv]
* lardly: routines for visualization of the data within a web-app or embedded into a Jupyter notebook [dep: larlite, larcv, lardly (external), python3 (external)]
* larflow: reconstruction algorithms built around the creation of 3D spacepoints directly from 2D images [dep: larlite, larcv, pytorch (external)]

Machine learning algorithms are also in use by this library. These include the following:

* LArMatch: 3D spacepoint reconstruction from 2D wireplane images along with keypoint detection
* SSNet: 2D pixel-labeling
* Sparse-MaskRCNN: (to be added)

## Building

Go to the [Wiki](https://github.com/NuTufts/icdl/wiki) to find info on building and other tutorials.

### On Tufts Machines

### On the Tufts Cluster

