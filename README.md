# Purpose
this is a Matlab implementation of Functional Faces: Groupwise Dense Correspondence using Functional Maps by Chao Zhang, William A.P. Smith, etc

# Dependency
the code is developed based on some existing codes, libraries, and data, they include:
1. pairwise func. map implementation developed by M. Ovsjanikov as to the paper Functional maps: a flexible representation of maps between shapes.
2. manopt: a Matlab toolbox for optimization on manifolds
3. the sample code uses three meshes from a high quality facial dataset:
reference: G. Stratou, A. Ghosh, P. Debevec, and L. Morency. Effect of illumination on automatic expression recognition: a novel 3D relightable facial database. In Proc. Face and Gesture, pages 611–618, 2011

# Usage
1. Download manopt if you were not using it and set up path in the code (see demo_ff.m)
2. Choose LB dimension and feature functions (see getMapOfTwoFaceShapes.m)
3. Run demo_ff.m

# Note
1. to compute the point-wise correspondence from functional maps, approximate nearest neighbor search might need to be compiled per to the operating system

# Reference
@inproceedings{zhang2016functional,
  title={Functional Faces: Groupwise Dense Correspondence using Functional Maps},
  author={Zhang, Chao and Smith, William AP and Dessein, Arnaud and Pears, Nick and Dai, Hang},
  booktitle={Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition},
  pages={5033--5041},
  year={2016}
}

