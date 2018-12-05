# QMVC #

## Overview ##

Reference implementation of the research paper:

[Mean value coordinates for quad cages in 3D]](https://www.telecom-paristech.fr/~boubek/papers/QMVC)
Jean-Marc Thiery, Pooran Memari and Tamy Boubekeur 
ACM Transactions on Graphics - Proc. SIGGRAPH Asia 2018 

This repository contains:
* QMVC/coordinates: a set methods to compute popular space coordinates, including:
	* Green Coordinates (GC)
	* Mean Value Coordinates (MVC)
	* Spherical Mean Value Coordinates (SMVC)
	* maximum Entropy coordinates (MEC) 
	* our Quad Mean Value Coordinates (QMVC)
* QMVC/viewer: an interactive 3D viewer allowing to:
	* compute the projection of a mesh within the space coordinates of a cage
	* deform the the mesh using the cage


Copyright(C) 2018
Jean-Marc Thiery, Pooran Memari and Tamy Boubekeur

## Release Notes ##

### v1.0 ###

Dry push of the paper's reference implementation. Future releases will be cleaner and easier to deploy. 

## Building the viewer ##

Viewer's dependencies:
- Qt5.5.1 or more recent
- libQGLViewer-2.6.1

Adapt the content of the file viewer/CageManip.pro to your setting (path to libGLViewer)

To be completed...

-------------------------------------------------------------------------------------------

At the beginning of CageManipInterface.h
you will see the following lines:
#define ALLOW_TRI_MVC
#define ALLOW_QMVC
#define ALLOW_QMVC_MEC // DO NOT COMPUTE QMVC_MEC IF YOU DO NOT COMPUTE QMVC !!!!
#define ALLOW_SMVC
#define ALLOW_GC

simply comment the corresponding lines if you do not want to compute the coordinates.

-------------------------------------------------------------------------------------------

See and example of computation of QMVC in viewer/CageManipInterface.h, in the function: computeQMVCCoordinates()


## Authors

* [**Jean-Marc Thiery**](https://www.telecom-paristech.fr/~thiery/) 
* [**Tamy Boubekeur**](https://www.telecom-paristech.fr/~boubek)

See also the list of [contributors](https://github.com/superboubek/QMVC/contributors) who participated in this project.

## Citation

Please cite the following paper in case you are using this code:
>**Mean value coordinates for quad cages in 3D** *Jean-Marc theiry, Pooran Memari and Tamy Boubekeur.* ACM Transactions on Graphics (Proc. SIGGRAPh Asia 2018), art. 29, 2018

## License

This project is licensed under a GNU-GPL license - see the [LICENSE.txt](LICENSE.txt) file for details.