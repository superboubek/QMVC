# QMVC Data #

This archive contains example data for the paper:

[Mean value coordinates for quad cages in 3D]](https://www.telecom-paristech.fr/~boubek/papers/QMVC)
Jean-Marc Thiery, Pooran Memari and Tamy Boubekeur 
ACM Transactions on Graphics - Proc. SIGGRAPH Asia 2018 

## Data description

For each 3D model <MODEL>, we provide the following files:
* <MODEL>.obj: the input high resolution mesh
* <MODEL>\_Cage.obj: the input cage in rest pose
* <MODEL>\_Cage\_Deformed.obj: the deformed cage, used to generate the paper's figures.

To reproduce the figures of the paper, proceed as follow:
* download, build and run [QMVCViewer](https://www.telecom-paristech.fr/~boubek/papers/QMVC)
* click on the "Open Mesh" button and select <MODEL>.obj
* click on the "Open Binding Cage" button and select <MODEL>\_Cage.obj (note: this can take few seconds as several coordinates may be computed, depending on the compiling options of QMVCViewer)
* click on the "Open Deformed Cage" button and select <MODEL>\_Cage\_Deformed.obj

Once loaded, the deformed cage is instantly used as the current pose for the binding cage and the mesh is deformed. 

## Models ##

* Beast
  * Mesh: 28388 vertices, 56772 faces
  * Cage: 152 vertices, 300 faces
* Bench
  * Mesh: 65430 vertices, 650406 faces
  * Cage: 24 vertices, 22 faces
* Cactus
  * Mesh: 98820 vertices, 98304 faces
  * Cage: 36 vertices, 34 faces
* FireHydrant
  * Mesh: 39028 vertices, 76444 faces
  * Cage: 48 vertices, 46 faces
* ManHead
  * Mesh: 24658 vertices, 24610 faces
  * Cage: 12 vertices, 10 faces
* Ogre
  * Mesh: 26155 vertices, 52306 faces
  * Cage: 117 vertices, 258 faces
* SpikyBox
  * Mesh: 60802 vertices, 60800 faces
  * Cage: 8 vertices, 6 faces
* WireSphere
  * Mesh: 48964 vertices, 49152 faces
  * Cage: 8 vertices, 6 faces

## Authors ##

* Jean-Marc Thiery
* Tamy Boubekeur

## Acknowledgements ##

We thank Kiaran Ritchie for the Ogre model. The beast model is courtesy of Autodesk. The ManHead model is courtesy Pixologic.  

## License ##

This data set is release under the [Creative Commons Attribution-NonCommercial 3.0 license (CC BY-NC 3.0)](https://creativecommons.org/licenses/by-nc/3.0/fr/deed.en)

Copyright(C) 2018
Jean-Marc Thiery and Tamy Boubekeur