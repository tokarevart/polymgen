# Polymgen
Polymgen is an Apache 2.0 licensed library for generating tetrahedral mesh of a polycrystal for subsequent use in solving boundary-value problems by the finite element method (FEM).
# Example:
```c++
// Number of cubes along the corresponding axis.
size_t nX = 3, nY = 3, nZ = 3;
polygen::PolyStruct polyStruct = generateCuboidsPolycrystal(nX, nY, nZ);
pmg::Polycrystal polycr(polyStruct);

double preferredTetrahedronEdgeLength = 0.45;
pmg::PolycrMesh* mesh = polycr.generateMesh(preferredTetrahedronEdgeLength);
...
delete mesh;
```
![Cube triangulation](https://github.com/Tokarevart/polycr-mesh-generator/blob/master/images/polymesh_1.png)
# License
Copyright © 2018-2019 Tokarev Artem. All rights reserved.

Licensed under the [MIT License](./LICENSE).
