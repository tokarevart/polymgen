# polycr-mesh-generator
# Example:
```c++
// Number of cubes along the corresponding axis.
size_t nX = 2, nY = 5, nZ = 9;
polygen::PolyStruct polyStruct = generateCuboidsPolycrystal(nX, nY, nZ);
pmg::Polycrystal polycr(polyStruct);

double preferredTetrahedronEdgeLength = 0.2;
pmg::PolycrMesh* mesh = polycr.generateMesh(preferredTetrahedronEdgeLength);
...
delete mesh;
```
![Cube triangulation](https://github.com/Tokarevart/polycr-mesh-generator/blob/master/images/shell_triang_3d_7.png)
