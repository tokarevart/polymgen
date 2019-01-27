# PolycrystalTriangulation
# Example:
```c++
PolyStruct polyStruct;
...
Polycrystal3 polycr("polyStruct.txt");
// Or
Polycrystal3 polycr(polyStruct);

double preferredTetrahedronEdgeLength = 0.2;
PolycrMesh* mesh = polycr.generateMesh(preferredTetrahedronEdgeLength);
...
delete mesh;
```
![Cube triangulation](https://github.com/Tokarevart/PolycrystalTriangulation/blob/master/images/shell_triang_3d_7.png)
