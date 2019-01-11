# PolycrystalTriangulation
# Example:
```c++
CrysesShell crysesShell;
...
Polycrystal3 polycr("crysesShell.txt");
// Or
Polycrystal3 polycr(crysesShell);

double preferredTetrahedronEdgeLength = 0.2;
PolycrMesh* mesh = polycr.generateMesh(preferredTetrahedronEdgeLength);
...
delete mesh;
```
![Cube triangulation](https://github.com/Tokarevart/PolycrystalTriangulation/blob/master/images/shell_triang_3d_7.png)
