# PolycrystalTriangulation
# Example:
```c++
CrysesShell cryses_shell;
...
Polycrystal3 polycr("cryses_shell.txt");
// Or
Polycrystal3 polycr(cryses_shell);

double preferred_tetrahedron_edge_length = 0.2;
PolycrMesh* mesh = polycr.TriangulatePolycrystal(preferred_tetrahedron_edge_length);
...
delete mesh;
```
![Cube triangulation](https://github.com/Tokarevart/PolycrystalTriangulation/blob/master/images/shell_triang_3d_7.png)
