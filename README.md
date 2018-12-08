# PolycrystalTriangulation
# Example:
```c++
CrystallitesShell crystallites_shell;
...
Polycrystal3 polycr("crystallites_shell.txt");
// Or
Polycrystal3 polycr(crystallites_shell);

double preferred_tetrahedron_edge_length = 0.2;
PolycrystalTriangulation* triangulation = polycr.TriangulatePolycrystal(preferred_tetrahedron_edge_length);
...
delete triangulation;
```
![Cube triangulation](https://github.com/Tokarevart/PolycrystalTriangulation/blob/master/images/shell_triang_3d_7.png)
