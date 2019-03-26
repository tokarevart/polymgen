# Polymgen
Polymgen is a library for generating tetrahedral mesh of many polyhedrons. 
It is created with the purpose of use in solving boundary-value problems for polycrystal by the finite element method.
# Example:
```c++
// Number of cubes along the corresponding axis.
size_t nX = 2, nY = 4, nZ = 8;
psg::PolyShell shell = psg::generateCuboids(nX, nY, nZ);
pmg::PolyhedralSet polyhedr(shell);

// Preferred tetrahedron edge length.
double preferredLength = 0.45;
polyhedr.generateMesh(preferredLength);
polyhedr.output(pmg::PolyhedralSet::FileType::WavefrontObj); // Or ...FileType::LsDynaKeyword
```
![Cube mesh](https://github.com/Tokarevart/polymgen/blob/master/images/polymesh_3.png)
![Cube mesh](https://github.com/Tokarevart/polymgen/blob/master/images/polymesh_2.png)
![Cube mesh](https://github.com/Tokarevart/polymgen/blob/master/images/polymesh_1.png)
# License
Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.

Licensed under the [MIT License](/LICENSE).
