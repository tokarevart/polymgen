# PolyMGen
PolyMGen is a library for parallel generating tetrahedral mesh of many polyhedrons. 
It is created with the purpose of use in solving boundary-value problems for polycrystal by the finite element method.
## Example:
```c++
// Number of cubes along the corresponding axis.
size_t nX = 2, nY = 4, nZ = 8;
psg::PolyShell shell = psg::generateCuboids(nX, nY, nZ);
pmg::PolyhedralSet polyhset(shell);

// Preferred tetrahedron edge length.
double preferredLength = 0.45;
polyhset.tetrahedralize(preferredLength);
polyhset.output(pmg::FileType::WavefrontObj); // Or pmg::FileType::LsDynaKeyword
```
![Cube mesh](https://github.com/Tokarevart/polymgen/blob/master/images/polymesh_3.png)
![Cube mesh](https://github.com/Tokarevart/polymgen/blob/master/images/polymesh_2.png)
![Cube mesh](https://github.com/Tokarevart/polymgen/blob/master/images/polymesh_1.png)
## License
Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.

Licensed under the [MIT License](/LICENSE).
