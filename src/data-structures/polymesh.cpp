// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "data-structures/polymesh.h"

using pmg::PolyMesh;




PolyMesh::~PolyMesh()
{
    delete[] nodesPositions;
    delete[] tetrs;
    delete[] nCrysesTetrs;
}
