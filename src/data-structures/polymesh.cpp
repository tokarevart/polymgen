// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Contacts: <tokarev28.art@gmail.com>
// Licensed under the MIT License.

#include "data-structures/polymesh.h"

using pmg::PolyMesh;




PolyMesh::~PolyMesh()
{
    delete[] nodesPositions;
    delete[] tetrs;
    delete[] nCrysesTetrs;
}
