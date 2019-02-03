#include "PolyGen.h"


std::unique_ptr<PolyStruct> polygen::generateCuboidsPolycrystal(size_t nX, size_t nY, size_t nZ, double dX, double dY, double dZ)
{
    const size_t w = nX;
    const size_t t = nY;
    const size_t h = nZ;

    std::unique_ptr<PolyStruct> shell = std::make_unique<PolyStruct>();
    shell->nNodes  = (w + 1) * (t + 1) * (h + 1);
    shell->nFacets = 2 * (w * h * (t + 1)
                        + t * h * (w + 1)
                        + w * t * (h + 1));
    shell->nCryses = w * t * h;

    shell->nCrysesFacets = new size_t[shell->nCryses];
    for (size_t i = 0; i < shell->nCryses; i++)
        shell->nCrysesFacets[i] = 12;
    
    shell->nodesPositions = new double[3 * shell->nNodes];
    for (size_t i = 0; i < shell->nNodes; i++)
    {
        shell->nodesPositions[3 * i]     = dX * (i % (w + 1));
        shell->nodesPositions[3 * i + 1] = dY * ((i % ((w + 1) * (t + 1))) / (w + 1));
        shell->nodesPositions[3 * i + 2] = dZ * (i / ((w + 1) * (t + 1)));
    }

    shell->facets = new size_t[3 * shell->nFacets];
    for (size_t i = 0; i < 2 * t * h * (w + 1); i += 2)
    {
        shell->facets[3 * i]     = i / 2 + (w + 1) * ((i / 2) / ((w + 1) * t));
        shell->facets[3 * i + 1] = shell->facets[3 * i] + (w + 1) * (t + 1);
        shell->facets[3 * i + 2] = shell->facets[3 * i + 1] + w + 1;

        shell->facets[3 * (i + 1)]     = shell->facets[3 * i];
        shell->facets[3 * (i + 1) + 1] = shell->facets[3 * (i + 1)]     + w + 1;
        shell->facets[3 * (i + 1) + 2] = shell->facets[3 * (i + 1) + 1] + (w + 1) * (t + 1);
    }
    for (size_t j = 0, j0 = 2 * t * h * (w + 1); j < 2 * w * h * (t + 1); j += 2)
    {
        shell->facets[3 * (j0 + j)]     = j / 2 + (j / 2) / w;
        shell->facets[3 * (j0 + j) + 1] = shell->facets[3 * (j0 + j)] + 1;
        shell->facets[3 * (j0 + j) + 2] = shell->facets[3 * (j0 + j) + 1] + (w + 1) * (t + 1);

        shell->facets[3 * (j0 + j + 1)]     = shell->facets[3 * (j0 + j)];
        shell->facets[3 * (j0 + j + 1) + 1] = shell->facets[3 * (j0 + j + 1)]     + (w + 1) * (t + 1);
        shell->facets[3 * (j0 + j + 1) + 2] = shell->facets[3 * (j0 + j + 1) + 1] + 1;
    }
    for (size_t k = 0, k0 = 2 * h * (t * (w + 1) + w * (t + 1)); k < 2 * w * t * (h + 1); k += 2)
    {
        shell->facets[3 * (k0 + k)]     = k / 2 + (k / 2) / w + (w + 1) * ((k / 2) / (w * t));
        shell->facets[3 * (k0 + k) + 1] = shell->facets[3 * (k0 + k)] + 1;
        shell->facets[3 * (k0 + k) + 2] = shell->facets[3 * (k0 + k) + 1] + w;

        shell->facets[3 * (k0 + k + 1)]     = shell->facets[3 * (k0 + k) + 1];
        shell->facets[3 * (k0 + k + 1) + 1] = shell->facets[3 * (k0 + k + 1)]     + w;
        shell->facets[3 * (k0 + k + 1) + 2] = shell->facets[3 * (k0 + k + 1) + 1] + 1;
    }

    shell->cryses = new size_t[12 * shell->nCryses];
    for (size_t i = 0; i < shell->nCryses; i++)
    {
        shell->cryses[12 * i]     = 2 * (i + i / w);
        shell->cryses[12 * i + 1] = shell->cryses[12 * i] + 1;

        shell->cryses[12 * i + 2] = shell->cryses[12 * i + 1] + 1;
        shell->cryses[12 * i + 3] = shell->cryses[12 * i + 2] + 1;

        shell->cryses[12 * i + 4] = 2 * (h * t * (w + 1) + i + w * (i / (w * t)));
        shell->cryses[12 * i + 5] = shell->cryses[12 * i + 4] + 1;

        shell->cryses[12 * i + 6] = shell->cryses[12 * i + 4] + 2 * w;
        shell->cryses[12 * i + 7] = shell->cryses[12 * i + 6] + 1;

        shell->cryses[12 * i + 8] = 2 * (h * (w * (t + 1) + t * (w + 1)) + i);
        shell->cryses[12 * i + 9] = shell->cryses[12 * i + 8] + 1;

        shell->cryses[12 * i + 10] = shell->cryses[12 * i + 8]  + 2 * w * t;
        shell->cryses[12 * i + 11] = shell->cryses[12 * i + 10] + 1;
    }

    return shell;
}
