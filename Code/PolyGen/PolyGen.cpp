#include "PolyGen.h"


PolyStruct* polygen::generateCuboidsPolycrystal(size_t nX, size_t nY, size_t nZ, double dX, double dY, double dZ)
{
    const size_t w = nX;
    const size_t t = nY;
    const size_t h = nZ;

    PolyStruct* shell = new PolyStruct;
    shell->nNodes  = (w + 1ull) * (t + 1ull) * (h + 1ull);
    shell->nFacets = 2ull * (w * h * (t + 1ull)
                           + t * h * (w + 1ull)
                           + w * t * (h + 1ull));
    shell->nCryses = w * t * h;

    shell->nCrysesFacets = new size_t[shell->nCryses];
    for (size_t i = 0ull; i < shell->nCryses; i++)
        shell->nCrysesFacets[i] = 12ull;
    
    shell->nodesPositions = new double[3ull * shell->nNodes];
    for (size_t i = 0ull; i < shell->nNodes; i++)
    {
        shell->nodesPositions[3ull * i]        = dX * (i % (w + 1ull));
        shell->nodesPositions[3ull * i + 1ull] = dY * ((i % ((w + 1ull) * (t + 1ull))) / (w + 1ull));
        shell->nodesPositions[3ull * i + 2ull] = dZ * (i / ((w + 1ull) * (t + 1ull)));
    }

    shell->facets = new size_t[3ull * shell->nFacets];
    for (size_t i = 0ull; i < 2ull * t * h * (w + 1ull); i += 2ull)
    {
        shell->facets[3ull * i]        = i / 2ull + (w + 1ull) * ((i / 2ull) / ((w + 1ull) * t));
        shell->facets[3ull * i + 1ull] = shell->facets[3ull * i] + (w + 1ull) * (t + 1ull);
        shell->facets[3ull * i + 2ull] = shell->facets[3ull * i + 1ull] + w + 1ull;

        shell->facets[3ull * (i + 1ull)]        = shell->facets[3ull * i];
        shell->facets[3ull * (i + 1ull) + 1ull] = shell->facets[3ull * (i + 1ull)]        + w + 1ull;
        shell->facets[3ull * (i + 1ull) + 2ull] = shell->facets[3ull * (i + 1ull) + 1ull] + (w + 1ull) * (t + 1ull);
    }
    for (size_t j = 0ull, j0 = 2ull * t * h * (w + 1ull); j < 2ull * w * h * (t + 1ull); j += 2ull)
    {
        shell->facets[3ull * (j0 + j)]        = j / 2ull + (j / 2ull) / w;
        shell->facets[3ull * (j0 + j) + 1ull] = shell->facets[3ull * (j0 + j)] + 1ull;
        shell->facets[3ull * (j0 + j) + 2ull] = shell->facets[3ull * (j0 + j) + 1ull] + (w + 1ull) * (t + 1ull);

        shell->facets[3ull * (j0 + j + 1ull)]        = shell->facets[3ull * (j0 + j)];
        shell->facets[3ull * (j0 + j + 1ull) + 1ull] = shell->facets[3ull * (j0 + j + 1ull)]        + (w + 1ull) * (t + 1ull);
        shell->facets[3ull * (j0 + j + 1ull) + 2ull] = shell->facets[3ull * (j0 + j + 1ull) + 1ull] + 1ull;
    }
    for (size_t k = 0ull, k0 = 2ull * h * (t * (w + 1ull) + w * (t + 1ull)); k < 2ull * w * t * (h + 1ull); k += 2ull)
    {
        shell->facets[3ull * (k0 + k)]        = k / 2ull + (k / 2ull) / w + (w + 1ull) * ((k / 2ull) / (w * t));
        shell->facets[3ull * (k0 + k) + 1ull] = shell->facets[3ull * (k0 + k)] + 1ull;
        shell->facets[3ull * (k0 + k) + 2ull] = shell->facets[3ull * (k0 + k) + 1ull] + w;

        shell->facets[3ull * (k0 + k + 1ull)]        = shell->facets[3ull * (k0 + k) + 1ull];
        shell->facets[3ull * (k0 + k + 1ull) + 1ull] = shell->facets[3ull * (k0 + k + 1ull)]        + w;
        shell->facets[3ull * (k0 + k + 1ull) + 2ull] = shell->facets[3ull * (k0 + k + 1ull) + 1ull] + 1ull;
    }

    shell->cryses = new size_t[12ull * shell->nCryses];
    for (size_t i = 0ull; i < shell->nCryses; i++)
    {
        shell->cryses[12ull * i]        = 2ull * (i + i / w);
        shell->cryses[12ull * i + 1ull] = shell->cryses[12ull * i] + 1ull;

        shell->cryses[12ull * i + 2ull] = shell->cryses[12ull * i + 1ull] + 1ull;
        shell->cryses[12ull * i + 3ull] = shell->cryses[12ull * i + 2ull] + 1ull;

        shell->cryses[12ull * i + 4ull] = 2ull * (h * t * (w + 1ull) + i + w * (i / (w * t)));
        shell->cryses[12ull * i + 5ull] = shell->cryses[12ull * i + 4ull] + 1ull;

        shell->cryses[12ull * i + 6ull] = shell->cryses[12ull * i + 4ull] + 2ull * w;
        shell->cryses[12ull * i + 7ull] = shell->cryses[12ull * i + 6ull] + 1ull;

        shell->cryses[12ull * i + 8ull] = 2ull * (h * (w * (t + 1ull) + t * (w + 1ull)) + i);
        shell->cryses[12ull * i + 9ull] = shell->cryses[12ull * i + 8ull] + 1ull;

        shell->cryses[12ull * i + 10ull] = shell->cryses[12ull * i + 8ull]  + 2ull * w * t;
        shell->cryses[12ull * i + 11ull] = shell->cryses[12ull * i + 10ull] + 1ull;
    }

    return shell;
}
