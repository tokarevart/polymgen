#include "polygen/polygen.h"

using polygen::PolyStruct;




PolyStruct polygen::generateCuboidsPolycrystal(size_t nX, size_t nY, size_t nZ, double dX, double dY, double dZ) noexcept
{
    const size_t w = nX;
    const size_t t = nY;
    const size_t h = nZ;

    PolyStruct shell;
    shell.nodes.resize((w + 1) * (t + 1) * (h + 1));
    shell.facets.resize(2 * (w * h * (t + 1)
                        + t * h * (w + 1)
                        + w * t * (h + 1)));
    shell.cryses.resize(w * t * h);

    for (size_t i = 0; i < shell.cryses.size(); i++)
        shell.cryses[i].assign(12, 0);

    for (size_t i = 0; i < shell.nodes.size(); i++)
    {
        shell.nodes[i][0] = dX * (i % (w + 1));
        shell.nodes[i][1] = dY * ((i % ((w + 1) * (t + 1))) / (w + 1));
        shell.nodes[i][2] = dZ * (i / ((w + 1) * (t + 1)));
    }

    for (size_t i = 0; i < 2 * t * h * (w + 1); i += 2)
    {
        shell.facets[i][0] = i / 2 + (w + 1) * ((i / 2) / ((w + 1) * t));
        shell.facets[i][1] = shell.facets[i][0] + (w + 1) * (t + 1);
        shell.facets[i][2] = shell.facets[i][1] + w + 1;

        shell.facets[i + 1][0] = shell.facets[i    ][0];
        shell.facets[i + 1][1] = shell.facets[i + 1][0] + w + 1;
        shell.facets[i + 1][2] = shell.facets[i + 1][1] + (w + 1) * (t + 1);
    }
    for (size_t j = 0, j0 = 2 * t * h * (w + 1); j < 2 * w * h * (t + 1); j += 2)
    {
        shell.facets[j0 + j][0]     = j / 2 + (j / 2) / w;
        shell.facets[j0 + j][1] = shell.facets[j0 + j][0] + 1;
        shell.facets[j0 + j][2] = shell.facets[j0 + j][1] + (w + 1) * (t + 1);

        shell.facets[j0 + j + 1][0] = shell.facets[j0 + j    ][0];
        shell.facets[j0 + j + 1][1] = shell.facets[j0 + j + 1][0] + (w + 1) * (t + 1);
        shell.facets[j0 + j + 1][2] = shell.facets[j0 + j + 1][1] + 1;
    }
    for (size_t k = 0, k0 = 2 * h * (t * (w + 1) + w * (t + 1)); k < 2 * w * t * (h + 1); k += 2)
    {
        shell.facets[k0 + k][0] = k / 2 + (k / 2) / w + (w + 1) * ((k / 2) / (w * t));
        shell.facets[k0 + k][1] = shell.facets[k0 + k][0] + 1;
        shell.facets[k0 + k][2] = shell.facets[k0 + k][1] + w;

        shell.facets[k0 + k + 1][0] = shell.facets[k0 + k    ][1];
        shell.facets[k0 + k + 1][1] = shell.facets[k0 + k + 1][0]     + w;
        shell.facets[k0 + k + 1][2] = shell.facets[k0 + k + 1][1] + 1;
    }

    for (size_t i = 0; i < shell.cryses.size(); i++)
    {
        shell.cryses[i][0] = 2 * (i + i / w);
        shell.cryses[i][1] = shell.cryses[i][0] + 1;

        shell.cryses[i][2] = shell.cryses[i][1] + 1;
        shell.cryses[i][3] = shell.cryses[i][2] + 1;

        shell.cryses[i][4] = 2 * (h * t * (w + 1) + i + w * (i / (w * t)));
        shell.cryses[i][5] = shell.cryses[i][4] + 1;

        shell.cryses[i][6] = shell.cryses[i][4] + 2 * w;
        shell.cryses[i][7] = shell.cryses[i][6] + 1;

        shell.cryses[i][8] = 2 * (h * (w * (t + 1) + t * (w + 1)) + i);
        shell.cryses[i][9] = shell.cryses[i][8] + 1;

        shell.cryses[i][10] = shell.cryses[i][8]  + 2 * w * t;
        shell.cryses[i][11] = shell.cryses[i][10] + 1;
    }

    return shell;
}
