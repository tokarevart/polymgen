// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "polygen/polygen.h"




psg::PolyStruct psg::generateCuboids(size_t nX, size_t nY, size_t nZ, real_t dX, real_t dY, real_t dZ) noexcept
{
    const size_t w = nX;
    const size_t t = nY;
    const size_t h = nZ;

    PolyStruct shell;
    shell.verts.resize((w + 1) * (t + 1) * (h + 1));
    shell.faces.resize(2 * (w * h * (t + 1)
                        + t * h * (w + 1)
                        + w * t * (h + 1)));
    shell.polyhedrons.resize(w * t * h);

    for (size_t i = 0; i < shell.polyhedrons.size(); i++)
        shell.polyhedrons[i].assign(12, 0);

    for (size_t i = 0; i < shell.verts.size(); i++)
    {
        shell.verts[i][0] = dX * (i % (w + 1));
        shell.verts[i][1] = dY * ((i % ((w + 1) * (t + 1))) / (w + 1));
        shell.verts[i][2] = dZ * (i / ((w + 1) * (t + 1)));
    }

    for (size_t i = 0; i < 2 * t * h * (w + 1); i += 2)
    {
        shell.faces[i][0] = i / 2 + (w + 1) * ((i / 2) / ((w + 1) * t));
        shell.faces[i][1] = shell.faces[i][0] + (w + 1) * (t + 1);
        shell.faces[i][2] = shell.faces[i][1] + w + 1;

        shell.faces[i + 1][0] = shell.faces[i    ][0];
        shell.faces[i + 1][1] = shell.faces[i + 1][0] + w + 1;
        shell.faces[i + 1][2] = shell.faces[i + 1][1] + (w + 1) * (t + 1);
    }
    for (size_t j = 0, j0 = 2 * t * h * (w + 1); j < 2 * w * h * (t + 1); j += 2)
    {
        shell.faces[j0 + j][0] = j / 2 + (j / 2) / w;
        shell.faces[j0 + j][1] = shell.faces[j0 + j][0] + 1;
        shell.faces[j0 + j][2] = shell.faces[j0 + j][1] + (w + 1) * (t + 1);

        shell.faces[j0 + j + 1][0] = shell.faces[j0 + j    ][0];
        shell.faces[j0 + j + 1][1] = shell.faces[j0 + j + 1][0] + (w + 1) * (t + 1);
        shell.faces[j0 + j + 1][2] = shell.faces[j0 + j + 1][1] + 1;
    }
    for (size_t k = 0, k0 = 2 * h * (t * (w + 1) + w * (t + 1)); k < 2 * w * t * (h + 1); k += 2)
    {
        shell.faces[k0 + k][0] = k / 2 + (k / 2) / w + (w + 1) * ((k / 2) / (w * t));
        shell.faces[k0 + k][1] = shell.faces[k0 + k][0] + 1;
        shell.faces[k0 + k][2] = shell.faces[k0 + k][1] + w;

        shell.faces[k0 + k + 1][0] = shell.faces[k0 + k    ][1];
        shell.faces[k0 + k + 1][1] = shell.faces[k0 + k + 1][0]     + w;
        shell.faces[k0 + k + 1][2] = shell.faces[k0 + k + 1][1] + 1;
    }

    for (size_t i = 0; i < shell.polyhedrons.size(); i++)
    {
        shell.polyhedrons[i][0] = 2 * (i + i / w);
        shell.polyhedrons[i][1] = shell.polyhedrons[i][0] + 1;

        shell.polyhedrons[i][2] = shell.polyhedrons[i][1] + 1;
        shell.polyhedrons[i][3] = shell.polyhedrons[i][2] + 1;

        shell.polyhedrons[i][4] = 2 * (h * t * (w + 1) + i + w * (i / (w * t)));
        shell.polyhedrons[i][5] = shell.polyhedrons[i][4] + 1;

        shell.polyhedrons[i][6] = shell.polyhedrons[i][4] + 2 * w;
        shell.polyhedrons[i][7] = shell.polyhedrons[i][6] + 1;

        shell.polyhedrons[i][8] = 2 * (h * (w * (t + 1) + t * (w + 1)) + i);
        shell.polyhedrons[i][9] = shell.polyhedrons[i][8] + 1;

        shell.polyhedrons[i][10] = shell.polyhedrons[i][8]  + 2 * w * t;
        shell.polyhedrons[i][11] = shell.polyhedrons[i][10] + 1;
    }

    return shell;
}
