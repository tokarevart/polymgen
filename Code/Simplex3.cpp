#include "Simplex3.h"


double Simplex3::computeVolume() const
{
	const Vec3& v0 = (*vertexes[0])->getPosition();
	return 0.16666666666666666 * abs(Vec3::mixedProduct(
		(*vertexes[1])->getPosition() - v0, 
		(*vertexes[2])->getPosition() - v0, 
		(*vertexes[3])->getPosition() - v0));
}

double Simplex3::computeQuality() const
{
	double prods[4];
	for (int i = 0; i < 4; i++)
	{
		prods[i] = 1.0;
		for (int j = 0; j < 4; j++)
			if (j != i)
				prods[i] *= (**vertexes[j] - **vertexes[i]).magnitude();
	}
	double max_prod = std::max({ prods[0], prods[1], prods[2], prods[3] });

	return 8.48528137423857 * computeVolume() / max_prod;
}

Simplex3::Simplex3() : unique_ptr_helper<Simplex3>(this) {}

Simplex3::Simplex3(Vertex3 &vert0, Vertex3 &vert1, Vertex3 &vert2, Vertex3 &vert3) : unique_ptr_helper<Simplex3>(this)
{
	vertexes[0] = vert0.getPtrToUPtr();
	vertexes[1] = vert1.getPtrToUPtr();
	vertexes[2] = vert2.getPtrToUPtr();
	vertexes[3] = vert3.getPtrToUPtr();
}

Simplex3::~Simplex3() {}