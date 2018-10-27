#include "Simplex3.h"


double Simplex3::Volume() const
{
	const Vector3& v0 = (*vertexes[0])->GetPosition();
	return 0.16666666666666666 * abs(Vector3::MixedProduct(
		(*vertexes[1])->GetPosition() - v0, 
		(*vertexes[2])->GetPosition() - v0, 
		(*vertexes[3])->GetPosition() - v0));
}

double Simplex3::Quality() const
{
	double prods[4];
	for (int i = 0; i < 4; i++)
	{
		prods[i] = 1.0;
		for (int j = 0; j < 4; j++)
			if (j != i)
				prods[i] *= (**vertexes[j] - **vertexes[i]).Magnitude();
	}
	double max_prod = std::max({ prods[0], prods[1], prods[2], prods[3] });

	return 8.48528137423857 * Volume() / max_prod;
}

Simplex3::Simplex3() : unique_ptr_helper<Simplex3>(this) {}

Simplex3::Simplex3(Vertex3 &vert0, Vertex3 &vert1, Vertex3 &vert2, Vertex3 &vert3) : unique_ptr_helper<Simplex3>(this)
{
	vertexes[0] = vert0.GetPtrToUniquePtr();
	vertexes[1] = vert1.GetPtrToUniquePtr();
	vertexes[2] = vert2.GetPtrToUniquePtr();
	vertexes[3] = vert3.GetPtrToUniquePtr();
}

Simplex3::~Simplex3() {}