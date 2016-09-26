
#ifndef TETRAHEDRON_H
#define TETRAHEDRON_H

#include <G3D/G3DAll.h>
class Tetrahedron{

	public:
		String toString();
		Tetrahedron(int, Array<int>);
		Tetrahedron();

		Array<int> indexes;
		Array<int> neighbors;
		int id;
};

#endif