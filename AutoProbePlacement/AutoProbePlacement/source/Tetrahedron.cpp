#include "Tetrahedron.h"

Tetrahedron::Tetrahedron(){

}

Tetrahedron::Tetrahedron(int _id, Array<int> _indexes){
	id = _id;
	indexes = _indexes;

}

String Tetrahedron::toString(){
	String s = "";
	char str[100];
    sprintf(str, " Tetrahedron # %i: %i %i %i %i ",id, indexes[0], indexes[1], indexes[2], indexes[3]);
    s.append(str);
    s.append(" Neighbors : ");
    for (int i = 0; i < neighbors.size(); i++){
    	sprintf(str,  " %d ",neighbors[i]);
    	s.append(str);
    }
    return s;
}