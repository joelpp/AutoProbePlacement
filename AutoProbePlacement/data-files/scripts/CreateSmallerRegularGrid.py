import os;
import sys;
import shutil;

def allodd(v):
	return (v[0] % 2 == 1) and (v[1] % 2 == 1) and (v[2] % 2 == 1);

def is_integer_vector(v):
	a = v[0].is_integer();
	b = v[1].is_integer();
	c = v[2].is_integer();
	# print ((v[0], v[0].is_integer()));
	# print ((v[1], v[1].is_integer()));
	# print ((v[2], v[2].is_integer()));

	return (a) and (b) and (c);
	
def LineToVectorF(line):
	vector = [float(x) for x in line.split()]
	# print(vector);
	return (vector);
	
newstep = 2;
OGlist = open("../ProbeStructures/RegularGrid_step05/probeList.txt",'r');
newList = open("../ProbeStructures/RegularGrid_step"+repr(newstep)+"/probeList.txt",'r+');

OGProbePath = "../ProbeStructures/RegularGrid_step05/coefficients/Coefficients_"
NewProbePath = "../ProbeStructures/RegularGrid_step"+repr(newstep)+"/coefficients/Coefficients_"

i = 0;
for probeNumber in xrange(4851):
	line = OGlist.readline();
	lineVF = LineToVectorF(line);
	if (is_integer_vector(lineVF)) and allodd(lineVF):
		newList.write(line);
		shutil.copy(OGProbePath+repr(probeNumber)+".txt", NewProbePath+repr(i)+".txt");
		
		i += 1;
		