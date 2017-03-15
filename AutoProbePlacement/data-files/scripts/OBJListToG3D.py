import sys;
import string

inputFile = open(sys.argv[1], 'r');

outputFile = open("out.txt", 'w');

arr = inputFile.readlines();


i = 0;
for line in arr:

	outputFile.write('model' + repr(i) + ' = ArticulatedModel::Specification { filename = "objs/mesh' + "%03d"%i + '.obj" };\n' );
	# outputFile.write('modelEntity' + repr(i) + ' = VisibleEntity { \nmodel = "model' + repr(i) + '";\n visible = true; \n}; \n' );

	i += 1;