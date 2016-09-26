

file = open("C:/temp/temp.txt", 'w');

for i in xrange(100):
	file.write(repr((float(i)+1) * 0.001) + "\n");
	
file.close();