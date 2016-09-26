

file = open("C:/temp/errorlog.txt", 'r');
file2 = open("C:/temp/errorlog2.txt", 'r');

lines = file.readlines();
lines2 = file2.readlines();

outfile = open("C:/temp/diff.txt", 'w');

for i in xrange(len(lines)):
	outfile.write(repr((float(lines[i]) - float(lines2[i])) ** 2) + "\n");

file.close();
file2.close();
outfile.close();