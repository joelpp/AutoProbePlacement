

FilePath = "../temp/samplePositionsSponza.txt"

File = open(FilePath);
i = 0;
points = [];
for line in File:
	# if (i % 8 == 1):
		# print line;
	
	if (i % 9 == 7):
		points.append(line);
	i = i + 1;
print "Number of read samples:"+repr(i);
Output = open("../temp/colorsOnly.txt", 'r+')
for index in xrange(len(points)):
	Output.write(points[index]+"\n")