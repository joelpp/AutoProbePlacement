
root = "C:/temp/"
filePaths = [root + "samplesRGB_ref.txt", root + "samplesRGB.txt"];
files = [open(x, 'r') for x in filePaths];
errorLogFile = open(root + "errorlog.txt",'a');

lines = [x.readlines() for x in files]
error = 0;
for i in xrange(len(lines[0])):
	error += (float(lines[0][i]) - float(lines[1][i])) ** 2;
	
print(error);

errorLogFile.write(repr(error) + '\n');