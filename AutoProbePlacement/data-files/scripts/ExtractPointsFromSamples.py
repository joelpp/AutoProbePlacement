
SceneName = "crytek-sponza";
SampleSetName = "95000_Samples_0";
FilePath = "../Scenes/" + SceneName + "/SampleSets/" + SampleSetName + "/SamplePositions.txt"

PositionLine = 1
NumberOfLinesPerSample = 9;

File = open(FilePath);
i = 0;
points = [];
for line in File:

	if (i % NumberOfLinesPerSample == PositionLine):
		points.append(line);
		
	i = i + 1;
	
print "Number of read samples:"+repr(i / NumberOfLinesPerSample);

Output = open("../Scenes/" + SceneName + "/SampleSets/" + SampleSetName + "/PointsOnly.txt", 'w')

for index in xrange(len(points)):
	Output.write(points[index])