import PIL
import sys;
from os import listdir
from os.path import isfile, join
import re
import os;
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

path = "C:/temp/screensbackups/" + sys.argv[1] + "/";
allthefiles = [f for f in listdir(path) if (isfile(join(path, f))) and (f[-4:] == ".jpg")];
allthefiles.sort(key=natural_keys);
for i in xrange(len(allthefiles)):
	stringnumber = allthefiles[i][:-4];
	print(allthefiles[i] + "\n");
	
	number = atoi(stringnumber);
	
	if not (number == i):
		# print(path + allthefiles[i]);
		os.rename(path + allthefiles[i], path + repr(i) + ".jpg");