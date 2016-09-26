import sys
import os
import LabelScreenshots
import shutil

save = (len(sys.argv) > 1) and (sys.argv[1] == "save");
src = "C:/temp/CurrentOptimization/screens/";
dest = "C:/temp/screensbackups/"

if ((not save) and (raw_input('Do you really want to reset without saving? ENTER will reset without saving, anything else will abort.') != "")):
	print("Operation cancelled.");
	exit(1);

if (save):
	
	LabelScreenshots.label();
	
	if (len(sys.argv) == 3):
		next = sys.argv[2];
	else:
		next = repr(len(os.walk(dest).next()[1]));
	
	directory = dest + next;
	if not os.path.exists(directory):
		os.makedirs(directory);
	print("Saving to " + directory);	
		
	src_files = os.listdir(src)
	for file_name in src_files:
		full_file_name = os.path.join(src, file_name)
		if (os.path.isfile(full_file_name)):
			shutil.copy(full_file_name, directory);
	try:
		shutil.copy("C:/temp/CurrentOptimization/dplog.txt", directory);
	except IOError:
		pass
		
	try:
		shutil.copy("C:/temp/CurrentOptimization/InitialConditions.txt", directory);
	except IOError:
		pass

	try:
		shutil.copy("C:/temp/CurrentOptimization/errorlog.txt", directory);
	except IOError:
		pass


try:
	open("C:/temp/CurrentOptimization/errorlog.txt", 'w').close();
except IOError:
	pass

[ os.remove(src + f) for f in os.listdir("C:/temp/CurrentOptimization/screens/") ]

try:
	os.remove("C:/temp/CurrentOptimization/dplog.txt");
except WindowsError:
	pass
	
try:
	os.remove("C:/temp/CurrentOptimization/InitialConditions.txt");
except WindowsError:
	pass

print("Reset successful.");