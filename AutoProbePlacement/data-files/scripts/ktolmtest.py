from helper import *

for i in xrange(30):
	print("i = "+repr(i));
	print("new = "+repr(kToLM_new(i)));
	print("old = "+repr(kToLM(i)));
	print();