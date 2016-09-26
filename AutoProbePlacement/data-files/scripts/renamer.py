import os;
# path = "/Users/joelpp/Documents/Maitrise/Images/Mitsuba/RedGreen/renamed/SphericalCamera_BetterScaleObj_";
# path2 = "/Users/joelpp/Documents/Maitrise/Code/IrradianceMapping/coefficients/RedGreen/renamed/QuadratureCoefficientsSloan_";

# num = [9,4,9];
# intIndex = 0;
# probeInfo = open("/Users/joelpp/Documents/Maitrise/code/combined/redgreenprobes.txt","r+");
# png = ".png";
# txt = ".txt";
# p = ".p";
# for i in range(num[0]):
# 	for j in range(num[1]):
# 		for k in range(num[2]):
# 			index = repr(i)+"_"+repr(j)+"_"+repr(k);
# 			# print(path+index+png);
# 			# os.rename(path+index+png, path+repr(intIndex)+png);
# 			# os.rename(path+index+txt, path2+repr(intIndex)+txt);
# 			# os.rename(path+index+p, path2+repr(intIndex)+p);
# 			probeInfo.write(repr(i-4) + " " +repr(j+1) + " " +repr(k-4)+"\n");
# 			intIndex = intIndex + 1;

path = "/Users/joelpp/Documents/Maitrise/Code/experiments/0/"

for i in xrange(80):
	# os.rename(path+"coefficients/SphericalCamera_BetterScaleObj_"+repr(i)+".p", path+"coefficients/Coefficients_"+repr(i)+".p")
	# os.rename(path+"coefficients/SphericalCamera_BetterScaleObj_"+repr(i)+".txt", path+"coefficients/Coefficients_"+repr(i)+".txt")
	os.rename(path+"probes/BetterScene_"+repr(i)+".png", path+"probes/Probe_"+repr(i)+".png")