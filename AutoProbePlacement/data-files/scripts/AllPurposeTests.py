 # -*- coding: utf-8 -*-
from sympy.mpmath import *
from scipy import misc

from helper import *
import numpy
import math
import cPickle
import sys
from Tkinter import *
import ImageTk

class App:

    def __init__(self, master):

        frame = Frame(master)
        frame.pack()

        self.button = Button(
            frame, text="QUIT", fg="red", command=frame.quit
            )
        self.button.pack(side=LEFT)

        self.hi_there = Button(frame, text="Hello", command=self.say_hi)
        self.hi_there.pack(side=LEFT);

        self.myImage = ImageTk.PhotoImage(file = '../experiments/57/misc/Probe_1.png');



        self.w = Canvas(frame, width=256, height=128,bd=1,relief='sunken');
        self.w.create_image(128,64,image = self.myImage);
        self.w.bind("<Button-1>", self.onCanvasClicked);
        self.w.pack();

        self.v = StringVar();
        self.v.set("Click to choose coordinates");

        self.coordLabel = Label(master, textvariable=self.v);
        self.coordLabel.pack();
    def onCanvasClicked(self, event):
    	canvas = event.widget
    	j = canvas.canvasx(event.x)
    	i = canvas.canvasy(event.y)
    	(u,v) = ijToUV_yup(i,j, 128, 256);
    	(phi, theta) = equiUVtoPT_yup(u,v);
    	(x,y,z) = sphericalToCartesian_yup(1,phi, theta);



    	self.v.set("Clicked at (i,j)="+repr((i,j))+"; (u,v)="+repr((u,v))+";\n (phi,theta)=("+repr((phi,theta))+"); (x,y,z)="+repr((x,y,z)));

    def say_hi(self):
        print "hi there, everyone!"

root = Tk()

app = App(root)
root.tkraise();
root.geometry("900x200")
root.mainloop()
root.destroy() # optional; see description below

height = 128;
width = 256;



maxBand = 9;
numCoeffs = 0;
for i in range(maxBand+1):
	numCoeffs = numCoeffs + (2*i + 1);

# shArray = numpy.zeros((numCoeffs, height,width));
# iArray = numpy.zeros((height,width));
# jArray = numpy.zeros((height,width));
# uArray = numpy.zeros((height,width));
# vArray = numpy.zeros((height,width));
# xArray = numpy.zeros((height,width));
# yArray = numpy.zeros((height,width));
# zArray = numpy.zeros((height,width));
allPurpose = numpy.zeros((height,width,3));
# longitudeArray = numpy.zeros((height,width));
# latitudeArray = numpy.zeros((height,width));
epsilon = 0.01;
total = height * width * numCoeffs;

perten = math.floor(total / 10);
percent = math.floor(total / 100);
progress = 0;

# allPurpose[height/4][width/2][0] = 255;
# allPurpose[height/4][width/2][1] = 255;
# allPurpose[height/4][width/2][2] = 255;

# allPurpose[height/3][width/2][0] = 255;
# allPurpose[height/3][width/2][1] = 0;
# allPurpose[height/3][width/2][2] = 0;

# allPurpose[2*height/3][width/3][0] = 0;
# allPurpose[2*height/3][width/3][1] = 255;
# allPurpose[2*height/3][width/3][2] = 0;

# allPurpose[2*height/3][2*width/3][0] = 0;
# allPurpose[2*height/3][2*width/3][1] = 0;
# allPurpose[2*height/3][2*width/3][2] = 255;
epsilon = 0.02
for i in range(height):
	print(repr(i)+"/"+repr(height)+ " complete")
	for j in range(width):

		(u,v) = ijToUV_yup(i,j, height, width);
		if ( abs((u*u + v*v) - 0.5) < epsilon ):
			(p,t) = uvToIJ_yup(u,v,width,height);
			if (u < 0):
				allPurpose[p][t][0]= 255.0;
				continue;
			if (u > 0):
				allPurpose[p][t][1]= 255.0;
			else:
				allPurpose[p][t][0]= 255.0;
				allPurpose[p][t][1]= 255.0;
				allPurpose[p][t][2]= 255.0;
		# (phi, theta) = equiUVtoPT_yup(u,v);
		# (x,y,z) = sphericalToCartesian_yup(1,phi, theta);


		pass;

		# if (i == 128):
		# 	print("i="+repr(i)+" j="+repr(j)+" u="+repr(u)+" v="+repr(v)+" phi="+repr(phi)+" theta="+repr(theta)+" x="+repr(x)+" y="+repr(y)+" z="+repr(z));

		# iArray[i][j] = float(i)/float(height);
		# jArray[i][j] = float(j)/float(width);
		# uArray[i][j] = (u+1)/2; #as u is between -1 and 1
		# vArray[i][j] = (v+1.0)/2.0; #as u is between -1 and 1
		# xArray[i][j] = x; #as u is between -1 and 1
		# yArray[i][j] = y; #as u is between -1 and 1
		# zArray[i][j] = z; #as u is between -1 and 1
		
		# if (phi == 0):
		# 	longitudeArray[i][j] = 1.0;
		# if (abs(phi-(pi/4)) < epsilon):
		# 	longitudeArray[i][j] = 1.0;
		# if (abs(phi-(pi/2)) < epsilon):
		# 	longitudeArray[i][j] = 1.0;
		# if (abs(phi-(3*pi/2)) < epsilon):
		# 	longitudeArray[i][j] = 1.0;
		# if (abs(phi+(pi/4)) < epsilon):
		# 	longitudeArray[i][j] = 1.0;
		# if (abs(phi+(pi/2)) < epsilon):
		# 	longitudeArray[i][j] = 1.0;
		# if (abs(phi+(3*pi/2)) < epsilon):
		# 	longitudeArray[i][j] = 1.0;
		#  #as u is between -1 and 1
		# if (abs(theta-(pi/4)) < epsilon):
		# 	latitudeArray[i][j] = 1.0;
		# if (abs(theta-(3*pi/4)) < epsilon):
		# 	latitudeArray[i][j] = 1.0;	
		# if (theta < epsilon):
		# 	latitudeArray[i][j] = 1.0;	
		# for k in range(numCoeffs):
		# 	shArray[k][i][j] = SHxyz_yup(kToLM(k), x,y,z);

# for i in xrange(numCoeffs):
# 	cPickle.dump(shArray[i], open('/Users/joelpp/Documents/Maitrise/Code/IrradianceMapping/SHvalues_yup/128x64/values_'+repr(i)+'.p', 'wb')) 

misc.imsave("/Users/joelpp/Documents/Maitrise/Code/Validation/allPurpose.png",allPurpose);
# misc.imsave("/Users/joelpp/Documents/Maitrise/Code/Validation/iArray.png",iArray);
# misc.imsave("/Users/joelpp/Documents/Maitrise/Code/Validation/jArray.png",jArray);
# misc.imsave("/Users/joelpp/Documents/Maitrise/Code/Validation/uArray.png",uArray);
# misc.imsave("/Users/joelpp/Documents/Maitrise/Code/Validation/vArray.png",vArray);
# misc.imsave("/Users/joelpp/Documents/Maitrise/Code/Validation/xArray.png",xArray);
# misc.imsave("/Users/joelpp/Documents/Maitrise/Code/Validation/yArray.png",yArray);
# misc.imsave("/Users/joelpp/Documents/Maitrise/Code/Validation/zArray.png",zArray);
# misc.imsave("/Users/joelpp/Documents/Maitrise/Code/Validation/latitudeArray.png",latitudeArray);
# misc.imsave("/Users/joelpp/Documents/Maitrise/Code/Validation/longitudeArray.png",longitudeArray);
