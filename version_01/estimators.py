"""
Copyright 2017, 2017 nowaythatsok

This file is part of GNU_lapser 0.1.

    GNU_lapser is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GNU_lapser is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GNU_lapser.  If not, see <http://www.gnu.org/licenses/>.
"""


#this files contains the IMAGE EXPOSURE estimators for "timelapser.py" and is a dependece of it
#
#I only included the best one...

#-------------------------------SWITCHES----------------------------------------

#none yet


#-------------------------------Official INCLUDES-------------------------------
import os				#to access files
from PIL import Image			#to open JPEGs
import numpy as np

#-------------------------------Custom   INCLUDES-------------------------------

import lookupTables as lT

#-------------------------------Function DEFINITIONS----------------------------


def fittingEstimator(inDir, file_name, EV_calc_local, EV_calc_global, LUT):

	EV_calc_local[:]	= [0.0]	#list of statistics
	EV_calc_global[:]	= []	

	
	#lattice parameter
	R=10

	#-------Opening image1
	#separate core of name, append jpg extension, open, BW
	name	= (file_name[0].split(".")[0] + ".jpg")
	img1 	= Image.open("./"+inDir+"/"+name)
	img1	= np.array( img1.convert("L") )

	#dimensions of image 
	N,M	= img1.shape


	#fill the first plot list
	cntr	= -1
	x1	= np.array( [0.0]*(len(range(0,N,R))*len(range(0,M,R))) )
	x2	= np.array( [0.0]*(len(range(0,N,R))*len(range(0,M,R))) )
	for i in range(0,N,R):
		for j in range(0,M,R):
			cntr	+= 1
			x1[cntr]	= LUT( img1[i,j] )
			

	for name in file_name[1:]:
	#-------Opening image 2
		#separate core of name, append jpg extension, open, BW
		name	= (name.split(".")[0] + ".jpg")
		img2 	= Image.open("./"+inDir+"/"+name)
		img2	= np.array( img2.convert('L') )
	
			

		cntr	= -1
		
		for i in range(0,N,R):
			for j in range(0,M,R):
			
				cntr	+= 1
				x2[cntr]	= LUT( img2[i,j] )
		
		#now compare x1  and x2
		#fitting x+b is equivalent to calculating 
		#print(np.average( x1-x2 ))
		
		xx1	= []
		xx2	= []

		for (x,y) in zip(x1,x2):
			if x>=1 and x<=7 and y>=1 and y<=7:
				xx1.append(x)
				xx2.append(y)
		if len(xx1) == 0 or len(xx2) == 0:
			print("WARNING!!") 
			print(x2) 
		EV_calc_local.append( np.average( np.array(xx2)-np.array(xx1) ) + EV_calc_local[-1] )

		#import matplotlib.pyplot as plt
		#plt.figure(1)
		#plt.plot(x1,x2,"x",xx1,xx2,"o")
		#plt.show()

		#change images
		img1[:]	= img2[:]
		x1[:]	= x2[:]
	EV_calc_global[:]	= EV_calc_local[:]















