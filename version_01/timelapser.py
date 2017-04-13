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


#this is the MAIN componenet of the GNU_lapser project
#use python 3 to run this script
#see the INSTRUCTIONS at ....




#the program will use the jpg 
#			IMAGES and xmp
#			SIDECAR FILES
#					located in the folder "files" in the script's directory (mod at SWITCHES)
#have fun timelapsing :)


#-------------------------------SWITCHES----------------------------------------

	#input folder
inDir	= "./files" #path

	#output folderfor the new XMPs
outDir	= "./outFiles" #path

	#visualization
v_read_data				= 1 #bool
v_set_differences			= 1 #bool
v_changes				= 1 #bool

v_plot_exp				= 1 #bool
v_plot_reference        		= 1 #bool

	#benchmarking outputs
b_print					= 1 #bool

    #threshold for local adjustments
EV_threshold            = .06 #EV

    #level of enviromental EXP preservation 
    #how many stops of liminance change is alowed in the course of the video
    #for a sunset time lapse this will for examle even everything out for 0 and darken by 1 stop by the end for value 1
    #in case EXP REFERENCE FRAMES are used this variable has no effect
exp_preserv             = 2.0  #EV
    
    #Are we making use of the reference frames? 
ref_on                  = 1 #bool
exp_ref                 = 0 #bool (not ready yet, besides I found out a nicer solution, so let it be 0)

    #the sigma of smoothening the transitions of the reference frames in SECONDS
smoothening             = 1 #positive float (you can make it 0)
    
    
#-------------------------------Official INCLUDES-------------------------------
import os				#to access files
import matplotlib.pyplot as plt		#to display modifications
import time				#for benchmarking
import re				#for using regular expressions
import math				#for EV calculation (log)
import numpy as np			#for fitting polynomials and evaluating them

#-------------------------------My INCLUDES-------------------------------------

import lookupTables as lT
import estimators as est
from xmp_io import *

#-------------------------------Function DEFINITIONS----------------------------



	#now sort results according to date taken
def timeOrdData(file_name,
			exposure_time,
			fstop_value,
			iso_value,
			date_taken,
			is_1_star):
	from operator import itemgetter
	

	nameXdate	=	zip(file_name[:],date_taken[:])
	file_name[:]	=	[x for (x,y) in sorted(nameXdate, key=itemgetter(1))]
	#print(file_name[25:35])

	expoXdate	=	zip(exposure_time,date_taken)
	exposure_time[:]=	[x for (x,y) in sorted(expoXdate, key=itemgetter(1))]
	#print(expo_time)

	fstopXdate	=	zip(fstop_value,date_taken)
	fstop_value[:]	=	[x for (x,y) in sorted(fstopXdate, key=itemgetter(1))]
	#print(fstop_value)

	isoXdate	=	zip(iso_value,date_taken)
	iso_value[:]	=	[x for (x,y) in sorted(isoXdate, key=itemgetter(1))]
	#print(iso_value)

	starXdate	=	zip(is_1_star,date_taken)
	is_1_star[:]	=	[x for (x,y) in sorted(starXdate, key=itemgetter(1))]
	#print(iso_value)

		

	#extend the values of the reference frames through interpolation 
def interpolateParameters(parameters, is_1_star, smoothening):
    N_images    = len(parameters[0])
    
    
    #in case only one reference frame is set 
    	#we just copy this value for all parameters
    if sum(is_1_star) == 1:
        ind	= is_1_star.index(1)
        for v in parameters:
            v[:]	= [v[ind]]*N_images 
             
             
    #if there are multile reference frames
        #we first perform a linear interpolation
        #then one can smoothen this with a gauss curve
    if sum(is_1_star) > 1:
        i1	= 0
        ind	= is_1_star.index(1)


        #add values on the ends of the interval
        for v in parameters:
            v[0]	= v[ind]

        r_is_1_star	= is_1_star[:]
        r_is_1_star.reverse()
        ind	= r_is_1_star.index(1)
	
        for v in parameters:
            v[-1]	= v[-ind-1]
            is_1_star[-1] = 1


        #now interpolate!
        for i2 in range(1,N_image):
            if is_1_star[i2] == 1:
                
                for v in parameters:
                    a	= float(v[i2] - v[i1]) / float(i2-i1)     #the tangent of the 2 points 
                    
                    for j in range(1,i2-i1):                    #between the two points it is linear  
                        v[i1+j]	= v[i1] + j*a 
                #step to the next interval     
                i1	= i2
                
    #smoothening 
    if 25*smoothening>=2:
        frames	= float(25*smoothening) #1/2 width in frames
        w	= int(3*frames)            #the size of the convolution window is 3*relaxation distance
        
        #lets pad the data
        for v in parameters:
        	a 	= v[1] - v[0]
        	pre	= list(np.arange(-w,0)*a + v[0])
        	a 	= v[-2] - v[-1]
        	post	= list(np.arange(0,w)*a + v[-1])
        	v[:]	= pre + v + post
                  
         #create convolving Gaussian
        window	= 1.0/np.sqrt(2*np.pi)/frames*np.exp(-np.arange(-w,w,1)**2/2/frames**2)
         
         #smoothen!!
        for v in parameters:
            v[:]	= np.convolve(v,window,'same')[w:-w]
            
        
        
def prepareRefExp(rExp_comp, EV_calc_local, is_1_star):
    
    for i in range(len(is_1_star)):
        if is_1_star[i] == 1:
            rExp_comp[i] -= EV_calc_local[i]
    
    
            
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>MAIN<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<í

#start timer
start	= time.time()

file_name	= [] #list of files we are manipulating
exposure_time	= [] #list of integration times in seconds
fstop_value	= [] #list of normed F-stops (apertures)
iso_value	= [] #list of sensor gain in ISO units
date_taken	= [] #list of dates in seconds within 24 hours (see extractSidecarData function)
is_1_star	= [] #list of one star evaluation statuses --> master frames

#this finction returns the abowe 
extractSidecarData(inDir, file_name, exposure_time, fstop_value, iso_value, date_taken,	is_1_star)

#ordering the results
timeOrdData(file_name, exposure_time, fstop_value, iso_value, date_taken, is_1_star)

#the number of images
N_image 	= len(file_name)

#point test
#print(file_name)

if v_read_data == 1:
	T1	= max(date_taken)-min(date_taken)
	T2	= len(date_taken)
	print("==================================================OOO==================================================")
	print("The time difference between the first and last images: {0:10.2f} sec".format(T1))
	print("At 25fps this is {0:5.1f} mins, 30 fps: {1:5.1f} mins (30).".format((T2/25.0/60.0),  (T2/30.0/60.0)) )
	print("At 25fps this gives a speedup of {0:5.1f}x.".format( T1/(T2/25.)) )
	#print("From now we use 25fps for the calculations.")
	print("We are using {0:2d} master frames.".format(sum(is_1_star)) )


if b_print == 1:
	print("************")
	print("XMP read: {0:5.1f}".format(time.time()-start) )
	prev = time.time()
	print("************")
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////#SET EXP
EV_set		= [] #list of Exposure Values calculated based on the camera settings (arbitrary units, mainly due to the aperture)
luminanace_set	= [] #list of set luminance values (other arbitrary units)

#calculate EV based on just obtained data
for i in range(0,len(file_name)):
	luminanace_set.append(exposure_time[i]*iso_value[i]/100/(fstop_value[i])**2)
	EV_set.append(-math.log(luminanace_set[-1],2)) #EV 
#point test
#print(EV_set)

if b_print == 1:
	print("************")
	print("Calc from SET: {0:5.1f}".format(time.time()-prev) )
	prev = time.time()
	print("************")

#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////#CALC EXP
#now comes the ricky part where we calculate the real exposure based on the low resolution JPGs. 

EV_calc_local	=[] #list of calculated exposure values according to one of the methods LOCAL version 	---- for minimum perception
EV_calc_global	=[] #same but GLOBAL version								---- for protecting highlights

#the different methods are equally likely to work in my opinion, the decision is made based on human perception   -- released version only contains the best

est.fittingEstimator(inDir, file_name, EV_calc_local, EV_calc_global, lT.percentile_61_3ColorLookUp)


if b_print == 1:
	print("************")
	print("Calc from IMAGES: {0:5.1f}".format(time.time()-prev) )
	prev = time.time()
	print("************")
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////#CORRECTIONS
#now that we have the EV we have to calculate the EXPOSURE DIFFERENCE which has to be applied on the images

local_e_EV		= EV_calc_local[:]	#list of local expected exposures 
global_e_EV		= EV_calc_global[:]	#list of global expected exposures 
calc_jump_flag		= [0]*N_image		#if the estimated exposure change is greater than THRESHOLD then set a FLAG

l_set_comp		= []		#local comensations based on SET values (EXIF data)
g_set_comp		= [0]*N_image	#same but global
set_jump_flag		= [0]*N_image	#if there is a change in the SETTINGS set a FLAG

#calculate compensation based on EXIF
for i in range(0,N_image-1):
	Delta_EV	= EV_set[i] - EV_set[i+1]
	#if settings changed
	if Delta_EV>0:
		set_jump_flag[i+1]	= 1

	g_set_comp[i+1]	= g_set_comp[i]+Delta_EV

if v_set_differences == 1:
	print("==================================================111==================================================")
	print("The exposure was changed {0:2d} times during the shoot through the settings.".format(sum(set_jump_flag)) )

#calculate local compensation based on JPG and apply it 
#the local corrections have the same effect on the 		LOCAL EXPECTED EV 
#								GLOBAL EXPECTED EV!!

counter	= 0	#how many local corrections have been performed? 
carry	= 0	#local cumulant correction 

comp	= [0]*N_image   #list of ^^

for i in range(0,N_image-1):
	Delta_EV	= EV_calc_local[i] - EV_calc_local[i+1] #for the local corrections a measure with greater local correlation is advised, this is MACRO

	#local corrections are only performed if the difference is great enough (>0.1, or EV_treshold)
	if abs(Delta_EV)>EV_threshold:
		#step with indicators
		counter			+= 1
		calc_jump_flag[i+1]	= 1
		#calc cumulant	
		carry			+= Delta_EV

		if v_changes == 1:
			if counter == 1:
				print("==================================================222==================================================")
			print("Change No. {0:3d} The difference of the indicators:".format(counter) )
			print("T1: {0:3.1f} sec; T2: {1:3.1f} sec.            Delta SET: {2:1.3f}; Delta JPG {3:1.3f}".format(exposure_time[i],exposure_time[i+1],EV_set[i] - EV_set[i+1], -Delta_EV) )
			#print("Delta SET: {0:1.1f}; Delta JPG {1:1.1f}".format(EV_set[i] - EV_set[i+1], Delta_EV) )
		
	#add to lists
	comp[i+1]		=  carry
	local_e_EV[i+1]		+= carry
	global_e_EV[i+1]	+= carry

if v_changes == 1:
	tmp	= [abs(EV_set[i] - EV_set[i+1]+ EV_calc_local[i] - EV_calc_local[i+1]) for i in range(0,N_image-1)]
	print("====")
	print("The maximal deviation is: {0:1.2f}, and the average absolute deviation is: {1:1.2f}".format(np.max(tmp), np.average(tmp) ) )	
	print("The maximal difference in enviromental exposure value was {0:1.2f} EV.".format(max(global_e_EV)-min(global_e_EV)))	

#now that local adjustments are taken care of it is time to make global adjustments, so that the real (environmental) exposure can be mimicked
#this is done by substracting a trend line from the timeline of the EVs. 
#This is actually a poliomial in our case of order 3. 

#fit polynomial on the time series. here the GLOBAL index is used to PRESERVE HIGHLIGHTS
z	= np.polyfit( range(1,N_image+1), global_e_EV,3)
poly	= np.poly1d(z)
fitted	= poly(range(1,N_image+1))

#set the strength of the "light change over time" effect
if exp_ref == 0: 
    finetuner   = fitted*(exp_preserv*1.0/max(abs(fitted)))

if ref_on == 1: 
    #define values to be found
    white_balance	= [0]*N_image
    tint		    = [0]*N_image
                
    rExp_comp   = [0]*N_image
    contrast	= [0]*N_image
    
    whites		= [0]*N_image
    blacks		= [0]*N_image
    shadows		= [0]*N_image
    highlights	= [0]*N_image
                 
    clarity		= [0]*N_image
    vibrance	= [0]*N_image
    saturation	= [0]*N_image
                 
    #compacify to pass for the reader function
    parameters  = [white_balance, tint, rExp_comp, contrast, whites, blacks, shadows, highlights, clarity, vibrance, saturation]
    
    #get reference data
    extractReferenceData(inDir, file_name, is_1_star, 
                         parameters)
    
    #the exposure needs a minimal preprocessing 
    rExp_comp_old = rExp_comp[:]
    prepareRefExp(rExp_comp, EV_calc_local, is_1_star)

    #do interpollation
    interpolateParameters(parameters, is_1_star, smoothening)

    if exp_ref==1:
        #compesnate in case exp_ref is used
        finetuner   = rExp_comp


Exposure_Compensation	= -np.array( EV_calc_local ) + finetuner  
#print(Exposure_Compensation)


if b_print == 1:
	print("************")
	print("Calculating Corrections: {0:5.1f}".format(time.time()-prev) )
	prev = time.time()
	print("************")

#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////#OUTPUT

#writeOutputIntoXMP1(inDir, outDir, file_name, Exposure_Compensation, calc_jump_flag)
writeOutputIntoXMP2(inDir, outDir, file_name, Exposure_Compensation, parameters, calc_jump_flag, is_1_star)


if b_print == 1:
	print("************")
	print("Writng Sidecars: {0:5.1f}".format(time.time()-prev) )
	prev = time.time()
	print("************")



if b_print == 1:
	print("/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\")
	print("All in all, the run took: {0:5.1f} seconds.".format(time.time()-start) )
	print("/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\")



#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////#PLOTS

if v_plot_exp				== 1:
    plt.figure(1)
    
    plt.subplot(311)
    plt.ylabel("EV")
    plt.title("Exposure curves")
    l1, = plt.plot( range(1,N_image+1), EV_set, label='Set EXP')
    l15, = plt.plot(range(1,N_image+1), EV_calc_global, label='Calculated EXP')
    l2, = plt.plot(range(1,N_image+1), global_e_EV, label='Calculated environmental EXP')
    l3, = plt.plot(range(1,N_image+1), fitted, label="poly" )
    l4, = plt.plot(range(1,N_image+1), comp, label="comp" )
    plt.legend(handles=[l1, l15, l2, l3, l4], loc="upper right")
    
    plt.subplot(312)
    plt.ylabel("EV")
    l5, = plt.plot( range(1,N_image+1), np.array(global_e_EV)-fitted, label='global_e_EV-poly')
    plt.legend(handles=[l5], loc="upper right")
    
    plt.subplot(313)
    plt.ylabel("EV")
    plt.xlabel("Imange number #n")
    l6, = plt.plot( range(1,N_image+1), Exposure_Compensation, label='EXP comp')
    plt.legend(handles=[l6], loc="upper right")
    

if v_plot_reference == 1 and ref_on == 1:
    plt.figure(2)
    
    plt.subplot(311)
    plt.title("Reference curves")
    plt.ylabel("Tempereature [K]")
    l1, = plt.plot( range(1,N_image+1), white_balance, label='Temp')
    plt.legend(handles=[l1], loc="upper right")

    plt.subplot(312)
    plt.ylabel("Appr. unit")
    l2, = plt.plot(range(1,N_image+1), tint, label='Tint')
    l3, = plt.plot(range(1,N_image+1), whites, label="Whites" )
    l4, = plt.plot(range(1,N_image+1), blacks, label="Blacks" )
    l5, = plt.plot(range(1,N_image+1), highlights, label="Highlights" )
    l6, = plt.plot(range(1,N_image+1), shadows, label="Shadows" )
    plt.legend(handles=[l2, l3, l4, l5, l6], loc="upper right")

    plt.subplot(313)
    plt.ylabel("Contrast [perc.]")
    plt.xlabel("Imange number #n")
    l7, = plt.plot(range(1,N_image+1), contrast, label="Contrast" )
    plt.legend(handles=[l7], loc="upper right")
    
    plt.figure(3)
    
    plt.subplot(211)
    plt.title("Reference curves 2")
    plt.ylabel("Appr. unit")
    l1, = plt.plot( range(1,N_image+1), vibrance, label='Vibrance')
    l2, = plt.plot( range(1,N_image+1), saturation, label='Saturation')
    plt.legend(handles=[l1, l2], loc="upper right")
    
    plt.subplot(212)
    plt.ylabel("Appr. unit")
    plt.xlabel("Imange number #n")    
    l3, = plt.plot( range(1,N_image+1), clarity, label='Clarity')
    plt.legend(handles=[l3], loc="upper right")
    
    if exp_ref == 1:
         plt.figure(4)
         plt.ylabel("EV")
         plt.xlabel("Imange number #n")
         l1, = plt.plot( range(1,N_image+1), Exposure_Compensation, label='EXP comp')
         l2, = plt.plot( range(1,N_image+1), -np.array( EV_calc_local ), label='Calculated compensation')
         
         l3, = plt.plot( [1+i for i in range(N_image) if is_1_star[i] == 1 ], [rExp_comp_old[i]  for i in range(N_image) if is_1_star[i] == 1 ], "x",label='Reference')
         plt.legend(handles=[l1, l2, l3], loc="upper right")
    
if (  v_plot_reference == 1 and ref_on == 1 ) or  (v_plot_exp == 1) :
    plt.show()








