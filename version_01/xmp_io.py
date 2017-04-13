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


#this file is responsible for the xmp INPUT/OUTPU operations 
#																for Lr 2012 process version
#it produces files that Lr accepts, but it knows that these are externally modified.
#the code is pretty much self explanatory, however a bit barbaric.


#-------------------------------Official INCLUDES-------------------------------
import os				#to access files
import matplotlib.pyplot as plt		#to display modifications
import time				#for benchmarking
import re				#for using regular expressions
import math				#for EV calculation (log)
import numpy as np			#for fitting polynomials and evaluating them

#-------------------------------Function DEFINITIONS----------------------------

#read sidecar information for calculations and return them 
def extractSidecarData(inDir, 
			file_name,
			exposure_time,
			fstop_value,
			iso_value,
			date_taken,
			is_1_star):

	#list filenames from inDir
	listDIR = os.listdir(inDir) 

	#we iterate through the filenames in the specified folder	
	for name in listDIR:
		#check if it is an original sidecar file (xmp)
		if name.endswith(".xmp") and not name.endswith("tmp.xmp"):

			#open the sidecar file, read data 
			f = open(inDir+"/"+name, 'r')
			xmp_dat	= f.read()
		
			file_name.append(name)


		#-------find the data needed using regular expressions
			strings = re.findall(r"exif:ExposureTime.*$", xmp_dat, flags=re.MULTILINE)	#integration time in seconds
			#point test
			#print(strings[0])
			value	=	strings[0].split('"')[1]
			value1	=	value.split("/")[0]
			value2	=	value.split("/")[1]
			exposure_time.append(float(value1)/float(value2))
			#point test
			#print(exposure_time[-1])

			

			strings = re.findall(r"exif:FNumber.*$", xmp_dat, flags=re.MULTILINE) 		#fstop is the normalized aperture
			#point test
			#print(strings[0])
			value	=	strings[0].split('"')[1]
			value1	=	value.split("/")[0]
			value2	=	value.split("/")[1]
			fstop_value.append(float(value1)/float(value2))
			#point test
			#print(exposure_time[-1])f.seek(0)

		

			strings = re.findall(r"<rdf:li>.*$", xmp_dat, flags=re.MULTILINE)		#sensor gain in ISO units
			#point test
			#print(strings[0])
			value	=	strings[0].split('>')[1].split('<')[0]
			iso_value.append(float(value))
			#point test
			#print(exposure_time[-1])

		
		

			strings = re.findall(r"exif:DateTimeOriginal.*$", xmp_dat, flags=re.MULTILINE)	#time the picture was taken in seconds ON GIVEN DAY !!!
			#point test
			#print(strings)
			value	=	strings[0].split('"')[1]
			hundr	=	value.split(".")[-1]
			leftovr	= 	value.split(".")[-2]
			sec	 =	leftovr.split(":")[-1]
			minute	=	leftovr.split(":")[-2]
			hour	=	(leftovr.split(":")[-3]).split("T")[-1]
			date_taken.append(float(hundr)/100+float(sec)+float(minute)*60+float(hour)*3600)
			#point test
			#print(date_taken[-1])




			strings = re.findall(r"xmp:Rating.*$", xmp_dat, flags=re.MULTILINE)		#is there a one star sign on this image (if yes it is going to be a reference image)
			#print(strings)
			if len(strings) > 0:
				value	=	strings[0].split('"')[1]
				if value == '1':
					is_1_star.append(1)
				else:
					is_1_star.append(0)
		
			else:
				is_1_star.append(0)
		#-------data obtained
		
			#close file 
			f.close()

#make output sidecar files
def writeOutputIntoXMP2(inDir, outDir, file_name, Exposure_Compensation, parameters, calc_jump_flag, is_1_star):
	#if there is no output folder, create it!
	if not os.path.exists("./"+outDir):
		os.makedirs("./"+outDir)

	#parameters
	white_balance	= parameters[0]
	tint		= parameters[1]
	rExp_comp	= parameters[2]
	contrast	= parameters[3]
	whites		= parameters[4]
	blacks		= parameters[5]
	shadows		= parameters[6]
	highlights	= parameters[7]
	clarity		= parameters[8]
	vibrance	= parameters[9]
	saturation	= parameters[10]

	#make an output for each file on the list
	
	#counter of flies
	i=-1

	for name in file_name:
		#in theory all are XMPs, but let's check
		if name.endswith(".xmp"):
			i+=1
			#open the original XMP file, most of it will stay unchanged 
			f = open("./"+inDir+"/"+name, 'r')
			#create the output file 
			k = open("./"+outDir+"/"+name.split(".")[0]+".xmp",'w+')

			
			#flags to help following where we are in the XMP file

			red_label_flag		= 0	#the image is labeled RED
			green_label_flag	= 0	#the image is labeled GREEN

			last_crs_line		= 0	#Attributes are over and the LAST EDIT LINE reached. If so and none of our EDIT data has been output yet now is the time.

			exp_flag		= 0	#exposure already printed

			wb_flag			= 0	#flag overwritten to CUSTOM
			temp_flag		= 0	#actual temp vale set
			tint_flag		= 0	#tint value set

			contrast_flag	= 0	#contrast calue set

			clar_flag	   = 0
			vibr_flag	   = 0
			satu_flag	   = 0
			
			label_veto	   = 0 #we wont include the original color labels... label duplication messes things up
			
			p_flags		= [exp_flag, wb_flag, temp_flag, tint_flag, contrast_flag, clar_flag, vibr_flag, satu_flag]

			#iterate through all the lines, and insert edit data where appropriate
			for line in f.read().split("\n"):
				

				label_veto		= 0
				if "xmp:Label=\"Red\"" in line:
					red_label_flag		= 1
					label_veto		= 1
				if "xmp:Label=\"Green\"" in line:
					green_label_flag	= 1
					label_veto		= 1
				if "crs:RawFileName=" in line:
					last_crs_line	= 1
	
				if "crs:Exposure2012" in line:
					exp_flag	= 1
					if Exposure_Compensation[i]<0:
						char	= ""
					else:
						char	= "+"
					line	= "   crs:Exposure2012=\"" + char + str(int(Exposure_Compensation[i]*100)/100.0) + "\""
					if calc_jump_flag[i] == 1 and red_label_flag == 0:
						k.write("   xmp:Label=\"Red\"\n")
						
					if i % 10 == 0 and green_label_flag == 0 and calc_jump_flag[i] == 0:
						k.write("   xmp:Label=\"Green\"\n")
				if sum(is_1_star)>0:#--------------------------------------in case we have master images we have other metadata to sync too!!
					if "crs:Temperature" in line:
						temp_flag	= 1
						line	= "   crs:Temperature=\"" + str(int(white_balance[i])) + "\""
					if "crs:WhiteBalance=" in line:
						wb_flag		= 1
						line	= "   crs:WhiteBalance=\"Custom\""
					if "crs:Tint" in line:
						tint_flag	= 1
						if tint[i]<0:
							char	= ""
						else: 
							char	= "+"
						line	= "   crs:Tint=\"" + char + str(int(tint[i])) + "\""


					if "crs:Contrast2012=" in line:
						contrast_flag	= 1
						if contrast[i]<0:
							char	= ""
						else:
							char	= "+"
						line	= "   crs:Contrast2012=\"" + char + str(int(contrast[i])) + "\""


					if "crs:Whites2012=" in line:
						if whites[i]<0:
							char	= ""
						else:
							char	= "+"
						line	= "   crs:Whites2012=\"" + char + str(int(whites[i])) + "\""
					if "crs:Highlights2012=" in line:
						if highlights[i]<0:
							char	= ""
						else:
							char	= "+"
						line	= "   crs:Highlights2012=\"" + char + str(int(highlights[i])) + "\""
					if "crs:Blacks2012=" in line:
						if blacks[i]<0:
							char	= ""
						else:
							char	= "+"
						line	= "   crs:Blacks2012=\"" + char + str(int(blacks[i])) + "\""
					if "crs:Shadows2012=" in line:
						if shadows[i]<0:
							char	= ""
						else:
							char	= "+"
						line	= "   crs:Shadows2012=\"" + char + str(int(shadows[i])) + "\""


					#"""
					if "crs:Clarity2012=" in line:
						clar_flag	= 1
						if clarity[i]<0:
							char	= ""
						else:
							char	= "+"
						line	= "   crs:Clarity2012=\"" + char + str(int(clarity[i])) + "\""
					if "crs:Vibrance=" in line:
						vibr_flag	= 1
						if vibrance[i]<0:
							char	= ""
						else:
							char	= "+"
						line	= "   crs:Vibrance=\"" + char + str(int(vibrance[i])) + "\""
					if "crs:Saturation=" in line:
						satu_flag	= 1
						if saturation[i]<0:
							char	= ""
						else:
							char	= "+"
						line	= "   crs:Saturation=\"" + char + str(int(saturation[i])) + "\""
					#"""
				#===========================#
				if last_crs_line == 1 and sum(is_1_star) == 0 and exp_flag == 0: # if for some reason there was not an EXP and ONLY exp is to be written
					last_crs_line	= 0
					exp_flag	= 1
					if Exposure_Compensation[i]<0:
						char	= ""
					else:
						char	= "+"
					k.write("   crs:Exposure2012=\"" + char + str(int(Exposure_Compensation[i]*100)/100.0) + "\"\n")
					if calc_jump_flag[i] == 1 and red_label_flag == 0:
						k.write("   xmp:Label=\"Red\"\n")
						
					if i % 10 == 0 and green_label_flag == 0 and calc_jump_flag[i] == 0:
						k.write("   xmp:Label=\"Green\"\n")
				if last_crs_line == 1 and sum(is_1_star) > 0 and (len(p_flags)-sum(p_flags))>0: # if for some reason some param entries were missing and we have MASTER images
					last_crs_line	= 0
					if exp_flag == 0:
						exp_flag	= 1
						if Exposure_Compensation[i]<0:
							char	= ""
						else:
							char	= "+"
						k.write("   crs:Exposure2012=\"" + char + str(int(Exposure_Compensation[i]*100)/100.0) + "\"\n")
						if calc_jump_flag[i] == 1 and red_label_flag == 0:
							k.write("   xmp:Label=\"Red\"\n")
							
						if i % 10 == 0 and green_label_flag == 0:
							k.write("   xmp:Label=\"Green\"\n")

						if whites[i]<0:
							char	= ""
						else:
							char	= "+"
						k.write("   crs:Whites2012=\"" + char + str(int(whites[i])) + "\"\n")

						if highlights[i]<0:
							char	= ""
						else:
							char	= "+"
						k.write("   crs:Highlights2012=\"" + char + str(int(highlights[i])) + "\"\n")

						if blacks[i]<0:
							char	= ""
						else:
							char	= "+"
						k.write("   crs:Blacks2012=\"" + char + str(int(blacks[i])) + "\"\n")

						if shadows[i]<0:
							char	= ""
						else:
							char	= "+"
						k.write("   crs:Shadows2012=\"" + char + str(int(shadows[i])) + "\"\n")
					if wb_flag == 0:
						wb_flag		= 1
						k.write("   crs:WhiteBalance=\"Custom\"\n")
					if temp_flag == 0:
						temp_flag	= 1
						k.write("   crs:Temperature=\"" + str(int(white_balance[i])) + "\"\n")
					if tint_flag == 0:
						tint_flag	= 1
						if tint[i]<0:
							char	= ""
						else: 
							char	= "+"
						k.write("   crs:Tint=\"" + char + str(int(tint[i])) + "\"\n")
					if contrast_flag == 0:
						contrast_flag	= 1
						if contrast[i]<0:
							char	= ""
						else:
							char	= "+"
						k.write("   crs:Contrast2012=\"" + char + str(int(contrast[i])) + "\"\n")
					#"""
					if clar_flag == 0:
						clar_flag	= 1
						if clarity[i]<0:
							char	= ""
						else:
							char	= "+"
						k.write("   crs:Clarity2012=\"" + char + str(int(clarity[i])) + "\"\n")
					if vibr_flag == 0:
						vibr_flag	= 1
						if vibrance[i]<0:
							char	= ""
						else:
							char	= "+"
						k.write("   crs:Vibrance=\"" + char + str(int(vibrance[i])) + "\"\n")
					if satu_flag == 0:
						satu_flag	= 1
						if saturation[i]<0:
							char	= ""
						else:
							char	= "+"
						k.write("   crs:Saturation=\"" + char + str(int(saturation[i])) + "\"\n")
					#"""
				#================#
				
				if label_veto == 0:			
					k.write( line +"\n")				
			f.close()
			k.close()

def extractReferenceData(inDir, file_name, is_1_star, parameters):
	N_image	 = len(file_name)

	#parameters
	white_balance	= parameters[0]
	tint		= parameters[1]
	rExp_comp	= parameters[2]
	contrast	= parameters[3]
	whites		= parameters[4]
	blacks		= parameters[5]
	shadows		= parameters[6]
	highlights	= parameters[7]
	clarity		= parameters[8]
	vibrance	= parameters[9]
	saturation	= parameters[10]


	for i in range(N_image):
		if is_1_star[i] == 1: #if it is a reference frame
		
			#open XMP
			f = open(inDir+"/"+file_name[i], 'r')
			xmp_dat	= f.read()
			f.close()
			
			#get reference values and add them to the lists
			#ZERO is not set for all values!!!!!!
			#whitebalance and tint are EXCEPTIONS
			
			
			strings = re.findall(r"crs:Temperature.*$", xmp_dat, flags=re.MULTILINE)
			if len(strings) > 0:
				value	=	strings[0].split('"')[1]
				white_balance[i]	= (float(value))
			else:
				white_balance[i]	= 5001  

	
			strings = re.findall(r"crs:Tint.*$", xmp_dat, flags=re.MULTILINE)
			#print(strings[0])
			if len(strings) > 0:
				value	=	strings[0].split('"')[1]
				tint[i]	= (float(value))
			else:
				tint[i]	= 10
		   
			
			
	
			strings = re.findall(r"crs:Exposure2012=.*$", xmp_dat, flags=re.MULTILINE)
			#print(strings[0])
			if len(strings) > 0:
				value	=	strings[0].split('"')[1]
				rExp_comp[i]	= (float(value))
			else:
				rExp_comp[i]	= 0
	
			strings = re.findall(r"crs:Contrast2012=.*$", xmp_dat, flags=re.MULTILINE)
			#print(strings[0])
			if len(strings) > 0:
				value	=	strings[0].split('"')[1]
				contrast[i]	= (float(value))
			else:
				contrast[i]	= 0
	
	
	
	
			strings = re.findall(r"crs:Whites2012.*$", xmp_dat, flags=re.MULTILINE)
			#print(strings[0])
			if len(strings) > 0:
				value	=	strings[0].split('"')[1]
				whites[i]	= (float(value))
			else:
				whites[i]	= 0
	
	
			strings = re.findall(r"crs:Blacks2012.*$", xmp_dat, flags=re.MULTILINE)
			#print(strings[0])
			if len(strings) > 0:
				value	=	strings[0].split('"')[1]
				blacks[i]	= (float(value))
			else:
				blacks[i]	= 0
	
	
			strings = re.findall(r"crs:Shadows2012.*$", xmp_dat, flags=re.MULTILINE)
			#print(strings[0])
			if len(strings) > 0:
				value	=	strings[0].split('"')[1]
				shadows[i]	= (float(value))
			else:
				shadows[i]	= 0
	
	
			strings = re.findall(r"crs:Highlights2012.*$", xmp_dat, flags=re.MULTILINE)
			#print(strings[0])
			if len(strings) > 0:
				value	=	strings[0].split('"')[1]
				highlights[i]	= (float(value))
			else:
				highlights[i]	= 0


			strings = re.findall(r"crs:Clarity2012.*$", xmp_dat, flags=re.MULTILINE)
			#print(strings[0])
			if len(strings) > 0:
				value	=	strings[0].split('"')[1]
				clarity[i]	= (float(value))
			else:
				clarity[i]	= 0


			strings = re.findall(r"crs:Vibrance=.*$", xmp_dat, flags=re.MULTILINE)
			#print(strings[0])
			if len(strings) > 0:
				value	=	strings[0].split('"')[1]
				vibrance[i]	= (float(value))
			else:
				vibrance[i]	= 0


			strings = re.findall(r"crs:Saturation=.*$", xmp_dat, flags=re.MULTILINE)
			#print(strings[0])
			if len(strings) > 0:
				value	=	strings[0].split('"')[1]
				saturation[i]	= (float(value))
			else:
				saturation[i]	= 0

