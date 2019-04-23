#!/usr/bin/env python

import math, os, sys, random
try:
	from optparse import OptionParser
except:
	from optik import OptionParser

def main():
	(inpdb,apix,ysize,alt,az,phi,output,long_coord) =  parse_command_line()
	a=open(inpdb,"r")
	pdb_data=a.readlines()
	o=open(output,"w")
	PI=3.14159265359
	angle1=-1.0*az
	angle2=180.0-alt
	angle3=-1.0*phi
	angle1=angle1*PI/180.0
	angle2=angle2*PI/180.0
	angle3=angle3*PI/180.0
	wobble=0
	occ=1
	element=[]
	ele_x=[]
	ele_y=[]
	ele_z=[]
	count_ele=0
	tmp=""
	ele_num=0
	abandon_ele=0
	x_min=0.0
	x_max=0.0
	x_aver=0.0
	y_min=0.0
	y_max=0.0
	y_aver=0.0
	z_min=0.0
	z_max=0.0
	z_aver=0.0
	x=0.0
	y=0.0
	z=0.0
	for i in range(0,len(pdb_data)):
		if(pdb_data[i].split()):
			if(str(pdb_data[i].split()[0])=="ATOM"):
				tmp=str(pdb_data[i].split()[2])[0]
				if(tmp=="C"):
					ele_num=6
				elif(tmp=="N"):
					ele_num=7
				elif(tmp=="O"):
					ele_num=8
				elif(tmp=="P"):
					ele_num=15
				elif(tmp=="S"):
					ele_num=16
				else:
					ele_num=0
				if(ele_num!=0):
					element.append([])
					ele_x.append([])
					ele_y.append([])
					ele_z.append([])
					element[count_ele]=ele_num
					if(long_coord==0):
						x=float(pdb_data[i].split()[6])
						y=float(pdb_data[i].split()[7])
						z=float(pdb_data[i].split()[8])
					else:
						x=float(pdb_data[i][30:38])
						y=float(pdb_data[i][38:46])
						z=float(pdb_data[i][46:54])
					x_aver+=x
					y_aver+=y
					z_aver+=z
									
					ele_x[count_ele]=x
					ele_y[count_ele]=y
					ele_z[count_ele]=z
					count_ele+=1
				else:
					abandon_ele+=1
				if(count_ele%1000000==0 and count_ele!=0):
					print("Prog has read "+str(count_ele)+" atoms.")
	print ("Finish reading pdb file, processing now.")
	print("Total number of "+str(count_ele)+" atoms, "+str(abandon_ele)+" atoms were dropped.")
	x_aver=x_aver/float(count_ele)
	y_aver=y_aver/float(count_ele)
	z_aver=z_aver/float(count_ele)

	a.close()

	xoff=x_aver
	yoff=y_aver
	zoff=z_aver
	cr = math.cos(angle1)
	sr = math.sin(angle1)
	ct = math.cos(angle2)
	st = math.sin(angle2)
	ctt = math.cos(angle3)
	stt = math.sin(angle3)
	
	for i in range (0,count_ele):
		ele_x[i] = ele_x[i] - xoff
		ele_y[i] = ele_y[i] - yoff
		ele_z[i] = ele_z[i] - zoff
		
		x0 = ele_x[i]
		y0 = ele_y[i]
		ele_x[i] = cr*x0 - sr*y0
		ele_y[i] = sr*x0 + cr*y0

		y0 = ele_y[i]
		z0 = ele_z[i]
		ele_y[i] = ct*y0 + st*z0
		ele_z[i] = -st*y0 + ct*z0
		
		x0 = ele_x[i]
		y0 = ele_y[i]
		ele_x[i] = ctt*x0 - stt*y0
		ele_y[i] = stt*x0 + ctt*y0
		
		if(ele_x[i]<x_min):
			x_min=ele_x[i]
		if(ele_y[i]<y_min):
			y_min=ele_y[i]
		if(ele_z[i]<z_min):
			z_min=ele_z[i]
		if(ele_x[i]>x_max):
			x_max=ele_x[i]
		if(ele_y[i]>y_max):
			y_max=ele_y[i]
		if(ele_z[i]>z_max):
			z_max=ele_z[i]
		if(i%1000000==0 and i!=0):
					print("Prog has rotate "+str(i)+" atoms.")	
	x1=math.fabs(x_min)
	x2=math.fabs(x_max)
	y1=math.fabs(y_min)
	y2=math.fabs(y_max)
	z1=math.fabs(z_min)
	z2=math.fabs(z_max)
	print("x1="+str(x1)+" x2="+str(x2)+" y1="+str(y1)+" y2="+str(y2)+" z1="+str(z1)+" z2="+str(z2))
#	print("x_min="+str(x_min)+" x_max="+str(x_max)+" y_min="+str(y_min)+" y_max="+str(y_max)+" z_min="+str(z_min)+" z_max="+str(z_max))
	print("x_aver="+str(x_aver)+" y_aver="+str(y_aver)+" z_aver="+str(z_aver))
	if(x1<=x2):
		x1=x2
	if(y1<=y2):
		y1=y2
	if(z1<=z2):
		z1=z2
	ax=ysize*apix
	cz=math.fabs(z_max-z_min)
	if(1.25*x1>=ax):
		ax=1.25*x1
	if(1.25*y1>=ax):
		ax=1.25*y1
	if(ax>ysize*apix):
		print ("The given apix is too small, change apix to "+str(ax/ysize))
	else:
		ax=ysize*apix
	print("Coord X range is "+str(x_min)+" to "+str(x_max))
	print("Coord Y range is "+str(y_min)+" to "+str(y_max))
	print("Coord Z range is "+str(z_min)+" to "+str(z_max))			
	xoff=ax/2
	yoff=ax/2
	zoff=math.fabs(z_min)
	o.write("pdb2xyz translation of "+str(inpdb)+"\n")
	o.write("%16.2f"%ax+"%16.2f"%ax+"%16.2f"%cz+"\n")
	
	for i in range(0,count_ele):
		x0=ele_x[i]+xoff
		y0=ele_y[i]+yoff
		z0=ele_z[i]+zoff
		o.write("%5d"%element[i]+"%14.3f"%x0+"%14.3f"%y0+"%14.3f"%z0+"%14.3f"%occ+"%14.3f"%wobble+"\n")
	o.write("-1\n")
	o.close()
	
def parse_command_line():
	usage="%prog <input pdb> <apix> <ySize> <EMAN alt> <az> <phi> <output> <long coord>"
	parser = OptionParser(usage=usage, version="%1")
	
	if len(sys.argv)<9: 
		print "<input pdb> <apix> <ySize> <EMAN alt> <az> <phi> <output> <long coord>"
		sys.exit(-1)
	
	(options, args)=parser.parse_args()
	
	inpdb = args[0]
	apix = float(args[1])
	ysize= float(args[2])
	alt=float(args[3])
	az=float(args[4])
	phi=float(args[5])
	output = args[6]
	long_coord=int(args[7])
	return (inpdb,apix,ysize,alt,az,phi,output,long_coord)

if __name__== "__main__":
	main()


			
