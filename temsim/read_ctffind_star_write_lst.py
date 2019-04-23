#!/usr/bin/env python

import math, os, sys
try:
	from optparse import OptionParser
except:
	from optik import OptionParser

def main():
	(starfile,metadataline,ny,output) =  parse_command_line()
	f=open(starfile,"r")
	starline=f.readlines()
	s=open(output,"w")
	s.write("#LST\n")
	b=open("eulerfile.txt","r")
	eufline=b.readlines()
	
	for i in range(4,metadataline):
		tmp=str(starline[i].split()[0])
		if(tmp=="_rlnMicrographName"):
			MN=int(starline[i].split('#')[1])
		if(tmp=="_rlnDefocusU"):
			dfu_n=int(starline[i].split('#')[1])
		if(tmp=="_rlnDefocusV"):
			dfv_n=int(starline[i].split('#')[1])
		if(tmp=="_rlnDefocusAngle"):
			dfa_n=int(starline[i].split('#')[1])	
		if(tmp=="_rlnMagnification"):
			mag_n=int(starline[i].split('#')[1])
		if(tmp=="_rlnDetectorPixelSize"):
			dect_n=int(starline[i].split('#')[1])
		if(tmp=="_rlnVoltage"):
			vol_n=int(starline[i].split('#')[1])
		if(tmp=="_rlnSphericalAberration"):
			cs_n=int(starline[i].split('#')[1])
		if(tmp=="_rlnAmplitudeContrast"):
			ampcont_n=int(starline[i].split('#')[1])
#	print(MN)
#	print(dfu_n)
#	print(dfv_n)
#	print(dfa_n)
#	print(mag_n)
#	print(dect_n)
#	print(cs_n)
#	print(ampcont_n)
#	print(vol_n)
	for i in range(metadataline,len(starline)):
		temp="0\t"
		if(starline[i].split()):
			name_number=int((starline[i].split()[MN-1]).split('_')[5].split('-')[0])-1
			euler1=float(eufline[name_number].split()[0])
			euler2=float(eufline[name_number].split()[1])
			euler3=float(eufline[name_number].split()[2])
			old_df=float(((starline[i].split()[MN-1]).split('.mrc')[0]).split('_')[7])
			temp=temp+str(starline[i].split()[MN-1])+"\t"
			dfv=float(starline[i].split()[dfv_n-1])
			dfu=float(starline[i].split()[dfu_n-1])
			angle=float(starline[i].split()[dfa_n-1])
			defocus = (dfu+dfv)/20000.
			dfdiff = math.fabs(dfu-dfv)/20000.
			apix=float(starline[i].split()[dect_n-1])/float(starline[i].split()[mag_n-1])*10000.0
			voltage=float(starline[i].split()[vol_n-1])
			cs=float(starline[i].split()[cs_n-1])
			ampc=float(starline[i].split()[ampcont_n-1])
			center=ny/2
			if dfu>dfv:
				dfang = math.fmod(angle+360.+90., 360.)
			else:
				dfang=angle
			temp=temp+"defocus="+str(defocus)+"\tdfdiff="+str(dfdiff)+"\tdfang="+str(dfang)+"\t"
			temp=temp+"apix="+str(apix)+"\tvoltage="+str(voltage)+"\tcs="+str(cs)+"\tampcont="+str(ampc)+"\t"
			defocuschange=defocus-old_df/10000.0
			temp=temp+"euler="+str(euler1)+","+str(euler2)+","+str(euler3)+"\tcenter="+str(center)+","+str(center-1)+"\tdefocusChange="+str(defocuschange)+"\n"
			s.write(temp)

	f.close()	
	s.close()
	
def parse_command_line():
	usage="%prog <starfile> <metadata line +4> <ny> <output>"
	parser = OptionParser(usage=usage, version="%1")
	
	if len(sys.argv)<5: 
		print "<starfile> <metadata line +4> <ny> <output>"
		sys.exit(-1)
	
	(options, args)=parser.parse_args()
	starfile=args[0]
	metadataline=int(args[1])
	ny=float(args[2])
	output=args[3]
	return (starfile,metadataline,ny,output)

if __name__== "__main__":
	main()


			
