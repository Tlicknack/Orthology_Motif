import os
import shutil
import re
seq = list()

file_folder = "/N/dc2/scratch/tlicknac/All_Aurelias_Upstream/"   		#Path upstream sequences

seq = os.listdir(file_folder)                         				#seq will have list of all file names

for i in seq:                                                           	#iterate through list of files
	real_i = file_folder + i
        filepath = "/N/dc2/scratch/tlicknac/All_Aurelias_MEME_Results/" +     #Location of output
        try:
            	comm1 = "mkdir " + "/N/dc2/scratch/tlicknac/All_Aurelias_MEME_Results/"                           	
        except:
               	pass
        os.system(comm1)
        comm2 = "/N/u/tlicknac/Carbonate/software/meme/install/bin/meme " + real_i + " -dna -oc " + filepath + " -nostatus -time 18000 -maxsize 60000 -nmotifs 5 -minw 5 -maxw 18"   #run my meme
        os.system(comm2)

