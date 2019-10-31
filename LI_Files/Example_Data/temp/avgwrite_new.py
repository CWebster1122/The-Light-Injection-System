#****************************************************************************************************************************************************
#This code is meant to process the temperature data for the LI system, and output the corresponding data to file "Temp_Output", which can then be used in the LI_Analysis.C macro. The logic of this
#is heavily reliant on the excel structure in which the data was originally taken, so if any part of that process was changed, this code will break pretty easily, but it is generaly an 
#easy/quick fix if you know what the procedure changed from and to.
#*****************************************************************************************************************************************************

#!/usr/bin/env python
import sys
import os
import time as t
start_time = t.mktime((2016,1,1,0,0,0,0,1,-1)) #time stamp of January 1, 2016
#print("starting time: " + str(t.localtime(start_time)))
dir = os.listdir('.') #obtain list of all files in current directory
dir.remove(os.path.basename(sys.argv[0])) #remove this file's name from the list
output_file=open("Temp_Output", "a")
dir.remove("Temp_Output")
dir.sort()
print(str(dir))
time_array = []
for i in range(len(dir)): #loop through all files in the directory
    original_temp_data = open(dir[i], 'rU') #open the i_th file in the directory list
    line_number = 1
    for line in original_temp_data:
        if(line_number == 1):
            line_number = line_number + 1
            continue
        line_split = list(map(str.strip, line.split("\t")))
        #print line_split
        print((line_split[1])+(line_split[2]))
        file_time = t.strptime(line_split[1]+line_split[2],'%m/%d/%y%H:%M:%S')
        #print(file_time)
        time_elapsed = (t.mktime(file_time) - start_time)
        #print(str(t.mktime(file_time)) + " " + str((t.mktime(file_time) - start_time)))
        time_array.append(time_elapsed)

        #print("Readout: " + str(t.localtime(time_elapsed)))
        #print("Readout without subtraction: " + str(t.localtime(t.mktime(file_time))))
        line_number = line_number + 1
        #print("\n \n")
        output_file.write(str(time_elapsed) + " " + str(line_split[3]) + " " + str(line_split[5]) + " " + str(line_split[7]) + " " + str(line_split[9]) + "\n")
