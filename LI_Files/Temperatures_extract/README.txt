Using "temp_write.py" and "mult_temp_write.py"

The files "temp_write.py" and "mult_temp_write.py" are macros which read-in tab-delimited .txt files of temperature data and write the unix timestamps and temp data to "Temp_Output".
Originally "temp_write.py" was written to accomodate the timestamps acquired through the old "doit_with_date.sh" bash script. So that necessary portion is included in the files but commented out.

temp_write.py:
	It is pretty straightforward how to use this. Just run it with a tab-delimited .txt file and it will produce "Temp_Output". 
	* Make sure to check the format of your columns in the excel file before running this on the .txt file. *

mult_temp_write.py:
	This file is used if you want to acquire temperature data from multiple files and write to Temp_Output. 
	* There must be no other files in the directory but this one and the .txt files for it to work. *
	* Make sure to check the format of your columns in the excel files before running this on the .txt files. *
