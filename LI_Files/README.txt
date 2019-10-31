The Three Key LI Data Analysis Files: Fitter.C, Run_Fitter.sh, and LI_Analysis.C
I will be thorough here so that I leave (hopefully) very few stones unturned.




1) Fitter.C	This macro reads-in .root data files, and then creates fits and extracts parameters from said files. 
		It then writes these parameters to 3 parameter files: "Parameters1.dat", "Parameters2.dat", and "Parameters3.dat".
		Parameters1.dat contains parameters extracted from the complete fits of Bi-207 and Am-241.
		Parameters2.dat and Parameters3.dat contain parameters extracted from restricted-range gaussian fits surrounding the peaks of Bi-207 and Am-241.
		This macro only writes parameters associated with fits that meet a defined chi^2/NDF threshold of < 8. 
		Otherwise, it will write a print statement to the file, e.g. 'ERROR', as an erasure indicator for Run_Fitter.sh.	
				

		* Make sure to check the name of the .root file you read-in on line 40 of Fitter.C *
		* Make sure to check Bool_t Write = 0 on line 20 when running on a single .root file (unless you want to write to the parameter files).
 		



2) Run_Fitter.sh	This is a bash script used to run Fitter.C on every .root file in the current directory. It also appends a timestamp at the end of each line in the parameter files.
			When it sees the string 'ERROR' written by Fitter.C, it will then delete the line (timestamp and all).		
					
			* Make sure the name of the input file in Fitter.C is "temp.root" when running Run_Fitter.sh. *
			* Make sure that Bool_t Write = 1 on line 20 of Fitter.C before running this bash script. *
			



4) LI_Analysis.C	This is our end-game analysis file which reads-in our parameter files, and then produces calculations, plots, and of course a "Money Plot".
			Specifically, its inputs are Parameters1.dat or Parameters2.dat or Parameters3.dat and Temp_Output. 
			It calculates Bi_pred/Bi_t, and the resulting time evolution plot is known as The Money Plot.    


			* Make sure to check what Parameter file you read-in on line 96. *


			** If you get a segmentation fault:

				1) If you DO NOT have temperature data, make sure Bool_t TempInt = 0.
				2) If you DO have temperature data, make sure Temp_Output is in the current directory.
 				3) There may have been an issue with Run_Fitter.sh. In the past, it deleted lines of bad fits, but would still append the timestamp to the file. 
				This should not happen anymore since I fixed the issue, but if there are 2 timestamps at the end of a line, this will cause a segmentation fault. 
			   



P.S. I added the new doit_with_date.sh script to this directory. This version appends the correct timestamp to the beginning of .root data files.
