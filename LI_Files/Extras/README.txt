 "Extra" files directory

This directory contains files which have failed temperature correction attempts, and also evidence of a time asynchronization issue centered around the old "doit_with_date.sh" bash script.
I included these files as an attempt to provide extra information/groundwork in case you decide to rethink the temperature correction. Also to draw attention to the necessity of correct timestamps.

1) Extra temperature corrections

	The file "LI_Analysis_with_NormCorr.C" is an almost exact copy of "LI_Analysis.C" except it has extra lines for a normalized temperature correction method. 
	This method attempted to normalize all the mean peak and temp vectors and subtract out the trends found in the temperature data. 

	The file "AvgTemp.C" is a file which attempts to use statistical methodology to correct for temperature. 
	Specifically, it treats the source/leds vectors ~ temperature vectors as bivariate normal distributions, and attempts to use a conditional expectation E[Y|X] formula for temp correction.
	This formula is supposed to account for the correlation between the two vectors. 
	It turned out that even if this method is valid, I did not have all the information necessary to completely carry this correction out as intended (No temperature std dev per sample).
	Then I realized if I averaged temperature values together, then I would have all the information I needed for the formula (though I of course was skeptical of its validity).
	It turned out, this averaging of temperatures was a fundamentally flawed method from the start. 
	I am glad I was finally able to have someone actually investigate this method I was concerned was fundamentally flawed, thanks Ramon.
	Though I still am not convinced this method (without averaging a bunch of things of course) would not work given I had all the necessary information.


2) The "Not_time_synched_Parameters" directory 

	In this directory I have the "Parameters1.dat" file filled with data associated with unix timestamps dating ~2018. I got these from adding a specific # of secs to the old timestamps. 
	
	The parameters found in the "LI_Files" directory are associated with unix timestamps dating ~1972.
	The reason for the ~1972 timestamps is because the "doit_with_date.sh" was hardcoded to get timestamps which turned out to be Current date - 2016.  
	This gave me problems when I wanted to do temperature analysis as a function of days and time, since I was not 100% sure the days/times were approximately correct.
	However, after trying to correct all the times, I found that the data timestamps and the temperature timestamps seemed out of sync. 

	* I would suggest using the corrected version of "doit_with_date.sh" to get the correct starting times of the runs (assuming you are still using the camac_daq executable). * 
