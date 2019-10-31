#include "TPad.h"
#include <TROOT.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TTree.h>
#include <TMath.h>
#include <TSpectrum.h>


//*******//
// Flags //
//*******//

Bool_t Debug = 0; // Toggle outpt of debugging comments
Bool_t TempInt = 1; // Toggle integration of temperature readings (must be set to 0 if no temperature data available)
Bool_t DrawLines = 1; // Toggle drawing of various plot lines
Bool_t SavePlots = 0; // Toggle saving a png of plots
Bool_t DrawPlots = 1; // Toggle drawing of plots


//**************************//
// Variable initializations //
//**************************//

Double_t average(Double_t numbers[], int numb);

// Parameter file vectors
Double_t hmev_mean, hmev_mean_err, hmev_sigma, hmev_sigma_err, hmev_chi2, hmev_chi2_NDF,
mev_mean, mev_mean_err, mev_sigma, mev_sigma_err, mev_chi2, mev_chi2_NDF,
led_calo_mean, led_calo_mean_err, led_calo_sigma, led_calo_sigma_err, led_calo_chi2, led_calo_chi2_NDF,  // First calo led
am_mean, am_mean_err, am_sigma, am_sigma_err, am_chi2, am_chi2_NDF, // First am source
led_ref_mean, led_ref_mean_err, led_ref_sigma, led_ref_sigma_err, led_ref_chi2, led_ref_chi2_NDF, // First ref led
am2_mean, am2_mean_err, am2_sigma, am2_sigma_err, // Used if a 2nd reference OM is present
led2_ref_mean, led2_ref_mean_err, led2_ref_sigma, led2_ref_sigma_err, // Used if a 2nd reference OM is present
time_data;

// Temperature file variables
Double_t time_temp, T_pb, T_db, T_lab, T_lab2;

// Parameter file vectors
std::vector<double> HMeV_Mean, HMeV_Mean_Err, HMeV_Sigma, HMeV_Sigma_Err, HMeV_chi2, HMeV_chi2_NDF,
MeV_Mean, MeV_Mean_Err, MeV_Sigma, MeV_Sigma_Err, MeV_chi2, MeV_chi2_NDF,
LED_Calo_Mean, LED_Calo_Mean_Err, LED_Calo_Sigma, LED_Calo_Sigma_Err, LED_Calo_chi2, LED_Calo_chi2_NDF, // First calo led
Am_Mean, Am_Mean_Err, Am_Sigma, Am_Sigma_Err, Am_chi2, Am_chi2_NDF, // First am source
LED_Ref_Mean, LED_Ref_Mean_Err, LED_Ref_Sigma, LED_Ref_Sigma_Err, LED_Ref_chi2, LED_Ref_chi2_NDF, // First ref led
Am2_Mean, Am2_Mean_Err, Am2_Sigma, Am2_Sigma_Err, // Used if a 2nd reference OM is present
LED2_Ref_Mean, LED2_Ref_Mean_Err, LED2_Ref_Sigma, LED2_Ref_Sigma_Err, // Used if a 2nd reference OM is present
Time_Data;

// Temperature file vectors
std::vector<double> Time_Temp, T_PB, T_DB, T_LAB, T_LAB2;

// Final time series vectors
std::vector<double> Bi_vs_Time, Bi_vs_Time_Err, Bi_vs_Time_chi2_NDF,
Bi_Pred_vs_Time, Bi_Pred_vs_Time_Err,
Am_vs_Time, Am_vs_Time_Err, Am_vs_Time_chi2_NDF,
LED_Calo_vs_Time, LED_Calo_vs_Time_Err, LED_Calo_vs_Time_chi2_NDF,
LED_Ref_vs_Time, LED_Ref_vs_Time_Err, LED_Ref_vs_Time_chi2_NDF,
Ratio_vs_Time, Ratio_vs_Time_Err,

// Temperature time series vectors
T_PB_vs_Time, T_DB_vs_Time, T_Lab_vs_Time, T_Lab2_vs_Time,

// Various ratios time series vectors for troubleshooting
Bi_Ratio_vs_Time, Am_Ratio_vs_Time, LED_Calo_Ratio_vs_Time, LED_Ref_Ratio_vs_Time,
Bi_Ratio_vs_Time_Err, Am_Ratio_vs_Time_Err, LED_Calo_Ratio_vs_Time_Err, LED_Ref_Ratio_vs_Time_Err,
Sources_Ratio_vs_Time, Sources_Ratio_vs_Time_Err, LEDs_Ratio_vs_Time, LEDs_Ratio_vs_Time_Err, 
Bi_Calo_Ratio_vs_Time, Bi_Calo_Ratio_vs_Time_Err, Am_Ref_Ratio_vs_Time, Am_Ref_Ratio_vs_Time_Err, 
All_Ratio_vs_Time, All_Ratio_vs_Time_Err, 
MeV_HMeV_Ratio_vs_Time, MeV_HMeV_Ratio_vs_Time_Err;


//*********************//
//*********************//
//    MAIN FUNCTION    //     
//*********************//
//*********************//

int LI_Analysis_with_NormCorr()
{


  //**********************************//
  // Read all parameters into vectors //
  //**********************************//

  // Read in file of all fit parameters
  ifstream in;
  in.open("Parameters1.dat");
  while(in
    >> hmev_mean >> hmev_mean_err >> hmev_sigma >> hmev_sigma_err
    >> mev_mean >> mev_mean_err >> mev_sigma >> mev_sigma_err >> mev_chi2 >> mev_chi2_NDF // chi2 and chi2/NDF are from full Bi fit
    >> led_calo_mean >> led_calo_mean_err >> led_calo_sigma >> led_calo_sigma_err >> led_calo_chi2 >> led_calo_chi2_NDF
    >> am_mean >> am_mean_err >> am_sigma >> am_sigma_err >> am_chi2 >> am_chi2_NDF
    >> led_ref_mean >> led_ref_mean_err >> led_ref_sigma >> led_ref_sigma_err >> led_ref_chi2 >> led_ref_chi2_NDF
    >> time_data)
  {
    // Calo OM Peaks
    HMeV_Mean.push_back(hmev_mean);
    HMeV_Mean_Err.push_back(hmev_mean_err);
    HMeV_Sigma.push_back(hmev_sigma);
    HMeV_Sigma_Err.push_back(hmev_sigma_err);

    MeV_Mean.push_back(mev_mean);
    MeV_Mean_Err.push_back(mev_mean_err);
    MeV_Sigma.push_back(mev_sigma);
    MeV_Sigma_Err.push_back(mev_sigma_err);
    MeV_chi2.push_back(mev_chi2);
    MeV_chi2_NDF.push_back(mev_chi2_NDF);

    LED_Calo_Mean.push_back(led_calo_mean);
    LED_Calo_Mean_Err.push_back(led_calo_mean_err);
    LED_Calo_Sigma.push_back(led_calo_sigma);
    LED_Calo_Sigma_Err.push_back(led_calo_sigma_err);
    LED_Calo_chi2.push_back(led_calo_chi2);
    LED_Calo_chi2_NDF.push_back(led_calo_chi2_NDF);
    
    // Ref OM Peaks
    Am_Mean.push_back(am_mean);
    Am_Mean_Err.push_back(am_mean_err);
    Am_Sigma.push_back(am_sigma);
    Am_Sigma_Err.push_back(am_sigma_err);
    Am_chi2.push_back(am_chi2);
    Am_chi2_NDF.push_back(am_chi2_NDF);

    LED_Ref_Mean.push_back(led_ref_mean);
    LED_Ref_Mean_Err.push_back(led_ref_mean_err);
    LED_Ref_Sigma.push_back(led_ref_sigma);
    LED_Ref_Sigma_Err.push_back(led_ref_sigma_err);
    LED_Ref_chi2.push_back(led_ref_chi2);
    LED_Ref_chi2_NDF.push_back(led_ref_chi2_NDF);
    // Time stamp for run
    Time_Data.push_back(time_data);
  }

  in.close();
  const int nentries = Time_Data.size();
  printf("Parameters Vector Size = %d\n",nentries);

  // Read in file of temperature values (if it exists)
  if(TempInt){
    ifstream in2;
    in2.open("Temp_Output");
    while(in2 >> time_temp >> T_lab >> T_db >> T_pb >> T_lab2){
      Time_Temp.push_back(time_temp); 
      T_PB.push_back(T_pb);
      T_DB.push_back(T_db);
      T_LAB.push_back(T_lab);
      T_LAB2.push_back(T_lab2);
    }
    in2.close();
  }  

  // Calculate Avg chi2/NDF for Sources and LEDs
  Double_t MeV_chi2_NDF_Avg = TMath::Mean(MeV_chi2_NDF.begin(), MeV_chi2_NDF.end());
  Double_t Am_chi2_NDF_Avg = TMath::Mean(Am_chi2_NDF.begin(), Am_chi2_NDF.end());
  Double_t LED_Calo_chi2_NDF_Avg = TMath::Mean(LED_Calo_chi2_NDF.begin(), LED_Calo_chi2_NDF.end());
  Double_t LED_Ref_chi2_NDF_Avg = TMath::Mean(LED_Ref_chi2_NDF.begin(), LED_Ref_chi2_NDF.end());

  cout << "\nMeV_chi2_NDF_Avg: " << MeV_chi2_NDF_Avg << endl;
  cout << "Am_chi2_NDF_Avg: " << Am_chi2_NDF_Avg << endl;
  cout << "LED_Calo_chi2_NDF_Avg: " << LED_Calo_chi2_NDF_Avg << endl;
  cout << "LED_Ref_chi2_NDF_Avg: " << LED_Ref_chi2_NDF_Avg << "\n" << endl;
  
  const int nentries2 = Time_Temp.size();

  printf("Temperature Vector Size = %d\n",nentries2);
  
  
  //*********************************//
  // Perform prediction calculations //
  //*********************************//
  
  // Start values (i.e. calibration run values) 
  Double_t Bi_0 = MeV_Mean[0];
  Double_t Bi_0_err = MeV_Mean_Err[0];
  Double_t Bi2_0 = HMeV_Mean[0];
  Double_t Bi2_0_err = HMeV_Mean_Err[0];
  Double_t LED_Calo_0 = LED_Calo_Mean[0];
  Double_t LED_Calo_0_err = LED_Calo_Mean_Err[0];
  Double_t LED_Ref_0 = LED_Ref_Mean[0];
  Double_t LED_Ref_0_err = LED_Ref_Mean_Err[0];
  Double_t Am_0 = Am_Mean[0];
  Double_t Am_0_err = Am_Mean_Err[0];

  // Later values at time t
  Double_t Bi_t, Bi_t_err, Bi_t_chi2_NDF, LED_Calo_t, LED_Calo_t_err, LED_Calo_t_chi2_NDF, LED_Ref_t, LED_Ref_t_err, LED_Ref_t_chi2_NDF, Am_t, Am_t_err, Am_t_chi2_NDF, LED_Ref_Pred, LED_Ref_Pred_err, LED_Calo_Pred, LED_Calo_Pred_err, Bi_Pred, Bi_Pred_err, Ratio, Ratio_err;
  Double_t Bi_last = Bi_0;
  Double_t LED_Calo_last = LED_Calo_0;
  int counter1 = 0, counter2 = 0, rejections = 0;
  for(int i=0; i<nentries; i++){ 
    Bi_t           = MeV_Mean[i];
    Bi_t_err       = MeV_Mean_Err[i];
    Bi_t_chi2_NDF  = MeV_chi2_NDF[i];
    LED_Calo_t     = LED_Calo_Mean[i];
    LED_Calo_t_err = LED_Calo_Mean_Err[i];
    LED_Calo_t_chi2_NDF = LED_Calo_chi2_NDF[i];
    LED_Ref_t      = LED_Ref_Mean[i];
    LED_Ref_t_err  = LED_Ref_Mean_Err[i];
    LED_Ref_t_chi2_NDF = LED_Ref_chi2_NDF[i];
    Am_t           = Am_Mean[i];
    Am_t_err       = Am_Mean_Err[i];
    Am_t_chi2_NDF  = Am_chi2_NDF[i];
    
    // Skip runs with large Bi peak deviations unaccounted for in LED peak
    if((TMath::Abs(Bi_t - Bi_last) > 10*Bi_t_err) &&
    !(TMath::Abs(LED_Calo_t - LED_Calo_last) > 10*LED_Calo_t_err)){
      if(Debug){ 
	printf("Chris: Skipping Run %d (timestamp: %f): Bi mean value = %f while LED mean value = %f\n",i, Time_Data[i], Bi_t, LED_Calo_t);
      }
      rejections++;
      continue;
    }
   
    // Update comparison values 
    Bi_last = Bi_t; 
    LED_Calo_last = LED_Calo_t;

    // Debugging outputs
    if(Debug) {printf("i = %d: Time = %f Bi mean = %f LED mean = %f\n", i, Time_Data[i], Bi_t, LED_Calo_t);}    

    // Use Am to predict Ref LED (correct for drifts)
    LED_Ref_Pred = LED_Ref_0 * (Am_t/Am_0);
    LED_Ref_Pred_err = sqrt(pow(LED_Ref_0_err*Am_t/Am_0,2) 
			    + pow(LED_Ref_0*Am_t_err/Am_0,2) 
			    + pow(LED_Ref_0*Am_t*Am_0_err/(Am_0*Am_0),2));
   
    // Use Ref LED to predict Calo LED
    LED_Calo_Pred = LED_Calo_t * (LED_Ref_Pred/LED_Ref_t);
    LED_Calo_Pred_err = sqrt(pow(LED_Calo_t_err*LED_Ref_Pred/LED_Ref_t,2) 
			     + pow(LED_Calo_t*LED_Ref_Pred_err/LED_Ref_t,2) 
			     + pow(LED_Calo_t*LED_Ref_Pred*LED_Ref_t_err/(LED_Ref_t*LED_Ref_t),2));
    
    // Use Calo LED to predict Bi
    Bi_Pred = Bi_0 * (LED_Calo_Pred/LED_Calo_0);  
    Bi_Pred_err = sqrt(pow(Bi_0_err*LED_Calo_Pred/LED_Calo_0,2) 
		       + pow(Bi_0*LED_Calo_Pred_err/LED_Calo_0,2) 
		       + pow(Bi_0*LED_Calo_Pred*LED_Calo_0_err/(LED_Calo_0*LED_Calo_0),2));

    // Calculate the ratio of predicted to measured Bi
    Ratio = Bi_Pred/Bi_t;
    Ratio_err = sqrt(pow(Bi_Pred_err/Bi_t,2) 
		     + pow(Bi_Pred*Bi_t_err/(Bi_t*Bi_t),2));

    // Cross check with full calculation 
    Double_t test = (Bi_0*LED_Calo_t/LED_Calo_0*LED_Ref_0/LED_Ref_t*Am_t/Am_0);

    // Fill Time Series Vectors
    Bi_vs_Time.push_back(Bi_t);
    Bi_vs_Time_Err.push_back(Bi_t_err);
    Bi_vs_Time_chi2_NDF.push_back(Bi_t_chi2_NDF);
    LED_Calo_vs_Time.push_back(LED_Calo_t);
    LED_Calo_vs_Time_Err.push_back(LED_Calo_t_err);
    LED_Calo_vs_Time_chi2_NDF.push_back(LED_Calo_t_chi2_NDF);
    LED_Ref_vs_Time.push_back(LED_Ref_t);
    LED_Ref_vs_Time_Err.push_back(LED_Ref_t_err);
    LED_Ref_vs_Time_chi2_NDF.push_back(LED_Ref_t_chi2_NDF);
    Am_vs_Time.push_back(Am_t);
    Am_vs_Time_Err.push_back(Am_t_err);
    Am_vs_Time_chi2_NDF.push_back(Am_t_chi2_NDF);

    Bi_Pred_vs_Time.push_back(Bi_Pred);
    Bi_Pred_vs_Time_Err.push_back(Bi_Pred_err);
    Ratio_vs_Time.push_back(Ratio);
    Ratio_vs_Time_Err.push_back(Ratio_err);

    // Fill Ratio Time Series Vectors for troubleshooting
    Bi_Ratio_vs_Time.push_back(Bi_0/Bi_t);
    Bi_Ratio_vs_Time_Err.push_back(sqrt(pow(Bi_0_err/Bi_t,2) 
					+ pow(Bi_0*Bi_t_err/(Bi_t*Bi_t),2)));

    Am_Ratio_vs_Time.push_back(Am_t/Am_0);
    Am_Ratio_vs_Time_Err.push_back(sqrt(pow(Am_t_err/Am_0,2) 
					+ pow(Am_t*Am_0_err/(Am_0*Am_0),2)));

    LED_Calo_Ratio_vs_Time.push_back(LED_Calo_t/LED_Calo_0);
    LED_Calo_Ratio_vs_Time_Err.push_back(sqrt(pow(LED_Calo_t_err/LED_Calo_0,2) 
					      + pow(LED_Calo_t*LED_Calo_0_err/(LED_Calo_0*LED_Calo_0),2)));

    LED_Ref_Ratio_vs_Time.push_back(LED_Ref_0/LED_Ref_t);
    LED_Ref_Ratio_vs_Time_Err.push_back(sqrt(pow(LED_Ref_0_err/LED_Ref_t,2) 
					     + pow(LED_Ref_0*LED_Ref_t_err/(LED_Ref_t*LED_Ref_t),2)));
    
 
    //*********************************//
    // Assign temperatures to each run //
    //*********************************//

    // Search through temperature vectors for temp values at closest matching time values
    if(TempInt){
      int j=0;
      
      // Use if you want the avg temperature for time elapsed between data points 
      /*Double_t T_PB_vs_Time_avg = 0;
      Double_t T_DB_vs_Time_avg = 0;
      Double_t T_Lab_vs_Time_avg = 0;
      Double_t T_Lab2_vs_Time_avg = 0;*/

      while(j < nentries2){
	if(TMath::Abs(Time_Data[i] - Time_Temp[j]) <= 600){ // Temperature readings taken every 600 seconds
	  if(Debug) {printf("Temp at closest time found:   index = %d   |   time = %f   |   time diff = %f mins   |   temp = %f\n",i,Time_Temp[j],(Time_Data[i]-Time_Temp[j])/60,T_PB[j]);} // Print statement to check how close temp and parameter times are
	  T_PB_vs_Time.push_back(T_PB[j]);
	  T_DB_vs_Time.push_back(T_DB[j]);
	  T_Lab_vs_Time.push_back(T_LAB[j]);
	  T_Lab2_vs_Time.push_back(T_LAB2[j]);
         
          // Use if you want the avg temperature for time elapsed between data points 
          /*T_PB_vs_Time_avg = (T_PB[j]+T_PB[j-1]+T_PB[j-2]+T_PB[j-3])/4;
          T_DB_vs_Time_avg = (T_DB[j]+T_DB[j-1]+T_DB[j-2]+T_DB[j-3])/4;
          T_Lab_vs_Time_avg = (T_LAB[j]+T_LAB[j-1]+T_LAB[j-2]+T_LAB[j-3])/4;
          T_Lab2_vs_Time_avg = (T_LAB2[j]+T_LAB2[j-1]+T_LAB2[j-2]+T_LAB2[j-3])/4;
     
          T_DB_sigma.push_back(sqrt((pow(T_DB[j]-T_DB_vs_Time_avg,2)
				     + pow(T_DB[j-1]-T_DB_vs_Time_avg,2)
				     + pow(T_DB[j-2]-T_DB_vs_Time_avg,2)
				     + pow(T_DB[j-3]-T_DB_vs_Time_avg,2))/4));
          T_Lab2_sigma.push_back(sqrt((pow(T_LAB2[j]-T_Lab2_vs_Time_avg,2)
				       + pow(T_LAB2[j-1]-T_Lab2_vs_Time_avg,2)
				       + pow(T_LAB2[j-2]-T_Lab2_vs_Time_avg,2)
				       + pow(T_LAB2[j-3]-T_Lab2_vs_Time_avg,2))/4));     

          T_PB_vs_Time.push_back(T_PB_vs_Time_avg);
          T_DB_vs_Time.push_back(T_DB_vs_Time_avg);
          T_Lab_vs_Time.push_back(T_Lab_vs_Time_avg);
          T_Lab2_vs_Time.push_back(T_Lab2_vs_Time_avg);
     
          cout << " T_DB_avg: " << T_DB_vs_Time_avg << endl;
     
          T_PB_vs_Time_avg = 0; 
          T_DB_vs_Time_avg = 0; 
          T_Lab_vs_Time_avg = 0; 
          T_Lab2_vs_Time_avg = 0;*/

	  counter2++;
	  break;
	  }
	else{
	  j++;
	}
      }
   
    }
    counter1++;  
  }

  const int size = Ratio_vs_Time.size();

  // Find min and max for DB and Lab Temperature
  Double_t T_DB_Max, T_DB_Min, T_Lab2_Max, T_Lab2_Min;
  if(TempInt){
    T_DB_Max = *std::max_element(T_DB_vs_Time.begin(),T_DB_vs_Time.end());
    T_DB_Min = *std::min_element(T_DB_vs_Time.begin(),T_DB_vs_Time.end());
    T_Lab2_Max = *std::max_element(T_Lab2_vs_Time.begin(),T_Lab2_vs_Time.end());
    T_Lab2_Min = *std::min_element(T_Lab2_vs_Time.begin(),T_Lab2_vs_Time.end());
  
    cout << "\nThe Max DB Temperature is: " << T_DB_Max << endl;
    cout << "The Min DB Temperature is: " << T_DB_Min << endl;
    cout << "The Max Lab Temperature is: " << T_Lab2_Max << endl;
    cout << "The Min Lab Temperature is: " << T_Lab2_Min << "\n" << endl;
  
    if(Debug){printf("The vector size is T_PB_vs_Time: %lu and T_DB_vs_Time: %lu and T_Lab_vs_Time: %lu and T_Lab2_vs_Time: %lu\n", T_PB_vs_Time.size(), T_DB_vs_Time.size(), T_Lab_vs_Time.size(), T_Lab2_vs_Time.size());}
  
  // Sanity Check for Debugging
    if(Debug){
      printf("Final Parameter Vector Size = %d\n",size);
      printf("Final Temperature Vector Size = %lu\n",T_DB_vs_Time.size());
      printf("counter 1 = %d  |  counter 2 = %d\n",counter1, counter2);
    }
  }
  
  // Find min and max for chi2/NDF to set histogram ranges
  Double_t Bi_vs_Time_chi2_NDF_Max = *std::max_element(Bi_vs_Time_chi2_NDF.begin(),Bi_vs_Time_chi2_NDF.end());
  Double_t Bi_vs_Time_chi2_NDF_Min = *std::min_element(Bi_vs_Time_chi2_NDF.begin(),Bi_vs_Time_chi2_NDF.end());
  Double_t LED_Calo_vs_Time_chi2_NDF_Max = *std::max_element(LED_Calo_vs_Time_chi2_NDF.begin(),LED_Calo_vs_Time_chi2_NDF.end());
  Double_t LED_Calo_vs_Time_chi2_NDF_Min = *std::min_element(LED_Calo_vs_Time_chi2_NDF.begin(),LED_Calo_vs_Time_chi2_NDF.end());
  Double_t Am_vs_Time_chi2_NDF_Max = *std::max_element(Am_vs_Time_chi2_NDF.begin(),Am_vs_Time_chi2_NDF.end());
  Double_t Am_vs_Time_chi2_NDF_Min = *std::min_element(Am_vs_Time_chi2_NDF.begin(),Am_vs_Time_chi2_NDF.end());
  Double_t LED_Ref_vs_Time_chi2_NDF_Max = *std::max_element(LED_Ref_vs_Time_chi2_NDF.begin(),LED_Ref_vs_Time_chi2_NDF.end());
  Double_t LED_Ref_vs_Time_chi2_NDF_Min = *std::min_element(LED_Ref_vs_Time_chi2_NDF.begin(),LED_Ref_vs_Time_chi2_NDF.end());
  
  // Define histograms to characterize Chi^2/NDF spreads
  TH1* Bi_chi2_NDF_Spread = new TH1D("Bi Chi^2/NDF Spread", "Bi Chi^2/NDF Spread", 300, Bi_vs_Time_chi2_NDF_Min, Bi_vs_Time_chi2_NDF_Max); 
  TH1* LED_Calo_chi2_NDF_Spread = new TH1D("LED_Calo Chi^2/NDF Spread", "LED_Calo Chi^2/NDF Spread", 300, LED_Calo_vs_Time_chi2_NDF_Min, LED_Calo_vs_Time_chi2_NDF_Max);
  TH1* Am_chi2_NDF_Spread = new TH1D("Am Chi^2/NDF Spread", "Am Chi^2/NDF Spread", 300, Am_vs_Time_chi2_NDF_Min, Am_vs_Time_chi2_NDF_Max);
  TH1* LED_Ref_chi2_NDF_Spread = new TH1D("LED_Ref Chi^2/NDF Spread", "LED_Ref Chi^2/NDF Spread", 300, LED_Ref_vs_Time_chi2_NDF_Min, LED_Ref_vs_Time_chi2_NDF_Max);

  for(int i=0; i<size; i++){
    Bi_chi2_NDF_Spread->Fill(Bi_vs_Time_chi2_NDF[i]);
    LED_Calo_chi2_NDF_Spread->Fill(LED_Calo_vs_Time_chi2_NDF[i]);
    Am_chi2_NDF_Spread->Fill(Am_vs_Time_chi2_NDF[i]); 
    LED_Ref_chi2_NDF_Spread->Fill(LED_Ref_vs_Time_chi2_NDF[i]);
  }

  //********************************************************************************************//
  // Comparison of Source and LED Signals vs Time (to investigate Bi_Pred:Bi Ratio deliniation) //
  //********************************************************************************************//

  for(int i=0; i<size; i++){

    // Fill Sources and LEDs Ratio Vectors  
    Sources_Ratio_vs_Time.push_back(Bi_Ratio_vs_Time[i]*Am_Ratio_vs_Time[i]);
    Sources_Ratio_vs_Time_Err.push_back(sqrt(pow(Bi_Ratio_vs_Time_Err[i]*Am_Ratio_vs_Time[i],2) 
				   	     + pow(Bi_Ratio_vs_Time[i]*Am_Ratio_vs_Time_Err[i],2)));

    LEDs_Ratio_vs_Time.push_back(LED_Calo_Ratio_vs_Time[i]*LED_Ref_Ratio_vs_Time[i]);
    LEDs_Ratio_vs_Time_Err.push_back(sqrt(pow(LED_Calo_Ratio_vs_Time_Err[i]*LED_Ref_Ratio_vs_Time[i],2) 
					  + pow(LED_Calo_Ratio_vs_Time[i]*LED_Ref_Ratio_vs_Time_Err[i],2)));

    // Fill MeV:HMeV Ratio Vectors
    MeV_HMeV_Ratio_vs_Time.push_back(MeV_Mean[i]/HMeV_Mean[i]);
    MeV_HMeV_Ratio_vs_Time_Err.push_back(sqrt(pow(MeV_Mean_Err[i]/HMeV_Mean[i],2) 
					      + pow(MeV_Mean[i]*HMeV_Mean_Err[i]/(HMeV_Mean[i]*HMeV_Mean[i]),2))); 

    // Fill Bi:Calo_LED and Am:Ref_LED Ratio Vectors	
    Bi_Calo_Ratio_vs_Time.push_back(Bi_Ratio_vs_Time[i]*LED_Calo_Ratio_vs_Time[i]);
    Bi_Calo_Ratio_vs_Time_Err.push_back(sqrt(pow(Bi_Ratio_vs_Time_Err[i]*LED_Calo_Ratio_vs_Time[i],2) 
					     + pow(Bi_Ratio_vs_Time[i]*LED_Calo_Ratio_vs_Time_Err[i],2)));
    Am_Ref_Ratio_vs_Time.push_back(Am_Ratio_vs_Time[i]*LED_Ref_Ratio_vs_Time[i]);
    Am_Ref_Ratio_vs_Time_Err.push_back(sqrt(pow(Am_Ratio_vs_Time_Err[i]*LED_Ref_Ratio_vs_Time[i],2) 
					    + pow(Am_Ratio_vs_Time[i]*LED_Ref_Ratio_vs_Time_Err[i],2)));
  } 

    // Fill Full Ratio Vectors
  for(int i=0; i<size; i++){
    All_Ratio_vs_Time.push_back(Sources_Ratio_vs_Time[i]/LEDs_Ratio_vs_Time[i]);
    All_Ratio_vs_Time_Err.push_back(sqrt(pow(Sources_Ratio_vs_Time_Err[i]/LEDs_Ratio_vs_Time[i],2) 
					 + pow(Sources_Ratio_vs_Time[i]*LEDs_Ratio_vs_Time_Err[i]/(LEDs_Ratio_vs_Time[i]*LEDs_Ratio_vs_Time[i]),2)));
  }


  //*************************************//
  // Build final arrays used for TGraphs //
  //*************************************//
   
  // Fill Time, Sources, LEDs, and Ratio Arrays
  Double_t Time[size], Time_Err[size],
  Bi_vs_Time_Array[size], Bi_vs_Time_Err_Array[size], Bi_vs_Time_chi2_NDF_Array[size],
  LED_Calo_vs_Time_Array[size], LED_Calo_vs_Time_Err_Array[size], LED_Calo_vs_Time_chi2_NDF_Array[size],
  Am_vs_Time_Array[size], Am_vs_Time_Err_Array[size], Am_vs_Time_chi2_NDF_Array[size],
  LED_Ref_vs_Time_Array[size], LED_Ref_vs_Time_Err_Array[size], LED_Ref_vs_Time_chi2_NDF_Array[size],
  Bi_Pred_vs_Time_Array[size], Bi_Pred_vs_Time_Err_Array[size],
  Ratio_vs_Time_Array[size], Ratio_vs_Time_Err_Array[size],
    
  // Ratio Arrays for Testing 
  Bi_Ratio_vs_Time_Array[size], Bi_Ratio_vs_Time_Err_Array[size],
  Am_Ratio_vs_Time_Array[size], Am_Ratio_vs_Time_Err_Array[size],
  LED_Ref_Ratio_vs_Time_Array[size], LED_Ref_Ratio_vs_Time_Err_Array[size],
  LED_Calo_Ratio_vs_Time_Array[size], LED_Calo_Ratio_vs_Time_Err_Array[size],
  Sources_Ratio_vs_Time_Array[size], Sources_Ratio_vs_Time_Err_Array[size],
  LEDs_Ratio_vs_Time_Array[size], LEDs_Ratio_vs_Time_Err_Array[size],
  Bi_Calo_Ratio_vs_Time_Array[size], Bi_Calo_Ratio_vs_Time_Err_Array[size],
  Am_Ref_Ratio_vs_Time_Array[size], Am_Ref_Ratio_vs_Time_Err_Array[size],
  All_Ratio_vs_Time_Array[size], All_Ratio_vs_Time_Err_Array[size], MeV_HMeV_Ratio_vs_Time_Array[size], MeV_HMeV_Ratio_vs_Time_Err_Array[size],
   
  // Temperature correlation arrays (all use the same error array with constant value 0.1 deg F)
  T_PB_vs_Time_Array[size], T_DB_vs_Time_Array[size],
  T_Lab_vs_Time_Array[size], T_Lab2_vs_Time_Array[size],
  T_Err_vs_Time_Array[size];
 
  // Loop to fill arrays from vectors
  for(int i=0;i<size;i++){

    // Sources and LEDs Arrays
    Bi_vs_Time_Array[i] = Bi_vs_Time[i];
    Bi_vs_Time_Err_Array[i] = Bi_vs_Time_Err[i];
    Bi_vs_Time_chi2_NDF_Array[i] = Bi_vs_Time_chi2_NDF[i];
    LED_Calo_vs_Time_Array[i] = LED_Calo_vs_Time[i];
    LED_Calo_vs_Time_Err_Array[i] = LED_Calo_vs_Time_Err[i];
    LED_Calo_vs_Time_chi2_NDF_Array[i] = LED_Calo_vs_Time_chi2_NDF[i];
    Am_vs_Time_Array[i] = Am_vs_Time[i];
    Am_vs_Time_Err_Array[i] = Am_vs_Time_Err[i];
    Am_vs_Time_chi2_NDF_Array[i] = Am_vs_Time_chi2_NDF[i];
    LED_Ref_vs_Time_Array[i] = LED_Ref_vs_Time[i];
    LED_Ref_vs_Time_Err_Array[i] = LED_Ref_vs_Time_Err[i];
    LED_Ref_vs_Time_chi2_NDF_Array[i] = LED_Ref_vs_Time_chi2_NDF[i];
    
    Time[i] = (Time_Data[i]-Time_Data[0])/86400.; // Set time as days since the first run
    Time_Err[i] = 0.5*30./1440.;  // Error = 1/2 * run duration in days (~30 min run = 1/48 days)
    
    // Money Plot Arrays
    Bi_Pred_vs_Time_Array[i] = Bi_Pred_vs_Time[i];
    Bi_Pred_vs_Time_Err_Array[i] = Bi_Pred_vs_Time_Err[i];
    Ratio_vs_Time_Array[i] = Ratio_vs_Time[i];
    Ratio_vs_Time_Err_Array[i] = Ratio_vs_Time_Err[i];

    // Arrays for testing ratios to investigate Bi_Pred:Bi Ratio deliniation
    Bi_Ratio_vs_Time_Array[i] = Bi_Ratio_vs_Time[i];
    Bi_Ratio_vs_Time_Err_Array[i] = Bi_Ratio_vs_Time_Err[i];
    Am_Ratio_vs_Time_Array[i] = Am_Ratio_vs_Time[i];
    Am_Ratio_vs_Time_Err_Array[i] = Am_Ratio_vs_Time_Err[i];
    LED_Calo_Ratio_vs_Time_Array[i] = LED_Calo_Ratio_vs_Time[i];
    LED_Calo_Ratio_vs_Time_Err_Array[i] = LED_Calo_Ratio_vs_Time_Err[i];
    LED_Ref_Ratio_vs_Time_Array[i] = LED_Ref_Ratio_vs_Time[i];
    LED_Ref_Ratio_vs_Time_Err_Array[i] = LED_Ref_Ratio_vs_Time_Err[i];

    Sources_Ratio_vs_Time_Array[i] = Sources_Ratio_vs_Time[i];
    Sources_Ratio_vs_Time_Err_Array[i] = Sources_Ratio_vs_Time_Err[i];
    LEDs_Ratio_vs_Time_Array[i] = LEDs_Ratio_vs_Time[i];
    LEDs_Ratio_vs_Time_Err_Array[i] = LEDs_Ratio_vs_Time_Err[i];
    All_Ratio_vs_Time_Array[i] = All_Ratio_vs_Time[i];
    All_Ratio_vs_Time_Err_Array[i] = All_Ratio_vs_Time_Err[i];
    Bi_Calo_Ratio_vs_Time_Array[i] = Bi_Calo_Ratio_vs_Time[i];
    Bi_Calo_Ratio_vs_Time_Err_Array[i] = Bi_Calo_Ratio_vs_Time_Err[i];
    Am_Ref_Ratio_vs_Time_Array[i] = Am_Ref_Ratio_vs_Time[i];
    Am_Ref_Ratio_vs_Time_Err_Array[i] = Am_Ref_Ratio_vs_Time_Err[i];
    
    // Fill MeV:HMeV Ratio Arrays to check consistancy of Bi peaks
    MeV_HMeV_Ratio_vs_Time_Array[i] = MeV_HMeV_Ratio_vs_Time[i];
    MeV_HMeV_Ratio_vs_Time_Err_Array[i] = MeV_HMeV_Ratio_vs_Time_Err[i];
  
    // Fill Temperature Arrays
    if(TempInt){    
      T_PB_vs_Time_Array[i] = T_PB_vs_Time[i];
      T_DB_vs_Time_Array[i] = T_DB_vs_Time[i];
      T_Lab_vs_Time_Array[i] = T_Lab_vs_Time[i];
      T_Lab2_vs_Time_Array[i] = T_Lab2_vs_Time[i];
      T_Err_vs_Time_Array[i] = 0.1; // Constant temperature logger error (0.1 deg F)
    }
  }


  //************************************************************//
  // Normalization Method for Temp Corr (Does Not Seem to Work) //
  //************************************************************//
    
  std::vector<double> RatioTemp_Lab_vs_Time, NormLab, NormRatio, RecNormTemp_Lab, Norm2Lab, NormDB, RecNormTemp_DB, Norm2DB, RecCorrRatio, FinalTemp_Lab, FinalTemp_DB, NormBi, NormAm, NormCalo, NormRef, Corr_Bi, Corr_Am, Corr_Calo, Corr_Ref, CorrRatio_Lab, PreCorrRatio_DB, CorrRatio_DB, CorrRatio_DB_Err, NormBi_Err, NormAm_Err, NormCalo_Err, NormRef_Err, Ratio_Line_Diff, Ratio_Line_Diff_Err;
  Double_t RatioTemp_Lab_vs_Time_Array[size], RecNormTemp_Lab_vs_Time_Array[size], RecNormTemp_DB_vs_Time_Array[size], NormRatio_vs_Time_Array[size], RatioNormTemp_Max, RecNormTemp_Lab_Max, RecNormTemp_DB_Max, CorrRatio_Lab_vs_Time_Array[size], CorrRatio_DB_vs_Time_Array[size], CorrRatio_DB_vs_Time_Err_Array[size], CorrRatio_vs_Time_Array[size], Norm_Bi_Array[size], Norm_Am_Array[size], Norm_LED_Calo_Array[size], Norm_LED_Ref_Array[size], Norm_Bi_Err_Array[size], Norm_Am_Err_Array[size], Norm_LED_Calo_Err_Array[size], Norm_LED_Ref_Err_Array[size], Ratio_Line_Diff_Array[size], Ratio_Line_Diff_Err_Array[size];
    
  /*Double_t Ratio_Max = *std::max_element(Ratio_vs_Time.begin(),Ratio_vs_Time.end());
  Double_t Bi_Max = *std::max_element(Bi_vs_Time.begin(),Bi_vs_Time.end());
  Double_t Am_Max = *std::max_element(Am_vs_Time.begin(),Am_vs_Time.end());
  Double_t Calo_Max = *std::max_element(LED_Calo_vs_Time.begin(),LED_Calo_vs_Time.end());
  Double_t Ref_Max = *std::max_element(LED_Ref_vs_Time.begin(),LED_Ref_vs_Time.end());
  Double_t Bi_Max_Err = *std::max_element(Bi_vs_Time_Err.begin(),Bi_vs_Time_Err.end());
  Double_t Am_Max_Err = *std::max_element(Am_vs_Time_Err.begin(),Am_vs_Time_Err.end()); 
  Double_t Calo_Max_Err = *std::max_element(LED_Calo_vs_Time_Err.begin(),LED_Calo_vs_Time_Err.end());
  Double_t Ref_Max_Err = *std::max_element(LED_Ref_vs_Time_Err.begin(),LED_Ref_vs_Time_Err.end());  

  for(int i=0; i<size; i++){
    NormLab.push_back(T_Lab2_vs_Time[i]/T_Lab2_Max);
    NormDB.push_back(T_DB_vs_Time[i]/T_DB_Max);
    NormRatio.push_back(Ratio_vs_Time[i]/Ratio_Max);
    NormBi.push_back(Bi_vs_Time[i]/Bi_Max);
    NormAm.push_back(Am_vs_Time[i]/Am_Max);
    NormCalo.push_back(LED_Calo_vs_Time[i]/Calo_Max);
    NormRef.push_back(LED_Ref_vs_Time[i]/Ref_Max);
    NormBi_Err.push_back(sqrt(pow(Bi_vs_Time_Err[i]/Bi_Max,2)+pow(Bi_vs_Time[i]*Bi_Max_Err/(Bi_Max*Bi_Max),2)));
    NormAm_Err.push_back(sqrt(pow(Am_vs_Time_Err[i]/Am_Max,2)+pow(Am_vs_Time[i]*Am_Max_Err/(Am_Max*Am_Max),2)));
    NormCalo_Err.push_back(sqrt(pow(LED_Calo_vs_Time_Err[i]/Calo_Max,2)+pow(LED_Calo_vs_Time[i]*Calo_Max_Err/(Calo_Max*Calo_Max),2)));
    NormRef_Err.push_back(sqrt(pow(LED_Ref_vs_Time_Err[i]/Ref_Max,2)+pow(LED_Ref_vs_Time[i]*Ref_Max_Err/(Ref_Max*Ref_Max),2)));
  }
	
  for(int i=0; i<size; i++){
    RecNormTemp_Lab.push_back((-1)/NormLab[i]);
    RecNormTemp_DB.push_back((-1)/NormDB[i]);  
  }

  RecNormTemp_Lab_Max = *std::max_element(RecNormTemp_Lab.begin(),RecNormTemp_Lab.end());
  RecNormTemp_DB_Max = *std::max_element(RecNormTemp_DB.begin(),RecNormTemp_DB.end());
	
  for(int i=0; i<size; i++){
    Norm2Lab.push_back(RecNormTemp_Lab[i]/RecNormTemp_Lab_Max);
    Norm2DB.push_back(RecNormTemp_DB[i]/RecNormTemp_DB_Max);
  }

  Double_t Norm2Lab_Max = *std::max_element(Norm2Lab.begin(),Norm2Lab.end());
  Double_t Norm2DB_Max = *std::max_element(Norm2DB.begin(),Norm2DB.end());
	
  for(int i=0; i<size; i++){
    FinalTemp_Lab.push_back(Norm2Lab[i]/Norm2Lab_Max);
    FinalTemp_DB.push_back(Norm2DB[i]/Norm2DB_Max);
    //cout << FinalTemp_Lab[i] << endl;
  }
	
  for(int i=0; i<size; i++){
    RatioTemp_Lab_vs_Time.push_back(FinalTemp_Lab[i]/NormRatio[i]);
  }

  for(int i=0; i<size; i++){
    CorrRatio_Lab.push_back((NormRatio[i]+((FinalTemp_Lab[i]-FinalTemp_Lab[1]))));
    //CorrRatio_DB.push_back((NormRatio[i]+((FinalTemp_DB[i]-FinalTemp_DB[0]))));
    Corr_Bi.push_back((NormBi[i]+((FinalTemp_DB[i]-FinalTemp_DB[0]))));
    Corr_Am.push_back((NormAm[i]+((FinalTemp_DB[i]-FinalTemp_DB[0]))));
    Corr_Calo.push_back((NormCalo[i]+((FinalTemp_Lab[i]-FinalTemp_Lab[0]))));
    Corr_Ref.push_back((NormRef[i]+((FinalTemp_Lab[i]-FinalTemp_Lab[0]))));
  }
  
  for(int i=0; i<size; i++){
    Norm_Bi_Array[i] = NormBi[i];
    Norm_Am_Array[i] = NormAm[i];
    Norm_LED_Calo_Array[i] = NormCalo[i];
    Norm_LED_Ref_Array[i] = NormRef[i];
    Norm_Bi_Err_Array[i] = NormBi_Err[i];
    Norm_Am_Err_Array[i] = NormAm_Err[i];
    Norm_LED_Calo_Err_Array[i] = NormCalo_Err[i];
    Norm_LED_Ref_Err_Array[i] = NormRef_Err[i];
  }
  */	
  
		
  //*******************//
  // Ratio vs Temp Fit //
  //*******************//
  
  Double_t Ratio_Temp_r, Ratio_Temp_yMean, Ratio_Temp_slope, Ratio_Temp_slope_err, Ratio_Temp_const, Ratio_Temp_const_err;
  TGraphErrors *Ratio_vs_Temp_Graph;
  if(TempInt){
    Ratio_vs_Temp_Graph = new TGraphErrors(size, T_DB_vs_Time_Array, Ratio_vs_Time_Array, T_Err_vs_Time_Array,Ratio_vs_Time_Err_Array);
    TF1 *Ratio_Temp = new TF1("Ratio_Temp","[0]*x+[1]",T_DB_Min,T_DB_Max); 
    Ratio_Temp->SetParameter(0,-0.33); 
    //Ratio_Temp->SetParLimits(0, -0.37, -0.25);
    Ratio_Temp->SetParameter(1,1.01);
    Ratio_Temp->SetLineColor(kRed);
    Ratio_Temp->SetLineWidth(4);
    Ratio_vs_Temp_Graph->Fit("Ratio_Temp", "RSB");
    Ratio_Temp_r = Ratio_vs_Temp_Graph->GetCorrelationFactor();
    Ratio_Temp_yMean = Ratio_vs_Temp_Graph->GetMean(2);
    Ratio_Temp_slope = Ratio_Temp->GetParameter(0);
    Ratio_Temp_slope_err = Ratio_Temp->GetParError(0);
    Ratio_Temp_const = Ratio_Temp->GetParameter(1);
    Ratio_Temp_const_err = Ratio_Temp->GetParError(1);
  }

  //***********************************************************//
  // Normalization Temp Corr Method Arrays and Temp Correction //
  //***********************************************************//
 
  /*for(int i=0; i<size; i++){
    RatioTemp_Lab_vs_Time_Array[i] = RatioTemp_Lab_vs_Time[i];
    //RecNormTemp_Lab_vs_Time_Array[i] = FinalTemp_Lab[i];
    //NormRatio_vs_Time_Array[i] = NormRatio[i];
    //CorrRatio_Lab_vs_Time_Array[i] = CorrRatio_Lab[i];	
  }
	
  for(int i=0; i<size; i++){ 
    CorrRatio_DB.push_back(Ratio_vs_Time_Array[i]-((Ratio_Temp_slope*20*(T_DB_vs_Time_Array[i]-T_DB_vs_Time_Array[0]))+0.01*Temp_DB_Line_Diff[i]));
    CorrRatio_DB_Err.push_back(sqrt(pow(Ratio_vs_Time_Err_Array[i]-(Ratio_Temp_slope*20*(T_DB_vs_Time_Array[i]-T_DB_vs_Time_Array[0])),2)
				    + pow(Ratio_vs_Time_Array[i]-(Ratio_Temp_slope_err*20*(T_DB_vs_Time_Array[i]-T_DB_vs_Time_Array[0])),2)
				    + pow(Ratio_vs_Time_Array[i]-(Ratio_Temp_slope*20*(T_Err_vs_Time_Array[i]-T_DB_vs_Time_Array[0])),2)
				    + pow(Ratio_vs_Time_Array[i]-(Ratio_Temp_slope*20*(T_DB_vs_Time_Array[i]-T_Err_vs_Time_Array[0])),2)));
    //CorrRatio_DB.push_back(Ratio_vs_Time_Array[i]+(0.003*(T_DB_vs_Time_Array[i]-T_DB_vs_Time_Array[0])));	
  }    

  for(int i=0; i<size; i++){
    CorrRatio_DB_vs_Time_Array[i] = CorrRatio_DB[i];
    CorrRatio_DB_vs_Time_Err_Array[i] = CorrRatio_DB_Err[i];
  }*/


  //*******************************//
  // Creating TGraphs for Plotting //
  //****************************** //
  
  // Basic time evolution of different peaks
  TGraphErrors *Bi_vs_Time_Graph = new TGraphErrors(size,Time,Bi_vs_Time_Array,Time_Err,Bi_vs_Time_Err_Array);
  TGraphErrors *LED_Calo_vs_Time_Graph = new TGraphErrors(size,Time,LED_Calo_vs_Time_Array,Time_Err,LED_Calo_vs_Time_Err_Array);
  TGraphErrors *LED_Ref_vs_Time_Graph = new TGraphErrors(size,Time,LED_Ref_vs_Time_Array,Time_Err,LED_Ref_vs_Time_Err_Array);
  TGraphErrors *Am_vs_Time_Graph = new TGraphErrors(size,Time,Am_vs_Time_Array,Time_Err,Am_vs_Time_Err_Array);

  // Time evolution for the chi2/NDF of different fits from Fitter.C
  TGraph *Bi_vs_Time_chi2_NDF_Graph = new TGraph(size,Time,Bi_vs_Time_chi2_NDF_Array);
  TGraph *LED_Calo_vs_Time_chi2_NDF_Graph = new TGraph(size,Time,LED_Calo_vs_Time_chi2_NDF_Array);
  TGraph *LED_Ref_vs_Time_chi2_NDF_Graph = new TGraph(size,Time,LED_Ref_vs_Time_chi2_NDF_Array);
  TGraph *Am_vs_Time_chi2_NDF_Graph = new TGraph(size,Time,Am_vs_Time_chi2_NDF_Array); 
  
  // Money plots: time evolution of prediction and of ratio between predicted and measured values
  TGraphErrors *Bi_Pred_vs_Time_Graph = new TGraphErrors(size,Time,Bi_Pred_vs_Time_Array,Time_Err,Bi_Pred_vs_Time_Err_Array);
  TGraphErrors *Ratio_vs_Time_Graph = new TGraphErrors(size, Time, Ratio_vs_Time_Array,Time_Err,Ratio_vs_Time_Err_Array);

  // Check Individual Ratios vs Time
  TGraphErrors *Bi_Ratio_vs_Time_Graph = new TGraphErrors(size,Time,Bi_Ratio_vs_Time_Array,Time_Err,Bi_Ratio_vs_Time_Err_Array);
  TGraphErrors *Am_Ratio_vs_Time_Graph = new TGraphErrors(size,Time,Am_Ratio_vs_Time_Array,Time_Err,Am_Ratio_vs_Time_Err_Array);
  TGraphErrors *LED_Calo_Ratio_vs_Time_Graph = new TGraphErrors(size,Time,LED_Calo_Ratio_vs_Time_Array,Time_Err,LED_Calo_Ratio_vs_Time_Err_Array);
  TGraphErrors *LED_Ref_Ratio_vs_Time_Graph = new TGraphErrors(size,Time,LED_Ref_Ratio_vs_Time_Array,Time_Err,LED_Ref_Ratio_vs_Time_Err_Array);
  
  // Check Various Ratios
  TGraphErrors *Sources_Ratio_vs_Time_Graph = new TGraphErrors(size,Time,Sources_Ratio_vs_Time_Array,Time_Err,Sources_Ratio_vs_Time_Err_Array);
  TGraphErrors *LEDs_Ratio_vs_Time_Graph = new TGraphErrors(size,Time,LEDs_Ratio_vs_Time_Array,Time_Err,LEDs_Ratio_vs_Time_Err_Array);
  TGraphErrors *Bi_Calo_Ratio_vs_Time_Graph = new TGraphErrors(size,Time,Bi_Calo_Ratio_vs_Time_Array,Time_Err,Bi_Calo_Ratio_vs_Time_Err_Array);
  TGraphErrors *Am_Ref_Ratio_vs_Time_Graph = new TGraphErrors(size,Time,Am_Ref_Ratio_vs_Time_Array,Time_Err,Am_Ref_Ratio_vs_Time_Err_Array);
  TGraphErrors *All_Ratio_vs_Time_Graph = new TGraphErrors(size,Time,All_Ratio_vs_Time_Array,Time_Err,All_Ratio_vs_Time_Err_Array);
  TGraphErrors *MeV_HMeV_Ratio_vs_Time_Graph = new TGraphErrors(size,Time,MeV_HMeV_Ratio_vs_Time_Array,Time_Err,MeV_HMeV_Ratio_vs_Time_Err_Array);

  
  // Normalized Signals vs Time
  /*TGraphErrors *Norm_Bi_vs_Time_Graph = new TGraphErrors(size,Time,Norm_Bi_Array,Time_Err,Norm_Bi_Err_Array);
  TGraphErrors *Norm_Am_vs_Time_Graph = new TGraphErrors(size,Time,Norm_Am_Array,Time_Err,Norm_Am_Err_Array);
  TGraphErrors *Norm_LED_Calo_vs_Time_Graph = new TGraphErrors(size,Time,Norm_LED_Calo_Array,Time_Err,Norm_LED_Calo_Err_Array);
  TGraphErrors *Norm_LED_Ref_vs_Time_Graph = new TGraphErrors(size,Time,Norm_LED_Ref_Array,Time_Err,Norm_LED_Ref_Err_Array);*/

  TGraphErrors *T_PB_vs_Time_Graph, *T_DB_vs_Time_Graph, *T_Lab_vs_Time_Graph, *T_Lab2_vs_Time_Graph,
  *Bi_vs_Temp_Graph, *Am_vs_Temp_Graph, *LED_Ref_vs_Temp_Graph, *LED_Calo_vs_Temp_Graph;

  if(TempInt){
    // Temperature over time 
    T_PB_vs_Time_Graph = new TGraphErrors(size,Time,T_PB_vs_Time_Array,Time_Err,T_Err_vs_Time_Array);
    T_DB_vs_Time_Graph = new TGraphErrors(size,Time,T_DB_vs_Time_Array,Time_Err,T_Err_vs_Time_Array);
    T_Lab_vs_Time_Graph = new TGraphErrors(size,Time,T_Lab_vs_Time_Array,Time_Err,T_Err_vs_Time_Array);
    T_Lab2_vs_Time_Graph = new TGraphErrors(size,Time,T_Lab2_vs_Time_Array,Time_Err,T_Err_vs_Time_Array);
  
    // Temperature correlations (more can be added)
    Bi_vs_Temp_Graph = new TGraphErrors(size,T_DB_vs_Time_Array,Bi_vs_Time_Array,T_Err_vs_Time_Array,Bi_vs_Time_Err_Array);
    Am_vs_Temp_Graph = new TGraphErrors(size,T_DB_vs_Time_Array,Am_vs_Time_Array,T_Err_vs_Time_Array,Am_vs_Time_Err_Array);
    LED_Ref_vs_Temp_Graph = new TGraphErrors(size,T_Lab2_vs_Time_Array,LED_Ref_vs_Time_Array,T_Err_vs_Time_Array,LED_Ref_vs_Time_Err_Array);
    LED_Calo_vs_Temp_Graph = new TGraphErrors(size, T_Lab2_vs_Time_Array, LED_Calo_vs_Time_Array, T_Err_vs_Time_Array, LED_Calo_vs_Time_Err_Array);
  }

  /*Double_t CorrAm_vs_Time_Err_Array[size], CorrBi_vs_Time_Err_Array[size], CorrCalo_vs_Time_Err_Array[size], CorrRef_vs_Time_Err_Array[size];
std::vector<double> Corr_Am_Err, Corr_Bi_Err, Corr_Calo_Err, Corr_Ref_Err;*/


  //*****************************//
  // Ratio and Temp vs Time Fits //
  //*****************************//

  Double_t Ratio_Time_r, Ratio_Time_yMean, Ratio_Time_slope, Ratio_Time_slope_err, Ratio_Time_const, Ratio_Time_const_err,
  	   Temp_DB_Time_r, Temp_DB_Time_yMean, Temp_DB_Time_slope, Temp_DB_Time_slope_err, Temp_DB_Time_const, Temp_DB_Time_const_err,
	   Temp_Lab_Time_r, Temp_Lab_Time_yMean, Temp_Lab_Time_slope, Temp_Lab_Time_slope_err, Temp_Lab_Time_const, Temp_Lab_Time_const_err;

  // Ratio vs Time Fit
  TF1 *Ratio_Time = new TF1("Ratio_Time","[0]*x+[1]",0,Time[size-1]);
  Ratio_Time->SetParameter(0,-0.33);
  //Ratio_Time->SetParLimits(0, -0.37, -0.25);
  Ratio_Time->SetParameter(1,1.01);
  Ratio_Time->SetLineColor(kRed);
  Ratio_Time->SetLineWidth(4);
  Ratio_vs_Time_Graph->Fit("Ratio_Time", "RSB");
  Ratio_Time_r = Ratio_vs_Time_Graph->GetCorrelationFactor();
  Ratio_Time_yMean = Ratio_vs_Time_Graph->GetMean(2);
  Ratio_Time_slope = Ratio_Time->GetParameter(0);
  Ratio_Time_slope_err = Ratio_Time->GetParError(0);
  Ratio_Time_const = Ratio_Time->GetParameter(1);
  Ratio_Time_const_err = Ratio_Time->GetParError(1);

  if(TempInt){
    // Temp DB vs Time Fit
    TF1 *Temp_DB_Time = new TF1("Temp_DB_Time","[0]*x+[1]",0,Time[size-1]);
    Temp_DB_Time->SetParameter(0,-0.33);
    //Temp_DB_Time->SetParLimits(0, -0.37, -0.25);
    Temp_DB_Time->SetParameter(1,1.01);
    Temp_DB_Time->SetLineColor(kRed);
    Temp_DB_Time->SetLineWidth(4);
    T_DB_vs_Time_Graph->Fit("Temp_DB_Time", "RSB");
    Temp_DB_Time_r = T_DB_vs_Time_Graph->GetCorrelationFactor();
    Temp_DB_Time_yMean = T_DB_vs_Time_Graph->GetMean(2);
    Temp_DB_Time_slope = Temp_DB_Time->GetParameter(0);
    Temp_DB_Time_slope_err = Temp_DB_Time->GetParError(0);
    Temp_DB_Time_const = Temp_DB_Time->GetParameter(1);
    Temp_DB_Time_const_err = Temp_DB_Time->GetParError(1);

    // Temp Lab vs Time Fit
    TF1 *Temp_Lab_Time = new TF1("Temp_Lab_Time","[0]*x+[1]",0,Time[size-1]);
    Temp_Lab_Time->SetParameter(0,-0.33);
    //Temp_Lab_Time->SetParLimits(0, -0.37, -0.25);
    Temp_Lab_Time->SetParameter(1,1.01);
    Temp_Lab_Time->SetLineColor(kRed);
    Temp_Lab_Time->SetLineWidth(4);
    T_Lab2_vs_Time_Graph->Fit("Temp_Lab_Time", "RSB");
    Temp_Lab_Time_r = T_Lab2_vs_Time_Graph->GetCorrelationFactor();
    Temp_Lab_Time_yMean = T_Lab2_vs_Time_Graph->GetMean(2);
    Temp_Lab_Time_slope = Temp_Lab_Time->GetParameter(0);
    Temp_Lab_Time_slope_err = Temp_Lab_Time->GetParError(0);
    Temp_Lab_Time_const = Temp_Lab_Time->GetParameter(1);
    Temp_Lab_Time_const_err = Temp_Lab_Time->GetParError(1);
  }


  //**************************//
  // Ratio and Temp Residuals //
  //**************************//

  std::vector<double> Temp_DB_Line_Diff, Temp_DB_Line_Diff_Err, Temp_Lab_Line_Diff, Temp_Lab_Line_Diff_Err, PreCorrRatio_DB_Err, Term_PreCorr_Err, Term_Corr_Err;
  Double_t Temp_DB_Line_Diff_Array[size], Temp_DB_Line_Diff_Err_Array[size], Temp_Lab_Line_Diff_Array[size], Temp_Lab_Line_Diff_Err_Array[size], PreCorrRatio_DB_vs_Time_Array[size], PreCorrRatio_DB_vs_Time_Err_Array[size];
  
  for(int i=0; i<size; i++){ 

    // Ratio vs Time Residuals Vectors
    Ratio_Line_Diff.push_back(Ratio_vs_Time_Array[i]-(Ratio_Time_slope*Time[i]+Ratio_Time_const));
    Ratio_Line_Diff_Err.push_back(Ratio_vs_Time_Err_Array[i]);
   
    if(TempInt){ 
      // Temp_DB vs Time Residuals Vectors		
      Temp_DB_Line_Diff.push_back(T_DB_vs_Time_Array[i]-(Temp_DB_Time_slope*Time[i]+Temp_DB_Time_const));
      Temp_DB_Line_Diff_Err.push_back(T_Err_vs_Time_Array[i]);

      // Temp_Lab vs Time Residuals Vectors
      Temp_Lab_Line_Diff.push_back(T_Lab2_vs_Time_Array[i]-(Temp_Lab_Time_slope*Time[i]+Temp_Lab_Time_const));
      Temp_Lab_Line_Diff_Err.push_back(T_Err_vs_Time_Array[i]);
    }
  }

  // Fill Residuals Arrays
  for(int i=0; i<size; i++){
    Ratio_Line_Diff_Array[i] = Ratio_Line_Diff[i];
    Ratio_Line_Diff_Err_Array[i] = Ratio_Line_Diff_Err[i];
	
    if(TempInt){	
      Temp_DB_Line_Diff_Array[i] = Temp_DB_Line_Diff[i];
      Temp_DB_Line_Diff_Err_Array[i] = Temp_DB_Line_Diff_Err[i];

      Temp_Lab_Line_Diff_Array[i] = Temp_Lab_Line_Diff[i];
      Temp_Lab_Line_Diff_Err_Array[i] = Temp_Lab_Line_Diff_Err[i];
    }
  }
  
  // TGraphs of Residuals
  TGraphErrors *Ratio_Line_Diff_Graph, *Temp_DB_Line_Diff_Graph, *Temp_Lab_Line_Diff_Graph;
  
  if(TempInt){
    Ratio_Line_Diff_Graph = new TGraphErrors(size,Time,Ratio_Line_Diff_Array,Time_Err,Ratio_Line_Diff_Err_Array);
    Temp_DB_Line_Diff_Graph = new TGraphErrors(size,Time,Temp_DB_Line_Diff_Array,Time_Err,Temp_DB_Line_Diff_Err_Array);
    Temp_Lab_Line_Diff_Graph = new TGraphErrors(size,Time,Temp_Lab_Line_Diff_Array,Time_Err,Temp_Lab_Line_Diff_Err_Array);
  }

  //******************************************************************************//
  // Simple Temperature Correction (Current Primary Correction though INCOMPLETE) //
  //******************************************************************************//		
 
  float precorr_c = 19; // constant to change
  float corr_c = 1.248; // constant to change

  if(TempInt){
    // First Temperature Correction Vector (to address obvious shape "signatures")
    for(int i=0; i<size; i++){
      //PreCorrRatio_DB.push_back(Ratio_vs_Time_Array[i]-(Ratio_Temp_slope*5.4*(T_DB_vs_Time_Array[i]-T_DB_vs_Time_Array[0]))); // Optimized for smallest std dev of 2nd correction
      PreCorrRatio_DB.push_back(Ratio_vs_Time_Array[i]-(Ratio_Temp_slope * precorr_c * (T_DB_vs_Time_Array[i]-T_DB_vs_Time_Array[0])));
      //PreCorrRatio_DB.push_back(Ratio_vs_Time_Array[i]);
      //PreCorrRatio_DB.push_back(Ratio_vs_Time_Array[i]-(Ratio_Temp_slope*41*(T_DB_vs_Time_Array[i]-T_DB_vs_Time_Array[0])));
      //PreCorrRatio_DB.push_back(Ratio_vs_Time_Array[i]-(Ratio_Temp_slope*5*(T_DB_vs_Time_Array[i]-T_DB_vs_Time_Array[0])));
      Term_PreCorr_Err.push_back(sqrt(pow(Ratio_Temp_slope_err*(T_DB_vs_Time_Array[i]-T_DB_vs_Time_Array[0]),2)
				    + pow(Ratio_Temp_slope*T_Err_vs_Time_Array[i],2))); // second "term" in PreCorrRatio_Err calc
    }
  
    // First Temperature Correction Error Vector
    for(int i=0; i<size; i++){
      PreCorrRatio_DB_Err.push_back(sqrt(precorr_c * (pow(Ratio_vs_Time_Err[i],2) + pow(Term_PreCorr_Err[i],2))));
      //PreCorrRatio_DB_Err.push_back(Ratio_vs_Time_Err[i]);
    }

    // Second Temperature Correction Vector (to address negative linear trend)
    for(int i=0; i<size; i++){
      //CorrRatio_DB.push_back(PreCorrRatio_DB[i]-((Time[i]-Time[0])*Ratio_Time_slope*1));
      CorrRatio_DB.push_back(PreCorrRatio_DB[i]-((Time[i]-Time[0])*Ratio_Time_slope * corr_c)); // Optimized for smallest linear trend
      //CorrRatio_DB.push_back(PreCorrRatio_DB[i]-((Time[i]-Time[0])*Ratio_Time_slope*1.7));
      Term_Corr_Err.push_back(sqrt(pow(Ratio_Time_slope_err*Time[i],2) + pow(Ratio_Time_slope*Time_Err[i],2))); // second "term" in CorrRatio_Err calc
      //CorrRatio_DB_Err.push_back(Ratio_vs_Time_Err[i]);
    }

    // Second Temperature Correction Error Vector
    for(int i=0; i<size; i++){
      CorrRatio_DB_Err.push_back(sqrt(corr_c * (pow(PreCorrRatio_DB_Err[i],2) + pow(Term_Corr_Err[i],2))));
      //cout << " CorrRatio_DB[i]: " << CorrRatio_DB[i] << " Ratio_vs_Time[i]: " << Ratio_vs_Time[i] << endl; 
    }
  }

  // Evaluate mins and maxes of ratios to restrict histogram ranges
  Double_t Ratio_Max, Ratio_Min, CorrRatio_Max, CorrRatio_Min, PreCorrRatio_Max, PreCorrRatio_Min;
  Ratio_Max = *std::max_element(Ratio_vs_Time.begin(),Ratio_vs_Time.end());
  Ratio_Min = *std::min_element(Ratio_vs_Time.begin(),Ratio_vs_Time.end());

  if(TempInt){
    PreCorrRatio_Max = *std::max_element(PreCorrRatio_DB.begin(),PreCorrRatio_DB.end());
    PreCorrRatio_Min = *std::min_element(PreCorrRatio_DB.begin(),PreCorrRatio_DB.end()); 
    CorrRatio_Max = *std::max_element(CorrRatio_DB.begin(),CorrRatio_DB.end());
    CorrRatio_Min = *std::min_element(CorrRatio_DB.begin(),CorrRatio_DB.end());
  }

  cout << "\nRatio_Max: " << Ratio_Max << " Ratio_Min: " << Ratio_Min << endl;
  
  if(TempInt){
    cout << "PreCorrRatio_Max: " << PreCorrRatio_Max << " PreCorrRatio_Min: " << PreCorrRatio_Min << " CorrRatio_Max: " << CorrRatio_Max << " CorrRatio_Min: " << CorrRatio_Min << endl;  
  } 
  
       
  // Histograms depicting the spread of each step along the ratio temp correction
  TH1 *Ratio_Spread, *PreCorrRatio_Spread, *CorrRatio_Spread;

  Ratio_Spread = new TH1D("Ratio Spread", "Ratio Spread", 300, 0.98, 1.01); //PreCorrRatio_Min
  if(TempInt){
    PreCorrRatio_Spread = new TH1D("1st Correction Spread", "1st Correction Spread", 300, 0.98, 1.01);
    CorrRatio_Spread = new TH1D("2nd Correction Spread", "2nd Correction Spread", 300, 0.98, 1.01);
  }
  
  // Fill histograms
  for(int i=0; i<size; i++){
    Ratio_Spread->Fill(Ratio_vs_Time[i]);
    if(TempInt){
      PreCorrRatio_Spread->Fill(PreCorrRatio_DB[i]); 
      CorrRatio_Spread->Fill(CorrRatio_DB[i]);
    }
  }
  
  // Print statement for checking constant values used in temp correction
  cout << " \nRatio_Time_slope: "<< Ratio_Time_slope <<" Temp_DB_Time_slope: " << Temp_DB_Time_slope << " Ratio_Temp_slope: " << Ratio_Temp_slope << "\n" <<  endl;
  
  // Fill temp corrected ratio arrays
  if(TempInt){
    for(int i=0; i<size; i++){
      PreCorrRatio_DB_vs_Time_Array[i] = PreCorrRatio_DB[i];
      PreCorrRatio_DB_vs_Time_Err_Array[i] = PreCorrRatio_DB_Err[i];
      CorrRatio_DB_vs_Time_Array[i] = CorrRatio_DB[i];
      CorrRatio_DB_vs_Time_Err_Array[i] = CorrRatio_DB_Err[i];
    }
  }

  /*if(TempInt){
    for(int i=0; i<size; i++){
  	if(T_DB_vs_Time_Array[i] == 73.7){
		cout << " Time at start of window is: " << Time[i] << "    " << Time_Data[i] << endl;
	}
	if(T_DB_vs_Time_Array[i] == 73.5){
		cout << " Time at end of window is: " << Time[i] << "    " << Time_Data[i] << endl;
	}
    }
  }*/
  

  //*************************************//
  // CorrRatio vs Time Fit and Residuals //
  //*************************************//

  TGraphErrors *RatioTemp_Lab_vs_Time_Graph, *NormRatio_vs_Time_Graph, *RecNormTemp_Lab_vs_Time_Graph, 
	       *CorrRatio_Lab_vs_Time_Graph, *CorrRatio_DB_vs_Time_Graph, *PreCorrRatio_DB_vs_Time_Graph;
  TF1 *CorrRatio_Time;
  Double_t CorrRatio_Time_r, CorrRatio_Time_yMean, CorrRatio_Time_slope, CorrRatio_Time_slope_err, CorrRatio_Time_const, CorrRatio_Time_const_err;

  if(TempInt){
    //RatioTemp_Lab_vs_Time_Graph = new TGraphErrors(size,Time,RatioTemp_Lab_vs_Time_Array,Time_Err,Ratio_vs_Time_Err_Array);
    //NormRatio_vs_Time_Graph = new TGraphErrors(size,Time,NormRatio_vs_Time_Array,Time_Err,Ratio_vs_Time_Err_Array);
    //RecNormTemp_Lab_vs_Time_Graph = new TGraphErrors(size,Time,RecNormTemp_Lab_vs_Time_Array,Time_Err,Ratio_vs_Time_Err_Array);
    CorrRatio_Lab_vs_Time_Graph = new TGraphErrors(size,Time,CorrRatio_Lab_vs_Time_Array,Time_Err,Ratio_vs_Time_Err_Array);
    CorrRatio_DB_vs_Time_Graph = new TGraphErrors(size,Time,CorrRatio_DB_vs_Time_Array,Time_Err,CorrRatio_DB_vs_Time_Err_Array);
    PreCorrRatio_DB_vs_Time_Graph = new TGraphErrors(size,Time,PreCorrRatio_DB_vs_Time_Array,Time_Err,PreCorrRatio_DB_vs_Time_Err_Array);

    // CorrRatio vs Time Fit
    CorrRatio_Time = new TF1("CorrRatio_Time","[0]*x+[1]",0,Time[size-1]);
    CorrRatio_Time->SetParameter(0,-0.33);
    //CorrRatio_Time->SetParLimits(0, -0.37, -0.25);
    CorrRatio_Time->SetParameter(1,1.01);
    CorrRatio_Time->SetLineColor(kRed);
    CorrRatio_Time->SetLineWidth(4);
    CorrRatio_DB_vs_Time_Graph->Fit("CorrRatio_Time", "RSB");
    CorrRatio_Time_r = CorrRatio_DB_vs_Time_Graph->GetCorrelationFactor();
    CorrRatio_Time_yMean = CorrRatio_DB_vs_Time_Graph->GetMean(2);
    CorrRatio_Time_slope = CorrRatio_Time->GetParameter(0);
    CorrRatio_Time_slope_err = CorrRatio_Time->GetParError(0);
    CorrRatio_Time_const = CorrRatio_Time->GetParameter(1);
    CorrRatio_Time_const_err = CorrRatio_Time->GetParError(1);
  }
  
  std::vector<double> CorrRatio_Line_Diff, CorrRatio_Line_Diff_Err, CorrRatio_Line_Ratio, CorrRatio_Line_Ratio_Err;
  Double_t CorrRatio_Line_Diff_Array[size], CorrRatio_Line_Diff_Err_Array[size], CorrRatio_Line_Ratio_Array[size], CorrRatio_Line_Ratio_Err_Array[size];
  TGraphErrors *CorrRatio_Line_Diff_Graph, *CorrRatio_Line_Ratio_Graph;

  if(TempInt){
    // Fill Vectors for CorrRatio vs Time Fit Residuals
    for(int i=0; i<size; i++){
      CorrRatio_Line_Diff.push_back(CorrRatio_DB_vs_Time_Array[i]-(CorrRatio_Time_slope*Time[i]+CorrRatio_Time_const));
      CorrRatio_Line_Diff_Err.push_back(CorrRatio_DB_vs_Time_Err_Array[i]);
      CorrRatio_Line_Ratio.push_back(CorrRatio_DB_vs_Time_Array[i]/(CorrRatio_Time_slope*Time[i]+CorrRatio_Time_const));
      CorrRatio_Line_Ratio_Err.push_back(CorrRatio_DB_vs_Time_Err_Array[i]);
    }    

    // Fill Arrays for CorrRatio vs Time Fit Residuals
    for(int i=0; i<size; i++){
      CorrRatio_Line_Diff_Array[i] = CorrRatio_Line_Diff[i];
      CorrRatio_Line_Diff_Err_Array[i] = CorrRatio_Line_Diff_Err[i];
      CorrRatio_Line_Ratio_Array[i] = CorrRatio_Line_Ratio[i];
      CorrRatio_Line_Ratio_Err_Array[i] = CorrRatio_Line_Ratio_Err[i];
    }    
  
    // Graphs for CorrRatio vs Time Fit Residuals
    CorrRatio_Line_Diff_Graph = new TGraphErrors(size,Time,CorrRatio_Line_Diff_Array,Time_Err,CorrRatio_Line_Diff_Err_Array);
    CorrRatio_Line_Ratio_Graph = new TGraphErrors(size,Time,CorrRatio_Line_Ratio_Array,Time_Err,CorrRatio_Line_Ratio_Err_Array);
  }


  //*************************************************************************************//
  // Temperature Corrected Ratio For Normalization Correction Method (Seems to Not Work) //
  //*************************************************************************************//
  
  //Double_t Bi_t, Bi_t_err, LED_Calo_t, LED_Calo_t_err, LED_Ref_t, LED_Ref_t_err, Am_t, Am_t_err, LED_Ref_Pred, LED_Ref_Pred_err, LED_Calo_Pred, LED_Calo_Pred_err, Bi_Pred, Bi_Pred_err, Ratio, Ratio_err;
  std::vector<double> Corr_Bi_Pred_vs_Time, Corr_Ratio_vs_Time, Bi_0_Bi_t, Am_t_Am_0, Ref_0_Ref_t, Calo_t_Calo_0, Corr_Bi_Pred_vs_Time_Err, Corr_Ratio_vs_Time_Err;
  Double_t Corr_Bi_Pred_vs_Time_Array[size], Corr_Ratio_vs_Time_Array[size], Corr_Ratio_vs_Time_Err_Array[size], Bi_0_Bi_t_Array[size], Am_t_Am_0_Array[size], Ref_0_Ref_t_Array[size], Calo_t_Calo_0_Array[size], Bi_to, LED_Calo_to, LED_Ref_to, Am_to, Bi_lasto, LED_Calo_lasto;
  /* 
  // Start values (i.e. calibration run values)
  Bi_0 = Corr_Bi[0];
  Bi_0_err = MeV_Mean_Err[0];  //Corr_Bi_Err[0];
  //Bi2_0 = HMeV_Mean[0];
  //Bi2_0_err = HMeV_Mean_Err[0];
  LED_Calo_0 = Corr_Calo[0];
  LED_Calo_0_err = LED_Mean_Err[0];  //Corr_Calo_Err[0];
  LED_Ref_0 = Corr_Ref[0];
  LED_Ref_0_err = LED2_Mean_Err[0];  //Corr_Ref_Err[0];
  Am_0 = Corr_Am[0];
  Am_0_err = Americium1_Mean_Err[0];
  //int corrsize;
  //corrsize = size-rejections;
  // Later values for LED and Am ONLY
  */

  
/*
  LED_Calo_last = LED_Calo_0;
  counter1 = 0; 
  counter2 = 0;
  for(int i=0; i<size; i++){
    Bi_to	   = Corr_Bi[0];  //MeV_Mean[i]
    Bi_t           = Corr_Bi[i];
    Bi_t_err       = MeV_Mean_Err[i];  //Corr_Bi_Err[i];
    LED_Calo_to    = Corr_Calo[0];  //LED_Mean[i];
    LED_Calo_t     = Corr_Calo[i];
    LED_Calo_t_err = LED_Mean_Err[i];  //Corr_Calo_Err[i];
    LED_Ref_to     = Corr_Ref[0];  //LED2_Mean[i];
    LED_Ref_t      = Corr_Ref[i];
    LED_Ref_t_err  = LED2_Mean_Err[i];  //Corr_Ref_Err[i];
    Am_to	   = Corr_Am[0];  //Americium1_Mean[i];
    Am_t           = Corr_Am[i];
    Am_t_err       = Americium1_Mean_Err[i];  //Corr_Am_Err[i];

    // Skip runs with large Bi peak deviations unaccounted for in LED peak
    if((TMath::Abs(Bi_to - Bi_lasto) > 10*Bi_t_err) &&
    !(TMath::Abs(LED_Calo_to - LED_Calo_lasto) > 10*LED_Calo_t_err)){
      if(Debug) {
        printf("Reese: Skipping Run %d (timestamp: %f): Bi mean value = %f while LED mean value = %f\n",i, Time_Data.at(i), Bi_t, LED_Calo_t);
        }
      continue;
    }
    
    else{
    //if(Am_t/Am_0 < 0.5 or Am_t/Am_0 < 2){
      continue;
    //}

  
    Bi_lasto = Bi_to;
    Bi_last = Bi_t; // Update comparison values
    LED_Calo_lasto = LED_Calo_to;
    //Am_t = Am_0;

    // Debugging outputs
    //if(Debug) printf("%f: Bi mean = %f LED mean  = %f", i, Time_Data.at(i), Bi_t, LED_Calo_t);    

    // Use Am to predict Ref LED (correct for drifts)
    LED_Ref_Pred = LED_Ref_0 * (Am_t/Am_0);
    LED_Ref_Pred_err = sqrt(pow(LED_Ref_0_err*Am_t/Am_0,2)
                            + pow(LED_Ref_0*Am_t_err/Am_0,2)
                            + pow(LED_Ref_0*Am_t*Am_0_err/(Am_0*Am_0),2));
    //if(Debug) printf("LED pred: %f*%f = %f\n", Am_t, Am_0, (Am_t/Am_0));
    // Use Ref LED to predict Calo LED
    LED_Calo_Pred = LED_Calo_t * (LED_Ref_Pred/LED_Ref_t);
    LED_Calo_Pred_err = sqrt(pow(LED_Calo_t_err*LED_Ref_Pred/LED_Ref_t,2)
                             + pow(LED_Calo_t*LED_Ref_Pred_err/LED_Ref_t,2)
                             + pow(LED_Calo_t*LED_Ref_Pred*LED_Ref_t_err/(LED_Ref_t*LED_Ref_t),2));
  
    // Use Calo LED to predict Bi
    Bi_Pred = Bi_0 * (LED_Calo_Pred/LED_Calo_0); //(Bi_0 * LED_Calo_t * LED_Ref_0 *Am_t)/(LED_Calo_0 * LED_Ref_t *Am_0); //just trying a more explicit implementation of the prediction algorithm, but same effect
    //if (Debug) printf("(%f*%f*%f*%f)/(%f*%f*%f)\n", Bi_0, LED_Calo_t, LED_Ref_0, Am_t, LED_Calo_0, LED_Ref_t, Am_0);
    Bi_Pred_err = sqrt(pow(Bi_0_err*LED_Calo_Pred/LED_Calo_0,2)
                       + pow(Bi_0*LED_Calo_Pred_err/LED_Calo_0,2)
                       + pow(Bi_0*LED_Calo_Pred*LED_Calo_0_err/(LED_Calo_0*LED_Calo_0),2));

    // Calculate the ratio of predicted to measured Bi
    Ratio = Bi_Pred/Bi_t;
    //cout << " Bi_Pred:  " << Bi_Pred << " Bi_t:  " << Bi_t << " Ratio:  " << Ratio << endl;
    Ratio_err = sqrt(pow(Bi_Pred_err/Bi_t,2) + pow(Bi_Pred*Bi_t_err/(Bi_t*Bi_t),2));

    // Cross check with full calculation 
    Double_t test = ((Bi_0/Bi_t)*(LED_Calo_t/LED_Calo_0)*(LED_Ref_0/LED_Ref_t)*(Am_t/Am_0));

    //Fill vectors
    Bi_vs_Time.push_back(Bi_t);
    Bi_vs_Time_Err.push_back(Bi_t_err);
    LED_Calo_vs_Time.push_back(LED_Calo_t);
    LED_Calo_vs_Time_Err.push_back(LED_Calo_t_err);
    LED_Ref_vs_Time.push_back(LED_Ref_t);
    LED_Ref_vs_Time_Err.push_back(LED_Ref_t_err);
    Am_vs_Time.push_back(Am_t);
    Am_vs_Time_Err.push_back(Am_t_err);
    
    //cout << "Am_t :  " << Am_t << "  at:  " << i << endl;  
    Corr_Bi_Pred_vs_Time.push_back(Bi_Pred);
    Corr_Bi_Pred_vs_Time_Err.push_back(Bi_Pred_err);
    Corr_Ratio_vs_Time.push_back(Ratio);
    Corr_Ratio_vs_Time_Err.push_back(Ratio_err);
    cout << "Corr Ratio_vs_Time_Err @" << i << ": " << Corr_Ratio_vs_Time_Err[i] << endl; 
    //cout << "Corr_Ratio_vs_Time[i]  " << Corr_Ratio_vs_Time[i] << endl; 
    
   
// Troubleshoot with Ratio Plots *****
     Bi_0_Bi_t.push_back(Bi_0/Bi_t); 
     Am_t_Am_0.push_back(Am_t/Am_0);
     Ref_0_Ref_t.push_back(LED_Ref_0/LED_Ref_t);
     Calo_t_Calo_0.push_back(LED_Calo_t/LED_Calo_0);

     }	               
  //cout << "  Bi_0_Bi_t:  " << Bi_0/Bi_t << "  Am_t_Am_0:  " << Am_t/Am_0 << "  Ref_0_Ref_t:  " << LED_Ref_0/LED_Ref_t << "  Calo_t_Calo_0:  " << LED_Calo_t/LED_Calo_0 << "  Ratio  " << test << endl;
	
      
  }
  
  //const int mysize = corrsize - rejections;
	cout << "Number of Analysis Rejections:  " << rejections << endl;
 */ 


  //******************//
  // Corrected Graphs //
  //******************//
/*
  // Corr_Am Graphs
  TGraphErrors *CorrAm_vs_Time_Graph = new TGraphErrors(size,Time,CorrAm_vs_Time_Array,Time_Err,CorrAm_vs_Time_Err_Array);
  TGraphErrors *CorrAm_vs_Temp_Graph = new TGraphErrors(size,T_DB_vs_Time_Array,CorrAm_vs_Time_Array,T_Err_vs_Time_Array, CorrAm_vs_Time_Err_Array);

  // Corr_Bi Graphs
  TGraphErrors *CorrBi_vs_Time_Graph = new TGraphErrors(size,Time,CorrBi_vs_Time_Array,Time_Err,CorrBi_vs_Time_Err_Array);
  TGraphErrors *CorrBi_vs_Temp_Graph = new TGraphErrors(size,T_DB_vs_Time_Array,CorrBi_vs_Time_Array,T_Err_vs_Time_Array, CorrBi_vs_Time_Err_Array);

  // Corr_Ref_LED Graphs
  TGraphErrors *CorrRef_vs_Time_Graph = new TGraphErrors(size,Time,CorrRef_vs_Time_Array,Time_Err, CorrRef_vs_Time_Err_Array);
  TGraphErrors *CorrRef_vs_Temp_Graph = new TGraphErrors(size,T_Lab2_vs_Time_Array,CorrRef_vs_Time_Array,T_Err_vs_Time_Array, CorrRef_vs_Time_Err_Array);

  // Corr_Calo_LED Graphs
  TGraphErrors *CorrCalo_vs_Time_Graph = new TGraphErrors(size,Time,CorrCalo_vs_Time_Array,Time_Err,CorrCalo_vs_Time_Err_Array);
  TGraphErrors *CorrCalo_vs_Temp_Graph = new TGraphErrors(size,T_Lab2_vs_Time_Array,CorrCalo_vs_Time_Array,T_Err_vs_Time_Array, CorrCalo_vs_Time_Err_Array);


  // Corrected Ratio Graphs
  TGraphErrors *Corr_Ratio_vs_Time_Graph = new TGraphErrors(size,Time,Corr_Ratio_vs_Time_Array,Time_Err,Corr_Ratio_vs_Time_Err_Array);
  TGraphErrors *Corr_Ratio_vs_Temp_Graph = new TGraphErrors(size,T_Lab2_vs_Time_Array,Corr_Ratio_vs_Time_Array,T_Err_vs_Time_Array, Corr_Ratio_vs_Time_Err_Array);

  //TGraph *Corr_Ratio_vs_Time_Graph = new TGraph(size,Time, Corr_Ratio_vs_Time_Array);
  //TGraph *Corr_Ratio_vs_Temp_Graph = new TGraph(size,T_DB_vs_Time_Array,Corr_Ratio_vs_Time_Array);
  
  for(int i=0; i<size; i++){
    cout << "Corrected Ratio @i:  " << Corr_Ratio_vs_Time_Array[i] << endl;
    cout << "Original Ratio @i:  "  << Ratio_vs_Time_Array[i] << endl;
  }*/



  //********************//
  // Time For Plotting! //
  //******************* //

  // Setting some style options
  gStyle->SetPadRightMargin(0.03); 
  gStyle->SetPadLeftMargin(0.06);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleFontSize(0.08);
  /*
  gStyle->SetPalette(1,0);
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(11111);
  */

  if(DrawPlots){
    // MeV:HMeV Ratio Plot
    TCanvas *c15 = new TCanvas("c15","c15", 50,50,2200,400);
    c15->cd();
    c15->SetGrid();
    MeV_HMeV_Ratio_vs_Time_Graph->SetTitle("MeV:HMeV Ratio over Time");
    MeV_HMeV_Ratio_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    MeV_HMeV_Ratio_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    MeV_HMeV_Ratio_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
    MeV_HMeV_Ratio_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.88);
    MeV_HMeV_Ratio_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
    MeV_HMeV_Ratio_vs_Time_Graph->GetYaxis()->SetTitle("Peaks Ratio MeV:HMeV");
    MeV_HMeV_Ratio_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
    MeV_HMeV_Ratio_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
    MeV_HMeV_Ratio_vs_Time_Graph->GetYaxis()->SetLabelSize(0.06);
    MeV_HMeV_Ratio_vs_Time_Graph->SetMarkerStyle(7);
    MeV_HMeV_Ratio_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
    MeV_HMeV_Ratio_vs_Time_Graph->SetMarkerColor(kRed);
    MeV_HMeV_Ratio_vs_Time_Graph->Draw("AP");
    if(SavePlots){
      c15->Print("Plots/MeVHMeVRatio_vs_Time.png");
    }
  
     
    if(TempInt){
    // Temp Corrected Ratio, Temp Corrected Ratio Fit Residuals, and a Ratio of Temp Corr Ratio over Fit Values
    TCanvas *c14 = new TCanvas("c14","c14", 10,10,2000,1200); 
    c14->Divide(1,3);
    c14->cd(1);
    c14->SetGrid();
    CorrRatio_DB_vs_Time_Graph->SetTitle("Corrected Ratio of Predicted to Measured ^{207}Bi Peak Over Time");
    CorrRatio_DB_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    CorrRatio_DB_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    CorrRatio_DB_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
    CorrRatio_DB_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.88);
    CorrRatio_DB_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
    CorrRatio_DB_vs_Time_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
    CorrRatio_DB_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
    CorrRatio_DB_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
    CorrRatio_DB_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
    CorrRatio_DB_vs_Time_Graph->GetYaxis()->SetRangeUser(0.975, 1.015); 
    CorrRatio_DB_vs_Time_Graph->SetMarkerStyle(7);
    CorrRatio_DB_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
    CorrRatio_DB_vs_Time_Graph->SetMarkerColor(kRed);
    CorrRatio_DB_vs_Time_Graph->Draw("AP");
    if(DrawLines){
      TLine *upper10 = new TLine(0,1.01,Time[size-1],1.01);
      upper10->SetLineColor(kRed);
      upper10->SetLineWidth(2);
      upper10->Draw(); 
      TLine *lower10 = new TLine(0,0.99,Time[size-1],0.99);
      lower10->SetLineColor(kRed);
      lower10->SetLineWidth(2);
      lower10->Draw();
    }

    c14->cd(2);
    c14->SetGrid();
    CorrRatio_Line_Diff_Graph->SetTitle("Corrected Ratio vs Time Fit Residuals");
    CorrRatio_Line_Diff_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    CorrRatio_Line_Diff_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    CorrRatio_Line_Diff_Graph->GetXaxis()->SetTitleSize(0.08);
    CorrRatio_Line_Diff_Graph->GetXaxis()->SetTitleOffset(0.88);
    CorrRatio_Line_Diff_Graph->GetXaxis()->SetLabelSize(0.08);
    CorrRatio_Line_Diff_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
    CorrRatio_Line_Diff_Graph->GetYaxis()->SetTitleSize(0.08);
    CorrRatio_Line_Diff_Graph->GetYaxis()->SetTitleOffset(0.4);
    CorrRatio_Line_Diff_Graph->GetYaxis()->SetLabelSize(0.08);
    CorrRatio_Line_Diff_Graph->GetYaxis()->SetRangeUser(-0.025,0.015);
    CorrRatio_Line_Diff_Graph->SetMarkerStyle(7);
    CorrRatio_Line_Diff_Graph->SetLineColorAlpha(kBlack, 0.6);
    CorrRatio_Line_Diff_Graph->SetMarkerColor(kRed);
    CorrRatio_Line_Diff_Graph->Draw("AP");
    if(DrawLines){
      TLine *upper9 = new TLine(0,0.01,Time[size-1],0.01);
      upper9->SetLineColor(kRed);
      upper9->SetLineWidth(2);
      upper9->Draw(); 
      TLine *lower9 = new TLine(0,-0.01,Time[size-1],-0.01);
      lower9->SetLineColor(kRed);
      lower9->SetLineWidth(2);
      lower9->Draw();
    }

    c14->cd(3);
    c14->SetGrid();
    CorrRatio_Line_Ratio_Graph->SetTitle("Ratio of Corrected Ratio:Fit Values vs Time");
    CorrRatio_Line_Ratio_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    CorrRatio_Line_Ratio_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    CorrRatio_Line_Ratio_Graph->GetXaxis()->SetTitleSize(0.08);
    CorrRatio_Line_Ratio_Graph->GetXaxis()->SetTitleOffset(0.88);
    CorrRatio_Line_Ratio_Graph->GetXaxis()->SetLabelSize(0.08);
    CorrRatio_Line_Ratio_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
    CorrRatio_Line_Ratio_Graph->GetYaxis()->SetTitleSize(0.08);
    CorrRatio_Line_Ratio_Graph->GetYaxis()->SetTitleOffset(0.4);
    CorrRatio_Line_Ratio_Graph->GetYaxis()->SetLabelSize(0.08);
    CorrRatio_Line_Ratio_Graph->GetYaxis()->SetRangeUser(0.975,1.015);
    CorrRatio_Line_Ratio_Graph->SetMarkerStyle(7);
    CorrRatio_Line_Ratio_Graph->SetLineColorAlpha(kBlack, 0.6);
    CorrRatio_Line_Ratio_Graph->SetMarkerColor(kRed);
    CorrRatio_Line_Ratio_Graph->Draw("AP");
    if(DrawLines){
      TLine *upper8 = new TLine(0,1.01,Time[size-1],1.01);
      upper8->SetLineColor(kRed);
      upper8->SetLineWidth(2);
      upper8->Draw();
      TLine *lower8 = new TLine(0,0.99,Time[size-1],0.99);
      lower8->SetLineColor(kRed);
      lower8->SetLineWidth(2);
      lower8->Draw();
    }
 
    if(SavePlots){ 
      c14->Print("Plots/CorrRatio_Line_Diff.png");
    }

    
    // vs Time Fit Residuals for Ratio and Temp plots
    TCanvas *c13 =new TCanvas("c13", "c13", 10,10,2000,1200);
    c13->Divide(1,3);
    c13->cd(1);
    c13->SetGrid();
    Ratio_Line_Diff_Graph->SetTitle("Ratio vs Time Fit Residuals");
    Ratio_Line_Diff_Graph->GetXaxis()->SetRangeUser(0, Time[size-1]);
    Ratio_Line_Diff_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    Ratio_Line_Diff_Graph->GetXaxis()->SetTitleSize(0.08);
    Ratio_Line_Diff_Graph->GetXaxis()->SetTitleOffset(0.88);
    Ratio_Line_Diff_Graph->GetXaxis()->SetLabelSize(0.08);
    Ratio_Line_Diff_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
    Ratio_Line_Diff_Graph->GetYaxis()->SetTitleSize(0.08);
    Ratio_Line_Diff_Graph->GetYaxis()->SetTitleOffset(0.4);
    Ratio_Line_Diff_Graph->GetYaxis()->SetLabelSize(0.08);
    Ratio_Line_Diff_Graph->SetMarkerStyle(7);
    Ratio_Line_Diff_Graph->SetLineColorAlpha(kBlack, 0.6);
    Ratio_Line_Diff_Graph->SetMarkerColor(kRed);
    Ratio_Line_Diff_Graph->Draw("AP");
    if(DrawLines){
      TLine *upper7 = new TLine(0,0.01,Time[size-1],0.01);
      upper7->SetLineColor(kRed);
      upper7->SetLineWidth(2);
      upper7->Draw();
      TLine *lower7 = new TLine(0,-0.01,Time[size-1],-0.01);
      lower7->SetLineColor(kRed);
      lower7->SetLineWidth(2);
      lower7->Draw();
    }    
  
    c13->cd(2);
    c13->SetGrid();
    Temp_DB_Line_Diff_Graph->SetTitle("Temp_DB vs Time Fit Residuals");
    Temp_DB_Line_Diff_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    Temp_DB_Line_Diff_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    Temp_DB_Line_Diff_Graph->GetXaxis()->SetTitleSize(0.08);
    Temp_DB_Line_Diff_Graph->GetXaxis()->SetTitleOffset(0.88);
    Temp_DB_Line_Diff_Graph->GetXaxis()->SetLabelSize(0.08);
    Temp_DB_Line_Diff_Graph->GetYaxis()->SetTitle("Temperature [Deg. F]");
    Temp_DB_Line_Diff_Graph->GetYaxis()->SetTitleSize(0.08);
    Temp_DB_Line_Diff_Graph->GetYaxis()->SetTitleOffset(0.4);
    Temp_DB_Line_Diff_Graph->GetYaxis()->SetLabelSize(0.08);
    Temp_DB_Line_Diff_Graph->SetMarkerStyle(7);
    Temp_DB_Line_Diff_Graph->SetLineColorAlpha(kBlack, 0.6);
    Temp_DB_Line_Diff_Graph->SetMarkerColor(kRed);
    Temp_DB_Line_Diff_Graph->Draw("AP");

    c13->cd(3);
    c13->SetGrid();
    Temp_Lab_Line_Diff_Graph->SetTitle("Temp_Lab vs Time Fit Residuals");
    Temp_Lab_Line_Diff_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    Temp_Lab_Line_Diff_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    Temp_Lab_Line_Diff_Graph->GetXaxis()->SetTitleSize(0.08);
    Temp_Lab_Line_Diff_Graph->GetXaxis()->SetTitleOffset(0.88);
    Temp_Lab_Line_Diff_Graph->GetXaxis()->SetLabelSize(0.08);
    Temp_Lab_Line_Diff_Graph->GetYaxis()->SetTitle("Temperature [Deg. F]");
    Temp_Lab_Line_Diff_Graph->GetYaxis()->SetTitleSize(0.08);
    Temp_Lab_Line_Diff_Graph->GetYaxis()->SetTitleOffset(0.4);
    Temp_Lab_Line_Diff_Graph->GetYaxis()->SetLabelSize(0.08);
    Temp_Lab_Line_Diff_Graph->SetMarkerStyle(7);
    Temp_Lab_Line_Diff_Graph->SetLineColorAlpha(kBlack, 0.6);
    Temp_Lab_Line_Diff_Graph->SetMarkerColor(kRed);
    Temp_Lab_Line_Diff_Graph->Draw("AP");
    if(SavePlots){
    c13->Print("Plots/Line_Diff_Plots.png");  
    }    
 

    /*TCanvas *c21 = new TCanvas("c21","c21", 50,50,2200,400);
    c21->cd();
    c21->SetGrid();
    Ratio_vs_Temp_Graph->SetTitle("Ratio vs Temp");
    Ratio_vs_Temp_Graph->GetXaxis()->SetRangeUser(70,80);
    Ratio_vs_Temp_Graph->GetXaxis()->SetTitle("Temperature");
    Ratio_vs_Temp_Graph->GetXaxis()->SetTitleSize(0.08);
    Ratio_vs_Temp_Graph->GetXaxis()->SetTitleOffset(0.8);
    Ratio_vs_Temp_Graph->GetXaxis()->SetLabelSize(0.08);
    Ratio_vs_Temp_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
    Ratio_vs_Temp_Graph->GetYaxis()->SetTitleSize(0.08);
    Ratio_vs_Temp_Graph->GetYaxis()->SetTitleOffset(0.4);
    Ratio_vs_Temp_Graph->GetYaxis()->SetLabelSize(0.06);
    Ratio_vs_Temp_Graph->SetMarkerStyle(7);
    Ratio_vs_Temp_Graph->SetLineColorAlpha(kBlack, 0.6);
    Ratio_vs_Temp_Graph->SetMarkerColor(kRed);
    Ratio_vs_Temp_Graph->Draw("AP");
    if(SavePlots){
      c21->Print("Plots/Ratio_vs_Temp.png");
    }*/
  

    // Histograms for progression of temp correcting ratio plot
    TCanvas *c12 =new TCanvas("c12", "c12", 10,10,2000,1200);
    c12->Divide(1,3);
    c12->cd(1);
    c12->SetGrid();
    Ratio_Spread->GetXaxis()->SetLabelSize(0.08);
    Ratio_Spread->GetYaxis()->SetLabelSize(0.06);
    Ratio_Spread->Draw();

    c12->cd(2);
    c12->SetGrid();
    PreCorrRatio_Spread->GetXaxis()->SetLabelSize(0.08);
    PreCorrRatio_Spread->GetYaxis()->SetLabelSize(0.06);
    PreCorrRatio_Spread->Draw();

    c12->cd(3);  
    c12->SetGrid();
    CorrRatio_Spread->GetXaxis()->SetLabelSize(0.08);
    CorrRatio_Spread->GetYaxis()->SetLabelSize(0.06);
    CorrRatio_Spread->Draw();  
    if(SavePlots){ 
      c12->Print("Plots/CorrRatioSpreads_vs_Time.png");  
    }


    // Progression of Temperature Correction of Ratio Plot
    TCanvas *c11 =new TCanvas("c11", "c11", 10,10,2000,1200);
    c11->Divide(1,3);
    c11->cd(1);
    c11->SetGrid();
    Ratio_vs_Time_Graph->SetTitle("Ratio of Predicted to Measured ^{207}Bi Peak Over Time");
    Ratio_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    Ratio_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    Ratio_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
    Ratio_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.88);
    Ratio_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
    Ratio_vs_Time_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
    Ratio_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
    Ratio_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
    Ratio_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
    Ratio_vs_Time_Graph->GetYaxis()->SetRangeUser(0.97,1.015);
    Ratio_vs_Time_Graph->SetMarkerStyle(7);
    Ratio_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
    Ratio_vs_Time_Graph->SetMarkerColor(kRed);
    Ratio_vs_Time_Graph->Draw("AP");
    if(DrawLines){
      TLine *upper6 = new TLine(0,1.01,Time[size-1],1.01);
      upper6->SetLineColor(kRed);
      upper6->SetLineWidth(2);
      upper6->Draw();
      TLine *lower6 = new TLine(0,0.99,Time[size-1],0.99);
      lower6->SetLineColor(kRed);
      lower6->SetLineWidth(2);
      lower6->Draw();
    }

    c11->cd(2);
    c11->SetGrid();
    PreCorrRatio_DB_vs_Time_Graph->SetTitle("Ratio vs Time After 1st Temperature Correction");
    PreCorrRatio_DB_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    PreCorrRatio_DB_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    PreCorrRatio_DB_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
    PreCorrRatio_DB_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.88);
    PreCorrRatio_DB_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
    PreCorrRatio_DB_vs_Time_Graph->GetYaxis()->SetRangeUser(0.97,1.015);
    PreCorrRatio_DB_vs_Time_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
    PreCorrRatio_DB_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
    PreCorrRatio_DB_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
    PreCorrRatio_DB_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
    PreCorrRatio_DB_vs_Time_Graph->SetMarkerStyle(7);
    PreCorrRatio_DB_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
    PreCorrRatio_DB_vs_Time_Graph->SetMarkerColor(kRed);
    PreCorrRatio_DB_vs_Time_Graph->Draw("AP");
    if(DrawLines){
      TLine *upper5 = new TLine(0,1.01,Time[size-1],1.01);
      upper5->SetLineColor(kRed);
      upper5->SetLineWidth(2);
      upper5->Draw();
      TLine *lower5 = new TLine(0,0.99,Time[size-1],0.99);
      lower5->SetLineColor(kRed);
      lower5->SetLineWidth(2);
      lower5->Draw();
    }

    c11->cd(3);
    c11->SetGrid();
    CorrRatio_DB_vs_Time_Graph->SetTitle("Ratio vs Time After 2nd Temperature Correction");
    CorrRatio_DB_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    CorrRatio_DB_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    CorrRatio_DB_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
    CorrRatio_DB_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.88); 
    CorrRatio_DB_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
    CorrRatio_DB_vs_Time_Graph->GetYaxis()->SetRangeUser(0.97,1.015);
    CorrRatio_DB_vs_Time_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
    CorrRatio_DB_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
    CorrRatio_DB_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
    CorrRatio_DB_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
    CorrRatio_DB_vs_Time_Graph->SetMarkerStyle(7);
    CorrRatio_DB_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
    CorrRatio_DB_vs_Time_Graph->SetMarkerColor(kRed);
    CorrRatio_DB_vs_Time_Graph->Draw("AP");
    if(DrawLines){
      TLine *upper4 = new TLine(0,1.01,Time[size-1],1.01);
      upper4->SetLineColor(kRed);
      upper4->SetLineWidth(2);
      upper4->Draw();
      TLine *lower4 = new TLine(0,0.99,Time[size-1],0.99);
      lower4->SetLineColor(kRed);
      lower4->SetLineWidth(2);
      lower4->Draw();
    }
    if(SavePlots){
      c11->Print("Plots/CorrRatios_vs_Time.png");
    }


      // These plots are comparing ratio vs time, corrected ratio vs time, and temperature
      TCanvas *c10 =new TCanvas("c10", "c10", 10,10,2000,1200);
      c10->Divide(1,3);
      c10->cd(1);
      c10->SetGrid();
      Ratio_vs_Time_Graph->SetTitle("Ratio of Predicted to Measured ^{207}Bi Peak Over Time");
      Ratio_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
      Ratio_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
      Ratio_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
      Ratio_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.88);
      Ratio_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
      Ratio_vs_Time_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
      Ratio_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
      Ratio_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
      Ratio_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
      Ratio_vs_Time_Graph->GetYaxis()->SetRangeUser(0.97,1.015);
      Ratio_vs_Time_Graph->SetMarkerStyle(7);
      Ratio_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
      Ratio_vs_Time_Graph->SetMarkerColor(kRed);
      Ratio_vs_Time_Graph->Draw("AP");
      if(DrawLines){
        TLine *upper3 = new TLine(0,1.01,Time[size-1],1.01);
        upper3->SetLineColor(kRed);
        upper3->SetLineWidth(2);
        upper3->Draw();
        TLine *lower3 = new TLine(0,0.99,Time[size-1],0.99);
        lower3->SetLineColor(kRed);
        lower3->SetLineWidth(2);
        lower3->Draw();
      }
    
      c10->cd(2);
      c10->SetGrid();
      CorrRatio_DB_vs_Time_Graph->SetTitle("Temperature Corrected Ratio vs Time");
      CorrRatio_DB_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
      CorrRatio_DB_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
      CorrRatio_DB_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
      CorrRatio_DB_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.88); 
      CorrRatio_DB_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
      CorrRatio_DB_vs_Time_Graph->GetYaxis()->SetRangeUser(0.97,1.015);
      CorrRatio_DB_vs_Time_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
      CorrRatio_DB_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
      CorrRatio_DB_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
      CorrRatio_DB_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
      CorrRatio_DB_vs_Time_Graph->SetMarkerStyle(7);
      CorrRatio_DB_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
      CorrRatio_DB_vs_Time_Graph->SetMarkerColor(kRed);
      CorrRatio_DB_vs_Time_Graph->Draw("AP");
      if(DrawLines){
        TLine *upper2 = new TLine(0,1.01,Time[size-1],1.01);
        upper2->SetLineColor(kRed);
        upper2->SetLineWidth(2);
        upper2->Draw();
        TLine *lower2 = new TLine(0,0.99,Time[size-1],0.99);
        lower2->SetLineColor(kRed);
        lower2->SetLineWidth(2);
        lower2->Draw();
      }
     
      c10->cd(3);
      c10->SetGrid();
      T_DB_vs_Time_Graph->SetTitle("Temperature Over Time in Dark Box");
      T_DB_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
      T_DB_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
      T_DB_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
      T_DB_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.88);
      T_DB_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
      T_DB_vs_Time_Graph->GetYaxis()->SetTitle("Temperature [Deg. F]");
      T_DB_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
      T_DB_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
      T_DB_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
      T_DB_vs_Time_Graph->SetMarkerStyle(7);
      T_DB_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
      T_DB_vs_Time_Graph->SetMarkerColor(kRed);
      T_DB_vs_Time_Graph->Draw("AP");
    
      if(SavePlots){
      c10->Print("Plots/CorrRatio&Ratio&Temp_vs_Time.png");
      }
    } 
    

    // Am Chi^2/NDF and LED_Ref Chi^2/NDF Histograms
    TCanvas *c9 =new TCanvas("c9", "c9", 10,10,1950,600);
    c9->Divide(1,2);
    c9->cd(1);
    c9->SetGrid();
    Am_chi2_NDF_Spread->GetXaxis()->SetLabelSize(0.08);
    Am_chi2_NDF_Spread->GetYaxis()->SetLabelSize(0.06);
    Am_chi2_NDF_Spread->Draw();

    c9->cd(2);
    c9->SetGrid();
    LED_Ref_chi2_NDF_Spread->GetXaxis()->SetLabelSize(0.08);
    LED_Ref_chi2_NDF_Spread->GetYaxis()->SetLabelSize(0.06);
    LED_Ref_chi2_NDF_Spread->Draw();
    
    if(SavePlots){
      c9->Print("Plots/Am&LED_Ref_chi2_NDF_Spreads.png");
    }
    

    // Bi Chi^2/NDF and LED_Calo Chi^2/NDF Histograms
    TCanvas *c8 =new TCanvas("c8", "c8", 10,10,1950,600);
    c8->Divide(1,2);
    c8->cd(1);
    c8->SetGrid();
    Bi_chi2_NDF_Spread->GetXaxis()->SetLabelSize(0.08);
    Bi_chi2_NDF_Spread->GetYaxis()->SetLabelSize(0.06);
    Bi_chi2_NDF_Spread->Draw();

    c8->cd(2);
    c8->SetGrid();
    LED_Calo_chi2_NDF_Spread->GetXaxis()->SetLabelSize(0.08);
    LED_Calo_chi2_NDF_Spread->GetYaxis()->SetLabelSize(0.06);
    LED_Calo_chi2_NDF_Spread->Draw();
    
    if(SavePlots){
      c8->Print("Plots/Bi&LED_Calo_chi2_NDF_Spreads.png");
    }
  
  
    // Am Chi^2 vs Time and LED_Ref Chi^2 vs Time Plots
    TCanvas *c7 = new TCanvas("c7","c7", 10,10,1950,600);
    c7->Divide(1,2);
    c7->cd(1);
    c7->SetGrid();
    Am_vs_Time_chi2_NDF_Graph->SetTitle("Am_Chi^2/NDF vs Time");
    Am_vs_Time_chi2_NDF_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    Am_vs_Time_chi2_NDF_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    Am_vs_Time_chi2_NDF_Graph->GetXaxis()->SetTitleSize(0.08);
    Am_vs_Time_chi2_NDF_Graph->GetXaxis()->SetTitleOffset(0.88);
    Am_vs_Time_chi2_NDF_Graph->GetXaxis()->SetLabelSize(0.08);
    Am_vs_Time_chi2_NDF_Graph->GetYaxis()->SetTitle("Am_Chi^2/NDF");
    Am_vs_Time_chi2_NDF_Graph->GetYaxis()->SetTitleSize(0.08);
    Am_vs_Time_chi2_NDF_Graph->GetYaxis()->SetTitleOffset(0.4);
    Am_vs_Time_chi2_NDF_Graph->GetYaxis()->SetLabelSize(0.08);
    //Am_vs_Time_chi2_NDF_Graph->GetYaxis()->SetRangeUser(.975,1.015);
    Am_vs_Time_chi2_NDF_Graph->SetMarkerStyle(7);
    Am_vs_Time_chi2_NDF_Graph->SetMarkerColor(kBlack);
    Am_vs_Time_chi2_NDF_Graph->Draw("AP");

    c7->cd(2);
    c7->SetGrid();
    LED_Ref_vs_Time_chi2_NDF_Graph->SetTitle("LED_Ref_Chi^2/NDF vs Time");
    LED_Ref_vs_Time_chi2_NDF_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    LED_Ref_vs_Time_chi2_NDF_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    LED_Ref_vs_Time_chi2_NDF_Graph->GetXaxis()->SetTitleSize(0.08);
    LED_Ref_vs_Time_chi2_NDF_Graph->GetXaxis()->SetTitleOffset(0.88);
    LED_Ref_vs_Time_chi2_NDF_Graph->GetXaxis()->SetLabelSize(0.08);
    LED_Ref_vs_Time_chi2_NDF_Graph->GetYaxis()->SetTitle("Chi^2/NDF");
    LED_Ref_vs_Time_chi2_NDF_Graph->GetYaxis()->SetTitleSize(0.08);
    LED_Ref_vs_Time_chi2_NDF_Graph->GetYaxis()->SetTitleOffset(0.4);
    LED_Ref_vs_Time_chi2_NDF_Graph->GetYaxis()->SetLabelSize(0.08);
    //LED_Ref_vs_Time_chi2_NDF_Graph->GetYaxis()->SetRangeUser(.975,1.015);
    LED_Ref_vs_Time_chi2_NDF_Graph->SetMarkerStyle(7);
    LED_Ref_vs_Time_chi2_NDF_Graph->SetMarkerColor(kBlack);
    LED_Ref_vs_Time_chi2_NDF_Graph->Draw("AP");
    if(SavePlots){
      c7->Print("Plots/Am&LED_Ref_chi2_NDF.png");
    }    


    // Bi Chi^2/NDF vs Time and LED_Calo Chi^2/NDF vs Time Plots
    TCanvas *c6 = new TCanvas("c6","c6", 10,10,1950,600);
    c6->Divide(1,2);
    c6->cd(1);
    c6->SetGrid();
    Bi_vs_Time_chi2_NDF_Graph->SetTitle("Bi_Chi^2/NDF vs Time");
    Bi_vs_Time_chi2_NDF_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    Bi_vs_Time_chi2_NDF_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    Bi_vs_Time_chi2_NDF_Graph->GetXaxis()->SetTitleSize(0.08);
    Bi_vs_Time_chi2_NDF_Graph->GetXaxis()->SetTitleOffset(0.88);
    Bi_vs_Time_chi2_NDF_Graph->GetXaxis()->SetLabelSize(0.08);
    Bi_vs_Time_chi2_NDF_Graph->GetYaxis()->SetTitle("Chi^2/NDF");
    Bi_vs_Time_chi2_NDF_Graph->GetYaxis()->SetTitleSize(0.08);
    Bi_vs_Time_chi2_NDF_Graph->GetYaxis()->SetTitleOffset(0.4);
    Bi_vs_Time_chi2_NDF_Graph->GetYaxis()->SetLabelSize(0.08);
    //Bi_vs_Time_chi2_NDF_Graph->GetYaxis()->SetRangeUser(.975,1.015);
    Bi_vs_Time_chi2_NDF_Graph->SetMarkerStyle(7);
    Bi_vs_Time_chi2_NDF_Graph->SetMarkerColor(kBlack);
    Bi_vs_Time_chi2_NDF_Graph->Draw("AP");

    c6->cd(2);
    c6->SetGrid();
    LED_Calo_vs_Time_chi2_NDF_Graph->SetTitle("LED_Calo_Chi^2/NDF vs Time");
    LED_Calo_vs_Time_chi2_NDF_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    LED_Calo_vs_Time_chi2_NDF_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    LED_Calo_vs_Time_chi2_NDF_Graph->GetXaxis()->SetTitleSize(0.08);
    LED_Calo_vs_Time_chi2_NDF_Graph->GetXaxis()->SetTitleOffset(0.88);
    LED_Calo_vs_Time_chi2_NDF_Graph->GetXaxis()->SetLabelSize(0.08);
    LED_Calo_vs_Time_chi2_NDF_Graph->GetYaxis()->SetTitle("Chi^2/NDF");
    LED_Calo_vs_Time_chi2_NDF_Graph->GetYaxis()->SetTitleSize(0.08);
    LED_Calo_vs_Time_chi2_NDF_Graph->GetYaxis()->SetTitleOffset(0.4);
    LED_Calo_vs_Time_chi2_NDF_Graph->GetYaxis()->SetLabelSize(0.08);
    //LED_Calo_vs_Time_chi2_NDF_Graph->GetYaxis()->SetRangeUser(.975,1.015);
    LED_Calo_vs_Time_chi2_NDF_Graph->SetMarkerStyle(7);
    LED_Calo_vs_Time_chi2_NDF_Graph->SetMarkerColor(kBlack);
    LED_Calo_vs_Time_chi2_NDF_Graph->Draw("AP");
    if(SavePlots){
      c6->Print("Plots/Bi&LED_Calo_chi2_NDF.png");
    }


    // LED_Calo vs Time and LED_Ref vs Time Plots
    TCanvas *c5 = new TCanvas("c5","c5", 10,10,1950,600);  
    c5->Divide(1,2);
    c5->cd(1);
    c5->SetGrid();
    LED_Calo_vs_Time_Graph->SetTitle("LED_Calo vs Time");
    LED_Calo_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    LED_Calo_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    LED_Calo_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
    LED_Calo_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.88);
    LED_Calo_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
    LED_Calo_vs_Time_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
    LED_Calo_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
    LED_Calo_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
    LED_Calo_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
    //LED_Calo_vs_Time_Graph->GetYaxis()->SetRangeUser(.975,1.015);
    LED_Calo_vs_Time_Graph->SetMarkerStyle(7);
    LED_Calo_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
    LED_Calo_vs_Time_Graph->SetMarkerColor(kRed);
    LED_Calo_vs_Time_Graph->Draw("AP");

    c5->cd(2);
    c5->SetGrid();
    LED_Ref_vs_Time_Graph->SetTitle("LED_Ref vs Time");
    LED_Ref_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    LED_Ref_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    LED_Ref_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
    LED_Ref_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.88);
    LED_Ref_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
    LED_Ref_vs_Time_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
    LED_Ref_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
    LED_Ref_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
    LED_Ref_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
    //LED_Ref_vs_Time_Graph->GetYaxis()->SetRangeUser(.975,1.015);
    LED_Ref_vs_Time_Graph->SetMarkerStyle(7);
    LED_Ref_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
    LED_Ref_vs_Time_Graph->SetMarkerColor(kRed);
    LED_Ref_vs_Time_Graph->Draw("AP");
    if(SavePlots){ 
      c5->Print("Plots/LEDs.png");
    }


    // Bi vs Time and Am vs Time Plots
    TCanvas *c4 = new TCanvas("c4","c4", 10,10,1950,600);  
    c4->Divide(1,2);
    c4->cd(1);
    c4->SetGrid();
    Bi_vs_Time_Graph->SetTitle("Bi vs Time");
    Bi_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    Bi_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    Bi_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
    Bi_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.88);
    Bi_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
    Bi_vs_Time_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
    Bi_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
    Bi_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
    Bi_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
    //Bi_vs_Time_Graph->GetYaxis()->SetRangeUser(.975,1.015);
    Bi_vs_Time_Graph->SetMarkerStyle(7);
    Bi_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
    Bi_vs_Time_Graph->SetMarkerColor(kRed);
    Bi_vs_Time_Graph->Draw("AP");

    c4->cd(2);
    c4->SetGrid();
    Am_vs_Time_Graph->SetTitle("Am vs Time");
    Am_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    Am_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    Am_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
    Am_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.88);
    Am_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
    Am_vs_Time_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
    Am_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
    Am_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
    Am_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
    //Am_vs_Time_Graph->GetYaxis()->SetRangeUser(.975,1.015);
    Am_vs_Time_Graph->SetMarkerStyle(7);
    Am_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
    Am_vs_Time_Graph->SetMarkerColor(kRed);
    Am_vs_Time_Graph->Draw("AP");
    if(SavePlots){ 
      c4->Print("Plots/Sources.png");
    }


    // Am vs Time and LED_Ref vs Time Plots
    TCanvas *c3 = new TCanvas("c3","c3", 10,10,1950,600);  
    c3->Divide(1,2);
    c3->cd(1);
    c3->SetGrid();
    Am_vs_Time_Graph->SetTitle("Am vs Time");
    Am_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    Am_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    Am_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
    Am_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.88);
    Am_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
    Am_vs_Time_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
    Am_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
    Am_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
    Am_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
    //Am_vs_Time_Graph->GetYaxis()->SetRangeUser(.975,1.015);
    Am_vs_Time_Graph->SetMarkerStyle(7);
    Am_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
    Am_vs_Time_Graph->SetMarkerColor(kRed);
    Am_vs_Time_Graph->Draw("AP");

    c3->cd(2);
    c3->SetGrid();
    LED_Ref_vs_Time_Graph->SetTitle("LED_Ref vs Time");
    LED_Ref_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    LED_Ref_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    LED_Ref_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
    LED_Ref_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.88);
    LED_Ref_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
    LED_Ref_vs_Time_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
    LED_Ref_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
    LED_Ref_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
    LED_Ref_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
    //LED_Ref_vs_Time_Graph->GetYaxis()->SetRangeUser(.975,1.015);
    LED_Ref_vs_Time_Graph->SetMarkerStyle(7);
    LED_Ref_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
    LED_Ref_vs_Time_Graph->SetMarkerColor(kRed);
    LED_Ref_vs_Time_Graph->Draw("AP");
    if(SavePlots){ 
      c3->Print("Plots/Am&LED_Ref.png");
    }


    // Bi_t vs Time and LED_Calo vs Time Plots 
    TCanvas *c2 = new TCanvas("c2","c2", 10,10,1950,600);  
    c2->Divide(1,2);
    c2->cd(1);
    c2->SetGrid();
    Bi_vs_Time_Graph->SetTitle("Bi vs Time");
    Bi_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    Bi_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    Bi_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
    Bi_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.88);
    Bi_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
    Bi_vs_Time_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
    Bi_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
    Bi_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
    Bi_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
    //Bi_vs_Time_Graph->GetYaxis()->SetRangeUser(.975,1.015);
    Bi_vs_Time_Graph->SetMarkerStyle(7);
    Bi_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
    Bi_vs_Time_Graph->SetMarkerColor(kRed);
    Bi_vs_Time_Graph->Draw("AP");

    c2->cd(2);
    c2->SetGrid();
    LED_Calo_vs_Time_Graph->SetTitle("LED_Calo vs Time");
    LED_Calo_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    LED_Calo_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    LED_Calo_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
    LED_Calo_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.88);
    LED_Calo_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
    LED_Calo_vs_Time_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
    LED_Calo_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
    LED_Calo_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
    LED_Calo_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
    //LED_Calo_vs_Time_Graph->GetYaxis()->SetRangeUser(.975,1.015);
    LED_Calo_vs_Time_Graph->SetMarkerStyle(7);
    LED_Calo_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
    LED_Calo_vs_Time_Graph->SetMarkerColor(kRed);
    LED_Calo_vs_Time_Graph->Draw("AP");
    if(SavePlots){ 
      c2->Print("Plots/Bi&LED_Calo.png");
    }


    // Money plot, Bi_Pred plot, and Bi_t Plot
    TCanvas *c1 =new TCanvas("c1", "c1", 10,10,2000,1200);
    c1->Divide(1,3);
    c1->cd(1);
    c1->SetGrid();
    Ratio_vs_Time_Graph->SetTitle("Ratio of Predicted to Measured ^{207}Bi Peak Over Time");
    Ratio_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    Ratio_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    Ratio_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
    Ratio_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.88);
    Ratio_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
    Ratio_vs_Time_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
    Ratio_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
    Ratio_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
    Ratio_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
    Ratio_vs_Time_Graph->GetYaxis()->SetRangeUser(0.97,1.015);
    Ratio_vs_Time_Graph->SetMarkerStyle(7);
    Ratio_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
    Ratio_vs_Time_Graph->SetMarkerColor(kRed);
    Ratio_vs_Time_Graph->Draw("AP");
    if(DrawLines){
      TLine *upper1 = new TLine(0,1.01,Time[size-1],1.01);
      upper1->SetLineColor(kRed);
      upper1->SetLineWidth(2);
      upper1->Draw();
      TLine *lower1 = new TLine(0,0.99,Time[size-1],0.99);
      lower1->SetLineColor(kRed);
      lower1->SetLineWidth(2);
      lower1->Draw();
    }    

    c1->cd(2);
    c1->SetGrid();
    Bi_Pred_vs_Time_Graph->SetTitle("Bi_Pred vs Time");
    Bi_Pred_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    Bi_Pred_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    Bi_Pred_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
    Bi_Pred_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.88);
    Bi_Pred_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
    //Bi_Pred_vs_Time_Graph->GetYaxis()->SetRangeUser(0.97,1.015);
    Bi_Pred_vs_Time_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
    Bi_Pred_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
    Bi_Pred_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
    Bi_Pred_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
    Bi_Pred_vs_Time_Graph->SetMarkerStyle(7);
    Bi_Pred_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
    Bi_Pred_vs_Time_Graph->SetMarkerColor(kRed);
    Bi_Pred_vs_Time_Graph->Draw("AP");

    c1->cd(3);
    c1->SetGrid();
    Bi_vs_Time_Graph->SetTitle("Bi vs Time");
    Bi_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[size-1]);
    Bi_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    Bi_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
    Bi_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.88); 
    Bi_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
    //Bi_vs_Time_Graph->GetYaxis()->SetRangeUser(0.97,1.015);
    Bi_vs_Time_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
    Bi_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
    Bi_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
    Bi_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
    Bi_vs_Time_Graph->SetMarkerStyle(7);
    Bi_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
    Bi_vs_Time_Graph->SetMarkerColor(kRed);
    Bi_vs_Time_Graph->Draw("AP");
    if(SavePlots){
      c1->Print("Plots/CorrRatios_vs_Time.png");
    }    
   
 
  }

  return(23);
}
