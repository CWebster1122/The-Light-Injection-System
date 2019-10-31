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

//***********************
// Flags
//***********************

Bool_t Debug = 1; // Toggle outpt of debugging comments
Bool_t TempInt = 1; // Toggle integration of temperature readings (must be set to 0 if no temperature data available)
Bool_t DrawLines = 1; // Toggle drawing of various plot lines
Bool_t DrawPlots = 1; // Toggle drawing of plots

//***********************
// Variable initializations
//***********************
Double_t average(Double_t numbers[], int numb);
// Parameter file variables
Double_t hmev_mean, hmev_mean_err, hmev_sigma, hmev_sigma_err, hmev_chi2, hmev_chi2_NDF,
  mev_mean, mev_mean_err, mev_sigma, mev_sigma_err, mev_chi2, mev_chi2_NDF,
  led_mean, led_mean_err, led_sigma, led_sigma_err, led_calo_chi2, led_calo_chi2_NDF,
  americium1_mean, americium1_mean_err, americium1_sigma, americium1_sigma_err, am_chi2, am_chi2_NDF,
  led2_mean, led2_mean_err, led2_sigma, led2_sigma_err, led_ref_chi2, led_ref_chi2_NDF, 
//  americium4_mean, americium4_mean_err, americium4_sigma, americium4_sigma_err, // Used if a 3rd reference OM is present
//  led3_mean, led3_mean_err, led3_sigma, led3_sigma_err, // Used if a 3rd reference OM is present
  time_data;

// Temperature file variables
Double_t time_temp, T_pb, T_db, T_lab, T_lab2;

// Parameter file vectors
std::vector<double> HMeV_Mean, HMeV_Mean_Err, HMeV_Sigma, HMeV_Sigma_Err, HMeV_chi2, HMeV_chi2_NDF,
  MeV_Mean, MeV_Mean_Err, MeV_Sigma, MeV_Sigma_Err, MeV_chi2, MeV_chi2_NDF,
  LED_Mean, LED_Mean_Err, LED_Sigma, LED_Sigma_Err, LED_Calo_chi2, LED_Calo_chi2_NDF,
  Americium1_Mean, Americium1_Mean_Err, Americium1_Sigma, Americium1_Sigma_Err, Am_chi2, Am_chi2_NDF,
  LED2_Mean, LED2_Mean_Err, LED2_Sigma, LED2_Sigma_Err, LED_Ref_chi2, LED_Ref_chi2_NDF,
  Americium4_Mean, Americium4_Mean_Err, Americium4_Sigma, Americium4_Sigma_Err, // Used if a 3rd reference OM is present
  LED3_Mean, LED3_Mean_Err, LED3_Sigma, LED3_Sigma_Err, // Used if a 3rd reference OM is present
  Time_Data;

// Temperature file vectors
std::vector<double> Time_Temp, T_PB, T_DB, T_LAB, T_LAB2;

// Final time series vectors
std::vector<double> Bi_vs_Time, Bi_vs_Time_Err,
  Bi_Pred_vs_Time, Bi_Pred_vs_Time_Err, 
  Am_vs_Time, Am_vs_Time_Err,
  LED_Calo_vs_Time, LED_Calo_vs_Time_Err,
  LED_Ref_vs_Time, LED_Ref_vs_Time_Err,
  Ratio_vs_Time, Ratio_vs_Time_Err,
  T_PB_vs_Time, T_DB_vs_Time, T_Lab_vs_Time, T_Lab2_vs_Time, T_DB_sigma, T_Lab2_sigma;

//***********************
//***********************
//    MAIN FUNCTION    //     
//***********************
//***********************
int AvgTemp()
{
  //***********************
  // Read all parameters into vectors
  //***********************

  // Read in file of all fit parameters
  ifstream in;
  in.open("Parameters1.dat");
  while(in
	>> hmev_mean >> hmev_mean_err >> hmev_sigma >> hmev_sigma_err
	>> mev_mean >> mev_mean_err >> mev_sigma >> mev_sigma_err >> mev_chi2 >> mev_chi2_NDF
	>> led_mean >> led_mean_err >> led_sigma >> led_sigma_err >> led_calo_chi2 >> led_calo_chi2_NDF
	>> americium1_mean >> americium1_mean_err >> americium1_sigma >> americium1_sigma_err >> am_chi2 >> am_chi2_NDF
	>> led2_mean >> led2_mean_err >> led2_sigma >> led2_sigma_err >> led_ref_chi2 >> led_ref_chi2_NDF
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

    LED_Mean.push_back(led_mean);
    LED_Mean_Err.push_back(led_mean_err);
    LED_Sigma.push_back(led_sigma);
    LED_Sigma_Err.push_back(led_sigma_err);

    // Ref OM Peaks
    Americium1_Mean.push_back(americium1_mean);
    Americium1_Mean_Err.push_back(americium1_mean_err);
    Americium1_Sigma.push_back(americium1_sigma);
    Americium1_Sigma_Err.push_back(americium1_sigma_err);

    LED2_Mean.push_back(led2_mean);
    LED2_Mean_Err.push_back(led2_mean_err);
    LED2_Sigma.push_back(led2_sigma);
    LED2_Sigma_Err.push_back(led2_sigma_err);

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
      Time_Temp.push_back(time_temp); // account for hour offset due to time change (was - 3600.)
      T_PB.push_back(T_pb);
      T_DB.push_back(T_db);
      T_LAB.push_back(T_lab);
      T_LAB2.push_back(T_lab2);
    }
    in2.close();
  }
  const int nentries2 = Time_Temp.size();
  printf("Temperature Vector Size = %d\n",nentries2);
  

  //***********************
  // Perform prediction calculations
  //***********************
  
  // Start values (i.e. calibration run values)
  Double_t Bi_0 = MeV_Mean[0];
  Double_t Bi_0_err = MeV_Mean_Err[0];
  Double_t Bi2_0 = HMeV_Mean[0];
  Double_t Bi2_0_err = HMeV_Mean_Err[0];
  Double_t LED_Calo_0 = LED_Mean[0];
  Double_t LED_Calo_0_err = LED_Mean_Err[0];
  Double_t LED_Ref_0 = LED2_Mean[0];
  Double_t LED_Ref_0_err = LED2_Mean_Err[0];
  Double_t Am_0 = Americium1_Mean[0];
  Double_t Am_0_err = Americium1_Mean_Err[0];

  // Later values for LED and Am ONLY
  Double_t Bi_t, Bi_t_err, LED_Calo_t, LED_Calo_t_err, LED_Ref_t, LED_Ref_t_err, Am_t, Am_t_err, LED_Ref_Pred, LED_Ref_Pred_err, LED_Calo_Pred, LED_Calo_Pred_err, Bi_Pred, Bi_Pred_err, Ratio, Ratio_err;
  Double_t Bi_last = Bi_0;
  Double_t LED_Calo_last = LED_Calo_0;
  int counter1 = 0, counter2 = 0, rejections = 0;
  for(int i=0; i<nentries; i++){
    Bi_t           = MeV_Mean[i];
    Bi_t_err       = MeV_Mean_Err[i];
    LED_Calo_t     = LED_Mean[i];
    LED_Calo_t_err = LED_Mean_Err[i];
    LED_Ref_t      = LED2_Mean[i];
    LED_Ref_t_err  = LED2_Mean_Err[i];
    Am_t           = Americium1_Mean[i];
    Am_t_err       = Americium1_Mean_Err[i];
    
    // Skip runs with large Bi peak deviations unaccounted for in LED peak
    if( (TMath::Abs(Bi_t - Bi_last) > 10*Bi_t_err) &&
	!(TMath::Abs(LED_Calo_t - LED_Calo_last) > 10*LED_Calo_t_err)){
      if(Debug) {
	printf("Reese: Skipping Run %d (timestamp: %f): Bi mean value = %f while LED mean value = %f\n",i, Time_Data.at(i), Bi_t, LED_Calo_t);
	}
      rejections++;
      continue;
    }
    
    Bi_last = Bi_t; // Update comparison values
    LED_Calo_last = LED_Calo_t;
    //Am_t = Am_0;
    // Debugging outputs
//    if(Debug) printf("%f: Bi mean = %f LED mean  = %f", i, Time_Data.at(i), Bi_t, LED_Calo_t);    

    // Use Am to predict Ref LED (correct for drifts)
    LED_Ref_Pred = LED_Ref_0 * (Am_t/Am_0);
    LED_Ref_Pred_err = sqrt(pow(LED_Ref_0_err*Am_t/Am_0,2) 
			    + pow(LED_Ref_0*Am_t_err/Am_0,2) 
			    + pow(LED_Ref_0*Am_t*Am_0_err/(Am_0*Am_0),2));
 //   if(Debug) printf("LED pred: %f*%f = %f\n", Am_t, Am_0, (Am_t/Am_0));
    // Use Ref LED to predict Calo LED
    LED_Calo_Pred = LED_Calo_t * (LED_Ref_Pred/LED_Ref_t);
    LED_Calo_Pred_err = sqrt(pow(LED_Calo_t_err*LED_Ref_Pred/LED_Ref_t,2) 
			     + pow(LED_Calo_t*LED_Ref_Pred_err/LED_Ref_t,2) 
			     + pow(LED_Calo_t*LED_Ref_Pred*LED_Ref_t_err/(LED_Ref_t*LED_Ref_t),2));
    
    // Use Calo LED to predict Bi
    Bi_Pred = Bi_0 * (LED_Calo_Pred/LED_Calo_0); /*(Bi_0 * LED_Calo_t * LED_Ref_0 *Am_t)/(LED_Calo_0 * LED_Ref_t *Am_0);*/ //just trying a more explicit implementation of the prediction algorithm, but same effect
    //if (Debug) printf("(%f*%f*%f*%f)/(%f*%f*%f)\n", Bi_0, LED_Calo_t, LED_Ref_0, Am_t, LED_Calo_0, LED_Ref_t, Am_0);
    Bi_Pred_err = sqrt(pow(Bi_0_err*LED_Calo_Pred/LED_Calo_0,2) 
		       + pow(Bi_0*LED_Calo_Pred_err/LED_Calo_0,2) 
		       + pow(Bi_0*LED_Calo_Pred*LED_Calo_0_err/(LED_Calo_0*LED_Calo_0),2));

    // Calculate the ratio of predicted to measured Bi
    Ratio = Bi_Pred/Bi_t;
    Ratio_err = sqrt(pow(Bi_Pred_err/Bi_t,2) + pow(Bi_Pred*Bi_t_err/(Bi_t*Bi_t),2));

    // Cross check with full calculation 
    Double_t test = (Bi_0*LED_Calo_t/LED_Calo_0*LED_Ref_0/LED_Ref_t*Am_t/Am_0);

    //Fill vectors
    Bi_vs_Time.push_back(Bi_t);
    Bi_vs_Time_Err.push_back(Bi_t_err);
    LED_Calo_vs_Time.push_back(LED_Calo_t);
    LED_Calo_vs_Time_Err.push_back(LED_Calo_t_err);
    LED_Ref_vs_Time.push_back(LED_Ref_t);
    LED_Ref_vs_Time_Err.push_back(LED_Ref_t_err);
    Am_vs_Time.push_back(Am_t);
    Am_vs_Time_Err.push_back(Am_t_err);

    Bi_Pred_vs_Time.push_back(Bi_Pred);
    Bi_Pred_vs_Time_Err.push_back(Bi_Pred_err);
    Ratio_vs_Time.push_back(Ratio);
    Ratio_vs_Time_Err.push_back(Ratio_err);

    //***********************
    // Assign temperatures to each run
    //***********************

    // Search through temperature vector for temp values at closest matching time values
    if(TempInt){
      int j=0;
      Double_t T_PB_vs_Time_avg = 0; 
      Double_t T_DB_vs_Time_avg = 0; 
      Double_t T_Lab_vs_Time_avg = 0; 
      Double_t T_Lab2_vs_Time_avg = 0;
 
      while(j < nentries2){
	if(TMath::Abs(Time_Data[i] - Time_Temp[j]) <= 600){ // Temperature readings taken every 600 seconds
	  printf("Temp at closest time found:   index = %d   |   time = %f   |   temp_time = %f   |   temp = %f\n",i,Time_Data[i],Time_Temp[j],T_DB[j]); // Print statement to check
	  T_PB_vs_Time_avg = (T_PB[j]+T_PB[j-1]+T_PB[j-2]+T_PB[j-3])/4;
	  T_DB_vs_Time_avg = (T_DB[j]+T_DB[j-1]+T_DB[j-2]+T_DB[j-3])/4;
	  T_Lab_vs_Time_avg = (T_LAB[j]+T_LAB[j-1]+T_LAB[j-2]+T_LAB[j-3])/4;
	  T_Lab2_vs_Time_avg = (T_LAB2[j]+T_LAB2[j-1]+T_LAB2[j-2]+T_LAB2[j-3])/4;
	  
	  T_DB_sigma.push_back(sqrt((pow(T_DB[j]-T_DB_vs_Time_avg,2)+pow(T_DB[j-1]-T_DB_vs_Time_avg,2)+pow(T_DB[j-2]-T_DB_vs_Time_avg,2)+pow(T_DB[j-3]-T_DB_vs_Time_avg,2))/4));
	  T_Lab2_sigma.push_back(sqrt((pow(T_LAB2[j]-T_Lab2_vs_Time_avg,2)+pow(T_LAB2[j-1]-T_Lab2_vs_Time_avg,2)+pow(T_LAB2[j-2]-T_Lab2_vs_Time_avg,2)+pow(T_LAB2[j-3]-T_Lab2_vs_Time_avg,2))/4));	

	  T_PB_vs_Time.push_back(T_PB_vs_Time_avg);
	  T_DB_vs_Time.push_back(T_DB_vs_Time_avg);
	  T_Lab_vs_Time.push_back(T_Lab_vs_Time_avg);
	  T_Lab2_vs_Time.push_back(T_Lab2_vs_Time_avg);
	  
	  cout << " T_DB_avg: " << T_DB_vs_Time_avg << endl;
	
	  T_PB_vs_Time_avg = 0;
	  T_DB_vs_Time_avg = 0;
	  T_Lab_vs_Time_avg = 0;
	  T_Lab2_vs_Time_avg = 0;
	  counter2++;
	  
	  if(Debug){
	  	cout << " T_DB[j-1]: " << T_DB[j-1] << " T_DB: " << T_DB_vs_Time[0] << " T_Lab2[j-2]: " << T_LAB2[j-2] << " T_Lab2: " << T_Lab2_vs_Time[0] << endl;
	  }
	  break;
	}else{
	  j++;
	}
      }
    //printf("The vector size is T_PB_vs_Time %lu and T_DB_vs_Time %lu and T_Lab_vs_Time %lu and T_Lab2_vs_Time %lu\n", T_PB_vs_Time.size(), T_DB_vs_Time.size(), T_Lab_vs_Time.size(), T_Lab2_vs_Time.size());
    
    }
    counter1++;  
  }
  printf("The vector size is T_PB_vs_Time %lu and T_DB_vs_Time %lu and T_Lab_vs_Time %lu and T_Lab2_vs_Time %lu\n", T_PB_vs_Time.size(), T_DB_vs_Time.size(), T_Lab_vs_Time.size(), T_Lab2_vs_Time.size());
  const int size = Ratio_vs_Time.size();
  if(Debug){
    printf("Final Parameter Vector Size = %d\n",size);
    printf("Final Temperature Vector Size = %lu\n",T_DB_vs_Time.size());
    printf("counter 1 = %d  |  counter 2 = %d\n",counter1, counter2);
  }
  const int nomatch = size-T_DB_vs_Time.size();
  const int newsize = size-nomatch;
  cout << "New Size  " << size << endl; 
  //***********************
  // Build final arrays used for TGraphs
  //***********************
   
  Double_t Time[newsize], Time_Err[newsize],
    Bi_vs_Time_Array[newsize], Bi_vs_Time_Err_Array[newsize],
    LED_Calo_vs_Time_Array[newsize], LED_Calo_vs_Time_Err_Array[newsize],
    Am_vs_Time_Array[newsize], Am_vs_Time_Err_Array[newsize], 
    LED_Ref_vs_Time_Array[newsize], LED_Ref_vs_Time_Err_Array[newsize], 
    Bi_Pred_vs_Time_Array[newsize], Bi_Pred_vs_Time_Err_Array[newsize],
    Ratio_vs_Time_Array[newsize], Ratio_vs_Time_Err_Array[newsize],
    // Temperature correlation arrays (all use the same error array with constant value)
    T_PB_vs_Time_Array[newsize], T_DB_vs_Time_Array[newsize],
    T_Lab_vs_Time_Array[newsize], T_Lab2_vs_Time_Array[newsize],
    T_Err_vs_Time_Array[newsize], x, y, T_Lab2_Min, T_Lab2_Max, T_DB_Min, T_DB_Max;
  // Loop to fill arrays from vectors
  x=0;
  y=0;
  T_Lab2_Min = 85;
  T_Lab2_Max = 0;
  T_DB_Min = 85;
  T_DB_Max = 0;
  for(int i=0;i<newsize;i++){
    Bi_vs_Time_Array[i] = Bi_vs_Time[i];
    Bi_vs_Time_Err_Array[i] = Bi_vs_Time_Err[i];
    LED_Calo_vs_Time_Array[i] = LED_Calo_vs_Time[i];
    LED_Calo_vs_Time_Err_Array[i] = LED_Calo_vs_Time_Err[i];
    LED_Ref_vs_Time_Array[i] = LED_Ref_vs_Time[i];
    LED_Ref_vs_Time_Err_Array[i] = LED_Ref_vs_Time_Err[i];
    Am_vs_Time_Array[i] = Am_vs_Time[i];
    Am_vs_Time_Err_Array[i] = Am_vs_Time_Err[i];
    Time[i] = (Time_Data[i]-Time_Data[0])/86400.; // Set time as days since the first run
    Time_Err[i] = 0.5*30./1440.;  // Error = 1/2 * run duration in days (~30 min run = 1/48 days)
    
    Bi_Pred_vs_Time_Array[i] = Bi_Pred_vs_Time[i];
    Bi_Pred_vs_Time_Err_Array[i] = Bi_Pred_vs_Time_Err[i];
    Ratio_vs_Time_Array[i] = Ratio_vs_Time[i];
    Ratio_vs_Time_Err_Array[i] = Ratio_vs_Time_Err[i];

    if(TempInt){
      T_PB_vs_Time_Array[i] = T_PB_vs_Time[i];
      T_DB_vs_Time_Array[i] = T_DB_vs_Time[i];
      T_Lab_vs_Time_Array[i] = T_Lab_vs_Time[i];
      T_Lab2_vs_Time_Array[i] = T_Lab2_vs_Time[i];
      T_Err_vs_Time_Array[i] = 0.1; // Constant temperature logger error (0.1 deg F)

        x = T_Lab2_vs_Time_Array[i]; 
        if(x>T_Lab2_Max and x<85){
        T_Lab2_Max=x;
	}
	if(x<T_Lab2_Min and x>1){
	T_Lab2_Min=x;
        }
	
	y = T_DB_vs_Time_Array[i]; 
        if(y>T_DB_Max){
        T_DB_Max=y;
        }    
        if(y<T_DB_Min and y>1){
        T_DB_Min=y;
      }
    }
  }
  cout << "The Max DB Temperature is: " << T_DB_Max << endl;
  cout << "The Min DB Temperature is: " << T_DB_Min << endl;
  cout << "The Max Lab Temperature is: " << T_Lab2_Max << endl;
  cout << "The Min Lab Temperature is: " << T_Lab2_Min << endl;

  
  //***********************
  // Creating graphs and plotting
  //***********************
  
  // Basic time evolution of different peaks
  TGraphErrors *Bi_vs_Time_Graph = new TGraphErrors(newsize,Time,Bi_vs_Time_Array,Time_Err,Bi_vs_Time_Err_Array);
  TGraphErrors *LED_Calo_vs_Time_Graph = new TGraphErrors(newsize,Time,LED_Calo_vs_Time_Array,Time_Err,LED_Calo_vs_Time_Err_Array);
  TGraphErrors *LED_Ref_vs_Time_Graph = new TGraphErrors(newsize,Time,LED_Ref_vs_Time_Array,Time_Err,LED_Ref_vs_Time_Err_Array);
  TGraphErrors *Am_vs_Time_Graph = new TGraphErrors(newsize,Time,Am_vs_Time_Array,Time_Err,Am_vs_Time_Err_Array);
  // Money plots: time evolution of prediction and of ratio between predicted and measured values
  TGraphErrors *Bi_Pred_vs_Time_Graph = new TGraphErrors(newsize,Time,Bi_Pred_vs_Time_Array,Time_Err,Bi_Pred_vs_Time_Err_Array);
  TGraphErrors *Ratio_vs_Time_Graph = new TGraphErrors(newsize, Time, Ratio_vs_Time_Array,Time_Err,Ratio_vs_Time_Err_Array);

    // Temperature over time 
    TGraphErrors *T_PB_vs_Time_Graph = new TGraphErrors(newsize,Time,T_PB_vs_Time_Array,Time_Err,T_Err_vs_Time_Array);
    TGraphErrors *T_DB_vs_Time_Graph = new TGraphErrors(newsize,Time,T_DB_vs_Time_Array,Time_Err,T_Err_vs_Time_Array);
    TGraphErrors *T_Lab_vs_Time_Graph = new TGraphErrors(newsize,Time,T_Lab_vs_Time_Array,Time_Err,T_Err_vs_Time_Array);
    TGraphErrors *T_Lab2_vs_Time_Graph = new TGraphErrors(newsize,Time,T_Lab2_vs_Time_Array,Time_Err,T_Err_vs_Time_Array);
   
    // Temperature correlations (more can be added)
    TGraphErrors *Bi_vs_Temp_Graph = new TGraphErrors(newsize,T_DB_vs_Time_Array,Bi_vs_Time_Array,T_Err_vs_Time_Array,Bi_vs_Time_Err_Array);
    TGraphErrors *Am_vs_Temp_Graph = new TGraphErrors(newsize,T_Lab2_vs_Time_Array,Am_vs_Time_Array,T_Err_vs_Time_Array,Am_vs_Time_Err_Array);
    TGraphErrors *LED_Ref_vs_Temp_Graph = new TGraphErrors(newsize,T_Lab2_vs_Time_Array,LED_Ref_vs_Time_Array,T_Err_vs_Time_Array,LED_Ref_vs_Time_Err_Array);
    TGraphErrors *LED_Calo_vs_Temp_Graph = new TGraphErrors(newsize, T_Lab2_vs_Time_Array, LED_Calo_vs_Time_Array, T_Err_vs_Time_Array, LED_Calo_vs_Time_Err_Array);


  
  /****** Fits for Temp Correction ******/

  //Double_t T_DB_Min = TMath::MinElement(newsize,T_DB_vs_Time_Graph->GetY());
  //cout << "Minimum DB Temp: " << T_DB_Min << endl;
  TF1 *Am_Temp = new TF1("Am_Temp","[0]*x+[1]",T_Lab2_Min,T_Lab2_Max);
  Am_Temp->SetParameter(0,-0.33);
  Am_Temp->SetParLimits(0, -0.37, -0.25);
  Am_Temp->SetParameter(1,126.);
  //Am_Temp->SetParLimits(1, 124.5, 126.5);
  Am_Temp->SetLineColor(kRed);
  Am_Temp->SetLineWidth(4);
  Am_vs_Temp_Graph->Fit("Am_Temp", "RSB");
  Double_t Am_Temp_slope = Am_Temp->GetParameter(0);
  Double_t Am_Temp_slope_err = Am_Temp->GetParError(0);
  Double_t Am_Temp_const = Am_Temp->GetParameter(1);
  Double_t Am_Temp_const_err = Am_Temp->GetParError(1);
  Double_t Am_Temp_chi2 = Am_Temp->GetChisquare();
  Double_t Am_Temp_NDF = Am_Temp->GetNDF();
  Double_t Am_Temp_r = Am_vs_Temp_Graph->GetCorrelationFactor(); 
  Double_t Am_Temp_xMean = Am_vs_Temp_Graph->GetMean(1); 
  Double_t Am_Temp_yMean = Am_vs_Temp_Graph->GetMean(2); 


  TF1 *Bi_Temp = new TF1("Bi_Temp","[0]*x+[1]",T_DB_Min,T_DB_Max);
  Bi_Temp->SetParameter(0,-0.27);
  Bi_Temp->SetParameter(1,151.);
  Bi_Temp->SetParLimits(0, -0.37, -0.25);
  //Bi_Temp->SetParLimits(1, 150., 151.5);
  Bi_Temp->SetLineColor(kRed);
  Bi_Temp->SetLineWidth(4);
  Bi_vs_Temp_Graph->Fit("Bi_Temp", "RS");
  Double_t Bi_Temp_slope = Bi_Temp->GetParameter(0);
  Double_t Bi_Temp_slope_err = Bi_Temp->GetParError(0);
  Double_t Bi_Temp_const = Bi_Temp->GetParameter(1);
  Double_t Bi_Temp_const_err = Bi_Temp->GetParError(1);
  Double_t Bi_Temp_chi2 = Bi_Temp->GetChisquare();
  Double_t Bi_Temp_NDF = Bi_Temp->GetNDF();
  Double_t Bi_Temp_r = Bi_vs_Temp_Graph->GetCorrelationFactor();
  Double_t Bi_Temp_xMean = Bi_vs_Temp_Graph->GetMean(1);
  Double_t Bi_Temp_yMean = Bi_vs_Temp_Graph->GetMean(2);
 
  TF1 *LED_Calo_Temp = new TF1("LED_Calo_Temp","[0]*x+[1]",T_Lab2_Min,T_Lab2_Max);
  LED_Calo_Temp->SetParameter(0,-3.5);
  LED_Calo_Temp->SetParLimits(0, -5., -3.);
  LED_Calo_Temp->SetParameter(1,500.); //320
  //LED_Calo_Temp->SetParLimits(1, 500., 520.);
  LED_Calo_Temp->SetLineColor(kRed);
  LED_Calo_Temp->SetLineWidth(4);
  LED_Calo_vs_Temp_Graph->Fit("LED_Calo_Temp", "RS");
  Double_t LED_Calo_Temp_slope = LED_Calo_Temp->GetParameter(0);
  Double_t LED_Calo_Temp_slope_err = LED_Calo_Temp->GetParError(0);
  Double_t LED_Calo_Temp_const = LED_Calo_Temp->GetParameter(1);
  Double_t LED_Calo_Temp_const_err = LED_Calo_Temp->GetParError(1);
  Double_t LED_Calo_Temp_chi2 = LED_Calo_Temp->GetChisquare();
  Double_t LED_Calo_Temp_NDF = LED_Calo_Temp->GetNDF();
  Double_t LED_Calo_Temp_r = LED_Calo_vs_Temp_Graph->GetCorrelationFactor(); 
  Double_t LED_Calo_Temp_xMean = LED_Calo_vs_Temp_Graph->GetMean(1);
  Double_t LED_Calo_Temp_yMean = LED_Calo_vs_Temp_Graph->GetMean(2); 


  TF1 *LED_Ref_Temp = new TF1("LED_Ref_Temp","[0]*x+[1]",T_Lab2_Min,T_Lab2_Max);
  LED_Ref_Temp->SetParameter(0,-2.5);
  LED_Ref_Temp->SetParameter(1,294.); //298
  LED_Ref_Temp->SetParLimits(0, -5., -2.);
  //LED_Ref_Temp->SetParLimits(1, 288., 298.);
  LED_Ref_Temp->SetLineColor(kRed);
  LED_Ref_Temp->SetLineWidth(4);
  LED_Ref_vs_Temp_Graph->Fit("LED_Ref_Temp", "RS");
  Double_t LED_Ref_Temp_slope = LED_Ref_Temp->GetParameter(0);
  Double_t LED_Ref_Temp_slope_err = LED_Ref_Temp->GetParError(0);
  Double_t LED_Ref_Temp_const = LED_Ref_Temp->GetParameter(1);
  Double_t LED_Ref_Temp_const_err = LED_Ref_Temp->GetParError(1);
  Double_t LED_Ref_Temp_chi2 = LED_Ref_Temp->GetChisquare();
  Double_t LED_Ref_Temp_NDF = LED_Ref_Temp->GetNDF();
  Double_t LED_Ref_Temp_r = LED_Ref_vs_Temp_Graph->GetCorrelationFactor();
  Double_t LED_Ref_Temp_xMean = LED_Ref_vs_Temp_Graph->GetMean(1);
  Double_t LED_Ref_Temp_yMean = LED_Ref_vs_Temp_Graph->GetMean(2);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);

  cout << " *DB*  Bi_Temp_chi2: " << Bi_Temp_chi2 << " Bi_Temp_NDF: " << Bi_Temp_NDF << " Bi_Temp_chi2/Bi_Temp_NDF: " << Bi_Temp_chi2/Bi_Temp_NDF << " Bi_Temp_r: " << Bi_Temp_r << endl;
  cout << " *Lab2*  Am_Temp_chi2: " << Am_Temp_chi2 << " Am_Temp_NDF: " << Am_Temp_NDF << " Am_Temp_chi2/Am_Temp_NDF: " << Am_Temp_chi2/Am_Temp_NDF << " Am_Temp_r: " << Am_Temp_r << endl;
  cout << " *Lab2*  LED_Calo_Temp_chi2: " << LED_Calo_Temp_chi2 << " LED_Calo_Temp_NDF: " << LED_Calo_Temp_NDF << " LED_Calo_Temp_chi2/LED_Calo_Temp_NDF: " << LED_Calo_Temp_chi2/LED_Calo_Temp_NDF << " LED_Calo_Temp_r: " << LED_Calo_Temp_r << endl;
  cout << " *Lab2*  LED_Ref_Temp_chi2: " << LED_Ref_Temp_chi2 << " LED_Ref_Temp_NDF: " << LED_Ref_Temp_NDF << " LED_Ref_Temp_chi2/LED_Ref_Temp_NDF: " << LED_Ref_Temp_chi2/LED_Ref_Temp_NDF << " LED_Ref_Temp_r: " << LED_Ref_Temp_r << endl;
/***** Now Assuming a Bivariate Normal Distribution...Temp Correction Time *****/

// E(Y|X) = E(Y) + rsig_Y*(X-E(X))/sig_X
// Compute E(Y) = E(Y|X) - rsig_Y*(X-E(X))/sig_X
//    
std::vector<double> Corr_Am, Corr_Bi, Corr_Calo, Corr_Ref, Conditional_Bi_Avg, Conditional_Am_Avg, Conditional_Calo_Avg, Conditional_Ref_Avg, Conditional_Bi_Error, Conditional_Am_Error, Conditional_Calo_Error, Conditional_Ref_Error;
Double_t CorrAm_vs_Time_Array[newsize], CorrBi_vs_Time_Array[newsize], CorrRef_vs_Time_Array[newsize], CorrCalo_vs_Time_Array[newsize], Am_sigma, Am_sum, T_sigma, T_sum, T2_sigma, T2_sum, Bi_sigma, Bi_sum, Calo_sigma, Calo_sum, Ref_sigma, Ref_sum, sum, sum_sigma, Conditional_Bi_Err, Conditional_Am_Err, Conditional_LED_Calo_Err, Conditional_LED_Ref_Err, Bi_Sum_Err, Am_Sum_Err, LED_Calo_Sum_Err, LED_Ref_Sum_Err, Conditional_Am_sum, Conditional_Bi_sum, Conditional_Calo_sum, Conditional_Ref_sum, Conditional_Bi_sigma, Conditional_Am_sigma, Conditional_Calo_sigma, Conditional_Ref_sigma, Conditional_Bi_AvgSum, Conditional_Am_AvgSum, Conditional_Calo_AvgSum, Conditional_Ref_AvgSum, Conditional_Bi_AvgMean, Conditional_Am_AvgMean, Conditional_Calo_AvgMean, Conditional_Ref_AvgMean, Bi_Temp_Cond_AvgMean, Am_Temp_Cond_AvgMean, Calo_Temp_Cond_AvgMean, Ref_Temp_Cond_AvgMean, sum2;
int counter3;

/**** First E(Y|X) must be approximated for every X=x *****/

        for(int i=0; i<newsize; i++){
                counter3 = 0; 
                sum = 0; 
		sum2 = 0;
                sum_sigma = 0; 
                for(int j=0; j<newsize; j++){
                        if(T_DB_vs_Time_Array[j] == T_DB_vs_Time_Array[i]){
                                sum = sum + Bi_vs_Time_Array[j]/pow(MeV_Sigma[j],2);
				sum2 = sum2 + 1/pow(MeV_Sigma[j],2);
                                sum_sigma = sum_sigma + 1/pow(Bi_vs_Time_Err_Array[j],2);
				//cout << "Bi_vs_Time_ErrSame: " << Bi_vs_Time_Err_Array[j] << endl;
                                counter3++;
                        }    
                        else{
                                continue;
                        }    
                }    

                if(counter3 != 0){
                        Conditional_Bi_Avg.push_back(sum/sum2);
                        Conditional_Bi_Error.push_back(sqrt(1/sum_sigma));
                                }    
                else{
                        Conditional_Bi_Avg.push_back(Bi_vs_Time_Array[i]);
                        Conditional_Bi_Error.push_back(Bi_vs_Time_Err_Array[i]);
			//cout << "Bi_vs_Time_ErrIndiv: " << Bi_vs_Time_Err_Array[i] << endl;     
                }    
        }    


        for(int i=0; i<newsize; i++){
                counter3 = 0; 
                sum = 0; 
		sum2 = 0;
                sum_sigma = 0; 
                for(int j=0; j<newsize; j++){
                        if(T_DB_vs_Time_Array[j] == T_DB_vs_Time_Array[i]){
                                //sum = sum + Am_vs_Time_Array[j]/pow(Am_vs_Time_Err_Array[j],2);
				sum = sum + Am_vs_Time_Array[j]/pow(Americium1_Sigma[j],2);
				//sum2 = sum2 + 1/pow(Am_vs_Time_Err_Array[j],2);
				sum2 = sum2 + 1/pow(Americium1_Sigma[j],2);
                                sum_sigma = sum_sigma + 1/pow(Am_vs_Time_Err_Array[j],2);
                                counter3++;
                        }    
                        else{
                                continue;
                        }    
                }    

                if(counter3 != 0){
                Conditional_Am_Avg.push_back(sum/sum2);
                Conditional_Am_Error.push_back(sqrt(1/sum_sigma));
                //cout << "Cond Am_Avg: " << Conditional_Am_Avg[i] << "  and error:  " << Conditional_Am_Error[i] << endl;
                }    
     
                else{     
                        Conditional_Am_Avg.push_back(Am_vs_Time_Array[i]);
                        Conditional_Am_Error.push_back(Am_vs_Time_Err_Array[i]);
                        //cout << "Cond Am_Avg: " << Conditional_Am_Avg[i] << endl;   
                }    
        }    

        for(int i=0; i<newsize; i++){
                counter3 = 0; 
                sum = 0; 
		sum2 = 0;
                sum_sigma = 0; 
                for(int j=0; j<newsize; j++){
                        if(T_Lab2_vs_Time_Array[j] == T_Lab2_vs_Time_Array[i]){
                                sum = sum + LED_Calo_vs_Time_Array[j]/pow(LED_Sigma[j],2);
				sum2 = sum2 + 1/pow(LED_Sigma[j],2);
                                sum_sigma = sum_sigma + 1/pow(LED_Calo_vs_Time_Err_Array[j],2);
                                counter3++;
                        }    
                        else{
                                continue;
                        }    
                }    

                if(counter3 != 0){
                Conditional_Calo_Avg.push_back(sum/sum2);
                Conditional_Calo_Error.push_back(sqrt(1/sum_sigma));
                //cout << "Cond_Calo_Avg: " << Conditional_Calo_Avg[i] << endl;
                }    
     
                else{     
                        Conditional_Calo_Avg.push_back(LED_Calo_vs_Time_Array[i]);
                        Conditional_Calo_Error.push_back(LED_Calo_vs_Time_Err_Array[i]);   
                }    
        }    


        for(int i=0; i<newsize; i++){
                counter3 = 0; 
                sum = 0; 
		sum2 = 0;
                sum_sigma = 0; 
                for(int j=0; j<newsize; j++){
                        if(T_Lab2_vs_Time_Array[j] == T_Lab2_vs_Time_Array[i]){
                                sum = sum + LED_Ref_vs_Time_Array[j]/pow(LED2_Sigma[j],2);
				sum2 = sum2 + 1/pow(LED2_Sigma[j],2);
                                sum_sigma = sum_sigma + 1/pow(LED_Ref_vs_Time_Err_Array[j],2);
                                counter3++;
                        }    
                        else{
                                continue;
                        }    
                }    

                if(counter3 != 0){
                Conditional_Ref_Avg.push_back(sum/sum2);
                Conditional_Ref_Error.push_back(sqrt(1/sum_sigma));
                //cout << "Cond_Ref_Avg: " << Conditional_LED_Ref_Avg[i] << endl;
                }    
     
                else{     
                        Conditional_Ref_Avg.push_back(LED_Ref_vs_Time_Array[i]);
                        Conditional_Ref_Error.push_back(LED_Ref_vs_Time_Err_Array[i]);   
                }    
        }    
     
	Conditional_Bi_AvgSum = 0;
	Conditional_Am_AvgSum = 0;
	Conditional_Calo_AvgSum = 0;
	Conditional_Ref_AvgSum = 0;
	Double_t Bi_Temp_Cond_AvgSum = 0;
	Double_t Am_Temp_Cond_AvgSum = 0;
	Double_t Calo_Temp_Cond_AvgSum = 0;
	Double_t Ref_Temp_Cond_AvgSum = 0;
	for(int i=0; i<newsize; i++){
		Conditional_Bi_AvgSum = Conditional_Bi_AvgSum + Conditional_Bi_Avg[i];
		Conditional_Am_AvgSum = Conditional_Am_AvgSum + Conditional_Am_Avg[i];	
		Conditional_Calo_AvgSum = Conditional_Calo_AvgSum + Conditional_Calo_Avg[i];
     		Conditional_Ref_AvgSum = Conditional_Ref_AvgSum + Conditional_Ref_Avg[i];
		Bi_Temp_Cond_AvgSum = Bi_Temp_Cond_AvgSum + T_DB_vs_Time_Array[i]*Conditional_Bi_Avg[i];
		Am_Temp_Cond_AvgSum = Am_Temp_Cond_AvgSum + T_DB_vs_Time_Array[i]*Conditional_Am_Avg[i];
		Calo_Temp_Cond_AvgSum = Calo_Temp_Cond_AvgSum + T_Lab2_vs_Time_Array[i]*Conditional_Calo_Avg[i];
		Ref_Temp_Cond_AvgSum = Ref_Temp_Cond_AvgSum + T_Lab2_vs_Time_Array[i]*Conditional_Ref_Avg[i];
	}
	Conditional_Bi_AvgMean = Conditional_Bi_AvgSum/newsize; 
        Conditional_Am_AvgMean = Conditional_Am_AvgSum/newsize; 
        Conditional_Calo_AvgMean = Conditional_Calo_AvgSum/newsize; 
        Conditional_Ref_AvgMean = Conditional_Ref_AvgSum/newsize;

		 	
	for(int i=0; i<newsize; i++){ 
	//cout <<"Bi_Sigma: " << MeV_Sigma[i]/MeV_Sigma[0] << " Am_Sigma: " << Americium1_Sigma[i]/Americium1_Sigma[0] << " Calo: " << LED_Sigma[i]/LED_Sigma[0] << " Ref: " << LED2_Sigma[i]/LED2_Sigma[0] << endl;
	} 
	/***** Need E[X*E[Y|X]] ******/   
	Bi_Temp_Cond_AvgMean = Bi_Temp_Cond_AvgSum/newsize;
	Am_Temp_Cond_AvgMean = Am_Temp_Cond_AvgSum/newsize;
	Calo_Temp_Cond_AvgMean = Calo_Temp_Cond_AvgSum/newsize;
	Ref_Temp_Cond_AvgMean = Ref_Temp_Cond_AvgSum/newsize;
	
	//cout << "ADC*Temp_AvgMean...Bi: " << Bi_Temp_Cond_AvgMean << " Am: " << Am_Temp_Cond_AvgMean << " Calo: " << Calo_Temp_Cond_AvgMean << " Ref: " << Ref_Temp_Cond_AvgMean << endl;
       /* Bi_Sum_Err = 0; 
        for(int i=0; i<newsize; i++){
                Bi_Sum_Err = Bi_Sum_Err + pow(Conditional_Bi_PreError[i],2);
                }    
        Conditional_Bi_Err = sqrt(Bi_Sum_Err/newsize);


        Am_Sum_Err = 0; 
        for(int i=0; i<newsize; i++){
                Am_Sum_Err = Am_Sum_Err + pow(Conditional_Am_PreError[i],2);
                }    
        Conditional_Am_Err = sqrt(Am_Sum_Err/newsize);


        LED_Calo_Sum_Err = 0; 
        for(int i=0; i<newsize; i++){
                LED_Calo_Sum_Err = LED_Calo_Sum_Err + pow(Conditional_LED_Calo_PreError[i],2);
                }    
        Conditional_LED_Calo_Err = sqrt(LED_Calo_Sum_Err/newsize);
     
     
        LED_Ref_Sum_Err = 0; 
        for(int i=0; i<newsize; i++){
                LED_Ref_Sum_Err = LED_Ref_Sum_Err + pow(Conditional_LED_Ref_PreError[i],2);
                }    
        Conditional_LED_Ref_Err = sqrt(LED_Ref_Sum_Err/newsize);
        */

/***** Compute std deviations *****/

  Double_t Am_r_err, Bi_r_err, Calo_r_err, Ref_r_err, T_sum2, T2_sum2, T_sigma2, T2_sigma2, Am_sum2, Bi_sum2, Calo_sum2, Ref_sum2, Am_sigma2, Bi_sigma2, Ref_sigma2, Calo_sigma2;
	Conditional_Am_sum = 0;
	Am_sum = 0;
	Am_sum2 = 0;
	for(int i=0; i<newsize; i++){
	Am_sum = Am_sum+pow(Am_vs_Time_Err_Array[i],2);
	Am_sum2 = Am_sum2+pow(Am_vs_Time_Array[i]-Am_Temp_yMean,2);
	Conditional_Am_sum = Conditional_Am_sum+pow((Conditional_Am_Avg[i]-Conditional_Am_AvgMean),2);
	}
	Am_sigma=sqrt(Am_sum/newsize);
	Am_sigma2=sqrt(Am_sum2/newsize);
	Conditional_Am_sigma=sqrt(Conditional_Am_sum/newsize);
	Am_r_err = sqrt((1-pow(Am_Temp_r,2))/(newsize-1));


	Conditional_Bi_sum = 0;
	Bi_sum = 0;
	Bi_sum2 = 0;
        for(int i=0; i<newsize; i++){
        Bi_sum = Bi_sum+pow(Bi_vs_Time_Err_Array[i],2);
	Bi_sum2 = Bi_sum2+pow(Bi_vs_Time_Array[i]-Bi_Temp_yMean,2);
	Conditional_Bi_sum = Conditional_Bi_sum+pow((Conditional_Bi_Avg[i]-Conditional_Bi_AvgMean),2);
        }   
        Bi_sigma=sqrt(Bi_sum/newsize);
	Bi_sigma2=sqrt(Bi_sum2/newsize);
	Conditional_Bi_sigma=sqrt(Conditional_Bi_sum/newsize);
	Bi_r_err = sqrt((1-pow(Bi_Temp_r,2))/(newsize-1));
	
	
	Conditional_Calo_sum = 0;
	Calo_sum = 0;
	Calo_sum2 = 0;
        for(int i=0; i<newsize; i++){
        Calo_sum = Calo_sum+pow(LED_Calo_vs_Time_Err_Array[i],2);
	Calo_sum2 = Calo_sum2+pow(LED_Calo_vs_Time_Array[i]-LED_Calo_Temp_yMean,2);
	Conditional_Calo_sum = Conditional_Calo_sum+pow((Conditional_Calo_Avg[i]-Conditional_Calo_AvgMean),2);
        }   
        Calo_sigma=sqrt(Calo_sum/newsize);
	Calo_sigma2=sqrt(Calo_sum2/newsize);
	Conditional_Calo_sigma=sqrt(Conditional_Calo_sum/newsize);
	Calo_r_err = sqrt((1-pow(LED_Calo_Temp_r,2))/(newsize-1));


	Conditional_Ref_sum = 0;
	Ref_sum = 0;
	Ref_sum2 = 0;
        for(int i=0; i<newsize; i++){
        Ref_sum = Ref_sum+pow(LED_Ref_vs_Time_Err_Array[i],2);
	Ref_sum2 = Ref_sum2+pow(LED_Ref_vs_Time_Array[i]-LED_Ref_Temp_yMean,2);
	Conditional_Ref_sum = Conditional_Ref_sum+pow((Conditional_Ref_Avg[i]-Conditional_Ref_AvgMean),2);
        }   
        Ref_sigma=sqrt(Ref_sum/newsize);
	Ref_sigma2=sqrt(Ref_sum2/newsize);
	Conditional_Ref_sigma=sqrt(Conditional_Ref_sum/(newsize));
	Ref_r_err = sqrt((1-pow(LED_Ref_Temp_r,2))/(newsize-1));

	
	T_sum2 = 0;
	T_sum = 0;
        for(int i=0; i<newsize; i++){
        T_sum = T_sum+pow(T_Err_vs_Time_Array[i],2);
	T_sum2 = T_sum2+pow(T_DB_vs_Time_Array[i]-Am_Temp_xMean,2);
        }   
        T_sigma=sqrt(T_sum/newsize);
	T_sigma2=sqrt(T_sum2/newsize);

	T2_sum2 = 0;
	T2_sum = 0;
        for(int i=0; i<newsize; i++){
        T2_sum = T2_sum+pow(T_Err_vs_Time_Array[i],2);
	T2_sum2 = T2_sum2+pow(T_Lab2_vs_Time_Array[i]-Am_Temp_xMean,2);
        }    
        T2_sigma=sqrt(T2_sum/newsize);
	T2_sigma2=sqrt(T2_sum2/newsize);
	
	Double_t Am_Cov_Sum, Bi_Cov_Sum, Calo_Cov_Sum, Ref_Cov_Sum, Am_Cov, Bi_Cov, Calo_Cov, Ref_Cov;
	Am_Cov_Sum = 0;
	Bi_Cov_Sum = 0;
	for(int i=0; i<newsize; i++){
		Am_Cov_Sum = Am_Cov_Sum + ((Conditional_Am_Avg[i]-Conditional_Am_AvgMean)*(T_DB_vs_Time_Array[i]-Am_Temp_xMean));
		Bi_Cov_Sum = Bi_Cov_Sum + ((Conditional_Bi_Avg[i]-Conditional_Bi_AvgMean)*(T_DB_vs_Time_Array[i]-Bi_Temp_xMean));
		Calo_Cov_Sum = Calo_Cov_Sum + ((Conditional_Calo_Avg[i]-Conditional_Calo_AvgMean)*(T_Lab2_vs_Time_Array[i]-LED_Calo_Temp_xMean));
		Ref_Cov_Sum = Ref_Cov_Sum + ((Conditional_Ref_Avg[i]-Conditional_Ref_AvgMean)*(T_Lab2_vs_Time_Array[i]-LED_Ref_Temp_xMean));
	}
	Am_Cov = Am_Cov_Sum/newsize;
	Bi_Cov = Bi_Cov_Sum/newsize;
	//cout << "Am_Cov: " << Am_Cov << endl;
	//cout << "Bi_Cov: " << Bi_Cov << endl;
	//cout << "Calo_Cov: " << Calo_Cov << endl;
	//cout << "Ref_Cov: " << Ref_Cov << endl;

/****** Fill Corrected Arrays ******/

Double_t CorrAm_vs_Time_Err_Array[newsize], CorrBi_vs_Time_Err_Array[newsize], CorrCalo_vs_Time_Err_Array[newsize], CorrRef_vs_Time_Err_Array[newsize], Am_Corr_Line_Ratio_Array[newsize], Bi_Corr_Line_Ratio_Array[newsize], Calo_Corr_Line_Ratio_Array[newsize], Ref_Corr_Line_Ratio_Array[newsize];
std::vector<double> Corr_Am_Err, Corr_Bi_Err, Corr_Calo_Err, Corr_Ref_Err, Bi_Temp_Line, Am_Temp_Line, Calo_Temp_Line, Ref_Temp_Line, Am_Corr_Line_Ratio, Bi_Corr_Line_Ratio, Calo_Corr_Line_Ratio, Ref_Corr_Line_Ratio;



// Am Corrected 
	for(int i=0; i<newsize; i++){
	//	Corr_Am.push_back((Am_Temp_slope*T_DB_vs_Time_Array[i]+Am_Temp_const)-(Am_Temp_r*Am_sigma*(T_DB_vs_Time_Array[i]-Am_Temp_xMean))/T_sigma);
		//Corr_Am.push_back(Am_vs_Time_Array[i]-(Am_Temp_r*Americium1_Sigma[i]*(T_DB_vs_Time_Array[i]-T_DB_vs_Time_Array[0]))/T_DB_sigma[i]); // Y|X
		Corr_Am.push_back(Conditional_Am_Avg[i]-(Am_Temp_r*Am_sigma2*(T_DB_vs_Time_Array[i]-Am_Temp_xMean))/T_sigma2); // E[Y|X]
		//Corr_Am_Err.push_back(sqrt(pow(Conditional_Am_sigma,2)+(pow(Am_sigma/T_sigma,2))*((pow(T_DB_vs_Time_Array[i]*Am_Temp_r,2))*(pow(Am_r_err/Am_Temp_r,2)+pow(T_sigma/T_DB_vs_Time_Array[i],2)))+pow(Am_Temp_xMean,2)*pow(Am_r_err,2)));
		//Corr_Am_Err.push_back(sqrt(pow(Conditional_Am_Error[i],2)+(pow(Am_Temp_r*Am_sigma/T_sigma,2))*(pow(T_Err_vs_Time_Array[i],2)+pow(T_sigma,2))));/*2*((Am_Temp_r*Am_sigma/T_sigma)*(Am_Temp_Cond_AvgMean-(Conditional_Am_AvgMean*Am_Temp_xMean)))));*/ // if r has no uncertainty
		Corr_Am_Err.push_back(sqrt(pow(Conditional_Am_Error[i],2)+(pow(Am_Temp_r*Am_sigma2/T_sigma2,2))*(pow(T_Err_vs_Time_Array[i],2)+pow(T_sigma,2)))); // Std Dev Formula: X-Xmean etc
		
		//Corr_Am_Err.push_back(pow(Conditional_Bi_sigma,2)+(pow(Bi_sigma/T_sigma,2))*((pow(T_DB_vs_Time_Array[i]*Bi_Temp_r,2))*(pow(Bi_r_err/Bi_Temp_r,2)+pow(T_sigma/T_DB_vs_Time_Array[i],2))+pow(Bi_Temp_xMean,2)*pow(Bi_r_err,2))); //Test where all have Corr_Bi_Error
		//cout << "Error Test CorrAm: " << Corr_Am_Err[i] << endl;
	}

	for(int i=0; i<newsize; i++){
		CorrAm_vs_Time_Array[i] = Corr_Am[i];
		CorrAm_vs_Time_Err_Array[i] = Corr_Am_Err[i];
	}

// Bi Corrected
	for(int i=0; i<newsize; i++){  // added -1
        //	Corr_Bi.push_back((Bi_Temp_slope*T_DB_vs_Time_Array[i]+Bi_Temp_const)-(Bi_Temp_r*Bi_sigma*(T_DB_vs_Time_Array[i]-Bi_Temp_xMean))/T_sigma);
		//Corr_Bi.push_back(Bi_vs_Time_Array[i]-(Bi_Temp_r*MeV_Sigma[i]*(T_DB_vs_Time_Array[i]-T_DB_vs_Time_Array[0]))/T_DB_sigma[i]); // Y|X 
		Corr_Bi.push_back(Conditional_Bi_Avg[i]-(Bi_Temp_r*Bi_sigma2*(T_DB_vs_Time_Array[i]-Bi_Temp_xMean))/T_sigma2); // E[Y|X]
		//Corr_Bi_Err.push_back(pow(Conditional_Bi_sigma,2)+(pow(Bi_sigma/T_sigma,2))*(pow(Bi_r_err/Bi_Temp_r,2)+pow(T_sigma/T_DB_vs_Time_Array[i],2)+pow(Bi_Temp_xMean,2)*pow(Bi_r_err,2)));
		//Corr_Bi_Err.push_back(sqrt(pow(Conditional_Bi_sigma,2)+(pow(Bi_sigma/T_sigma,2))*((pow(T_DB_vs_Time_Array[i]*Bi_Temp_r,2))*(pow(Bi_r_err/Bi_Temp_r,2)+pow(T_sigma/T_DB_vs_Time_Array[i],2))+pow(Bi_Temp_xMean,2)*pow(Bi_r_err,2))));
		//Corr_Bi_Err.push_back(sqrt(pow(Conditional_Bi_Error[i],2)+(pow(Bi_Temp_r*Bi_sigma/T_sigma,2))*(pow(T_Err_vs_Time_Array[i],2)+pow(T_sigma,2)))); // Std Error of Mean for all
		Corr_Bi_Err.push_back(sqrt(pow(Conditional_Bi_Error[i],2)+(pow(Bi_Temp_r*Bi_sigma2/T_sigma2,2))*(pow(T_Err_vs_Time_Array[i],2)+pow(T_sigma,2)))); // Std Err X-Xmean for formula 
		//cout << "Error Test CorrBi: " << Corr_Bi_Err[i] << endl;
        }

        for(int i=0; i<newsize; i++){
                CorrBi_vs_Time_Array[i] = Corr_Bi[i];
		CorrBi_vs_Time_Err_Array[i] = Corr_Bi_Err[i];
	}

// Ref LED Corrected
	for(int i=0; i<newsize; i++){
        //	Corr_Ref.push_back((LED_Ref_Temp_slope*T_Lab2_vs_Time_Array[i]+LED_Ref_Temp_const)-(LED_Ref_Temp_r*Ref_sigma*(T_Lab2_vs_Time_Array[i]-LED_Ref_Temp_xMean))/T2_sigma);
		//Corr_Ref.push_back(LED_Ref_vs_Time_Array[i]-(LED_Ref_Temp_r*LED2_Sigma[i]*(T_Lab2_vs_Time_Array[i]-T_Lab2_vs_Time_Array[0]))/T_Lab2_sigma[i]); // Y|X
		Corr_Ref.push_back(Conditional_Ref_Avg[i]-(LED_Ref_Temp_r*Ref_sigma2*(T_Lab2_vs_Time_Array[i]-LED_Ref_Temp_xMean))/T2_sigma2);  // E[Y|X]
		//Corr_Ref_Err.push_back(pow(Conditional_Ref_sigma,2)+(pow(Ref_sigma/T2_sigma,2))*(pow(Ref_r_err/LED_Ref_Temp_r,2)+pow(T_sigma/T_Lab2_vs_Time_Array[i],2)+pow(LED_Ref_Temp_xMean,2)*pow(Ref_r_err,2)));
		//Corr_Ref_Err.push_back(sqrt(pow(Conditional_Ref_sigma,2)+(pow(Ref_sigma/T2_sigma,2))*((pow(T_Lab2_vs_Time_Array[i]*LED_Ref_Temp_r,2))*((pow(Ref_r_err/LED_Ref_Temp_r,2)+pow(T_sigma/T_Lab2_vs_Time_Array[i],2)))+pow(LED_Ref_Temp_xMean,2)*pow(Ref_r_err,2))));
		//Corr_Ref_Err.push_back(sqrt(pow(Conditional_Ref_Error[i],2)+(pow(LED_Ref_Temp_r*Ref_sigma/T2_sigma,2))*(pow(T_Err_vs_Time_Array[i],2)+pow(T2_sigma,2)))); // Std Error of Mean
		Corr_Ref_Err.push_back(sqrt(pow(Conditional_Ref_Error[i],2)+(pow(LED_Ref_Temp_r*Ref_sigma2/T2_sigma2,2))*(pow(T_Err_vs_Time_Array[i],2)+pow(T2_sigma,2))));
		//Corr_Ref_Err.push_back(pow(Conditional_Bi_sigma,2)+(pow(Bi_sigma/T_sigma,2))*((pow(T_DB_vs_Time_Array[i]*Bi_Temp_r,2))*(pow(Bi_r_err/Bi_Temp_r,2)+pow(T_sigma/T_DB_vs_Time_Array[i],2))+pow(Bi_Temp_xMean,2)*pow(Bi_r_err,2))); //Bi_Error test
		//cout << "Error Test CorrRef: " << Corr_Ref_Err[i] << endl;
        }

        for(int i=0; i<newsize; i++){
                CorrRef_vs_Time_Array[i] = Corr_Ref[i];
		CorrRef_vs_Time_Err_Array[i] = Corr_Ref_Err[i];
        }

// Calo LED Corrected
	for(int i=0; i<newsize; i++){
        //	Corr_Calo.push_back((LED_Calo_Temp_slope*T_Lab2_vs_Time_Array[i]+LED_Calo_Temp_const)-(LED_Calo_Temp_r*Calo_sigma*(T_Lab2_vs_Time_Array[i]-LED_Calo_Temp_xMean))/T2_sigma);
		//Corr_Calo.push_back(LED_Calo_vs_Time_Array[i]-(LED_Calo_Temp_r*LED_Sigma[i]*(T_Lab2_vs_Time_Array[i]-T_Lab2_vs_Time_Array[0]))/T_Lab2_sigma[i]); // Y|X
		Corr_Calo.push_back(Conditional_Calo_Avg[i]-(LED_Calo_Temp_r*Calo_sigma2*(T_Lab2_vs_Time_Array[i]-LED_Calo_Temp_xMean))/T2_sigma2);  // E[Y|X]
		//Corr_Calo_Err.push_back(pow(Conditional_Calo_sigma,2)+(pow(Calo_sigma/T2_sigma,2))*(pow(Calo_r_err/LED_Calo_Temp_r,2)+pow(T_sigma/T_Lab2_vs_Time_Array[i],2)+pow(LED_Calo_Temp_xMean,2)*pow(Calo_r_err,2)));
		//Corr_Calo_Err.push_back(sqrt(pow(Conditional_Calo_sigma,2)+(pow(Calo_sigma/T2_sigma,2))*((pow(T_Lab2_vs_Time_Array[i]*LED_Calo_Temp_r,2))*(pow(Calo_r_err/LED_Calo_Temp_r,2)+pow(T_sigma/T_Lab2_vs_Time_Array[i],2))+pow(LED_Calo_Temp_xMean,2)*pow(Calo_r_err,2)))); //including amplitude term
		//Corr_Calo_Err.push_back(sqrt(pow(Conditional_Calo_Error[i],2)+(pow(LED_Calo_Temp_r*Calo_sigma/T2_sigma,2))*(pow(T_Err_vs_Time_Array[i],2)+pow(T2_sigma,2)))); // Std Error of Mean
		Corr_Calo_Err.push_back(sqrt(pow(Conditional_Calo_Error[i],2)+(pow(LED_Calo_Temp_r*Calo_sigma2/T2_sigma2,2))*(pow(T_Err_vs_Time_Array[i],2)+pow(T2_sigma,2))));
		//Corr_Calo_Err.push_back(pow(Conditional_Bi_sigma,2)+(pow(Bi_sigma/T_sigma,2))*((pow(T_DB_vs_Time_Array[i]*Bi_Temp_r,2))*(pow(Bi_r_err/Bi_Temp_r,2)+pow(T_sigma/T_DB_vs_Time_Array[i],2))+pow(Bi_Temp_xMean,2)*pow(Bi_r_err,2))); //Bi_Error test	
		//cout << "Error Test CorrCalo: " << Corr_Calo_Err[i] << "  Cond_Calo_sigma: " << Conditional_Calo_sigma << endl;
        }

	for(int i=0; i<newsize; i++){
                Am_Temp_Line.push_back(Am_Temp_slope*T_DB_vs_Time_Array[i]+Am_Temp_const);
                Bi_Temp_Line.push_back(Bi_Temp_slope*T_DB_vs_Time_Array[i]+Bi_Temp_const);
                Ref_Temp_Line.push_back(LED_Ref_Temp_slope*T_Lab2_vs_Time_Array[i]+LED_Ref_Temp_const);
                Calo_Temp_Line.push_back(LED_Calo_Temp_slope*T_Lab2_vs_Time_Array[i]+LED_Calo_Temp_const);
                }    

        for(int i=0; i<newsize; i++){
                Am_Corr_Line_Ratio.push_back(Conditional_Am_Avg[i]/Am_Temp_Line[i]);
                Bi_Corr_Line_Ratio.push_back(Conditional_Bi_Avg[i]/Bi_Temp_Line[i]);
                Calo_Corr_Line_Ratio.push_back(Conditional_Calo_Avg[i]/Calo_Temp_Line[i]);
                Ref_Corr_Line_Ratio.push_back(Conditional_Ref_Avg[i]/Ref_Temp_Line[i]);
                }    

        for(int i=0; i<newsize; i++){
                Am_Corr_Line_Ratio_Array[i] = Am_Corr_Line_Ratio[i];
                Bi_Corr_Line_Ratio_Array[i] = Bi_Corr_Line_Ratio[i];
                Calo_Corr_Line_Ratio_Array[i] = Calo_Corr_Line_Ratio[i];
                Ref_Corr_Line_Ratio_Array[i] = Ref_Corr_Line_Ratio[i];
                }


        for(int i=0; i<newsize; i++){
                CorrCalo_vs_Time_Array[i] = Corr_Calo[i];
		CorrCalo_vs_Time_Err_Array[i] = Corr_Calo_Err[i];
        }
	//cout << "Term by term : " << pow(Conditional_Calo_Error[0], 2) << " + " << (pow(LED_Calo_Temp_r*Calo_sigma/T2_sigma,2))*pow(T_Err_vs_Time_Array[0],2) << " - " << 2*((LED_Calo_Temp_r*Ref_sigma/T_sigma)*(Calo_Temp_Cond_AvgMean-(Conditional_Calo_AvgMean*LED_Calo_Temp_xMean))) << endl;
	//cout << "Covar term Split : " << 2*(LED_Calo_Temp_r*Calo_sigma/T2_sigma) << " and " << (Calo_Temp_Cond_AvgMean-(Conditional_Calo_AvgMean*LED_Calo_Temp_xMean)) << endl;
	//cout << "Covar term2 Split : " << Calo_Temp_Cond_AvgMean << " and " << Conditional_Calo_AvgMean*LED_Calo_Temp_xMean << endl;
	//cout << "NewCovar term : " << (-2*(LED_Calo_Temp_r*Ref_sigma/T_sigma)) << endl;
/******* Temperature Corrected Ratio ********/

// Start values (i.e. calibration run values)
  Bi_0 = Corr_Bi[0];
  Bi_0_err = Corr_Bi_Err[0];
  //Bi2_0 = HMeV_Mean[0];
  //Bi2_0_err = HMeV_Mean_Err[0];
  LED_Calo_0 = Corr_Calo[0];
  LED_Calo_0_err = Corr_Calo_Err[0];
  LED_Ref_0 = Corr_Ref[0];
  LED_Ref_0_err = Corr_Ref_Err[0];
  Am_0 = Corr_Am[0];
  //Am_0_err = Americium1_Mean_Err[0];
  //int corrnewsize;
  //corrnewsize = newsize-rejections;
  // Later values for LED and Am ONLY
  /*Double_t Bi_t, Bi_t_err, LED_Calo_t, LED_Calo_t_err, LED_Ref_t, LED_Ref_t_err, Am_t, Am_t_err, LED_Ref_Pred, LED_Ref_Pred_err, LED_Calo_Pred, LED_Calo_Pred_err, Bi_Pred, Bi_Pred_err, Ratio, Ratio_err;*/
  std::vector<double> Corr_Bi_Pred_vs_Time, Corr_Ratio_vs_Time, Bi_0_Bi_t, Am_t_Am_0, Ref_0_Ref_t, Calo_t_Calo_0, Corr_Bi_Pred_vs_Time_Err, Corr_Ratio_vs_Time_Err;
  Double_t Corr_Bi_Pred_vs_Time_Array[newsize], Corr_Ratio_vs_Time_Array[newsize], Corr_Ratio_vs_Time_Err_Array[newsize], Bi_0_Bi_t_Array[newsize], Am_t_Am_0_Array[newsize], Ref_0_Ref_t_Array[newsize], Calo_t_Calo_0_Array[newsize], Bi_to, LED_Calo_to, LED_Ref_to, Am_to, Bi_lasto, LED_Calo_lasto;
 
  Bi_last = Bi_0;
  LED_Calo_last = LED_Calo_0;
  counter1 = 0; 
  counter2 = 0;
  for(int i=0; i<newsize; i++){
    Bi_to	   = MeV_Mean[0];  //MeV_Mean[i]
    Bi_t           = Corr_Bi[i];
    Bi_t_err       = Corr_Bi_Err[i];  //MeV_Mean_Err[i];
    LED_Calo_to    = Corr_Calo[0];  //LED_Mean[i];
    LED_Calo_t     = Corr_Calo[i];
    LED_Calo_t_err = Corr_Calo_Err[i];  //LED_Mean_Err[i];
    LED_Ref_to     = Corr_Ref[0];  //LED2_Mean[i];
    LED_Ref_t      = Corr_Ref[i];
    LED_Ref_t_err  = Corr_Ref_Err[i];  //LED2_Mean_Err[i];
    Am_to	   = Corr_Am[0];  //Americium1_Mean[i];
    Am_t           = Corr_Am[i];
    Am_t_err       = Corr_Am_Err[i];  //Americium1_Mean_Err[i];

    // Skip runs with large Bi peak deviations unaccounted for in LED peak
    if( (TMath::Abs(Bi_to - Bi_lasto) > 10*Bi_t_err) &&
        !(TMath::Abs(LED_Calo_to - LED_Calo_lasto) > 10*LED_Calo_t_err)){
      if(Debug) {
        printf("Reese: Skipping Run %d (timestamp: %f): Bi mean value = %f while LED mean value = %f\n",i, Time_Data.at(i), Bi_t, LED_Calo_t);
        }
      continue;
    }
    
    else{
    /*if(Am_t/Am_0 < 0.5 or Am_t/Am_0 < 2){
	continue;
    }*/

  
    Bi_lasto = Bi_to;
    Bi_last = Bi_t; // Update comparison values
    LED_Calo_lasto = LED_Calo_to;
    //Am_t = Am_0;
    // Debugging outputs
//    if(Debug) printf("%f: Bi mean = %f LED mean  = %f", i, Time_Data.at(i), Bi_t, LED_Calo_t);    

    // Use Am to predict Ref LED (correct for drifts)
    LED_Ref_Pred = LED_Ref_0 * (Am_t/Am_0);
    LED_Ref_Pred_err = sqrt(pow(LED_Ref_0_err*Am_t/Am_0,2)
                            + pow(LED_Ref_0*Am_t_err/Am_0,2)
                            + pow(LED_Ref_0*Am_t*Am_0_err/(Am_0*Am_0),2));
 //   if(Debug) printf("LED pred: %f*%f = %f\n", Am_t, Am_0, (Am_t/Am_0));
    // Use Ref LED to predict Calo LED
    LED_Calo_Pred = LED_Calo_t * (LED_Ref_Pred/LED_Ref_t);
    LED_Calo_Pred_err = sqrt(pow(LED_Calo_t_err*LED_Ref_Pred/LED_Ref_t,2)
                             + pow(LED_Calo_t*LED_Ref_Pred_err/LED_Ref_t,2)
                             + pow(LED_Calo_t*LED_Ref_Pred*LED_Ref_t_err/(LED_Ref_t*LED_Ref_t),2));
  
    // Use Calo LED to predict Bi
    Bi_Pred = Bi_0 * (LED_Calo_Pred/LED_Calo_0); /*(Bi_0 * LED_Calo_t * LED_Ref_0 *Am_t)/(LED_Calo_0 * LED_Ref_t *Am_0);*/ //just trying a more explicit implementation of the prediction algorithm, but same effect
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
   /* Bi_vs_Time.push_back(Bi_t);
    Bi_vs_Time_Err.push_back(Bi_t_err);
    LED_Calo_vs_Time.push_back(LED_Calo_t);
    LED_Calo_vs_Time_Err.push_back(LED_Calo_t_err);
    LED_Ref_vs_Time.push_back(LED_Ref_t);
    LED_Ref_vs_Time_Err.push_back(LED_Ref_t_err);
    Am_vs_Time.push_back(Am_t);
    Am_vs_Time_Err.push_back(Am_t_err);*/
    
    //cout << "Am_t :  " << Am_t << "  at:  " << i << endl;  
    Corr_Bi_Pred_vs_Time.push_back(Bi_Pred);
    Corr_Bi_Pred_vs_Time_Err.push_back(Bi_Pred_err);
    Corr_Ratio_vs_Time.push_back(Ratio);
    Corr_Ratio_vs_Time_Err.push_back(Ratio_err);
    //cout << "Corr Ratio_vs_Time_Err @" << i << ": " << Corr_Ratio_vs_Time_Err[i] << endl; 
    //cout << "Corr_Ratio_vs_Time[i]  " << Corr_Ratio_vs_Time[i] << endl; 
    
   
/**** Troubleshoot with Ratio Plots *****/
     Bi_0_Bi_t.push_back(Bi_0/Bi_t); 
     Am_t_Am_0.push_back(Am_t/Am_0);
     Ref_0_Ref_t.push_back(LED_Ref_0/LED_Ref_t);
     Calo_t_Calo_0.push_back(LED_Calo_t/LED_Calo_0);

     }	               
//cout << "  Bi_0_Bi_t:  " << Bi_0/Bi_t << "  Am_t_Am_0:  " << Am_t/Am_0 << "  Ref_0_Ref_t:  " << LED_Ref_0/LED_Ref_t << "  Calo_t_Calo_0:  " << LED_Calo_t/LED_Calo_0 << "  Ratio  " << test << endl;
	
      
  }

  //const int mynewsize = corrnewsize - rejections;
	cout << "Number of Analysis Rejections:  " << rejections << endl;
	for(int i=0; i<newsize; i++){
		//cout << " Test " << endl;
		Corr_Ratio_vs_Time_Array[i] = Corr_Ratio_vs_Time[i];
		Corr_Ratio_vs_Time_Err_Array[i] = Corr_Ratio_vs_Time_Err[i];
		
		if(Corr_Ratio_vs_Time[i]<0.99){
			cout << "Bad Ratio @DB_Temp: " << T_DB_vs_Time_Array[i] << "  Bad Ratio @Lab_Temp: " << T_Lab2_vs_Time_Array[i] << endl;
		}
		//cout << " Corr_Ratio_vs_Time[i]: " << Corr_Ratio_vs_Time[i] << endl;
	//	cout << " Corr_Ratio_vs_Time_Array[i]:  " << Corr_Ratio_vs_Time_Array[i] << endl;
	}

	/*for(int i; i<newsize; i++){
		Bi_0_Bi_t.push_back(Bi_0/Bi_t) 
		Am_t_Am_0.push_back(Am_t/Am_0)
		Ref_0_Ref_t.push_back(LED_Ref_0/LED_Ref_t)
		Calo_t_Calo_0.push_back(LED_Calo_t/LED_Calo_0)
	}

	for(int i; i<newsize; i++){
		cout << "  Bi_0_Bi_t:  " << Bi_0_Bi_t[i] << "  Am_t_Am_0:  " << Am_t_Am_0[i] << "  Ref_0_Ref_t:  " << Ref_0_Ref_t[i] << "  Calo_t_Calo_0:  " << Calo_t_Calo_0 << endl;
	}*/




  
  /***** Corrected Graphs ******/
// Corr_Am Graphs
  TGraphErrors *CorrAm_vs_Time_Graph = new TGraphErrors(newsize,Time,CorrAm_vs_Time_Array,Time_Err,CorrAm_vs_Time_Err_Array);
  TGraphErrors *CorrAm_vs_Temp_Graph = new TGraphErrors(newsize,T_DB_vs_Time_Array,CorrAm_vs_Time_Array,T_Err_vs_Time_Array, CorrAm_vs_Time_Err_Array);

// Corr_Bi Graphs
  TGraphErrors *CorrBi_vs_Time_Graph = new TGraphErrors(newsize,Time,CorrBi_vs_Time_Array,Time_Err,CorrBi_vs_Time_Err_Array);
  TGraphErrors *CorrBi_vs_Temp_Graph = new TGraphErrors(newsize,T_DB_vs_Time_Array,CorrBi_vs_Time_Array,T_Err_vs_Time_Array, CorrBi_vs_Time_Err_Array);

// Corr_Ref_LED Graphs
  TGraphErrors *CorrRef_vs_Time_Graph = new TGraphErrors(newsize,Time,CorrRef_vs_Time_Array,Time_Err, CorrRef_vs_Time_Err_Array);
  TGraphErrors *CorrRef_vs_Temp_Graph = new TGraphErrors(newsize,T_Lab2_vs_Time_Array,CorrRef_vs_Time_Array,T_Err_vs_Time_Array, CorrRef_vs_Time_Err_Array);

// Corr_Calo_LED Graphs
  TGraphErrors *CorrCalo_vs_Time_Graph = new TGraphErrors(newsize,Time,CorrCalo_vs_Time_Array,Time_Err,CorrCalo_vs_Time_Err_Array);
  TGraphErrors *CorrCalo_vs_Temp_Graph = new TGraphErrors(newsize,T_Lab2_vs_Time_Array,CorrCalo_vs_Time_Array,T_Err_vs_Time_Array, CorrCalo_vs_Time_Err_Array);


// Corrected Signals vs Fit Line Ratio Plots
  TGraph *CorrBi_Line_Ratio_Graph = new TGraph(newsize,Time,Bi_Corr_Line_Ratio_Array);
  TGraph *CorrAm_Line_Ratio_Graph = new TGraph(newsize,Time,Am_Corr_Line_Ratio_Array);
  TGraph *CorrCalo_Line_Ratio_Graph = new TGraph(newsize,Time,Calo_Corr_Line_Ratio_Array);
  TGraph *CorrRef_Line_Ratio_Graph = new TGraph(newsize,Time,Ref_Corr_Line_Ratio_Array);


//!!!!! Corrected Ratio Graphs !!!!//
 TGraphErrors *Corr_Ratio_vs_Time_Graph = new TGraphErrors(newsize,Time,Corr_Ratio_vs_Time_Array,Time_Err,Corr_Ratio_vs_Time_Err_Array);
 TGraphErrors *Corr_Ratio_vs_Temp_Graph = new TGraphErrors(newsize,T_Lab2_vs_Time_Array,Corr_Ratio_vs_Time_Array,T_Err_vs_Time_Array, Corr_Ratio_vs_Time_Err_Array);

 //TGraph *Corr_Ratio_vs_Time_Graph = new TGraph(newsize,Time, Corr_Ratio_vs_Time_Array);
 //TGraph *Corr_Ratio_vs_Temp_Graph = new TGraph(newsize,T_DB_vs_Time_Array,Corr_Ratio_vs_Time_Array);
  
  /*for(int i=0; i<newsize; i++){
  	cout << "Corrected Ratio @i:  " << Corr_Ratio_vs_Time_Array[i] << endl;
  	cout << "Original Ratio @i:  "  << Ratio_vs_Time_Array[i] << endl;
  }*/
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
  TCanvas *c14 = new TCanvas("c14","c14", 10, 10, 1950,600);
  c14->Divide(1,2);
  c14->cd(1);
  c14->SetGrid();
  Corr_Ratio_vs_Time_Graph->SetTitle("Corrected Ratio vs Time");
  Corr_Ratio_vs_Time_Graph->GetXaxis()->SetRangeUser(0, Time[newsize-1]);
  Corr_Ratio_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  Corr_Ratio_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  Corr_Ratio_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  Corr_Ratio_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  Corr_Ratio_vs_Time_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
  Corr_Ratio_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  Corr_Ratio_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  Corr_Ratio_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  Corr_Ratio_vs_Time_Graph->SetMarkerStyle(7);
  Corr_Ratio_vs_Time_Graph->SetMarkerColor(kRed);
  Corr_Ratio_vs_Time_Graph->Draw("AP");

  c14->cd(2);
  c14->SetGrid();
  Corr_Ratio_vs_Temp_Graph->SetTitle("Corrected Ratio vs T_Lab");
  Corr_Ratio_vs_Temp_Graph->GetXaxis()->SetRangeUser(T_Lab2_Min - 0.2,T_Lab2_Max + 0.2);
  Corr_Ratio_vs_Temp_Graph->GetXaxis()->SetTitle("Temperature [Deg. F]");
  Corr_Ratio_vs_Temp_Graph->GetXaxis()->SetTitleSize(0.08);
  Corr_Ratio_vs_Temp_Graph->GetXaxis()->SetTitleOffset(0.9);
  Corr_Ratio_vs_Temp_Graph->GetXaxis()->SetLabelSize(0.08);
  Corr_Ratio_vs_Temp_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
  Corr_Ratio_vs_Temp_Graph->GetYaxis()->SetTitleSize(0.08);
  Corr_Ratio_vs_Temp_Graph->GetYaxis()->SetTitleOffset(0.4);
  Corr_Ratio_vs_Temp_Graph->GetYaxis()->SetLabelSize(0.08);
  Corr_Ratio_vs_Temp_Graph->SetMarkerStyle(7);
  Corr_Ratio_vs_Temp_Graph->SetMarkerColor(kRed);
  Corr_Ratio_vs_Temp_Graph->Draw("AP");
  c14->Print("Plots/Corr_Ratio_Time&Temp.png");  
 
 
  TCanvas *c10 = new TCanvas("c10","c10", 10,10,1950,600);
  c10->Divide(1,2);
  c10->cd(1);
  c10->SetGrid();
  CorrBi_vs_Time_Graph->SetTitle("Corrected Bi vs Time");
  CorrBi_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  CorrBi_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  CorrBi_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  CorrBi_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  CorrBi_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  CorrBi_vs_Time_Graph->GetYaxis()->SetTitle("Peak Positions [ADC]");
  CorrBi_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  CorrBi_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  CorrBi_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  CorrBi_vs_Time_Graph->SetMarkerStyle(7);
  CorrBi_vs_Time_Graph->SetMarkerColor(kRed);
  CorrBi_vs_Time_Graph->Draw("AP");

  c10->cd(2);
  c10->SetGrid();
  CorrBi_vs_Temp_Graph->SetTitle("Corrected Bi vs T_DB");
  CorrBi_vs_Temp_Graph->GetXaxis()->SetRangeUser(T_DB_Min-0.2,T_DB_Max+0.2);
  CorrBi_vs_Temp_Graph->GetXaxis()->SetTitle("Temperature [Deg. F]");
  CorrBi_vs_Temp_Graph->GetXaxis()->SetTitleSize(0.08);
  CorrBi_vs_Temp_Graph->GetXaxis()->SetTitleOffset(0.9);
  CorrBi_vs_Temp_Graph->GetXaxis()->SetLabelSize(0.08);
  CorrBi_vs_Temp_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
  CorrBi_vs_Temp_Graph->GetYaxis()->SetTitleSize(0.08);
  CorrBi_vs_Temp_Graph->GetYaxis()->SetTitleOffset(0.4);
  CorrBi_vs_Temp_Graph->GetYaxis()->SetLabelSize(0.08);
  CorrBi_vs_Temp_Graph->SetMarkerStyle(7);
  CorrBi_vs_Temp_Graph->SetMarkerColor(kRed);
  CorrBi_vs_Temp_Graph->Draw("AP");
  c10->Print("Plots/Corr_Bi_Time&Temp.png"); 
 

  TCanvas *c11 = new TCanvas("c11","c11", 10,10,1950,600);
  c11->Divide(1,2);
  c11->cd(1);
  c11->SetGrid();
  CorrAm_vs_Time_Graph->SetTitle("Corrected Am vs Time");
  CorrAm_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  CorrAm_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  CorrAm_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  CorrAm_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  CorrAm_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  CorrAm_vs_Time_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
  CorrAm_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  CorrAm_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  CorrAm_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  CorrAm_vs_Time_Graph->SetMarkerStyle(7);
  CorrAm_vs_Time_Graph->SetMarkerColor(kRed);
  CorrAm_vs_Time_Graph->Draw("AP");

  c11->cd(2);
  c11->SetGrid();
  CorrAm_vs_Temp_Graph->SetTitle("Corrected Am vs T_Lab");
  CorrAm_vs_Temp_Graph->GetXaxis()->SetRangeUser(T_Lab2_Min-0.2,T_Lab2_Max+0.2);
  CorrAm_vs_Temp_Graph->GetXaxis()->SetTitle("Temperature [Deg. F]");
  CorrAm_vs_Temp_Graph->GetXaxis()->SetTitleSize(0.08);
  CorrAm_vs_Temp_Graph->GetXaxis()->SetTitleOffset(0.9);
  CorrAm_vs_Temp_Graph->GetXaxis()->SetLabelSize(0.08);
  CorrAm_vs_Temp_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
  CorrAm_vs_Temp_Graph->GetYaxis()->SetTitleSize(0.08);
  CorrAm_vs_Temp_Graph->GetYaxis()->SetTitleOffset(0.4);
  CorrAm_vs_Temp_Graph->GetYaxis()->SetLabelSize(0.08);
  CorrAm_vs_Temp_Graph->SetMarkerStyle(7);
  CorrAm_vs_Temp_Graph->SetMarkerColor(kRed);
  CorrAm_vs_Temp_Graph->Draw("AP");
  c11->Print("Plots/Corr_Am_Time&Temp.png");


  TCanvas *c12 = new TCanvas("c12","c12", 10,10,1950,600);
  c12->Divide(1,2);
  c12->cd(1);
  c12->SetGrid();
  CorrCalo_vs_Time_Graph->SetTitle("Corrected Calo LED vs Time");
  CorrCalo_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  CorrCalo_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  CorrCalo_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  CorrCalo_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  CorrCalo_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  CorrCalo_vs_Time_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
  CorrCalo_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  CorrCalo_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  CorrCalo_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  CorrCalo_vs_Time_Graph->SetMarkerStyle(7);
  CorrCalo_vs_Time_Graph->SetMarkerColor(kRed);
  CorrCalo_vs_Time_Graph->Draw("AP");

  c12->cd(2);
  c12->SetGrid();
  CorrCalo_vs_Temp_Graph->SetTitle("Corrected Calo LED vs T_Lab");
  CorrCalo_vs_Temp_Graph->GetXaxis()->SetRangeUser(T_Lab2_Min-0.2,T_Lab2_Max+0.2);
  CorrCalo_vs_Temp_Graph->GetXaxis()->SetTitle("Temperature [Deg. F]");
  CorrCalo_vs_Temp_Graph->GetXaxis()->SetTitleSize(0.08);
  CorrCalo_vs_Temp_Graph->GetXaxis()->SetTitleOffset(0.9);
  CorrCalo_vs_Temp_Graph->GetXaxis()->SetLabelSize(0.08);
  CorrCalo_vs_Temp_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
  CorrCalo_vs_Temp_Graph->GetYaxis()->SetTitleSize(0.08);
  CorrCalo_vs_Temp_Graph->GetYaxis()->SetTitleOffset(0.4);
  CorrCalo_vs_Temp_Graph->GetYaxis()->SetLabelSize(0.08);
  CorrCalo_vs_Temp_Graph->SetMarkerStyle(7);
  CorrCalo_vs_Temp_Graph->SetMarkerColor(kRed);
  CorrCalo_vs_Temp_Graph->Draw("AP");
  c12->Print("Plots/Corr_Calo_Time&Temp.png");


  TCanvas *c13 = new TCanvas("c13","c13", 10,10,1950,600);
  c13->Divide(1,2);
  c13->cd(1);
  c13->SetGrid();
  CorrRef_vs_Time_Graph->SetTitle("Corrected Ref LED vs Time");
  CorrRef_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  CorrRef_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  CorrRef_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  CorrRef_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  CorrRef_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  CorrRef_vs_Time_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
  CorrRef_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  CorrRef_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  CorrRef_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  CorrRef_vs_Time_Graph->SetMarkerStyle(7);
  CorrRef_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  CorrRef_vs_Time_Graph->SetMarkerColor(kRed);
  CorrRef_vs_Time_Graph->Draw("AP");

  c13->cd(2);
  c13->SetGrid();
  CorrRef_vs_Temp_Graph->SetTitle("Corrected Ref LED vs T_Lab");
  CorrRef_vs_Temp_Graph->GetXaxis()->SetRangeUser(T_Lab2_Min-0.2,T_Lab2_Max+0.2);
  CorrRef_vs_Temp_Graph->GetXaxis()->SetTitle("Temperature [Deg. F]");
  CorrRef_vs_Temp_Graph->GetXaxis()->SetTitleSize(0.08);
  CorrRef_vs_Temp_Graph->GetXaxis()->SetTitleOffset(0.9);
  CorrRef_vs_Temp_Graph->GetXaxis()->SetLabelSize(0.08);
  CorrRef_vs_Temp_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
  CorrRef_vs_Temp_Graph->GetYaxis()->SetTitleSize(0.08);
  CorrRef_vs_Temp_Graph->GetYaxis()->SetTitleOffset(0.4);
  CorrRef_vs_Temp_Graph->GetYaxis()->SetLabelSize(0.08);
  CorrRef_vs_Temp_Graph->SetMarkerStyle(7);
  CorrRef_vs_Temp_Graph->SetLineColorAlpha(kBlack, 0.6);
  CorrRef_vs_Temp_Graph->SetMarkerColor(kRed);
  CorrRef_vs_Temp_Graph->Draw("AP");
  c13->Print("Plots/Corr_Ref_Time&Temp.png");


 
  // Plotting
  // Measured 207Bi peak over time
  TCanvas *c1 = new TCanvas("c1","c1", 50,50,2200,400);
  c1->cd();
  c1->SetGrid();
  Bi_vs_Time_Graph->SetTitle("Fitted ^{207}Bi Peak Position Over Time");
  Bi_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  Bi_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  Bi_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  Bi_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  Bi_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  Bi_vs_Time_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
  Bi_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  Bi_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  Bi_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  Bi_vs_Time_Graph->SetMarkerStyle(7);
  Bi_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  Bi_vs_Time_Graph->SetMarkerColor(kRed);
  Bi_vs_Time_Graph->Draw("AP");
  c1->Print("Plots/Bismuth_Measured.png");
  
  // Predicted 207Bi peak over time
  TCanvas *c2 = new TCanvas("c2","c2", 50,50,2200,400);
  c2->cd();
  c2->SetGrid();
  Bi_Pred_vs_Time_Graph->SetTitle("Predicted ^{207}Bi Peak Position Over Time");
  Bi_Pred_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  Bi_Pred_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  Bi_Pred_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  Bi_Pred_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  Bi_Pred_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  Bi_Pred_vs_Time_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
  Bi_Pred_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  Bi_Pred_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  Bi_Pred_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  Bi_Pred_vs_Time_Graph->SetMarkerStyle(7);
  Bi_Pred_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  Bi_Pred_vs_Time_Graph->SetMarkerColor(kRed);
  Bi_Pred_vs_Time_Graph->Draw("AP");
  c2->Print("Plots/Bismuth_Predicted.png");

  // Money Plot: Ratio of prediction to measurement over time
  TCanvas *c3 = new TCanvas("c3","c3", 10,10,1950,600);
  int count_check = 0;
  for (int x = 0; x < newsize; x++) {
        if(Ratio_vs_Time_Graph->Eval(Time[x]) < .9) { 
                count_check = count_check+1;
                cout << x << " " ; 
                cout << fixed;
                cout << Time_Data.at(x) << ": " << Ratio_vs_Time_Graph->Eval(Time[x]) << endl;
        }   
  }
  cout << count_check << " + " << newsize <<  endl;
  
  c3->Divide(1,2);
  c3->cd(1);
  c3->SetGrid();
  Ratio_vs_Time_Graph->SetTitle("Ratio of Predicted to Measured ^{207}Bi Peak Over Time");
  Ratio_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  Ratio_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  Ratio_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  Ratio_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  Ratio_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  Ratio_vs_Time_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
  Ratio_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  Ratio_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  Ratio_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  Ratio_vs_Time_Graph->GetYaxis()->SetRangeUser(.985,1.015);
  Ratio_vs_Time_Graph->SetMarkerStyle(7);
  Ratio_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  Ratio_vs_Time_Graph->SetMarkerColor(kRed);
  Ratio_vs_Time_Graph->Draw("AP");
  if(DrawLines){
    TLine *upper1 = new TLine(0,1.01,Time[newsize-1],1.01);
    upper1->SetLineColor(kRed);
    upper1->SetLineWidth(2);
    upper1->Draw();
    TLine *lower1 = new TLine(0,0.99,Time[newsize-1],0.99);
    lower1->SetLineColor(kRed);
    lower1->SetLineWidth(2);
    lower1->Draw();
  }
  c3->cd(2);
  c3->SetGrid();
  T_DB_vs_Time_Graph->SetTitle("Temperature Over Time in Dark Box");
  T_DB_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  T_DB_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  T_DB_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  T_DB_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  T_DB_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  T_DB_vs_Time_Graph->GetYaxis()->SetTitle("Temperature [Deg. F]");
  T_DB_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  T_DB_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  T_DB_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  T_DB_vs_Time_Graph->SetMarkerStyle(7);
  T_DB_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  T_DB_vs_Time_Graph->SetMarkerColor(kRed);
  T_DB_vs_Time_Graph->Draw("AP"); 
  c3->Print("Plots/Ratio&Temp.png");

  // Plotting the time evolution of various ref OM peaks alongside temperature evolution
  TCanvas *c4 = new TCanvas("c4","c4", 10,10,2000,1200);
  c4->Divide(1,3);
  c4->cd(1);
  c4->SetGrid();
  LED_Calo_vs_Time_Graph->SetTitle("Calo LED vs Time");
  LED_Calo_vs_Time_Graph->GetXaxis()->SetRangeUser(0, Time[newsize-1]);
  LED_Calo_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  LED_Calo_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  LED_Calo_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  LED_Calo_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  LED_Calo_vs_Time_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
  LED_Calo_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  LED_Calo_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  LED_Calo_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  LED_Calo_vs_Time_Graph->SetMarkerStyle(7);
  LED_Calo_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  LED_Calo_vs_Time_Graph->SetMarkerColor(kRed);
  LED_Calo_vs_Time_Graph->Draw("AP");
  
  c4->cd(2);
  c4->SetGrid();
  LED_Ref_vs_Time_Graph->SetTitle("Fitted Reference LED Peak Position Over Time");
  LED_Ref_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  LED_Ref_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  LED_Ref_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  LED_Ref_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  LED_Ref_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  LED_Ref_vs_Time_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
  LED_Ref_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  LED_Ref_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  LED_Ref_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  LED_Ref_vs_Time_Graph->SetMarkerStyle(7);
  LED_Ref_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  LED_Ref_vs_Time_Graph->SetMarkerColor(kRed);
  LED_Ref_vs_Time_Graph->Draw("AP");
  if(TempInt){
    c4->cd(3);
    c4->SetGrid();
    T_Lab2_vs_Time_Graph->SetTitle("Temperature Over Time in Lab");
    T_Lab2_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
    T_Lab2_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    T_Lab2_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
    T_Lab2_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
    T_Lab2_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
    T_Lab2_vs_Time_Graph->GetYaxis()->SetTitle("Temperature [Deg. F]");
    T_Lab2_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
    T_Lab2_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
    T_Lab2_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
    T_Lab2_vs_Time_Graph->SetMarkerStyle(7);
    T_Lab2_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
    T_Lab2_vs_Time_Graph->SetMarkerColor(kRed);
    T_Lab2_vs_Time_Graph->Draw("AP");
  }

  c4->Print("Plots/LEDsand_Temp.png");

  TCanvas *c5 =new TCanvas("c5", "c5", 50,50,2200,400);

  c5->cd();
  c5->SetGrid();
  Am_vs_Time_Graph->SetTitle("Fitted ^{241}Am Peak Position Over Time");
  Am_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  Am_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  Am_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  Am_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  Am_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  Am_vs_Time_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
  Am_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  Am_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  Am_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  Am_vs_Time_Graph->SetMarkerStyle(7);
  Am_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  Am_vs_Time_Graph->SetMarkerColor(kRed);
  Am_vs_Time_Graph->Draw("AP");
  c5->Print("Plots/Am_vs_Time.png");
 /*
  LED_Calo_vs_Time_Graph->SetTitle("Calo LED vs Time");
  LED_Calo_vs_Time_Graph->GetXaxis()->SetRangeUser(0, Time[newsize-1]);
  LED_Calo_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  LED_Calo_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  LED_Calo_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  LED_Calo_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  LED_Calo_vs_Time_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
  LED_Calo_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  LED_Calo_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  LED_Calo_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  LED_Calo_vs_Time_Graph->SetMarkerStyle(7);
  LED_Calo_vs_Time_Graph->SetMarkerColor(kRed);
  LED_Calo_vs_Time_Graph->Draw("AP");
  */  
  TCanvas *c6 =new TCanvas("c6", "c6", 10,10,2000,1200);
  c6->Divide(1,3);
  c6->cd(1);
  c6->SetGrid();
  LED_Calo_vs_Time_Graph->SetTitle("Calo LED vs Time");
  LED_Calo_vs_Time_Graph->GetXaxis()->SetRangeUser(0, Time[newsize-1]);
  LED_Calo_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  LED_Calo_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  LED_Calo_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  LED_Calo_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  LED_Calo_vs_Time_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
  LED_Calo_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  LED_Calo_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  LED_Calo_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  LED_Calo_vs_Time_Graph->SetMarkerStyle(7);
  LED_Calo_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  LED_Calo_vs_Time_Graph->SetMarkerColor(kRed);
  LED_Calo_vs_Time_Graph->Draw("AP");

  c6->cd(2);
  c6->SetGrid();
  LED_Ref_vs_Time_Graph->SetTitle("Fitted Reference LED Peak Position Over Time");
  LED_Ref_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  LED_Ref_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  LED_Ref_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  LED_Ref_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  LED_Ref_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  LED_Ref_vs_Time_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
  LED_Ref_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  LED_Ref_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  LED_Ref_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  LED_Ref_vs_Time_Graph->SetMarkerStyle(7);
  LED_Ref_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  LED_Ref_vs_Time_Graph->SetMarkerColor(kRed);
  LED_Ref_vs_Time_Graph->Draw("AP");

    c6->cd(3);
    c6->SetGrid();
    T_Lab2_vs_Time_Graph->SetTitle("Temperature Over Time in Lab");
    T_Lab2_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
    T_Lab2_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    T_Lab2_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
    T_Lab2_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
    T_Lab2_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
    T_Lab2_vs_Time_Graph->GetYaxis()->SetTitle("Temperature [Deg. F]");
    T_Lab2_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
    T_Lab2_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
    T_Lab2_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
    T_Lab2_vs_Time_Graph->SetMarkerStyle(7);
    T_Lab2_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
    T_Lab2_vs_Time_Graph->SetMarkerColor(kRed);
    T_Lab2_vs_Time_Graph->Draw("AP");
  c6->Print("Plots/LED_Ref&Calo&Temp_vs_Time.png");
  
  
  
  TCanvas *c8 =new TCanvas("c8", "c8", 10,10,1950,600); 
  c8->Divide(1,2);
  c8->cd(1);
  c8->SetGrid();
  LED_Calo_vs_Temp_Graph->SetTitle("Calo LED vs Lab Temperature");
  LED_Calo_vs_Temp_Graph->GetXaxis()->SetRangeUser(70,80);
  LED_Calo_vs_Temp_Graph->GetXaxis()->SetTitle("Temperature [Deg. F]");
  LED_Calo_vs_Temp_Graph->GetXaxis()->SetTitleSize(0.08);
  LED_Calo_vs_Temp_Graph->GetXaxis()->SetTitleOffset(0.9);
  LED_Calo_vs_Temp_Graph->GetXaxis()->SetLabelSize(0.08);
  LED_Calo_vs_Temp_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
  LED_Calo_vs_Temp_Graph->GetYaxis()->SetTitleSize(0.08);
  LED_Calo_vs_Temp_Graph->GetYaxis()->SetTitleOffset(0.4);
  LED_Calo_vs_Temp_Graph->GetYaxis()->SetLabelSize(0.08);
  LED_Calo_vs_Temp_Graph->SetMarkerStyle(7);
  LED_Calo_vs_Temp_Graph->SetLineColorAlpha(kBlack, 0.6);
  LED_Calo_vs_Temp_Graph->SetMarkerColorAlpha(kBlack, 0.6);
  LED_Calo_vs_Temp_Graph->Draw("AP");
  
  c8->cd(2);
  c8->SetGrid();
  LED_Ref_vs_Temp_Graph->SetTitle("Ref LED vs Lab Temperature");
  LED_Ref_vs_Temp_Graph->GetXaxis()->SetRangeUser(70,80);
  LED_Ref_vs_Temp_Graph->GetXaxis()->SetTitle("Temperature [Deg. F]");
  LED_Ref_vs_Temp_Graph->GetXaxis()->SetTitleSize(0.08);
  LED_Ref_vs_Temp_Graph->GetXaxis()->SetTitleOffset(0.9);
  LED_Ref_vs_Temp_Graph->GetXaxis()->SetLabelSize(0.08);
  LED_Ref_vs_Temp_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
  LED_Ref_vs_Temp_Graph->GetYaxis()->SetTitleSize(0.08);
  LED_Ref_vs_Temp_Graph->GetYaxis()->SetTitleOffset(0.4);
  LED_Ref_vs_Temp_Graph->GetYaxis()->SetLabelSize(0.08);
  LED_Ref_vs_Temp_Graph->SetMarkerStyle(7);
  LED_Ref_vs_Temp_Graph->SetLineColorAlpha(kBlack, 0.6);
  LED_Ref_vs_Temp_Graph->SetMarkerColorAlpha(kBlack, 0.6);
  LED_Ref_vs_Temp_Graph->Draw("AP");
  c8->Print("Plots/LED_Ref&Calo_vs_Temp.png");

  
  TCanvas *c9 =new TCanvas("c9", "c9", 10,10,1950,600);;
  c9->Divide(1,2);
  c9->cd(1);
  c9->SetGrid();
  Am_vs_Temp_Graph->SetTitle("Am vs DB Temperature");
  Am_vs_Temp_Graph->GetXaxis()->SetRangeUser(70,80);
  Am_vs_Temp_Graph->GetXaxis()->SetTitle("Temperature [Deg. F]");
  Am_vs_Temp_Graph->GetXaxis()->SetTitleSize(0.08);
  Am_vs_Temp_Graph->GetXaxis()->SetTitleOffset(0.9);
  Am_vs_Temp_Graph->GetXaxis()->SetLabelSize(0.08);
  Am_vs_Temp_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
  Am_vs_Temp_Graph->GetYaxis()->SetTitleSize(0.08);
  Am_vs_Temp_Graph->GetYaxis()->SetTitleOffset(0.4);
  Am_vs_Temp_Graph->GetYaxis()->SetLabelSize(0.08);
  Am_vs_Temp_Graph->SetMarkerStyle(7);
  Am_vs_Temp_Graph->SetLineColorAlpha(kBlack, 0.6);
  Am_vs_Temp_Graph->SetMarkerColorAlpha(kBlack, 0.6);
  Am_vs_Temp_Graph->Draw("AP");
  
  c9->cd(2);
  c9->SetGrid();
  Bi_vs_Temp_Graph->SetTitle("Bi vs DB Temperature");
  Bi_vs_Temp_Graph->GetXaxis()->SetRangeUser(70,80);
  Bi_vs_Temp_Graph->GetXaxis()->SetTitle("Temperature [Deg. F]");
  Bi_vs_Temp_Graph->GetXaxis()->SetTitleSize(0.08);
  Bi_vs_Temp_Graph->GetXaxis()->SetTitleOffset(0.9);
  Bi_vs_Temp_Graph->GetXaxis()->SetLabelSize(0.08);
  Bi_vs_Temp_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
  Bi_vs_Temp_Graph->GetYaxis()->SetTitleSize(0.08);
  Bi_vs_Temp_Graph->GetYaxis()->SetTitleOffset(0.4);
  Bi_vs_Temp_Graph->GetYaxis()->SetLabelSize(0.08);
  Bi_vs_Temp_Graph->SetMarkerStyle(7);
  Bi_vs_Temp_Graph->SetLineColorAlpha(kBlack, 0.6);
  Bi_vs_Temp_Graph->SetMarkerColorAlpha(kBlack, 0.6);
  Bi_vs_Temp_Graph->Draw("AP");
  c9->Print("Plots/Am&Bi_vs_Temp.png");
 

  TCanvas *c15 =new TCanvas("c15", "c15", 10,10,2000,1200);
  c15->Divide(1,3); 
  c15->cd(1);
  c15->SetGrid();
  Bi_vs_Time_Graph->SetTitle("Fitted ^{207}Bi Peak Position Over Time");
  Bi_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  Bi_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  Bi_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  Bi_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  Bi_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  Bi_vs_Time_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
  Bi_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  Bi_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  Bi_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  Bi_vs_Time_Graph->SetMarkerStyle(7);
  Bi_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  Bi_vs_Time_Graph->SetMarkerColor(kRed);
  Bi_vs_Time_Graph->Draw("AP");
 
  c15->cd(2);
  c15->SetGrid();
  Am_vs_Time_Graph->SetTitle("Fitted ^{241}Am Peak Position Over Time");
  Am_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  Am_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  Am_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  Am_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  Am_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  Am_vs_Time_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
  Am_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  Am_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  Am_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  Am_vs_Time_Graph->SetMarkerStyle(7);
  Am_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  Am_vs_Time_Graph->SetMarkerColor(kRed);
  Am_vs_Time_Graph->Draw("AP");  
  
  c15->cd(3);
  c15->SetGrid();
  T_DB_vs_Time_Graph->SetTitle("Temperature Over Time in Dark Box");
  T_DB_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  T_DB_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  T_DB_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  T_DB_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  T_DB_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  T_DB_vs_Time_Graph->GetYaxis()->SetTitle("Temperature [Deg. F]");
  T_DB_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  T_DB_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  T_DB_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  T_DB_vs_Time_Graph->SetMarkerStyle(7);
  T_DB_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  T_DB_vs_Time_Graph->SetMarkerColor(kRed);
  T_DB_vs_Time_Graph->Draw("AP");
  c15->Print("Plots/Bi&Am&Temp_vs_Time.png");


  TCanvas *c16=new TCanvas("c16", "c16", 10, 10, 1950,600); //2100 , 600
  c16->Divide(1,2); 
  c16->cd(1);
  c16->SetGrid();
  Ratio_vs_Time_Graph->SetTitle("Ratio of Predicted to Measured ^{207}Bi Peak Over Time");
  Ratio_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  Ratio_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  Ratio_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  Ratio_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  Ratio_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  Ratio_vs_Time_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
  Ratio_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  Ratio_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  Ratio_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  Ratio_vs_Time_Graph->GetYaxis()->SetRangeUser(.985,1.015);
  Ratio_vs_Time_Graph->SetMarkerStyle(7);
  Ratio_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  Ratio_vs_Time_Graph->SetMarkerColor(kRed);
  Ratio_vs_Time_Graph->Draw("AP");
  if(DrawLines){
    TLine *upper2 = new TLine(0,1.01,Time[newsize-1],1.01);
    upper2->SetLineColor(kRed);
    upper2->SetLineWidth(2);
    upper2->Draw();
    TLine *lower2 = new TLine(0,0.99,Time[newsize-1],0.99);
    lower2->SetLineColor(kRed);
    lower2->SetLineWidth(2);
    lower2->Draw();
  }
  
  c16->cd(2);
  c16->SetGrid();
  Corr_Ratio_vs_Time_Graph->SetTitle("Temperature Corrected Ratio vs Time");
  Corr_Ratio_vs_Time_Graph->GetXaxis()->SetRangeUser(0, Time[newsize-1]);
  Corr_Ratio_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  Corr_Ratio_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  Corr_Ratio_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  Corr_Ratio_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  Corr_Ratio_vs_Time_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
  Corr_Ratio_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  Corr_Ratio_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  Corr_Ratio_vs_Time_Graph->GetYaxis()->SetRangeUser(0.985, 1.015);
  Corr_Ratio_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  Corr_Ratio_vs_Time_Graph->SetMarkerStyle(7);
  Corr_Ratio_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  Corr_Ratio_vs_Time_Graph->SetMarkerColor(kRed);
  Corr_Ratio_vs_Time_Graph->Draw("AP"); 
  if(DrawLines){
  TLine *upper3 = new TLine(0,1.01,Time[newsize-1],1.01);
    upper3->SetLineColor(kRed);
    upper3->SetLineWidth(2);
    upper3->Draw();
  TLine *lower3 = new TLine(0,0.99,Time[newsize-1],0.99);
    lower3->SetLineColor(kRed);
    lower3->SetLineWidth(2);
    lower3->Draw();
  }
  c16->Print("Plots/Ratio&CorrRatio_vs_Time.png");
  

  
  TCanvas *c17=new TCanvas("c17", "c17", 10, 10, 1950,600); //2100 , 600
  c17->Divide(1,2); 
  c17->cd(1);
  c17->SetGrid();
  Corr_Ratio_vs_Time_Graph->SetTitle("Temperature Corrected Ratio vs Time");
  Corr_Ratio_vs_Time_Graph->GetXaxis()->SetRangeUser(0, Time[newsize-1]);
  Corr_Ratio_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  Corr_Ratio_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  Corr_Ratio_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  Corr_Ratio_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  Corr_Ratio_vs_Time_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
  Corr_Ratio_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  Corr_Ratio_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  Corr_Ratio_vs_Time_Graph->GetYaxis()->SetRangeUser(0.994, 1.004);
  Corr_Ratio_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  Corr_Ratio_vs_Time_Graph->SetMarkerStyle(7);
  Corr_Ratio_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  Corr_Ratio_vs_Time_Graph->SetMarkerColor(kRed);
  Corr_Ratio_vs_Time_Graph->Draw("AP");
  if(DrawLines){
  TLine *upper4 = new TLine(0,1.01,Time[newsize-1],1.01);
    upper4->SetLineColor(kRed);
    upper4->SetLineWidth(2);
    upper4->Draw();
  TLine *lower4 = new TLine(0,0.99,Time[newsize-1],0.99);
    lower4->SetLineColor(kRed);
    lower4->SetLineWidth(2);
    lower4->Draw();
  }
  
  c17->cd(2);
  c17->SetGrid();
  T_DB_vs_Time_Graph->SetTitle("Temperature Over Time");
  T_DB_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  T_DB_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  T_DB_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  T_DB_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  T_DB_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  T_DB_vs_Time_Graph->GetYaxis()->SetTitle("Temperature [Deg. F]");
  T_DB_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  T_DB_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  T_DB_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  T_DB_vs_Time_Graph->SetMarkerStyle(7);
  T_DB_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  T_DB_vs_Time_Graph->SetMarkerColor(kRed);
  T_DB_vs_Time_Graph->Draw("AP");
  c17->Print("Plots/CorrRatio&Temp_vs_Time.png");  

  
  TCanvas *c18 = new TCanvas("c18","c18", 10,10,2000,1200);
  c18->Divide(1,3);
  c18->cd(1);
  CorrCalo_vs_Time_Graph->SetTitle("Corrected Calorimeter LED Peak Position Over Time");
  CorrCalo_vs_Time_Graph->GetXaxis()->SetRangeUser(0, Time[newsize-1]);
  CorrCalo_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  CorrCalo_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  CorrCalo_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  CorrCalo_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  CorrCalo_vs_Time_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
  CorrCalo_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  CorrCalo_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  CorrCalo_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  CorrCalo_vs_Time_Graph->SetMarkerStyle(7);
  CorrCalo_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  CorrCalo_vs_Time_Graph->SetMarkerColor(kRed);
  CorrCalo_vs_Time_Graph->Draw("AP");

  c18->cd(2);
  c18->SetGrid();
  CorrRef_vs_Time_Graph->SetTitle("Corrected Reference LED Peak Position Over Time");
  CorrRef_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  CorrRef_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  CorrRef_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  CorrRef_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  CorrRef_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  CorrRef_vs_Time_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
  CorrRef_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  CorrRef_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  CorrRef_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  CorrRef_vs_Time_Graph->SetMarkerStyle(7);
  CorrRef_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  CorrRef_vs_Time_Graph->SetMarkerColor(kRed);
  CorrRef_vs_Time_Graph->Draw("AP");
  if(TempInt){
    c18->cd(3);
    c18->SetGrid();
    T_Lab2_vs_Time_Graph->SetTitle("Temperature Over Time in Lab");
    T_Lab2_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
    T_Lab2_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
    T_Lab2_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
    T_Lab2_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
    T_Lab2_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
    T_Lab2_vs_Time_Graph->GetYaxis()->SetTitle("Temperature [Deg. F]");
    T_Lab2_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
    T_Lab2_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
    T_Lab2_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
    T_Lab2_vs_Time_Graph->SetMarkerStyle(7);
    T_Lab2_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
    T_Lab2_vs_Time_Graph->SetMarkerColor(kRed);
    T_Lab2_vs_Time_Graph->Draw("AP");
  }

  c18->Print("Plots/CorrLEDsand_Temp.png");  

  TCanvas *c19 =new TCanvas("c19", "c19", 10,10,2000,1200);
  c19->Divide(1,3);
  c19->cd(1);
  c19->SetGrid();
  CorrBi_vs_Time_Graph->SetTitle("Corrected ^{207}Bi Peak Position Over Time");
  CorrBi_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  CorrBi_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  CorrBi_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  CorrBi_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  CorrBi_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  CorrBi_vs_Time_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
  CorrBi_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  CorrBi_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  CorrBi_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  CorrBi_vs_Time_Graph->SetMarkerStyle(7);
  CorrBi_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  CorrBi_vs_Time_Graph->SetMarkerColor(kRed);
  CorrBi_vs_Time_Graph->Draw("AP");

  c19->cd(2);
  c19->SetGrid();
  CorrAm_vs_Time_Graph->SetTitle("Corrected ^{241}Am Peak Position Over Time");
  CorrAm_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  CorrAm_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  CorrAm_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  CorrAm_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  CorrAm_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  CorrAm_vs_Time_Graph->GetYaxis()->SetTitle("Peak Position [ADC]");
  CorrAm_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  CorrAm_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  CorrAm_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  CorrAm_vs_Time_Graph->SetMarkerStyle(7);
  CorrAm_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  CorrAm_vs_Time_Graph->SetMarkerColor(kRed);
  CorrAm_vs_Time_Graph->Draw("AP");
  
  c19->cd(3);
  c19->SetGrid();
  T_DB_vs_Time_Graph->SetTitle("Temperature Over Time in Dark Box");
  T_DB_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  T_DB_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  T_DB_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  T_DB_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  T_DB_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  T_DB_vs_Time_Graph->GetYaxis()->SetTitle("Temperature [Deg. F]");
  T_DB_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  T_DB_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  T_DB_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  T_DB_vs_Time_Graph->SetMarkerStyle(7);
  T_DB_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  T_DB_vs_Time_Graph->SetMarkerColor(kRed);
  T_DB_vs_Time_Graph->Draw("AP");
  c19->Print("Plots/CorrBi&Am&Temp_vs_Time.png");


  TCanvas *c20 =new TCanvas("c20", "c20", 10,10,2000,1200);
  c20->Divide(1,3);
  c20->cd(1);
  c20->SetGrid();
  Ratio_vs_Time_Graph->SetTitle("Ratio of Predicted to Measured ^{207}Bi Peak Over Time");
  Ratio_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  Ratio_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  Ratio_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  Ratio_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  Ratio_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  Ratio_vs_Time_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
  Ratio_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  Ratio_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  Ratio_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  Ratio_vs_Time_Graph->SetMarkerStyle(7);
  Ratio_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  Ratio_vs_Time_Graph->SetMarkerColor(kRed);
  Ratio_vs_Time_Graph->Draw("AP");
  if(DrawLines){
  TLine *upper5 = new TLine(0,1.01,Time[newsize-1],1.01);
    upper5->SetLineColor(kRed);
    upper5->SetLineWidth(2);
    upper5->Draw();
  TLine *lower5 = new TLine(0,0.99,Time[newsize-1],0.99);
    lower5->SetLineColor(kRed);
    lower5->SetLineWidth(2);
    lower5->Draw();
  }

  
  c20->cd(2);
  c20->SetGrid();
  T_DB_vs_Time_Graph->SetTitle("Temperature Over Time in Dark Box");
  T_DB_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  T_DB_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  T_DB_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  T_DB_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9);
  T_DB_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  T_DB_vs_Time_Graph->GetYaxis()->SetTitle("Temperature [Deg. F]");
  T_DB_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  T_DB_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  T_DB_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  T_DB_vs_Time_Graph->SetMarkerStyle(7);
  T_DB_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  T_DB_vs_Time_Graph->SetMarkerColor(kRed);
  T_DB_vs_Time_Graph->Draw("AP");

  c20->cd(3);
  c20->SetGrid();
  Corr_Ratio_vs_Time_Graph->SetTitle("Temperature Corrected Ratio vs Time");
  Corr_Ratio_vs_Time_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  Corr_Ratio_vs_Time_Graph->GetXaxis()->SetTitle("Time Since Start [Days]");
  Corr_Ratio_vs_Time_Graph->GetXaxis()->SetTitleSize(0.08);
  Corr_Ratio_vs_Time_Graph->GetXaxis()->SetTitleOffset(0.9); // 0.8
  Corr_Ratio_vs_Time_Graph->GetXaxis()->SetLabelSize(0.08);
  Corr_Ratio_vs_Time_Graph->GetYaxis()->SetRangeUser(0.985,1.015);
  Corr_Ratio_vs_Time_Graph->GetYaxis()->SetTitle("MeV Peak Position [ADC]");
  Corr_Ratio_vs_Time_Graph->GetYaxis()->SetTitleSize(0.08);
  Corr_Ratio_vs_Time_Graph->GetYaxis()->SetTitleOffset(0.4);
  Corr_Ratio_vs_Time_Graph->GetYaxis()->SetLabelSize(0.08);
  Corr_Ratio_vs_Time_Graph->SetMarkerStyle(7);
  Corr_Ratio_vs_Time_Graph->SetLineColorAlpha(kBlack, 0.6);
  Corr_Ratio_vs_Time_Graph->SetMarkerColor(kRed);
  Corr_Ratio_vs_Time_Graph->Draw("AP");
  if(DrawLines){
  TLine *upper6 = new TLine(0,1.01,Time[newsize-1],1.01);
    upper6->SetLineColor(kRed);
    upper6->SetLineWidth(2);
    upper6->Draw();
  TLine *lower6 = new TLine(0,0.99,Time[newsize-1],0.99);
    lower6->SetLineColor(kRed);
    lower6->SetLineWidth(2);
    lower6->Draw();
  }
  c20->Print("Plots/CorrRatio&Ratio&Temp_vs_Time.png");

  
  TCanvas *c21 =new TCanvas("c21", "c21", 10,10,1950,600);;
  c21->Divide(1,2);
  c21->cd(1);
  c21->SetGrid();
  CorrBi_Line_Ratio_Graph->SetTitle("Conditional Bi Avg and Temp Fit Line Ratio vs Time");
  CorrBi_Line_Ratio_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  CorrBi_Line_Ratio_Graph->GetXaxis()->SetTitle("Time [Days]");
  CorrBi_Line_Ratio_Graph->GetXaxis()->SetTitleSize(0.08);
  CorrBi_Line_Ratio_Graph->GetXaxis()->SetTitleOffset(0.9);
  CorrBi_Line_Ratio_Graph->GetXaxis()->SetLabelSize(0.08);
  CorrBi_Line_Ratio_Graph->GetYaxis()->SetTitle("Peak Ratio [ADC/ADC]");
  CorrBi_Line_Ratio_Graph->GetYaxis()->SetTitleSize(0.08);
  CorrBi_Line_Ratio_Graph->GetYaxis()->SetTitleOffset(0.4);
  CorrBi_Line_Ratio_Graph->GetYaxis()->SetLabelSize(0.08);
  CorrBi_Line_Ratio_Graph->SetMarkerStyle(7);
  CorrBi_Line_Ratio_Graph->SetLineColorAlpha(kBlack, 0.6);
  CorrBi_Line_Ratio_Graph->SetMarkerColorAlpha(kBlack, 0.9);
  CorrBi_Line_Ratio_Graph->Draw("AP");
  
  c21->cd(2);
  c21->SetGrid();
  CorrAm_Line_Ratio_Graph->SetTitle("Conditional Am Avg and Temp Fit Line Ratio vs Time");
  CorrAm_Line_Ratio_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  CorrAm_Line_Ratio_Graph->GetXaxis()->SetTitle("Time [Days]");
  CorrAm_Line_Ratio_Graph->GetXaxis()->SetTitleSize(0.08);
  CorrAm_Line_Ratio_Graph->GetXaxis()->SetTitleOffset(0.9);
  CorrAm_Line_Ratio_Graph->GetXaxis()->SetLabelSize(0.08);
  CorrAm_Line_Ratio_Graph->GetYaxis()->SetTitle("Peak Ratio [ADC/ADC]");
  CorrAm_Line_Ratio_Graph->GetYaxis()->SetTitleSize(0.08);
  CorrAm_Line_Ratio_Graph->GetYaxis()->SetTitleOffset(0.4);
  CorrAm_Line_Ratio_Graph->GetYaxis()->SetLabelSize(0.08);
  CorrAm_Line_Ratio_Graph->SetMarkerStyle(7);
  CorrAm_Line_Ratio_Graph->SetLineColorAlpha(kBlack, 0.6);
  CorrAm_Line_Ratio_Graph->SetMarkerColorAlpha(kBlack, 0.9);
  CorrAm_Line_Ratio_Graph->Draw("AP");
  c21->Print("Plots/Bi&Am_Ratio.png");

  
  TCanvas *c22 =new TCanvas("c22", "c22", 10,10,1950,600);;
  c22->Divide(1,2);
  c22->cd(1);
  c22->SetGrid();
  CorrCalo_Line_Ratio_Graph->SetTitle("Conditional Calo LED Avg and Temp Fit Line Ratio vs Time");
  CorrCalo_Line_Ratio_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  CorrCalo_Line_Ratio_Graph->GetXaxis()->SetTitle("Time [Days]");
  CorrCalo_Line_Ratio_Graph->GetXaxis()->SetTitleSize(0.08);
  CorrCalo_Line_Ratio_Graph->GetXaxis()->SetTitleOffset(0.9);
  CorrCalo_Line_Ratio_Graph->GetXaxis()->SetLabelSize(0.08);
  CorrCalo_Line_Ratio_Graph->GetYaxis()->SetTitle("Peak Ratio [ADC/ADC]");
  CorrCalo_Line_Ratio_Graph->GetYaxis()->SetTitleSize(0.08);
  CorrCalo_Line_Ratio_Graph->GetYaxis()->SetTitleOffset(0.4);
  CorrCalo_Line_Ratio_Graph->GetYaxis()->SetLabelSize(0.08);
  CorrCalo_Line_Ratio_Graph->SetMarkerStyle(7);
  CorrCalo_Line_Ratio_Graph->SetLineColorAlpha(kBlack, 0.6);
  CorrCalo_Line_Ratio_Graph->SetMarkerColorAlpha(kBlack, 0.9);
  CorrCalo_Line_Ratio_Graph->Draw("AP");
  
  c22->cd(2);
  c22->SetGrid();
  CorrRef_Line_Ratio_Graph->SetTitle("Conditional Ref LED Avg and Temp Fit Line Ratio vs Time");
  CorrRef_Line_Ratio_Graph->GetXaxis()->SetRangeUser(0,Time[newsize-1]);
  CorrRef_Line_Ratio_Graph->GetXaxis()->SetTitle("Time [days]");
  CorrRef_Line_Ratio_Graph->GetXaxis()->SetTitleSize(0.08);
  CorrRef_Line_Ratio_Graph->GetXaxis()->SetTitleOffset(0.9);
  CorrRef_Line_Ratio_Graph->GetXaxis()->SetLabelSize(0.08);
  CorrRef_Line_Ratio_Graph->GetYaxis()->SetTitle("Peak Ratio [ADC/ADC]");
  CorrRef_Line_Ratio_Graph->GetYaxis()->SetTitleSize(0.08);
  CorrRef_Line_Ratio_Graph->GetYaxis()->SetTitleOffset(0.4);
  CorrRef_Line_Ratio_Graph->GetYaxis()->SetLabelSize(0.08);
  CorrRef_Line_Ratio_Graph->SetMarkerStyle(7);
  CorrRef_Line_Ratio_Graph->SetLineColorAlpha(kBlack, 0.6);
  CorrRef_Line_Ratio_Graph->SetMarkerColorAlpha(kBlack, 0.9);
  CorrRef_Line_Ratio_Graph->Draw("AP");
  c22->Print("Plots/Calo&Ref_Ratio.png");
  
  }  
  return(23);
}


Double_t average(Double_t numbers[], int numb) {
    Double_t sum = 0;
    for (int x = 0; x < numb; x++)
    {
        sum += numbers[x];
    }
    return sum/(Double_t)numb;
}



 
 
