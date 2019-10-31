#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TMath.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <typeinfo>

using namespace std;

Double_t threegaussandledfit(Double_t *x, Double_t *par);

int Fitter()
{
  Bool_t Write = 1; // Toggle writing to parameter files  *Note: Must be ON when running Run_Fitter.sh and OFF when checking individual root file*
  Bool_t SavePlots = 0; // Toggle saving of plots 
  Bool_t ExtraFits = 0; // Toggle canvases for extra fits

 
  //******************//
  // Inititalizations //
  //******************//

  TCanvas *c2, *c3, *c5, *c6;

  // Some canvases that will be used for debugging and checking results
  TCanvas *c1 = new TCanvas("c1","c1",50,50,600,400); // for Source_Calo and LED_Calo
  if(ExtraFits){c2 = new TCanvas("c2","c2",350,50,600,400);} // for hmev2 and mev2
  if(ExtraFits){c3 = new TCanvas("c3","c3",650,50,600,400);} // for hmev3 and mev3
  TCanvas *c4 = new TCanvas("c4","c4",50,50,600,400); // for Source_Ref and LED_Ref 
  if(ExtraFits){c5 = new TCanvas("c5","c5",350,450,600,400);} // for am2
  if(ExtraFits){c6 = new TCanvas("c6","c6",650,450,600,400);} // for am3

  // Pull info from TTrees from file "f" 
  TFile *f = new TFile("temp.root"); // "temp.root" when running Run_Fitter.sh
  TTree *adctree = (TTree *)f->Get("ADC_VALUES");
  const int NBRCH = 12;
  int adcruns[NBRCH] = {0};
  for(int i = 0; i < NBRCH; i++) {
    adctree->SetBranchAddress(Form("ADC_Runs_%d", i),&adcruns[i]);
  }
  // Create histograms
  TH1F *ADC4 = (TH1F*)f->Get("ADC_4");
  TH1F *ADC5 = (TH1F*)f->Get("ADC_5");
  TH1F *ADC4ped = new TH1F("ADC4ped","ADC4ped",2048,0,2047);
  TH1F *ADC5ped = new TH1F("ADC5ped","ADC5ped",2048,0,2047);
  TH1F* hPedSubADC4 = new TH1F("hPedSubADC4","",2048,0,2047);
  TH1F* hPedSubADC5 = new TH1F("hPedSubADC5","",2048,0,2047);
  int z = 0;

  // Find Pedestals 
  while(adctree->GetEntry(z)) {
    if(adcruns[5]>105 && adcruns[5]<700 && adcruns[4]<70) {
      ADC4ped->Fill(adcruns[4]);
    }
    if (adcruns[4]>75 && adcruns[4]<700 && adcruns[5]<100) {
      ADC5ped->Fill(adcruns[5]);
    }
    z++; 
  }
  ADC5ped->GetXaxis()->SetRangeUser(5, 150);
  Double_t mean4 = ADC4ped->GetMaximumBin();
  Double_t mean5 = ADC5ped->GetMaximumBin();
  std::cout << "ped adc4 = " << mean4 << " | ped adc5 = " << mean5 << std::endl;
  Double_t contents = 0;
  
  // Correct ADC data for pedestals
  for(Int_t i = 1; i < 2047; i++) {
    Double_t rawADC4 = ADC4->GetBinCenter(i);
    Double_t pedSubADC4 = rawADC4 - mean4;
    if (pedSubADC4>0) {
      Int_t pedSubBin = hPedSubADC4->FindBin(pedSubADC4);
      Double_t y1 = ADC4->GetBinContent(i);
      hPedSubADC4->SetBinContent(pedSubBin, y1);
      contents+=ADC4->GetBinContent(i);
    }
    Double_t rawADC5 = ADC5->GetBinCenter(i);
    Double_t pedSubADC5 = rawADC5-mean5;
    if (pedSubADC5>0){
        Int_t pedSubBin2 = hPedSubADC5->FindBin(pedSubADC5);
        Double_t y2 = ADC5->GetBinContent(i);
        hPedSubADC5->SetBinContent(pedSubBin2, y2);
    }
  }
  // Clone histograms for extra fits
  TH1D* hPedSubADC4_2 = (TH1D*)hPedSubADC4->Clone("hPedSubADC4_2");
  TH1D* hPedSubADC5_2 = (TH1D*)hPedSubADC5->Clone("hPedSubADC4_2");
  TH1D* hPedSubADC4_3 = (TH1D*)hPedSubADC4->Clone("hPedSubADC4_3");
  TH1D* hPedSubADC5_3 = (TH1D*)hPedSubADC5->Clone("hPedSubADC4_3");

  
  //***************************************//
  // Time to find peaks from Calo Spectrum //
  //***************************************//

  TSpectrum *sp4 = new TSpectrum();
  sp4->Search(hPedSubADC4, 4.0, "");
  int npeaks4 = sp4->GetNPeaks();
  cout << "Number of peaks found in ADC4: " << npeaks4 << endl;
  Double_t *peaks4;
  peaks4 = sp4->GetPositionX();
  sort(peaks4, peaks4 + npeaks4);
  for (int j = 0; j < npeaks4; j++) {
    cout << "Peak " << j+1 << " at " << peaks4[j] << endl; 
  }

  // Series of if statements to correctly identify peaks. These were determined via trial and error. 
  Float_t hmev_peak, mev_peak, led_peak;
  if (npeaks4 == 7) { //reading 2 double peaks
    if (abs(peaks4[1] - peaks4[2]) < 25 && abs(peaks4[5] - peaks4[6]) < 100 && !(abs(peaks4[3] - peaks4[4]) < 20)) { //reading double peak on LED & peak before hmev and we are sure that there is no double peak on hmev
      hmev_peak = peaks4[3];
      mev_peak = peaks4[4];
      led_peak = (peaks4[5]+peaks4[6])/2;
      cout << "7 peaks: 1 option" << endl;
    }
    else {
    hmev_peak = (peaks4[3] + peaks4[4])/2;
    mev_peak = peaks4[5];
    led_peak = peaks4[6];
    cout << "7 peaks: 2 option" << endl;
    } 
  }
  else if (npeaks4 == 6) {
    if (peaks4[4] > 400) { //reading double peak on LED
     hmev_peak = peaks4[2];
     mev_peak = peaks4[3];
     led_peak = (peaks4[4] + peaks4[5])/2;
     cout << "6 peaks: " << "first option" << endl;
    }
    else if (peaks4[1] - peaks4[2] < 20) { // reading double peak on peak before hmev
      hmev_peak = peaks4[3];
      mev_peak = peaks4[4];
      led_peak = peaks4[5];
      cout << "6 peaks: " << "second option" <<  endl;
   } 
  }
  else if (peaks4[0]<5 && peaks4[1]>40){
    hmev_peak = peaks4[1];
    mev_peak = peaks4[2];
    led_peak = peaks4[3];
    cout << "1" << endl;
  } 
  // If there are 5 peaks and both first peaks are too small then
  // [2], [3], and [4] must be the hmev, mev, and LED respectively
  else if (peaks4[0]<5 && peaks4[1]<40){
    hmev_peak = peaks4[2];
    mev_peak = peaks4[3];
    led_peak = peaks4[4];
    cout << "2" << endl;
  } 
  
  else if (peaks4[0]>5 && peaks4[0]<20){
    hmev_peak = peaks4[1];
    mev_peak = peaks4[2];
    led_peak = peaks4[3];
    cout << "3" << endl;
  } 
  
  else if (peaks4[0]>20 && peaks4[0]<30){
    hmev_peak = peaks4[1];
    mev_peak = peaks4[2];
    led_peak = peaks4[3];
    cout << "4" << endl;
  } 
  // If the first peak is greater than 30 ADC then it's probably the hmev peak then
  // [0], [1], and [2] must be the hmev, mev, and LED respectively
  else if (peaks4[0]>30){
    hmev_peak = peaks4[0];
    mev_peak = peaks4[1];
    led_peak = peaks4[2];
    cout << "5" << endl;
  }

  // Print out Mean Peak Values 
  cout << "***********************" << endl;
  cout <<" base peak = "<<peaks4[0]<<endl;
  cout <<" hmev peak = "<<hmev_peak<<endl;
  cout <<" mev peak = "<<mev_peak<<endl;
  cout <<" led peak = "<<led_peak<<endl;
  cout << "***********************" << endl;


  //*********************************************//  
  // Bismuth Fit which calls threegaussandledfit //
  //*********************************************//

  TF1 *gauss = new TF1("gauss", "gaus", hmev_peak-5, hmev_peak+5);
  hPedSubADC4->Fit("gauss", "RQS0");
  double bkg_scaler = gauss->Eval(hmev_peak);
  cout << "***********************" << endl;
  cout <<" Background scaler: "<<bkg_scaler<<endl;
  cout << "***********************" << endl; 
  TF1 *Source_Calo = new TF1("Source_Calo", threegaussandledfit, hmev_peak-26.,2.*mev_peak-32.4,13); // hmev_peak-30., mev_peak+40, 13
  TF1 *LED_Calo = new TF1("LED_Calo", "gaus", led_peak-22, led_peak+22);
  Double_t myParameters[13];
  myParameters[0]=contents/1e4;                   // normalisation start value
  myParameters[1]=mev_peak;                       // mean of 1MeV_ peak in ADC
  myParameters[2]=1.0;                            // width in some units??
  myParameters[3]=1.6*myParameters[0];            // relative normalisation of 2nd peak
  myParameters[4]=hmev_peak;                      // mean of 0.5MeV_ peak
  myParameters[5]=60.0;                           // bkg decay const (bkg: exp(-x[0]/[5])*[6]+[7])
  myParameters[6]=bkg_scaler;                     // bkg amplitude
  myParameters[7]=-5.0;                           // bkg const
  myParameters[8]=0.3*myParameters[0];            // relative normalisation of 3rd peak
  myParameters[9]=hmev_peak+(mev_peak-hmev_peak); // mean of 1.6MeV_ peak (relative to mev and hmev)
  myParameters[10]=500.0;                         // amplitude of LED peak if included in fit
  myParameters[11]=led_peak;                      // mean of LED peak if included in fit
  myParameters[12]=20.0;

  
  //************************************************//
  // Americium Fits and Reduced Range Gaussian Fits //
  //************************************************//

  // Gaussian Fit for Finding Mean peak of Am 
  TF1 *americium_prefit = new TF1("americium_prefit", "gaus", 50, 120); //80, 200
  hPedSubADC5->Fit("americium_prefit", "RQ0");
  hPedSubADC5->Draw();
  double ampeak = americium_prefit->GetMaximumX();
  cout << "***********************" << endl;
  cout <<" Americium peak estimate: "<<ampeak<<endl;
  cout << "***********************" << endl;
  

  TF1 *Source_Ref = new TF1("Source_Ref","[3]*TMath::Gaus(x,[0],[1])*TMath::Landau(x,[0],[2])+exp(-x/[4])*[5]",ampeak-80,ampeak+110); 
  Source_Ref->SetParameter(0,ampeak);   //mean
  Source_Ref->SetParameter(1,10.);       //sigma gauss
  Source_Ref->SetParameter(2,10.);       //sigma landau
  Source_Ref->SetParameter(3,100.);      //amplitude 100
  Source_Ref->SetParameter(4,1000.);     //exp decay constant  (bkg)
  Source_Ref->SetParameter(5,5.);        //exp amplitude       (bkg)
  //Source_Ref->SetParameter(6,0.);        //constant            (bkg)
  
  // Gaussian fit to find mean peak for ref LED
  TF1 *LED_Ref = new TF1("LED_Ref", "gaus", 250, 330);

  // Reduced Range Guassian Fits
  TF1 *hmev2 = new TF1("hmev2","gaus",hmev_peak-17.,hmev_peak+17.);
  TF1 *mev2 = new TF1("mev2","gaus",mev_peak-35.,mev_peak+30.);
  
  // Just a slightly narrower fit range
  TF1 *hmev3 = new TF1("hmev3","gaus",hmev_peak-20.,hmev_peak+15.);
  TF1 *mev3 = new TF1("mev3","gaus",mev_peak-25.,mev_peak+20.);
  
  // Repeat for 241Am_ peak
  TF1 *am2 = new TF1("am2","gaus",ampeak+25.,ampeak-25.);
  TF1 *am3 = new TF1("am3","gaus",ampeak+10.,ampeak-10.);


  //******************************************************************//
  // Draw Extra Gaussian Fits for Calo and Ref OM and Save Parameters //
  //******************************************************************//
  
  // Draw Gaussian fits for Calo Source and save parameters
  if(ExtraFits){c2->cd();}
  hmev2->SetLineColor(kRed+2);
  mev2->SetLineColor(kRed+2);
  hPedSubADC4_2->Fit("hmev2", "RQ");
  hPedSubADC4_2->Fit("mev2", "RQ+");
  hPedSubADC4_2->GetXaxis()->SetRangeUser(5, 200);
  if(ExtraFits){hPedSubADC4_2->Draw();}

  Double_t HMeV2_Mean = hmev2->GetParameter(1);
  Double_t HMeV2_MeanErr = hmev2->GetParError(1);
  Double_t HMeV2_Sigma = hmev2->GetParameter(2);
  Double_t HMeV2_SigmaErr = hmev2->GetParError(2);
  Double_t MeV2_Mean = mev2->GetParameter(1);
  Double_t MeV2_MeanErr = mev2->GetParError(1);
  Double_t MeV2_Sigma = mev2->GetParameter(2);
  Double_t MeV2_SigmaErr = mev2->GetParError(2);
  Double_t HMeV2_chi2 = hmev2->GetChisquare();
  Double_t HMeV2_NDF = hmev2->GetNDF();
  Double_t MeV2_chi2 = mev2->GetChisquare();
  Double_t MeV2_NDF = mev2->GetNDF(); 
   
  if(ExtraFits){c3->cd();}
  hmev3->SetLineColor(kRed+4);
  mev3->SetLineColor(kRed+4);
  hPedSubADC4_3->Fit("hmev3", "RQ");
  hPedSubADC4_3->Fit("mev3","RQ+");
  hPedSubADC4_3->GetXaxis()->SetRangeUser(5,200);
  if(ExtraFits){hPedSubADC4_3->Draw();}
  
  Double_t HMeV3_Mean = hmev3->GetParameter(1);
  Double_t HMeV3_MeanErr = hmev3->GetParError(1);
  Double_t HMeV3_Sigma = hmev3->GetParameter(2);
  Double_t HMeV3_SigmaErr = hmev3->GetParError(2);
  Double_t MeV3_Mean = mev3->GetParameter(1);
  Double_t MeV3_MeanErr = mev3->GetParError(1);
  Double_t MeV3_Sigma = mev3->GetParameter(2);
  Double_t MeV3_SigmaErr = mev3->GetParError(2);

  // Draw Gaussian Fits for Ref Source and Save Parameters
  if(ExtraFits){c5->cd();}
  am2->SetLineColor(kRed);
  hPedSubADC5_2->Fit("am2", "RQ");
  hPedSubADC5_2->GetXaxis()->SetRangeUser(5,200);
  if(ExtraFits){hPedSubADC5_2->Draw("same");}

  Double_t Am2_Mean = am2->GetParameter(1);
  Double_t Am2_MeanErr = am2->GetParError(1);
  Double_t Am2_Sigma = am2->GetParameter(2);
  Double_t Am2_SigmaErr = am2->GetParError(2);
  Double_t Am2_chi2 = am2->GetChisquare();
  Double_t Am2_NDF = am2->GetNDF();

  if(ExtraFits){c6->cd();}
  am3->SetLineColor(kRed+4);
  hPedSubADC5_3->Fit("am3","RQ");
  hPedSubADC5_3->GetXaxis()->SetRangeUser(5,200);
  if(ExtraFits){hPedSubADC5_3->Draw("same");}

  Double_t Am3_Mean = am3->GetParameter(1);
  Double_t Am3_MeanErr = am3->GetParError(1);
  Double_t Am3_Sigma = am3->GetParameter(2);
  Double_t Am3_SigmaErr = am3->GetParError(2);
 
 
  //***********************************************************//
  // Draw Full Calo and Reference PMT Fits and Save Parameters // Note: These are after the extra fits due to case of extra fits drawing instead of full
  //***********************************************************//
  
  // Draw full Calo PMT fits 
  c1->cd();
  Source_Calo->SetParameters(myParameters);
  Source_Calo->SetLineColor(kRed);
  LED_Calo->SetLineColor(kBlue);
  //hPedSubADC4->Draw();
  hPedSubADC4->Fit("Source_Calo", "RS");
  hPedSubADC4->Fit("LED_Calo", "RQ+");
  hPedSubADC4->GetXaxis()->SetRangeUser(5, 360);
  hPedSubADC4->Draw();
  if (SavePlots){
  c1->Print("Plots/Bi_LED_fits.png");
  }

  // Check Background Fit with Source Fit
  //hPedSubADC4->GetXaxis()->SetTitle("ADC");
  //hPedSubADC4->GetXaxis()->SetTitleSize(0.08);
  //gStyle->SetOptFit(1); 
  // Study of the Calo OM background fit function
  //TF1 *Bkg_Calo = new TF1("calo_background","([1]*exp(-x/[0]))+[2]",0,300);
  //Bkg_Calo->SetParameter(0,myParameters[5]); // Seed parameters (prior to fit)
  //Bkg_Calo->SetParameter(1,myParameters[6]);    
  //Bkg_Calo->SetParameter(2,myParameters[7]);
  /*Bkg_Calo->SetParameter(0,Source_Calo->GetParameter(5));  // Resultant parameters (after fit)
  Bkg_Calo->SetParameter(1,Source_Calo->GetParameter(6));
  Bkg_Calo->SetParameter(2,Source_Calo->GetParameter(7));    
  Bkg_Calo->SetLineColor(kMagenta);
  Bkg_Calo->Draw("same");
  */

  // Save Parameters for Calo PMT's Bi Fit
  Source_Calo->GetParameters(myParameters);
  Double_t HMeV_Mean = Source_Calo->GetParameter(4);
  Double_t HMeV_MeanErr = Source_Calo->GetParError(4);
  Double_t HMeV_Sigma = sqrt(HMeV_Mean);
  Double_t HMeV_SigmaErr = (HMeV_MeanErr)/(2.0*HMeV_Sigma);
  Double_t MeV_Mean = Source_Calo->GetParameter(1);
  Double_t MeV_MeanErr = Source_Calo->GetParError(1);
  Double_t MeV_Sigma = sqrt(MeV_Mean);
  Double_t MeV_SigmaErr = (MeV_MeanErr)/(2.0*MeV_Sigma);

  Double_t Source_Calo_chi2 = Source_Calo->GetChisquare();
  Double_t Source_Calo_NDF = Source_Calo->GetNDF();

  // Save Parameters for Calo PMT's LED Fit
  Double_t CaloLED_Mean = LED_Calo->GetParameter(1);
  Double_t CaloLED_MeanErr = LED_Calo->GetParError(1);
  Double_t CaloLED_Sigma = LED_Calo->GetParameter(2);
  Double_t CaloLED_SigmaErr = LED_Calo->GetParError(2);
  Double_t CaloLED_chi2 = LED_Calo->GetChisquare();
  Double_t CaloLED_NDF = LED_Calo->GetNDF();


  // Draw full Ref PMT fits  
  c4->cd();
  Source_Ref->SetLineColor(kRed);
  LED_Ref->SetLineColor(kBlue);
  //hPedSubADC5->Draw();
  hPedSubADC5->Fit("Source_Ref", "RS");
  hPedSubADC5->Fit("LED_Ref", "RQ+");
  hPedSubADC5->GetXaxis()->SetRangeUser(5,360);
  hPedSubADC5->Draw();
  if (SavePlots){
  c4->Print("Plots/Am_LED_fits.png");
  }

  // Check Ref PMT Fits with Background
 /* TF1 *Bkg_Ref = new TF1("Bkg_Ref","([1]*exp(-x/[0]))+[2]",80, 200);
  Bkg_Ref->SetParameter(0,Source_Ref->GetParameter(5));
  Bkg_Ref->SetParameter(1,Source_Ref->GetParameter(4));
  Bkg_Ref->SetParameter(2,Source_Ref->GetParameter(6));
  Bkg_Ref->SetLineColor(kMagenta);
  Bkg_Ref->Draw("same");
 */

  // Save Parameters for Ref PMT's Am Fit 
  Double_t Am_Mean = Source_Ref->GetParameter(0);
  Double_t Am_MeanErr = Source_Ref->GetParError(0);
  Double_t Am_Sigma = Source_Ref->GetParameter(1);
  Double_t Am_SigmaErr = Source_Ref->GetParError(1);
  Double_t Am_chi2 = Source_Ref->GetChisquare();
  Double_t Am_NDF = Source_Ref->GetNDF();

  // Save Parameters for Ref PMT's LED Fit
  Double_t RefLED_Mean = LED_Ref->GetParameter(1);
  Double_t RefLED_MeanErr = LED_Ref->GetParError(1);
  Double_t RefLED_Sigma = LED_Ref->GetParameter(2);
  Double_t RefLED_SigmaErr = LED_Ref->GetParError(2);
  Double_t RefLED_chi2 = LED_Ref->GetChisquare();
  Double_t RefLED_NDF = LED_Ref->GetNDF();


  //******************************************************************//
  // Print Relevant Fit Stats and Write Parameters to Parameter Files //
  //******************************************************************//
  
  printf("Source_Calo: chi2 = %f; NDF = %f; chi2/NDF = %f \n", Source_Calo_chi2, Source_Calo_NDF, (Source_Calo_chi2/Source_Calo_NDF));
  printf("HMeV: chi2 = %f; NDF = %f; chi2/NDF = %f \n", HMeV2_chi2, HMeV2_NDF, (HMeV2_chi2/HMeV2_NDF));
  printf("MeV: chi2 = %f; NDF = %f; chi2/NDF = %f \n", MeV2_chi2, MeV2_NDF, (MeV2_chi2/MeV2_NDF));
  printf("Calo LED: chi2 = %f; NDF = %f; chi2/NDF = %f \n", CaloLED_chi2, CaloLED_NDF, (CaloLED_chi2/CaloLED_NDF));
  printf("Am: chi2 = %f; NDF = %f; chi2/NDF = %f \n", Am_chi2, Am_NDF, (Am_chi2/Am_NDF));
  printf("Ref LED: chi2 = %f; NDF = %f; chi2/NDF = %f \n", RefLED_chi2, RefLED_NDF, (RefLED_chi2/RefLED_NDF));
  TString outfile1 = "Parameters1.dat";
  TString outfile2 = "Parameters2.dat";
  TString outfile3 = "Parameters3.dat";
  cout << "Printed to File: \n" 
	"HMev:" <<HMeV_Mean<<" "<<HMeV_MeanErr<<" "<<HMeV_Sigma<<" "<<HMeV_SigmaErr<<"\n" 
	"MeV:" <<MeV_Mean<<" "<<MeV_MeanErr<<" "<<MeV_Sigma<<" "<<MeV_SigmaErr<<" "<<Source_Calo_chi2<<" "<<Source_Calo_chi2/Source_Calo_NDF<<"\n" 
	"CaloLED:" <<CaloLED_Mean<<" "<<CaloLED_MeanErr<<" "<<CaloLED_Sigma<<" "<<CaloLED_SigmaErr<<" "<<CaloLED_chi2<<" "<<CaloLED_chi2/CaloLED_NDF<<"\n" 
	"Am:" <<Am_Mean<<" "<<Am_MeanErr<<" "<<Am_Sigma<<" "<<Am_SigmaErr<<" "<<Am_chi2<<" "<<Am_chi2/Am_NDF<<"\n" 
	"RefLED:" <<RefLED_Mean<<" "<<RefLED_MeanErr<<" "<<RefLED_Sigma<<" "<<RefLED_SigmaErr<<" "<<RefLED_chi2<<" "<<RefLED_chi2/RefLED_NDF<< endl;
	ofstream out1, out2, out3;

  // Ensure that only the parameters from fits which meet the chi2/NDF threshold are written to files  
  if(Write){
    if((Source_Calo_chi2/Source_Calo_NDF)<8 && (Am_chi2/Am_NDF)<8) {
      
      // Primary output file
      out1.open(outfile1,ios_base::app);
      out1<<
	HMeV_Mean<<" "<<HMeV_MeanErr<<" "<<HMeV_Sigma<<" "<<HMeV_SigmaErr<<" "<<
	MeV_Mean<<" "<<MeV_MeanErr<<" "<<MeV_Sigma<<" "<<MeV_SigmaErr<<" "<<Source_Calo_chi2<<" "<<Source_Calo_chi2/Source_Calo_NDF<<" "<<
	CaloLED_Mean<<" "<<CaloLED_MeanErr<<" "<<CaloLED_Sigma<<" "<<CaloLED_SigmaErr<<" "<<CaloLED_chi2<<" "<<CaloLED_chi2/CaloLED_NDF<<" "<<
	Am_Mean<<" "<<Am_MeanErr<<" "<<Am_Sigma<<" "<<Am_SigmaErr<<" "<<Am_chi2<<" "<<Am_chi2/Am_NDF<<" "<<
	RefLED_Mean<<" "<<RefLED_MeanErr<<" "<<RefLED_Sigma<<" "<<RefLED_SigmaErr<<" "<<RefLED_chi2<<" "<<RefLED_chi2/RefLED_NDF<<
      	endl;
      out1.close();
      cout << "***********************" << endl;
      cout << "    GOOD FIT!!!  :)  " << endl;
      cout << "***********************" << endl;
    }

    // Notify if there is a bad fit and output a string that will be removed from the Parameters1.dat file
    else {
      if((Source_Calo_chi2/Source_Calo_NDF)>8){
      cout << "***********************" << endl;
      cout << "    BAD BISMUTH FIT   :(   " << endl;                                                 
      cout << "***********************" << endl;
      out1.open(outfile1,ios_base::app);
      out1<<"      ERROR   --    BAD BISTHMUTH FIT   --    ERASE LINE "<<endl;
      out1.close();
      }
      if((Am_chi2/Am_NDF)>8){
      cout << "***********************" << endl;
      cout << "    BAD AMERICIUM_2 FIT   :(   " << endl;
      cout << "***********************" << endl;
      out1.open(outfile1,ios_base::app);
      out1<<"      ERROR   --    BAD AMERICIUM FIT   --    ERASE LINE "<<endl;
      out1.close();
      } 
    }
  }

  // Write first alternate fit parameters to Parameters2.dat
  if(Write){
    out2.open(outfile2,ios_base::app);
    out2<<
    HMeV2_Mean<<" "<<HMeV2_MeanErr<<" "<<HMeV2_Sigma<<" "<<HMeV2_SigmaErr<<" "<<
    MeV2_Mean<<" "<<MeV2_MeanErr<<" "<<MeV2_Sigma<<" "<<MeV2_SigmaErr<<" "<<
    CaloLED_Mean<<" "<<CaloLED_MeanErr<<" "<<CaloLED_Sigma<<" "<<CaloLED_SigmaErr<<" "<<
    Am2_Mean<<" "<<Am2_MeanErr<<" "<<Am2_Sigma<<" "<<Am2_SigmaErr<<" "<<
    RefLED_Mean<<" "<<RefLED_MeanErr<<" "<<RefLED_Sigma<<" "<<RefLED_SigmaErr<<
    endl;
    out2.close();

    // Write second alternate fit parameters to Parameters3.dat
    out3.open(outfile3,ios_base::app);
    out3<<
    HMeV3_Mean<<" "<<HMeV3_MeanErr<<" "<<HMeV3_Sigma<<" "<<HMeV3_SigmaErr<<" "<<
    MeV3_Mean<<" "<<MeV3_MeanErr<<" "<<MeV3_Sigma<<" "<<MeV3_SigmaErr<<" "<<
    CaloLED_Mean<<" "<<CaloLED_MeanErr<<" "<<CaloLED_Sigma<<" "<<CaloLED_SigmaErr<<" "<<
    Am3_Mean<<" "<<Am3_MeanErr<<" "<<Am3_Sigma<<" "<<Am3_SigmaErr<<" "<<
    RefLED_Mean<<" "<<RefLED_MeanErr<<" "<<RefLED_Sigma<<" "<<RefLED_SigmaErr<<
    endl;
   out3.close();
  }
  

return(0);  
}
 
/************* Source_Calo Fit Function **********/

Double_t threegaussandledfit(Double_t *x, Double_t *par)
{   
  //*** copy of parameters for reference ***
  /*
    myParameters[0]=ADC6->GetEntries()/1e4;         normalisation
    myParameters[1]=181.0;                          mean of 1MeV peak
    myParameters[2]=1.0;                            width
    myParameters[3]=4.0*myParameters[0];            relative norm of 2nd peak
    myParameters[4]=76.0;                           mean of 0.5MeV peak
    myParameters[8]=0.1*myParameters[0];            relative norm of 3rd peak
    myParameters[9]=250.0;                          mean of 1.6MeV peak
    myParameters[5]=150000.0;                       bkg decay const
    myParameters[6]=50.0;                           bkg amp
    myParameters[7]=-300.0;                         bkg const
  */
  double convert = par[1]/976.0;
  double deltae = par[2]*sqrt(par[1]);
  double deltae2 = par[2]*sqrt(par[4]);
  double deltae3 = par[2]*sqrt(par[9]);
  
  double george = par[0]*
    (7.11*TMath::Gaus(x[0],par[1],deltae) + //centered 
     1.84*TMath::Gaus(x[0],(par[1]+convert*72),deltae) + //offset  
     0.441*TMath::Gaus(x[0],(par[1]+convert*84),deltae)) 
    +
    par[3]*
    (1.548*TMath::Gaus(x[0],par[4],deltae2) + //braching ratios 
     0.429*TMath::Gaus(x[0],(par[4]+convert*73),deltae2)+
     0.1057*TMath::Gaus(x[0],(par[4]+convert*85),deltae2))
    +
    par[8]*
    (1.52*TMath::Gaus(x[0],par[9],deltae3) +
     0.438*TMath::Gaus(x[0],(par[9]+convert*73),deltae3)+
     0.140*TMath::Gaus(x[0],(par[9]+convert*85),deltae3))
    +
    //  par[10]*TMath::Gaus(x[0],par[11],par[12])
    //  par[5]/(x[0]+par[6])+par[7];
    exp(-x[0]/par[5])*par[6]+par[7];
  return george;
}
