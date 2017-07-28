
#include "pulse.h"

#include <cmath>
#include <fstream>
#include <iostream>

#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TH1F.h"
#include "TLeaf.h"
#include "TProfile.h"
#include "TTree.h"
#include "TVirtualFFT.h"

#define TIME_WINDOW 2048

typedef struct {
        int32_t detid;
        Int_t  nsamples;
        //UInt_t  nsamples;
        bb::daqint_t * data;
        Double_t ts;
} Event;

Event  _e;


void set_branches(TTree * t)
{
        t->SetBranchAddress("detid", &_e.detid);
        // handle variable size array
        // allocate maximum dimension (cf. TTreePlayer::MakeClass code)
        TLeaf * l = t->GetLeaf("nsamples", "nsamples");
        _e.data = (bb::daqint_t *)calloc(l->GetMaximum(), sizeof(bb::daqint_t));//
        t->SetBranchAddress("raw_pulse", _e.data);
        t->SetBranchAddress("nsamples", &_e.nsamples);
        t->SetBranchAddress("detid", &_e.detid);
        t->SetBranchAddress("time", &_e.ts);
        //return br;
/*         t->SetBranchAddress("raw_signal", _e.data);
        t->SetBranchAddress("nsamples", &_e.nsamples);
        t->SetBranchAddress("detid", &_e.detid);
       // t->SetBranchAddress("time", &_e.ts); */
}


//float signal_index (float t);

int main()
{
        TFile * fout = TFile::Open("histos.root", "recreate");
        TProfile * p_average_signal= new TProfile("p_average_signal","p_average_signal", TIME_WINDOW, 0., TIME_WINDOW);
        TProfile * p_average_noise= new TProfile("p_average_noise","p_average_noise", TIME_WINDOW, 0., TIME_WINDOW);
        TH1F     * h_ampl_raw               = new TH1F("h_ampl_raw", "h_ampl_raw", 1000, -5, 5);
        TH1F     * h_ampl                   = new TH1F("h_ampl", "h_ampl", 500, 0., 1000);
        TH1F     * h_ampl_light             = (TH1F*)h_ampl_raw->Clone("h_ampl_light");
        TH1F     * h_dif                   = (TH1F*)h_ampl_raw->Clone("h_dif");
        TH1F     * h_ampl_heater            = (TH1F*)h_ampl_raw->Clone("h_ampl_heater");
        TH1F     * h_ampl_keV               = new TH1F("h_ampl_keV", "h_ampl_keV", 1500, 0., 6000.);
        TH1F     * h_ampl_res               = new TH1F("h_ampl_res", "h_ampl_res", 1000, -10., 10.);
        TH1F     * h_ped                    = new TH1F("h_ped", "h_ped", 2<<13, 0., 2<<17);
        TH1F     * h_ped_raw                = new TH1F("h_ped_raw", "h_ped_raw", 1000, -1e5, 1e5);
        TH1F     * h_ped_rms                = new TH1F("h_ped_rms", "h_ped_rms", 100, 0., 100.);
        TH1F     * h_dt                     = new TH1F("h_dt", "h_dt", 5000, 0., 5000000.);
        //TH1      * h_fft                     = 0;

        TH1F     * h_time_s1                = new TH1F("h_time_s1", "h_time_s1", 20000, 0., 50.);
        TGraph   * g_decay_vs_rise          = new TGraph();
        g_decay_vs_rise->SetNameTitle("g_decay_vs_rise", "g_decay_vs_rise");
        gDirectory->Add(g_decay_vs_rise);

        TGraph * g_baseline_vs_ampl            = new TGraph();
        g_baseline_vs_ampl->SetNameTitle("g_baseline_vs_ampl", "g_baseline_vs_ampl");
        gDirectory->Add(g_baseline_vs_ampl);

        TGraph * g_baseline_vs_ampl_stab            = new TGraph();
        g_baseline_vs_ampl_stab->SetNameTitle("g_baseline_vs_ampl_stab", "g_baseline_vs_ampl_stab");
        gDirectory->Add(g_baseline_vs_ampl_stab);

         TGraph * g_shape_ampl            = new TGraph();
        g_shape_ampl->SetNameTitle("g_shape_ampl", "g_shape_ampl");
        gDirectory->Add(g_shape_ampl);
         TGraph * g_scatter            = new TGraph();
        g_scatter->SetNameTitle("g_scatter", "g_scatter");
        gDirectory->Add(g_scatter);

        TGraph * g_shape_ampl_int            = new TGraph();
        g_shape_ampl_int->SetNameTitle("g_shape_ampl_int", "g_shape_ampl_int");
        gDirectory->Add(g_shape_ampl_int);

        TGraph * g_ampl_vs_rise            = new TGraph();
        g_ampl_vs_rise->SetNameTitle("g_ampl_vs_rise", "g_ampl_vs_rise");
        gDirectory->Add(g_ampl_vs_rise);


        TGraph * g_ampl_vs_decay = new TGraph();
        g_ampl_vs_decay->SetNameTitle("g_ampl_vs_decay", "g_ampl_vs_decay");
        gDirectory->Add(g_ampl_vs_decay);


        TGraph * g_ampl_vs_rise_light            = new TGraph();
        g_ampl_vs_rise_light->SetNameTitle("g_ampl_vs_rise_light", "g_ampl_vs_rise_light");
        gDirectory->Add(g_ampl_vs_rise_light);

        TGraph * g_ampl_vs_decay_light = new TGraph();
        g_ampl_vs_decay_light->SetNameTitle("g_ampl_vs_decay_light", "g_ampl_vs_decay_light");
        gDirectory->Add(g_ampl_vs_decay_light);

        /* TGraph * g_ampl_all_vs_t        = new TGraph();
        g_ampl_all_vs_t->SetNameTitle("g_ampl_all_vs_t", "g_ampl_all_vs_t");
        gDirectory->Add(g_ampl_all_vs_t);

       */ // creation of an output is done

        TFile * fin = TFile::Open("lsm_with_light.root");
        TTree * t = (TTree*)fin->Get("ntu");

        fout->cd();
        set_branches(t);
        Long64_t nentries = t->GetEntries(); //get amount of samples
        Long64_t gcnt = 0, gcntl=0, gcnt1=0, gcnt2 = 0, gcnta = 0, gcnt_rd = 0, sc=0;
        UInt_t otmaxdaq = 1; //what is it?
        FILE * fp = fopen("pulses.dat", "w+");
        FILE * f_t = fopen("pulse_sample.dat", "r");
        FILE * fn = fopen("noise.dat", "w+");
        size_t ipulse = 0, totpulse = 1, i=0, k=1;
        std::ofstream ofs;
        ofs.open ("pappa.dat", std::ofstream::out);
        float model[TIME_WINDOW];
        float number;
        size_t ind_mean=0, count=0, is=0; Int_t size=2048;

        TVirtualFFT *fft_noise = TVirtualFFT::FFT(1, &size, "R2C K");
                Double_t for_ttf[size];
   		Double_t *re_noise = new Double_t[size];
                Double_t *im_noise = new Double_t[size];
                Double_t *re_noise_temp = new Double_t[size];
                Double_t *im_noise_temp = new Double_t[size];
                for (size_t is =0; is < size; ++is)
                 {re_noise[is]=0;
                  im_noise[is]=0;
                  re_noise_temp[is]=0;
                  im_noise_temp[is]=0; 
		for_ttf[is]=0;}
                 while( fscanf(f_t, "%f \n", &number) > 0 ) // parse %d followed by ','
                {
                model[i]= number; // instead of sum you could put your numbers in an array
                i++;
                }
                fclose(f_t);
        for (Long64_t ien = 0; ien < nentries; ++ien)
                {

                t->GetEntry(ien);
                 bb::pulse p(_e.nsamples);
                p.set_data(_e.data);
                
                 float M, fM;
                // select only one channel
                if ( _e.detid != 5 && _e.detid != 1005)
                {  
                           
               continue;}
                
                // in pulse.cc - remove the pedestal
                if ( ien % 2==0)
                {
                // signal analysis




               //=======================COLLECTING NOISE FROM PEDESTALS
		p.pre_process(500);
                 float * data = p.data();
                
                  
                 auto res = p.maximum(0, 512);
                size_t iM = res.first;                        // index of the max_ampl from the beginning of the window
                M = res.second;  
                float Min=p.minimum(0,512);
               // size_t is;                 
                if(count<1000 && ien % 2 == 0 && M<20 && Min>-20){
                   
                if (k>4){k=1;} //if (count<20){ std::cerr << "noise value at 256 "<< k << " "<< 1+size*(k-1)/4 <<"\n";}
                for(is=(1+size*(k-1)/4); is<=(size*k/4); ++is)
                 { 

                 //if (count<20){ for (size_t i=0; i<1; ++i) {std::cerr << "is value  "<< is <<"\n";}}

                  for_ttf[is]=data[is-512*(k-1)];
                //  im_noise[is]+=im_noise_temp[is-512*(k-1)];                
                 }
		++count;
                 ++k;
  
              if (k>=4){
                fft_noise->SetPoints(for_ttf);
   		fft_noise->Transform();
                fft_noise->GetPointsComplex(re_noise_temp,im_noise_temp);
        
                   for (size_t is =0; is < size; ++is)
                 {  p_average_noise->Fill(is, for_ttf[is]);
                             }
                for (size_t is =0; is < size; ++is)
                 {
                  re_noise_temp[i]=+re_noise_temp[i];
                  im_noise_temp[i]=+im_noise_temp[i];
               }}
		
}
		//=======================END OF COLLECTION


                
                float ped_raw = p.average(0,100);
                h_ped_raw->Fill(ped_raw);
                p.pre_process(500);
                float ped = p.average(0, 100); // simple average
                float ped_rms = p.rms(0, 100);
                h_ped->Fill(ped);
                h_ped_rms->Fill(ped_rms);
                res = p.maximum(100, p.n_samples());
                iM = res.first;                        // index of the max_ampl from the beginning of the window
                M = res.second;
                res = p.maximum_fitted(100, p.n_samples()); // maximum defined by fit
                float fiM = res.first; // index again?
                fM = res.second;   // Maximum position.q
                float start_index=p.pulse_start(0, p.n_samples());
                float P_surf=p.surface(1024, 1800);
                h_ampl_raw->Fill(M);
                float SI=p.shape(1024,p.n_samples());
                float ft_daq = _e.ts + fiM; // in seconds
                //trying to use pulses with good SI for new mean pulse:
                 data = p.data();
                 float * mean;
                         // pause until Enter key pressed
                    //   std::cerr << "ien" << ien << " amplitude " << fM << " dumped, press [Enter] to continue...\n";
                        //getchar(); // uncomment if you want to pause
                
//-------------------AVERAGE PULSE COLLECTING USING CHI-SQUARE-------------------
                float chi_square=0;
                   for (size_t is =1; is < p.n_samples()-1; ++is)
                { 
                chi_square+=(data[is]/M-model[is])*(data[is]/M-model[is])/(model[is]);
                //if (totpulse<2) std::cerr << "chi square: "<< chi_square <<"\n";
                }
              // std::cerr << "chi square: "<< chi_square <<"\n";
                if(chi_square<1.2 && chi_square>0.8) {
                for (size_t is =0; is < p.n_samples(); ++is)
                {
                p_average_signal->Fill(is, data[is]/M);
                }
                ind_mean++;
                } 
//-------------------END OF COLLECTION-------------------------------------------


/* //--------------------old stuff
                if (SI<1.405 && SI>1.395 && fM>5000 && fM<8000 )
                {
                 // p.pre_process(1000);
                for (size_t is =0; is < p.n_samples(); ++is)
                {
                 //mean[is]=data[is]/fM;
                p_average_signal->Fill(is, data[is]/M);
                }
                ind_mean++;
                }
*/


                 
             //---------------------------------------
                 /* if( M<30){
                for (size_t is =0; is < 2048; ++is)
                 {
                 p_average_noise->Fill(is, data[is]);
                 } ++count;
                 }
               */
                 //---------------Shape parameter using mean walue---------------
                   float s1 = 0, s2 = 0;
                  for (size_t in = 1024; in < p.n_samples(); ++in)
                  {
                      s1 +=data[in]*model[in];
                      s2 +=data[in];
                  }
                  float SI_adv=s1/s2;
                //-------------------------------------------------------------



                //---------------Shape parameter using integral---------------

                 float s_l=0, s_h=0, m_l=0,m_h=0;
                  for (size_t in = 1024; in < p.n_samples(); ++in)
                  {
                      s_l+=data[in]/fM;
                      m_l+=model[in];
                  }
                  for (size_t in = (1024+(p.n_samples()-1024)/4); in < p.n_samples(); ++in)
                  {
                    s_h+=data[in]/fM;
                    m_h+=model[in];
                  }
                  float SI_int=s_h*s_l/(m_l*m_h);
               //------------------------------------------------------------
                h_ampl_res->Fill((fM - M) / fM);
                //std::cerr << "--> " << iM << " " << M << " " << fiM << " " << fM << "\n";
               // float ft_daq = _e.ts + fiM; // in seconds

                float trise = p.rise_time_interpolated(iM, 0.05, fM, fiM) - p.rise_time_interpolated(iM, 0.95, fM, fiM);

                 //trise = p.rise_time(iM, 0.10) - p.rise_time(iM, 0.90);

                float tdecay = p.decay_time_interpolated(iM, 0.30, fM, fiM) - p.decay_time_interpolated(iM, 0.90, fM, fiM);

                // tdecay = p.decay_time(iM, 0.30) - p.decay_time(iM, 0.90);
                if(trise<100 && tdecay<300 && trise>0 && tdecay>0 )
                g_decay_vs_rise->SetPoint(gcnt_rd++, trise, tdecay);
                
              /*  fprintf(stdout, "%f %f %lu %u %u %lu %lu %f %f %f %f %f %f\n",
                        _e.ampl, M, iM, _e.t_max, _e.t_max_daq,
                        signal_rise_time(fata, iM, 0.05), signal_decay_time(fata, iM, 0.20),
                        trise,
                        tdecay,
                        fM, fiM, ped, ped_rms);*/
                
                 if (fM > 0 && fM < 12000 && SI_int>-5 && SI_int<5 )
                       g_shape_ampl_int->SetPoint(gcnt2++, fM, SI_adv);

               if (otmaxdaq && fM > 9500 && fM < 10500)
                {

                        h_ampl_heater->Fill(fM);
                        otmaxdaq = _e.ts;
                }
                 if (fM > 0 && fM < 12000 && SI>-4 && SI<4 )
                      g_shape_ampl->SetPoint(gcnt2++, fM, SI);
        	//if (fM > 9000 && fM < 11000) {
                        h_dif->Fill(ft_daq);
                        h_dt->Fill(_e.ts - otmaxdaq);
                //}
               //  g_ampl_all_vs_t->SetPoint(gcnta++, ft_daq / 3600., fM);
                if (fM > 1e5) h_time_s1->Fill(trise - 0.01 * fM / 1000.);
                //else if (_e.ampl > 101.2e+3 && _e.ampl < 101.4e+3) g_ampl_vs_decay->SetPoint(gcnt2++, ft_daq / 2000. / 3600., fM);
                if ( fM<12000 && fM>0 ){//&& ped_raw<-10000 && ped_raw>-30000){
                g_baseline_vs_ampl->SetPoint(gcnt1++, ped_raw, fM);}
                
               
               //-----------work with pulse area data-----------
               /* if (fM<12000 && fM>0 && SI<1.5 && SI>1.3 && P_surf>0 && P_surf<2700000){
                g_ampl_vs_decay->SetPoint(gcnta++, P_surf/100, fM);
                g_decay_vs_rise->SetPoint(gcnt_rd++, ped_raw, P_surf/1000);}
                  
                //======================
                if (fM<12000 && fM>0 && SI<1.5 && SI>1.3) 
                {
                float P_surf_f=p.surface_fitted(1024, 1800, fM);
                h_ampl->Fill(P_surf_f);
                }

               */

               //------------------------------------
               //STABILISATION WITH HEATER
               //------------------------------------
               float Ampl_ref=10000;
                float p0=9138.3, p1=-0.0452516;
                float Ampl_stab=fM*Ampl_ref / (p0+p1*ped_raw);
                 if (fM<12000 && fM>0 && ped_raw<-10000 && ped_raw>-30000){
                g_baseline_vs_ampl_stab->SetPoint(gcnt1++, ped_raw, Ampl_stab);}
               //---------------------------------------
               if (Ampl_stab<12000 &&SI<1.8 &&SI>1.1) h_ampl->Fill(Ampl_stab);
                float adc2keV = 4783. / 10468.;
		if(Ampl_stab<12000 ){
               
                h_ampl_keV->Fill(Ampl_stab * adc2keV);}

               // printf("%f, %f, %f, \n", trise, fM, fiM);
                //getchar();
               if (fM > 0 && fM < 12000 && trise>0 ){//&& tdecay<134 && tdecay>122){ 
g_ampl_vs_rise->SetPoint(gcnt++, fM, trise);}

              if (fM > 0 && fM < 12000) g_ampl_vs_decay->SetPoint(gcnta++, fM, tdecay);



//--------------------COMMENTED: write to file the pulse after fillter

       //  detailed check of pulses if conditions applies: after filtering
                if( ipulse<7 && SI<1.4 && SI>1.36 && fM>10200 && fM<10700 && trise>32 && trise<33){
                        p.filter(p.n_samples(),iM,M);              //wiener filter
                        ofs << "# pulse number:" << ipulse << "\n";
                        p.inspect(ofs);
                        ofs << "\n\n";
                        ofs.flush();
                        // pause until Enter key pressed
                       //std::cerr << "pulse " << ipulse << " of " << totpulse << " dumped, press [Enter] to continue...\n";
                        //getchar(); // uncomment if you want to pause
                         ++ipulse; 

                        continue; // do not perform the rest of the analysis
                              }


              }
        else
        {
          //----------Light signals----------------
                p.pre_process(500);
                auto res = p.maximum(100, p.n_samples());
                size_t iM_l = res.first;                        // index of the max_ampl from the beginning of the window
                float M_l = res.second;
                res = p.maximum_fitted(100, p.n_samples()); // maximum defined by fit
               float fM_l = res.second;
                float fiM_l = res.first; // index again?
                 float trise_light = p.rise_time_interpolated(iM_l, 0.10, fM_l, fiM_l) - p.rise_time_interpolated(iM_l, 0.90, fM_l, fiM_l);
                 //trise = p.rise_time(iM, 0.10) - p.rise_time(iM, 0.90);
                float tdecay_light = p.decay_time_interpolated(iM_l, 0.30, fM_l, fiM_l) - p.decay_time_interpolated(iM_l, 0.90, fM_l, fiM_l);
                if (tdecay_light<200 && fM_l<5000 && fM_l>0) g_ampl_vs_decay_light->SetPoint(gcntl++, fM_l, tdecay_light);

                h_ampl_light->Fill(fM_l);
              if (fM<12000 &&fM>0 && fM_l>0 && fM_l<1200)  g_scatter->SetPoint(sc++, fM, fM_l);
                  ++totpulse;

                }

        }


  //--------------Write mean pulse and noise in the file-----------
        double  model1;
        Double_t data_fft[size];
          for (size_t o = 0; o < TIME_WINDOW; ++o)
                {
         model1=p_average_signal->GetBinContent(o);
        fprintf(fp,"%f \n", model1); }
        // re_noise[o]=re_noise[o]/count;
         //im_noise[o]=im_noise[o]/count;}

          TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &size, "C2R K");
               fft_back->SetPointsComplex(re_noise_temp,im_noise_temp);
               fft_back->Transform();
               fft_back->GetPoints(data_fft); 
          for (size_t o = 0; o < size; ++o)
                {
         model1=data_fft[o]/count;
        fprintf(fn,"%f \n", model1); }

          std::cerr << "totpulse "<< totpulse <<"\n";

        ofs.close();
        if (_e.data) free(_e.data);

        // ---------------GAUS fit with ped-----------------------
 std::cerr << "mean collection: "<< ind_mean <<"\n";
 std::cerr << "noise collection: "<< count <<"\n";
        fout->Write();
        fout->Close();
        fclose(fp);
        fclose(fn);
        return 0;
}
