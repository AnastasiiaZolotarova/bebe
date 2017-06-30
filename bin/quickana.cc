
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

}

//float signal_index (float t);

int main()
{
        TFile * fout = TFile::Open("histos.root", "recreate");
        TProfile * p_average_noise= new TProfile("p_average_noise","p_average_noise", TIME_WINDOW, 0., TIME_WINDOW);
        TH1F     * h_ampl_heat               = new TH1F("h_ampl_heat", "h_ampl_heat", 1000, -5, 5);
        TH1F     * h_ampl_light             = (TH1F*)h_ampl_heat->Clone("h_ampl_light");
        TH1F     * h_ampl_res               = new TH1F("h_ampl_res", "h_ampl_res", 1000, -10., 10.);
        TH1F     * h_ped_rms                = new TH1F("h_ped_rms", "h_ped_rms", 100, 0., 100.);
 
        TGraph * g_baseline_vs_ampl            = new TGraph();
        g_baseline_vs_ampl->SetNameTitle("g_baseline_vs_ampl", "g_baseline_vs_ampl");
        gDirectory->Add(g_baseline_vs_ampl);

          TGraph * g_baseline_vs_time            = new TGraph();
        g_baseline_vs_ampl->SetNameTitle("g_baseline_vs_time", "g_baseline_vs_time");
        gDirectory->Add(g_baseline_vs_time);

         TGraph * g_scatter            = new TGraph();
        g_scatter->SetNameTitle("g_scatter", "g_scatter");
        gDirectory->Add(g_scatter);

       // creation of an output is done

        TFile * fin = TFile::Open("lsm_with_light.root");
        TTree * t = (TTree*)fin->Get("ntu");

        fout->cd();
        set_branches(t);
        Long64_t nentries = t->GetEntries(); //get amount of samples
        UInt_t otmaxdaq = 1; //what is it?
        FILE * fn = fopen("noise.dat", "w+");
        size_t ipulse = 0, totpulse = 1, i=0, k=1;
        Long64_t sc=0, gcnt1=0;
        std::ofstream ofs;
        ofs.open ("pappa.dat", std::ofstream::out);     
        size_t ind_mean=0, count=0, is=0; Int_t size=2048;
        float M, fM;
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
                
        for (Long64_t ien = 0; ien < nentries; ++ien)
                {

                t->GetEntry(ien);
                 bb::pulse p(_e.nsamples);
                p.set_data(_e.data);
                
                 
                // select only one channel
                if ( _e.detid != 5 && _e.detid != 1005) continue;
                
                // in pulse.cc - remove the pedestal
                if ( ien % 2==0)
                {
                // signal analysis

               //=======================COLLECTING NOISE FROM PEDESTALS
                float ped_raw = p.average(0,200);
                p.pre_process(500);
                 float * data = p.data();
                
                  
                 auto res = p.maximum(0, 512);
                size_t iM = res.first;                        // index of the max_ampl from the beginning of the window
                M = res.second;  
                float Min=p.minimum(0,512);
               // size_t is;                 
                if(count<1000 && ien % 2 == 0 && M<20 && Min>-20)
           {                   
                if (k>4)k=1; 
                for(is=(1+size*(k-1)/4); is<=(size*k/4); ++is)
                 {for_ttf[is]=data[is-512*(k-1)];   }         
                 
		++count;
                 ++k;
  
              if (k>=4)
                {
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
                 }
		}
           }
		//=======================END OF COLLECTION

                float ped = p.average(0, 100); // simple average
                float ped_rms = p.rms(0, 100);
                h_ped_rms->Fill(ped_rms);
                res = p.maximum(100, p.n_samples());
                iM = res.first;                        // index of the max_ampl from the beginning of the window
                M = res.second;
                res = p.maximum_fitted(100, p.n_samples()); // maximum defined by fit
                float fiM = res.first; // index again?
                fM = res.second;   // Maximum position.q
                h_ampl_heat->Fill(fM);
                 float ft_daq = _e.ts + fiM; // in seconds
                
               //------------------------------------------------------------
                h_ampl_res->Fill((fM - M) / fM);
              //  h_dt->Fill(_e.ts - otmaxdaq);
             
                g_baseline_vs_ampl->SetPoint(gcnt1++, ped_raw, fM);
                
                 g_baseline_vs_time->SetPoint(gcnt1++, _e.ts - otmaxdaq, fM);
                
                         ++ipulse; 

                        //continue; // do not perform the rest of the analysis
                              


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
                
                h_ampl_light->Fill(fM_l);
              if (fM>0 && fM_l>0 && fM<20000 && fM_l<20000) { g_scatter->SetPoint(sc++, fM, fM_l);}
                  ++totpulse;

        }
    
        }
          std::cerr << "totpulse "<< totpulse <<"\n";

        ofs.close();
        if (_e.data) free(_e.data);      
 std::cerr << "noise collection: "<< count <<"\n";
        fout->Write();
        fout->Close();
        fclose(fn);
        return 0;

}
