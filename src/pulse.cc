#include "pulse.h"
#include "TMath.h"
#include "TComplex.h"
#include "TF1.h"
#include "TGraph.h"
#include "TVirtualFFT.h"
#include <algorithm>
#include <limits>

bb::pulse::pulse(size_t nsamples, const daqint_t * data) :
        _nsamples(nsamples)
{
        _data = (real_t *)calloc(_nsamples, sizeof(real_t));
        if (data) set_data(data);
}

bb::pulse::~pulse()
{
        //if (_data) free(_data);
}

void bb::pulse::set_data(const daqint_t * data)
{
        for (size_t i = 0; i < _nsamples; ++i) {
                _data[i] = (real_t)data[i];

        }
}


//void bb::pulse::print(FILE * fd)
//{
//        for (size_t i = 0; i < _nsamples; ++i) {
//                fprintf(fd, "%u %d\n", i, _data[i]); // FIXME
//        }
//}


void bb::pulse::print_data(std::ostream & os, real_t max, real_t time)
{
        for (size_t i = 0; i < _nsamples; ++i) {
              os << _data[i] <<" "<<  "\n";  //os << i << " " << _data[i] << " " << max << "\n";
        }
}


void bb::pulse::inspect(std::ostream & os)
{
        os << "#        average: " << average(0, _nsamples) << "\n";
        os << "#  average[0:50]: " << average(0, std::min((size_t)50, _nsamples)) << "\n";
        os << "# average[0:100]: " << average(0, std::min((size_t)100, _nsamples)) << "\n";
        os << "#            rms: " << rms(0, _nsamples) << "\n";
        os << "#      rms[0:50]: " << rms(0, std::min((size_t)50, _nsamples)) << "\n";
        auto p = maximum(0, _nsamples);
        os << "#            max: " << p.second << " at " << p.first << "\n";
        p = maximum_fitted(0, _nsamples);
        os << "#     fitted max: " << p.second << " at " << p.first << "\n";
        print_data(os, p.second, p.first);
}


void bb::pulse::pre_process(size_t ped_samples)
{
        real_t ped = average(0, ped_samples);
        for (size_t i = 0; i < _nsamples; ++i) _data[i] = _data[i] - ped;
}


void bb::pulse::filter(size_t size, size_t iM, real_t M)
{
                  FILE * f_ms = fopen("pulse_sample.dat", "r");
                  FILE * f_ns = fopen("noise_sample.dat", "r");
                  FILE * f_win = fopen("tukey.dat", "r");
                  FILE * f_ham = fopen("hamming.dat", "r");
                  Double_t model[size]={0}, model_win[size]={0}, noise[size]={0}, re_wiener[size]={0}, data_fft[size*2]={0},win[size]={0}, ham[size]={0},tranf[size]={0},trans_mod[size*2]={0};
                   float number=0,freq=992.056; size_t i=0, s=0;
                TComplex c_noise[size/2]=0, k=0, transfer_func[size/2]=0, c_tf_mod[size]=0,c_wiener[size/2]=0, c_model[size/2]=0, c_model_win[size/2]=0, c_data[size]=0;
                 while( fscanf(f_win, "%f \n", &number) > 0 ) 
               {    
                win[i]=number; 
                 
                 i++;
                } i=0;number=0;
                while( fscanf(f_ham, "%f \n", &number) > 0 ) 
               {    
                ham[i]=number; 
                i++; 
                } i=0;number=0; 
                while( fscanf(f_ms, "%f \n", &number) > 0 ) // loading model and writing data for fft transform
                {    
                model[i]= number; 
                i++; 
                } i=0;number=0;
               while( fscanf(f_ns, "%f \n", &number) > 0 ) 
               {    
                noise[i]=number; 
                model[i]=model[i];
                i++; 
                } 
                fclose(f_ms);
                fclose(f_ns);
                fclose(f_win);
                fclose(f_ham);

 /*/---------------making diff/---------
                 for (s = 0; s < size-1; ++s)
                 {
                 model[i]=model[i+1]-model[i];
                 _data[i]=_data[i+1]-_data[i];  
                 } model[size]=0;_data[size]=0;
//------------------------------------*/
                for (i=0; i<size*2; ++i)              
                {
                if(i>=size/2 && i<size*2-size/2) {data_fft[i]=model[i-size/2]*win[i-size/2]*M;}// if (i>1020 && i<1023){std::cerr << "data value  "<< win[i]<<" number: "<< i <<"\n";}
                else data_fft[i]=0;
                }
             
                Int_t n_size = size, n_size2 = size*2;
                
   		TVirtualFFT *fft_model = TVirtualFFT::FFT(1, &n_size, "R2C K");
   		fft_model->SetPoints(model);
   		fft_model->Transform();
                Double_t *re_model = new Double_t[size/2];
                Double_t *im_model = new Double_t[size/2];
                fft_model->GetPointsComplex(re_model,im_model);
		for (s = 0; s < size/2; ++s) { c_model[s] = TComplex(re_model[s], im_model[s]);  c_noise[s] = TComplex(noise[s], 0); } 
   
   		TVirtualFFT *fft_data = TVirtualFFT::FFT(1, &n_size2, "R2C K");
   		fft_data->SetPoints(data_fft);
   		fft_data->Transform();
                Double_t *re_data = new Double_t[size];
                Double_t *im_data = new Double_t[size];
                fft_data->GetPointsComplex(re_data,im_data);
		for (s = 0; s < size; ++s) { c_data[s] = TComplex(re_data[s], im_data[s]);model_win[s]=model[s]*win[s]; } 


                TVirtualFFT *fft_model_win = TVirtualFFT::FFT(1, &n_size, "R2C K");
   		fft_model_win->SetPoints(model_win);
   		fft_model_win->Transform();
                Double_t *re_model_win = new Double_t[size/2];
                Double_t *im_model_win = new Double_t[size/2];
                Double_t *re_tr = new Double_t[size/2];
                Double_t *im_tr = new Double_t[size/2];
                fft_model_win->GetPointsComplex(re_model_win,im_model_win);
		for (s = 0; s < size/2; ++s) { c_model_win[s] = TComplex(re_model_win[s], im_model_win[s]); } 
      /*       // WIENER FILTER
               for (i=0; i<size/2; ++i)
                {
                c_wiener[i]=(TComplex::Abs(c_model[i])*TComplex::Abs(c_model[i]))/((TComplex::Abs(c_model[i])*TComplex::Abs(c_model[i])) + c_noise[i]);              
                     //    std::cerr << "k re value  "<< c_wiener[i]<<"\n";                          
                }
//GATTI MANFREDI FILTER */
                 k=0;
                for (i=0; i<size/2; ++i)
                {
              
                   k+=(TComplex::Abs(c_model_win[i]))*(TComplex::Abs(c_model_win[i]))/(noise[i]); ////
                   //   if (i>1023 && i<1026){std::cerr << "k re value  "<< c_noise[i].Re() <<"\n";}
                } k=size/(k*2); 
              //   std::cerr << "k re value  "<< k <<"\n";
                  for (i=0; i<size/2; ++i)
                {
                transfer_func[i]=k*(((TComplex::Conjugate(c_model[i]))/(noise[i]))*exp(-TComplex::I()*2*(TMath::Pi())*iM*i));
               // transfer_func[i]=c_wiener[i];
               //    std::cerr << "transf func re value  "<<transfer_func[i] <<"\n";
                 }                
               for (size_t s = 0; s < size/2; ++s) { 
               //if (s>0){std::cerr << "noise value  "<< transfer_func[s]<<"s number"<< s <<"\n";}
                re_tr[s]=transfer_func[s].Re();
                im_tr[s]=transfer_func[s].Im();
  
                 }  
                
               TVirtualFFT *fft_back_transf = TVirtualFFT::FFT(1, &n_size, "C2R K");
               fft_back_transf->SetPointsComplex(re_tr,im_tr);
               fft_back_transf->Transform();
               fft_back_transf->GetPoints(tranf);                          
                
                 //int n=size/2;

                for (i=0; i<size*2; ++i)
                
                {//if(i<size){ std::cerr << "data value  "<< tranf[i]<<" number: "<< i <<"\n";}
                if(i>=0 && i<size/2) {trans_mod[i]=win[i+size/2]*tranf[i]/size;}
                if(i>=size/2 && i<size*3/2) {trans_mod[i]=0;}
                if(i>=size*3/2 && i<size*2) {trans_mod[i]=win[i-size*3/2]*tranf[i-size]/size;}
                // if (i>1020 && i<1030){
               // std::cerr << "data value  "<< trans_mod[i]<<" number: "<< i <<"\n";//}
                 }
                TVirtualFFT *fft_transf = TVirtualFFT::FFT(1, &n_size2, "R2C K");
   		fft_transf->SetPoints(trans_mod);
   		fft_transf->Transform();
                Double_t *re_tm = new Double_t[size];
                Double_t *im_tm = new Double_t[size];
                fft_transf->GetPointsComplex(re_tm,im_tm);
		for (s = 0; s < size; ++s) { c_tf_mod[s] = TComplex(re_tm[s], im_tm[s]); }                 

                for (i=0; i<size; ++i)
                {c_data[i]=c_tf_mod[i]*c_data[i];//*c_data[i];
                re_data[i]=c_data[i].Re();//*c_tf_mod[i].Re()-c_data[i].Im()*c_tf_mod[i].Im();
                im_data[i]=c_data[i].Im();//*c_tf_mod[i].Im()+c_data[i].Im()*c_tf_mod[i].Re();
                
                }
               // printf("%f, \n", k);
               // getchar();
               TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &n_size2, "C2R K");
               fft_back->SetPointsComplex(re_data,im_data);
               fft_back->Transform();
               fft_back->GetPoints(data_fft);
                for (i=size/2; i<size*3/2; ++i)
                {_data[i-size/2]=data_fft[i]/size/2;}
                                              
              
     
}


bb::real_t bb::pulse::average(size_t start, size_t size)
{
        real_t m = 0;
        for (size_t i = start; i < size; ++i) {
                m += _data[i];
        }
        return (real_t)m / (real_t)size;
}


bb::real_t bb::pulse::surface(size_t start, size_t size)
{
        real_t surface=0;
        for (size_t i = start; i < size; ++i) {
       surface+=_data[i];      
        }
        return surface;
}
bb::real_t bb::pulse::surface_fitted(size_t start, size_t size, real_t M)
{
         real_t surface=0;
         TGraph * g1 = new TGraph();
        for (int i = start; i < size; ++i) g1->SetPoint(i, i-start, (real_t)_data[i]);
  
        TF1 *myfit = new TF1("myfit","([0]*exp([1]*x)+[2]*exp([3]*x)*[4])", 0, size-start);
        myfit->SetParameter(0, -1.9);
        myfit->SetParameter(1, -0.009);
        myfit->SetParameter(2, 2.77);
        myfit->SetParameter(3, -0.0016);
        myfit->SetParameter(4, M);
        g1->Fit("myfit", "Q", 0, size-start);
        TF1 * fited = g1->GetFunction("myfit");
        surface=fited->Integral(0,size-start);

        return surface;
}
float model(size_t x)
{
        return (-1.932*exp(-0.008964*x)+2.77*exp(-0.00166*x));
      
}

bb::real_t bb::pulse::shape(size_t start, size_t size)
{
        real_t s1 = 0;
        real_t s2 = 0;
        for (size_t i = start; i < size; ++i) {
                s1 += _data[i]*model(i-start);
                s2 += _data[i];
        }
        return (real_t)s1 / (real_t)s2;
}

bb::real_t bb::pulse::shape_adv(size_t start, size_t size, real_t* mode)
{
        real_t s1 = 0;
        real_t s2 = 0;
        for (size_t i = start; i < size; ++i) {
                s1 += _data[i]*mode[i];
                s2 += _data[i];
        }
        return (real_t)s1 / (real_t)s2;
}

bb::real_t bb::pulse::pulse_start(size_t start, size_t size)
{
        real_t hi = 0;
        size_t i = start;
        while (hi<10 && i<size){
           hi += _data[i]/size;
           i++;
        }
        return i;
}


bb::real_t bb::pulse::rms(size_t start, size_t size)
{
        real_t m = 0, mm = 0;
        for (size_t i = start; i < size; ++i) {
                m  += _data[i];
                mm += _data[i] * _data[i];
        }
        m /= (real_t)(size - start);
        return sqrt(mm / (real_t)(size - start) - m * m);
}

bb::real_t bb::pulse::minimum(size_t start, size_t size)
{
        real_t m = 0;
                for (size_t i = start; i < size; ++i) {
                //fprintf(stderr, "--> %d %f %f\n", i, m, _data[i]);
                if (_data[i] < m) {
                        m = _data[i];
                                        }
                //getchar();
        }
        return m;
}


std::pair<bb::real_t, bb::real_t> bb::pulse::maximum(size_t start, size_t size)
{
        real_t m = -std::numeric_limits<real_t>::max();
        //std::cerr << "--> m = " << m << "\n";
        //fprintf(stderr, "--> m = %lf\n", m);
        real_t im = 0;
        for (size_t i = start; i < size; ++i) {
                //fprintf(stderr, "--> %d %f %f\n", i, m, _data[i]);
                if (_data[i] > m) {
                        m = _data[i];
                        im = i;
                }
                //getchar();
        }
        return std::make_pair(im, m);
}


std::pair<bb::real_t, bb::real_t> bb::pulse::maximum_fitted(size_t start, size_t size)
{
        real_t m = -std::numeric_limits<real_t>::max();
        size_t im = 0;
        for (size_t i = start; i < size; ++i) {
                if (_data[i] > m) {
                        m = _data[i];
                        im = i;
                }
        }
        //char tmp[16];
        //sprintf(tmp, "pulse_%06lu", _ipulse);
        TGraph * g = new TGraph();
        //g->SetNameTitle(tmp, tmp);
        int npoints = 10;
        for (int i = 0; i < npoints; ++i) g->SetPoint(i, (real_t)im - npoints / 2 + i, (real_t)_data[im - npoints / 2 + i]);
        g->Fit("pol2", "Q");
        TF1 * f = g->GetFunction("pol2");
        real_t p0 = f->GetParameter(0);
        real_t p1 = f->GetParameter(1);
        real_t p2 = f->GetParameter(2);
        delete g;
        // return tmax, max
        return std::make_pair(-0.5 * p1 / p2, p0 - 0.25 * p1 * p1 / p2);
}


bb::real_t bb::pulse::decay_time(size_t imax, real_t fraction)
{
        real_t m = _data[imax];
        for (size_t i = imax; i < _nsamples; ++i) {
                if (_data[i] / m <= fraction) return i - imax;
        }
        return 0;
}


bb::real_t bb::pulse::rise_time(size_t imax, real_t fraction)
{
        real_t m = _data[imax];
        for (size_t i = imax; i >= 0; --i) {
                if (_data[i] / m <= fraction) return imax - i;
        }
        return 0;
}


bb::real_t bb::pulse::rise_time_interpolated(size_t imax, real_t fraction, real_t amplitude, real_t tmax)
{
        real_t m = _data[imax];
        if (amplitude) m = amplitude;
        for (size_t i = imax; i >= 0; --i) {
                if (_data[i] / m <= fraction) {
                        //return imax - i + (m * fraction - _data[i]) / (_data[i + 1] - _data[i]) * 1.;
                      // return tmax - (i  + (_data[i + 1] - m * fraction) / (_data[i + 1] - _data[i]) * 1.);
                         return tmax - i + (m * fraction - _data[i] ) / (_data[i + 1] - _data[i]);
                }
        }
        return -1;
}


bb::real_t bb::pulse::decay_time_interpolated(size_t imax, real_t fraction, real_t amplitude, real_t tmax)
{
        real_t m = _data[imax];
        if (amplitude) m = amplitude;
        for (size_t i = imax; i < _nsamples; ++i) {
                if (_data[i] / m <= fraction) {
                        return i - (m * fraction - _data[i]) / (_data[i - 1] - _data[i]) - tmax;
                }
        }
        return -1;
}
