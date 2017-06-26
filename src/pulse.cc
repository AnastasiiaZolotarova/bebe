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


void bb::pulse::filter(size_t size, size_t iM)
{
                  FILE * f_ms = fopen("pulse_sample.dat", "r");
                  FILE * f_ns = fopen("noise_sample.dat", "r");
                  Double_t model[size], noise[size], re_wiener[size], data_fft[size];
                   float number,freq=992.056, trg=1024; size_t i=0, TIME_WINDOW=2048;
                TComplex k=0;
                TComplex transfer_func[size], c_model[size], c_noise[size], c_data[size];
                while( fscanf(f_ms, "%f \n", &number) > 0 ) // loading model and writing data for fft transform
                {    
                model[i]= number; data_fft[i]=_data[i];
                i++; 
                } i=0;
               while( fscanf(f_ns, "%f \n", &number) > 0 ) 
               {    
                noise[i]=number; 
                i++; 
                }
                fclose(f_ms);
                fclose(f_ns);
                Int_t n_size = size;
                
   		TVirtualFFT *fft_model = TVirtualFFT::FFT(1, &n_size, "R2C K");
   		fft_model->SetPoints(model);
   		fft_model->Transform();
                Double_t *re_model = new Double_t[size];
                Double_t *im_model = new Double_t[size];
                fft_model->GetPointsComplex(re_model,im_model);
		for (size_t s = 0; s < size; ++s) { c_model[s] = TComplex(re_model[s], im_model[s]); } 
      
                  
   		TVirtualFFT *fft_noise = TVirtualFFT::FFT(1, &n_size, "R2C K");
   		fft_noise->SetPoints(noise);
   		fft_noise->Transform();
                Double_t *re_noise = new Double_t[size];
                Double_t *im_noise = new Double_t[size];
                fft_noise->GetPointsComplex(re_noise,im_noise);
		for (size_t s = 0; s < size; ++s) { c_noise[s] = TComplex(re_noise[s], im_noise[s]); } 
   
   		TVirtualFFT *fft_data = TVirtualFFT::FFT(1, &n_size, "R2C K");
   		fft_data->SetPoints(data_fft);
   		fft_data->Transform();
                Double_t *re_data = new Double_t[size];
                Double_t *im_data = new Double_t[size];
                fft_data->GetPointsComplex(re_data,im_data);
		for (size_t s = 0; s < size; ++s) { c_data[s] = TComplex(re_data[s], im_data[s]); } 

               // WIENER FILTER
               for (i=0; i<size; ++i)
                {
                re_wiener[i]=(TMath::Abs(re_model[i]*re_model[i]+im_model[i]*im_model[i]))/(TMath::Abs(re_model[i]*re_model[i]+im_model[i]*im_model[i])+TMath::Abs(re_noise[i]*re_noise[i]+im_noise[i]*im_noise[i]));              
              //  re_data[i]=re_wiener[i]*re_data[i];
               // im_data[i]=re_wiener[i]*im_data[i];
                                             
                }
//GATTI MANFREDI FILTER 
                for (i=0; i<size; ++i)
                {
              //   k+=TMath::Abs(re_model[i]*re_model[i]+im_model[i]*im_model[i])/(re_noise[i]+im_noise[i]);
                      k+=TComplex::Abs(c_model[i])/(c_noise[i]);
                }
                  for (i=0; i<size; ++i)
                {
                //transfer_func[i]=(re_model[i]-TComplex::I()*im_model[i])/(re_noise[i]+TComplex::I()*im_noise[i])*exp(-TComplex::I()*(TMath::Pi())*iM*freq);

                transfer_func[i]=k*((TComplex::Conjugate(c_model[i]))/(c_noise[i])*exp(-TComplex::I()*(TMath::Pi())*iM*freq));
                c_data[i]=transfer_func[i]*c_data[i];
                //re_data[i]=transfer_func[i]*re_data[i];
                //im_data[i]=transfer_func[i]*im_data[i];

               for (size_t s = 0; s < size; ++s) { 
               
                re_data[s]=c_data[s].Re();
                im_data[s]=c_data[s].Im();
  
                 }  
                //transfer_func                             
                }
               // printf("%f, \n", k);
               // getchar();
               TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &n_size, "C2R K");
               fft_back->SetPointsComplex(re_data,im_data);
               fft_back->Transform();
               fft_back->GetPoints(data_fft);
                for (i=0; i<size; ++i)
                {_data[i]=data_fft[i]/size;}
                                              
              
     
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
