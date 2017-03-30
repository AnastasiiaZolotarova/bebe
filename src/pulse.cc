#include "pulse.h"

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


void bb::pulse::print_data(std::ostream & os, real_t max)
{
        for (size_t i = 0; i < _nsamples; ++i) {
              os << _data[i] << " "  "\n";  //os << i << " " << _data[i] << " " << max << "\n";
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
        print_data(os, p.second);
}


void bb::pulse::pre_process(size_t ped_samples)
{
        real_t ped = average(0, ped_samples);
        for (size_t i = 0; i < _nsamples; ++i) _data[i] = _data[i] - ped;
}


void bb::pulse::filter(size_t size)
{
                  FILE * f_ms = fopen("pulse_sample.dat", "r");
                  FILE * f_ns = fopen("noise_sample.dat", "r");
                  Double_t model[size], noise[size], re_wiener[size], data_fft[size];
                   float number; size_t i=0;
                while( fscanf(f_ms, "%f \n", &number) > 0 ) // loading model and writing data for fft transform
                {    
                model[i]= number; data_fft[i]=_data[i];
                i++; 
                } i=0;
               while( fscanf(f_ns, "%f \n", &number) > 0 ) 
               {    
                noise[i]= number; 
                i++; 
                }
                fclose(f_ms);
                fclose(f_ns);
                Int_t n_size = size+1;
                
   		TVirtualFFT *fft_model = TVirtualFFT::FFT(1, &n_size, "R2C K");
   		fft_model->SetPoints(model);
   		fft_model->Transform();
                Double_t *re_model = new Double_t[size];
                Double_t *im_model = new Double_t[size];
                fft_model->GetPointsComplex(re_model,im_model); 
      
                  
   		TVirtualFFT *fft_noise = TVirtualFFT::FFT(1, &n_size, "R2C K");
   		fft_noise->SetPoints(noise);
   		fft_noise->Transform();
                Double_t *re_noise = new Double_t[size];
                Double_t *im_noise = new Double_t[size];
                fft_noise->GetPointsComplex(re_noise,im_noise);
   
   		TVirtualFFT *fft_data = TVirtualFFT::FFT(1, &n_size, "R2C K");
   		fft_data->SetPoints(data_fft);
   		fft_data->Transform();
                Double_t *re_data = new Double_t[size];
                Double_t *im_data = new Double_t[size];
                fft_data->GetPointsComplex(re_data,im_data);
                for (i=0; i<size; ++i)
                {
                re_wiener[i]=(TMath::Abs(re_model[i]*re_model[i]+im_model[i]*im_model[i]))/(TMath::Abs(re_model[i]*re_model[i]+im_model[i]*im_model[i])+TMath::Abs(re_noise[i]*re_noise[i]+im_noise[i]*im_noise[i]));              
                re_data[i]=re_wiener[i]*re_data[i];
               
                               
                }
               TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &n_size, "C2R K");
               fft_back->SetPointsComplex(re_data,im_data);
               fft_back->Transform();
               fft_back->GetPoints(data_fft);
                for (i=0; i<size; ++i)
                {_data[i]=data_fft[i];}
                                              
              
     
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
                        return tmax - (i + 1 + (_data[i + 1] - m * fraction) / (_data[i + 1] - _data[i]) * 1.);
                }
        }
        return -1;
}


bb::real_t bb::pulse::decay_time_interpolated(size_t imax, real_t fraction, real_t amplitude, real_t tmax)
{
        real_t m = _data[imax];
        if (amplitude) m = amplitude;
        for (size_t i = imax; i < _nsamples; ++i) {
                if (_data[i] / m < fraction) {
                        return i - (m * fraction - _data[i]) / (_data[i - 1] - _data[i]) - tmax;
                }
        }
        return -1;
}
