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



bb::real_t bb::pulse::average(size_t start, size_t size)
{
        real_t m = 0;
        for (size_t i = start; i < size; ++i) {
                m += _data[i];
        }
        return (real_t)m / (real_t)size;
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



