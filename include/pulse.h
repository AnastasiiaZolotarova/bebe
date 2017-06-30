#ifndef PULSE_H
#define PULSE_H
/*
 *        Class: pulse
 *  Description: handling of sample and basic properties
 *       Author: Federico Ferri, CEA/Saclay
 */
#include "types.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>

namespace bb {

        class pulse
        {
                public:
                        pulse(size_t nsamples, const daqint_t * data = 0);

                        ~pulse();

                        real_t * data() { return _data; }
                        //void print(FILE * fd);
                        void print_data(std::ostream & os, real_t max, real_t time);
                        void inspect(std::ostream & os);

                        void set_data(const daqint_t * data);
                        //void copyData(daqint_t * data) { }

                        void pre_process(size_t ped_samples);
                    
                        size_t n_samples() { return _nsamples; };

                        // pulse analysis
                        real_t average(size_t start, size_t size);
                        real_t rms(size_t start, size_t size);
                        real_t minimum(size_t start, size_t size);
                        std::pair<real_t, real_t> maximum(size_t start, size_t size);
                        std::pair<real_t, real_t> maximum_fitted(size_t start, size_t size);


                private:
                        size_t _nsamples;
                        real_t * _data;
                        
        };

}

#endif
