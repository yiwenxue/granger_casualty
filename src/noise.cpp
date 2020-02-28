#include <noise.h>
#include <string>
#include <iostream>

double*
pinkNoise(double exponent,
          int size)
{
    double *normal = (double *)malloc(sizeof(double)*size);
    std::default_random_engine generator;
    generator.seed(time(NULL));
    std::normal_distribution<double> distribution(0.0,10.0);
    for (int i = 0; i < size; ++i) {
        normal[i] = distribution(generator);
    }
    int complex_size = ((int)(size*1.0+0.5))/2+1;
    double freq = 0; 

    fftw_complex * middle = fftw_alloc_complex(complex_size);
    fftw_plan plan = fftw_plan_dft_r2c_1d (size, normal, middle, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    for (int i = 0; i < complex_size; ++i) {
        freq = (i+1);
        freq = pow(freq, exponent/2.0);
        /* std::cout << "freq["<< i <<"] = " << freq << std::endl; */
        middle[i][0] /= freq;
        middle[i][1] /= freq;
    }
    plan = fftw_plan_dft_c2r_1d(size, middle, normal, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    fftw_free(middle);

    freq = gsl_stats_mean(normal, 1, size);
    for (int i = 0; i < size; ++i) {
        normal[i] /= freq;
    }
    return normal;
}

int
GenpinkNoise(double *noise, 
             double exponent,
             int size)
{
    std::default_random_engine generator;
    generator.seed(clock());
    std::normal_distribution<double> distribution(0.0,10.0);
    for (int i = 0; i < size; ++i) {
        noise[i] = distribution(generator);
    }
    int complex_size = ((int)(size*1.0+0.5))/2+1;
    double freq = 0; 

    fftw_complex * middle = fftw_alloc_complex(complex_size);
    fftw_plan plan = fftw_plan_dft_r2c_1d (size, noise, middle, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    for (int i = 0; i < complex_size; ++i) {
        freq = (i+1);
        freq = pow(freq, exponent/2.0);
        /* std::cout << "freq["<< i <<"] = " << freq << std::endl; */
        middle[i][0] /= freq;
        middle[i][1] /= freq;
    }
    plan = fftw_plan_dft_c2r_1d(size, middle, noise, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    fftw_free(middle);

    freq = gsl_stats_mean(noise, 1, size);
    for (int i = 0; i < size; ++i) {
        noise[i] /= freq;
    }
    return 0;
}

int
test_noise(void )
{
    auto *x = new double[1000];
    auto *y = new double[1000];
    GenpinkNoise(y,2.0, 1000);
    for (int i = 0; i < 1000; ++i) {
        x[i] = i;
    }

    TCanvas *cr = new TCanvas("canvas1");
    cr->cd();
    TGraph *gr = new TGraph(1000,x,y);
    gr->SetTitle("Noise signal");
    gr->Draw("ALC");
    cr->Update();
    cr->Print("noise.pdf","Title:Noise");
    cr->Clear();
    
    int complex_size = ((int)(1000*1.0+0.5))/2+1;
    fftw_complex *spec = fftw_alloc_complex(complex_size); 
    fftw_plan plan = fftw_plan_dft_r2c_1d(1000, y, spec, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    for (int i = 0; i < complex_size; ++i) {
        y[i] = spec[i][0] * spec[i][0] + spec[i][1] * spec[i][1];
        y[i] = sqrt(y[i]);
    }
    delete gr; 
    gr = new TGraph(complex_size,x,y);
    gr->SetTitle("Power spectrum");
    gr->Draw("ALC");
    cr->Update();
    cr->Print("powerspec.pdf","Title:Powerspec");
    cr->Clear();
    delete cr;
    delete gr; 
    delete[] x;
    delete[] y;
    fftw_free(spec);
    
    return 0;
}
