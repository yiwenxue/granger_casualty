#ifndef CASUALTY_POLYFIT
#define CASUALTY_POLYFIT

#include <mathematics.h>

/// Mathematical function for least squre polynominal fit 
/// @brief T Any float-point type such as float, double or long double
class Polyfit{
    public:
        Polyfit(double *datax, double *datay, int _size, int _order);
        Polyfit(int _size, int _order);
        Polyfit();

        virtual ~Polyfit();
        double Func(double input);
        inline int get_degree(){
            return degree;
        }
        inline double get_rsq(){
            return coeff;
        }
        int Init(int new_size, int new_degree);
        int LoadData(double *datax, double *datay);

        void recal(double *datax, double *datay); // size of the new seq should have the same size as the setting one;
        int Solve();

        double benchmark(int runtimes);
    public:
        std::string func; 
        double *store;
    private:
        gsl_multifit_linear_workspace *ws;
        gsl_matrix *cov, *X;
        gsl_vector *y, *c;
        int size;
        int degree;
        double coeff;
        double meanx;
        double tss;
        int i, j;
};

/// Mathematical function for least squre polynominal fit 
/// @brief T Any float-point type such as float, double or long double
class Polyfit2{
    public:
        Polyfit2(double *datax, double *datay, int _size, int _order);
        Polyfit2(int _size, int _order);
        Polyfit2();

        virtual ~Polyfit2();
        double Func(double input);
        inline int get_degree(){
            return degree;
        }
        inline double get_rsq(){
            return coeff;
        }
        int Init(int new_size, int new_degree);
        int LoadData(double *datax, double *datay);

        void recal(double *datax, double *datay); // size of the new seq should have the same size as the setting one;
        int Solve();

        double benchmark(int runtimes);
    public:
        std::string func; 
        double * store;
    private:
        double * data;
        gsl_permutation *p;
        gsl_matrix *A;
        gsl_vector *x, *b;
        int size;
        int degree;
        double coeff;
        double meanx;
        double tss;
};






#endif
