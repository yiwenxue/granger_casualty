#ifndef CASUALTY_VAR
#define CASUALTY_VAR

#include <cstdio>
#include <cstdlib>

#include <mathematics.h>
#include <noise.h>

#include <polyfit.h>

typedef struct{
    double **y;
    int size;
}Ttime; 

class TAr{
    public:
        TAr();
        ~TAr();

    public:
        int Init();
        int LoadData();
        int Solve();

        int freemem();

    private:

};

class AR{
    public:
        AR();
        ~AR();

    public: 
        int Init(int _length, int _latency);
        int LoadData(double *_input, int _length);
        int SurGen(int _length, int _latency, double _noise);
        int Solve();
        int freemem();

    private:
        int length;
        int latency;
        int Dimen;
        int loy;
        int noy;

        double *data;
        double *errors;
        double *paramenter;
        double **y;

        gsl_matrix *A;
        gsl_vector *x, *b;
        gsl_permutation *p;
};


class TVar{
    public:
        TVar();
        ~TVar();

    public:
        /// @details Init() is the function to initialize this Object, it will free memories, and realloc new memories for new analysis 
        int Init(int _length, int _size, int _latency);
        /// @details LoadData() will load in the data after you initialized memories.
        int LoadData(double **input, int length, int size);
        /// surrogate data generate 
        int SurGen(int _length, int _size, int _latency);
        /// @details
        int Solve();
        /// @details 
        void *GetPointer(std::string options);
        /// @details 
        int freemem();

    private:
        /// @details This parameter determines the dimension of the matrix of parameters, for example, y is the vector y consists of 4 tracks, the dimension will be 4 (4x4).
        int dimension;
        /// @details latency determines the number of matrics will appear in the optimize. And it is also the range of influence of a value acts on the future.
        int latency;
        /// @details the length of input data
        int data_length;
        /// @details data stores the input arrays, data[raw][column], Or data[dimension][data_length].
        double **data; // allocated
        /// @details errors is the error term in Var model, same size as data, the square sum of every raw is the Sigma we want to calculate Casuality.
        double **errors; // allocated
        /// @details The square ssum of every raws of errors.
        double *sigma;
        /// @details constants the constant terms appear in the var model.
        double *constants; // allocated
        /// @details model stores the parameters we need to find. In this shape model[latency][dimension][dimension]
        double ***model; // allocated
        /// @details 
        // number of y elements y[noy]
        int noy;
        // length of y, y[noy][loy]
        int loy;
        // pointers to the data, this make up the x_i array for multi-linear;
        double **y; // allocated

        gsl_matrix *A;
        int Dimen;
        gsl_vector *x, *b;
        gsl_permutation *p;
};

int printMatrix(double **data, int raw, int column);

#endif
