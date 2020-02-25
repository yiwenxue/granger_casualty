#include <cstdio>
#include <cstdlib>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_statistics_double.h>
#include <mathematics.h>
#include <var.h>

TVar::TVar(){
    // The process should be 
    // 1. Init model( allocate memories & change size ).
    // 2. LoadData form outer array.
    // 3. Perform the fit for whole system.
    errors = NULL; data = NULL; model = NULL; A = NULL; x = NULL; b = NULL; p = NULL; y = NULL;
    constants = NULL; sigma = NULL;
}

TVar::~TVar(){
    freemem();
}

int 
TVar::freemem(void)
{
    if( errors != NULL )
    {
        for (int i = 0; i < dimension; ++i) {
            free(errors[i]);
        }
        free(errors);
    }

    if( data != NULL )
    {
        for (int i = 0; i < dimension; ++i) {
            free(data[i]);
        }
        free(data);
    }

    if( model != NULL )
    {
        for (int i = 0; i < latency; ++i) {
            for (int j = 0; j < dimension; ++j) {
                free(model[i][j]);
            }
            free(model[i]);
        }
        free(model);
    }
    if (A != NULL)
        gsl_matrix_free(A);
    if (x != NULL)
        gsl_vector_free(x);
    if (b != NULL)
        gsl_vector_free(b);
    if (p != NULL)
        gsl_permutation_free(p);
    if(sigma != NULL)
        free(sigma);
    if(constants != NULL)
        free(constants);
    if (y != NULL)
        free(y);

    return 0;
}

int 
TVar::Init(int _length,
           int _size, 
           int _latency)
{
    dimension = _size;
    data_length = _length;
    latency = _latency;
    /* Resize the memory to store data */
    // Don't need to check more
    // The logic is that: 
    // 1.with constructor function, set all pointers to nullptr,
    // 2.in the init function First check if the pointers are null,
    //      1.(null) allocate memories.
    //      2.(not null) free all memory, and reallocate memories.
    // Reason: you can change the size of data with this function, and don't need extra if and else.
    freemem();

    data = (double **)malloc(sizeof(double *) * dimension);
    errors = (double **)malloc(sizeof(double *) * dimension);
    for (int i = 0; i < dimension; ++i) {
        data[i] = (double *)malloc(sizeof(double ) * data_length);
        errors[i] = (double *)malloc(sizeof(double ) * data_length);
        for (int j = 0; j < data_length; ++j) {
            data[i][j] = 0.0;
            errors[i][j] = 0.0;
        }
    }

    sigma = (double *)malloc(sizeof(double) * dimension);
    constants = (double *)malloc(sizeof(double) * dimension);

    model = (double ***)malloc(sizeof(double **) * latency);
    for (int i = 0; i < latency; ++i) {
        model[i] = (double **)malloc(sizeof(double *) * dimension);
        for (int j = 0; j < dimension; ++j) {
            model[i][j] = (double *) malloc (sizeof(double) * dimension);
            for (int k = 0; k < dimension; ++k) {
                model[i][j][k] = 0.0;
            }
        }
    }

    Dimen = dimension * latency + 1;
    
    A = gsl_matrix_alloc(Dimen, Dimen);
    x = gsl_vector_alloc(Dimen);
    b = gsl_vector_alloc(Dimen);
    p = gsl_permutation_alloc(Dimen);
    noy = (latency + 1) * dimension;
    loy = data_length - latency;
    // Only pointers of arrays, not 2d array of double.
    y = (double **) malloc(sizeof(double *) * noy);

    /* First place */
    // data in y is: 
    // (latency=0, dimansion=1,2,3...)-(latency=1,......)-(latency=2,....)-(latency=3,.....)
    latency ++;
    for (int i = 0; i < latency; ++i) {
        for (int j = 0; j < dimension; ++j) {
            y[i*dimension+j] = data[j]+latency - i;
        }
    }
    latency --;
    return 0;
}

int TVar::LoadData(double **input, 
                   int length, 
                   int size)
{
    if ( length != data_length || dimension != size )
    {
        fprintf(stderr, "Size of data unmatched!\n");
        return -1;
    }

    for (int i = 0; i < dimension; ++i) 
    {
        for (int j = 0; j < data_length; ++j) {
            data[i][j] = input[i][j];
        }
    }
    return 0;
}

int
TVar::SurGen(int _length,
             int _dimension,
             int _latency)
{
    dimension = _dimension;
    data_length = _length;
    latency = _latency;
    Init(data_length, dimension, latency);
    for (int i = 0; i < dimension; ++i) {
        GenpinkNoise(data[i], 2.0, data_length);
    }
    for (int i = 0; i < dimension; ++i) {
        double mean=0;
        for (int j = 0; j < data_length; ++j) {
            mean += data[i][j] > 0 ? data[i][j]:-data[i][j];
        }
        mean /= data_length;
        for (int j = 0; j < data_length; ++j) {
            data[i][j] /= mean;
            data[i][j] += data[(i+2)%dimension][(j-2+data_length)%data_length];
        }
    }
    for (int i = 0; i < data_length; ++i) {
        printf("[%d]: ",i);
        for (int j = 0; j < dimension; ++j) {
            printf("%11.3lf,",data[j][i]);
        }
        printf("\n");
    }
    Solve();
    int Dimen = dimension * latency + 1;
    printf("The matrix:\n");
    for (int i = 0; i < Dimen; ++i) {
        for (int j = 0; j < Dimen; ++j) {
            printf("%11.3lf,",gsl_matrix_get(A, i, j));
        }
        printf("\n");
    }
    printf("The fit result:\n");
    for (int i = 0; i < latency; ++i) {
        printf("parameter[%d]:\n",i+1);
        for (int j = 0; j < dimension; ++j) {
            for (int k = 0; k < dimension; ++k) {
                printf("%11.3lf,", model[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }
    printf("All finished");
    return 0;
}

int
TVar::Solve()
{

    int s;
    // Fill the matrix with data.
    double temp;
    /* Second place */
    // Fill (0,0) with loy
    gsl_matrix_set(A, 0, 0, loy);
    // Fill (0,1~Dimen-1) with summition of y
    for (int j = 1; j < Dimen; ++j) {
        temp = 0.0;
        for (int k = 0; k < loy; ++k) {
            temp += y[j+dimension-1][k];
        }
        gsl_matrix_set(A, 0, j, temp);
    }
    // Fill the upper part 
    for (int i = 1; i < Dimen; ++i) {
        for (int j = i; j < Dimen; ++j) {
            temp = 0.0;
            for (int k = 0; k < loy; ++k) {
                temp += y[j+dimension-1][k] * y[i+dimension-1][k];
            }
            gsl_matrix_set(A, i, j, temp);
        }
    }
    // Fill 
    for (int i = 1; i < Dimen; ++i) {
        for (int j = 0; j < i; ++j) {
            gsl_matrix_set(A, i, j, gsl_matrix_get(A, j, i));
        }
    }
    std::cout << "Filled matrix" << std::endl;


    /* Third place */
    // Then Fill the vectors;
    for (int i = 0; i < dimension; ++i) {
        temp = 0.0; 
        // Fill first element in b(0);
        for (int k = 0; k < loy; ++k) {
            temp += y[i][k];
        }
        gsl_vector_set(b, 0, temp);
        // Fill the rest b(1~Dimen-1)
        for (int j = 1; j < Dimen; ++j) {
            temp = 0.0; 
            for (int k = 0; k < loy; ++k) {
                temp += y[j+dimension-1][k] * y[i][k];
            }
            gsl_vector_set(b, j, temp);
        }
        gsl_linalg_LU_decomp(A, p, &s);
        gsl_linalg_LU_solve(A, p, b, x);
        constants[i] = gsl_vector_get(x, 0);

        printf("%dth b\n",i);
        for (int j = 0;j  < Dimen; ++j) {
            printf("%10.3lf,",gsl_vector_get(b, j));
        }
        printf("\n");

        printf("%dth parameter\n",i);
        for (int j = 0;j  < Dimen; ++j) {
            printf("%10.3lf,",gsl_vector_get(x, j));
        }
        printf("\n");
        
        // Fill the result into matrix;
        for (int j = 0; j < latency ; ++j) {
            for (int k = 0; k < dimension; ++k) {
                // The shape of this matrix is (d x d)->(d x d)->(d x d)->....->(d x d)
                model[j][i][k] = gsl_vector_get(x, j*dimension+k+1);
            }
        }
    }
    return 0;
}

void *TVar::GetPointer(std::string options)
{
    return (void *)(data);
}

int printMatrix(double **data, 
                int raw, 
                int column)
{
    printf("A = {\n");
    for (int i = 0; i < column; ++i) {
        printf("[");
        for (int j = 0; j < raw; ++j) {
            if(j!=0)
                printf(",");
            printf("%8.3lf",data[j][i]);
        }
        printf("}");
        if(i != column - 1)
            printf(",\n");
    }
    printf("]\n");
    return 0;
}
