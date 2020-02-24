#include <polyfit.h>

Polyfit::Polyfit(double *datax, double *datay, int _size, int _degree)
    : size(_size), degree(_degree)
{
    store = NULL; X = NULL; y = NULL; cov = NULL; c = NULL; ws = NULL;
    Init(_size, _degree);
    LoadData(datax, datay);
    Solve();
}

Polyfit::Polyfit(int _size, int _order)
{
    store = NULL; X = NULL; y = NULL; cov = NULL; c = NULL; ws = NULL;
    Init(_size,_order);
}

Polyfit::Polyfit(){
    store = NULL; X = NULL; y = NULL; cov = NULL; c = NULL; ws = NULL;
    printf("Hello World\n");
}

Polyfit::~Polyfit()
{
    delete[] store;
    gsl_multifit_linear_free(ws);
    gsl_matrix_free(X);
    gsl_matrix_free(cov);
    gsl_vector_free(y);
    gsl_vector_free(c);
}

int
Polyfit::Init(int new_size, int new_degree)
{
    degree = new_degree;
    size = new_size;

    /* Alloc memories for gsl linear fit */
    if (store != NULL)
        free(store);
    store = (double *)malloc(sizeof(double) * degree);

    if (X != NULL)
        gsl_matrix_free(X);
    X = gsl_matrix_alloc(size, degree);
    if (y != NULL)
        gsl_vector_free(y);
    y = gsl_vector_alloc(size);
    if (c != NULL)
        gsl_vector_free(c);
    c = gsl_vector_alloc(degree);
    if (cov != NULL)
        gsl_matrix_free(cov);
    cov = gsl_matrix_alloc(degree, degree);
    if (ws != NULL)
        gsl_multifit_linear_free(ws);
    ws = gsl_multifit_linear_alloc(size, degree);

    return 0;
}

double
Polyfit::Func(double input)
{
    double temp = store[0];
    for (i = 1; i < degree; ++i) {
        temp += store[i] * pow(input,i);
    }
    return temp;
}

int
Polyfit::LoadData(double *datax,double *datay)
{
    meanx = gsl_stats_mean(datax,1,size); 

    /* This is a fucking terrible code And there is a way to optimize it
     * For example, don't put the original data in to the matrix, but the other versioin 
     * get  */

    for(i=0; i < size; i++) {
        gsl_matrix_set(X, i, 0, 1.0);
        for(j=1; j < degree; j++) {
            gsl_matrix_set(X, i, j, pow(datax[i]-meanx, j));
        }
        gsl_vector_set(y, i, datay[i]);
    }
    tss = gsl_stats_tss(datay,1,size);

    return 0;
}


int
Polyfit::Solve()
{
    // calculate
    gsl_multifit_linear(X, y, c, cov, &coeff, ws);

    coeff = 1.0 - coeff/tss;

    for(i=0; i < degree; i++){
        store[i] = gsl_vector_get(c, i);
    }
    for(i=0;i<degree;i++){
        for(int j=i+1;j<degree;j++){
            int order = gsl_sf_fact(j)/gsl_sf_fact(i)/gsl_sf_fact(j-i);
            store[i]+= store[j]*order*pow(-1*meanx,j-i);
        }
    }
    // this will give you a string to describe the regression function
    func = std::to_string(store[0]);
    for (i = 1; i < degree; ++i) {
        func += " + " + std::to_string(store[i]) + "*x**" + std::to_string(i);
    }
    return 0;
}

double
Polyfit::benchmark(int runtimes)
{
    double *x = new double[1000];
    double *y = new double[1000];
    double tempx;
    srand(time(NULL));
    for (int i = 0; i < 1000; ++i) {
        x[i] = tempx = i*0.01 + (rand()%100)/10000.0;
        y[i] = - 0.56*tempx*tempx + 0.37*tempx + 1.68 + (rand()%100)/10000.0 ; 
    }

    clock_t current, stop;
    current = clock();

    for (int i = 0; i < runtimes; ++i) {
        Init(1000, 3);
        LoadData(x, y);
        Solve();
    }
    stop = clock();
    free(x);
    free(y);
    return (double)(stop - current)/CLOCKS_PER_SEC/runtimes;
}


Polyfit2::Polyfit2(double *datax, double *datay, int _size, int _degree)
    : size(_size), degree(_degree)
{
    A = NULL; x = NULL; b = NULL; p = NULL; store = NULL; data = NULL;
    Init(_size, _degree);
    LoadData(datax, datay);
    Solve();
}

Polyfit2::Polyfit2(int _size, int _order)
{
    A = NULL; x = NULL; b = NULL; p = NULL; store = NULL; data = NULL;
    Init(_size,_order);
}

Polyfit2::Polyfit2(){
    A = NULL; x = NULL; b = NULL; p = NULL; store = NULL; data = NULL;
    printf("Hello World\n");
}

Polyfit2::~Polyfit2()
{
    free(store);
    free(data);
    gsl_permutation_free(p);
    gsl_matrix_free(A);
    gsl_vector_free(x);
    gsl_vector_free(b);
}

int
Polyfit2::Init(int new_size, int new_degree)
{
    degree = new_degree;
    size = new_size;

    /* Alloc memories for gsl linear fit */
    if (store != NULL)
        free(store);
    store = (double *)malloc(sizeof(double) * degree);
    if(data != NULL)
        free(data);
    data = (double *) malloc(sizeof(double) * (degree*2-1) );

    if (A != NULL)
        gsl_matrix_free(A);
    A = gsl_matrix_alloc(degree, degree);
    if (x != NULL)
        gsl_vector_free(x);
    x = gsl_vector_alloc(degree);
    if (b != NULL)
        gsl_vector_free(b);
    b = gsl_vector_alloc(degree);
    if (p != NULL)
        gsl_permutation_free(p);
    p = gsl_permutation_alloc(degree);

    return 0;
}

double
Polyfit2::Func(double input)
{
    double temp = store[0];
    for (int i = 1; i < degree; ++i) {
        temp += store[i] * pow(input,i);
    }
    return temp;
}

int
Polyfit2::LoadData(double *datax,double *datay)
{
    meanx = gsl_stats_mean(datax,1,size); 

    /* This is a fucking terrible code And there is a way to optimize it
     * For example, don't put the original data in to the matrix, but the other versioin 
     * get  */

    int data_size = degree*2-1;

    for (int i = 0; i < data_size; ++i) {
        data[i] = 0;
        for (int j = 0; j < size; j++) {
            data[i] += pow(datax[j]-meanx,i);
        }
    }

    for (int i = 0; i < degree; ++i) {
        for (int j = 0; j < degree; ++j) {
            gsl_matrix_set(A, i, j, data[i+j]);
        }
    }

    for (int i = 0; i < degree; ++i) {
        data[i] = 0.0;
        for (int j = 0; j < size; j++) {
            data[i] += datay[j] * pow(datax[j]-meanx,i);
        }
        gsl_vector_set(b, i, data[i]);
    }

    tss = gsl_stats_tss(datay,1,size);

    return 0;
}


int
Polyfit2::Solve()
{

    int s;
    gsl_linalg_LU_decomp(A, p, &s);
    gsl_linalg_LU_solve(A, p, b, x);

    for (int i = 0; i < degree; ++i) {
        store[i] = gsl_vector_get(x,i);
    }

    for(int i=0;i<degree;i++){
        for(int j=i+1;j<degree;j++){
            int order = gsl_sf_fact(j)/gsl_sf_fact(i)/gsl_sf_fact(j-i);
            store[i]+= store[j]*order*pow(-1*meanx,j-i);
        }
    }
    // this will give you a string to describe the regression function
    func = std::to_string(store[0]);
    for (int i = 1; i < degree; ++i) {
        func += " + " + std::to_string(store[i]) + "*x**" + std::to_string(i);
    }
    return 0;
}

double
Polyfit2::benchmark(int runtimes)
{
    double *x = new double[1000];
    double *y = new double[1000];
    double tempx;
    srand(time(NULL));
    for (int i = 0; i < 1000; ++i) {
        x[i] = tempx = i*0.01 ;
        y[i] = - 0.56*tempx*tempx + 0.37*tempx + 1.68 ; 
    }

    clock_t current, stop;
    current = clock();

    for (int i = 0; i < runtimes; ++i) {
        Init(1000, 3);
        LoadData(x, y);
        Solve();
    }
    stop = clock();
    free(x);
    free(y);
    return (double)(stop - current)/CLOCKS_PER_SEC/runtimes;
}
