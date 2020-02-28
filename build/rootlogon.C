#include <cstdio>

{
    printf("-----------------------------\n");
    printf("Project: Granger Casualty\n");
    printf("-----------------------------\n\n");
    printf("Loading Include Path...");
    gROOT->ProcessLine(".I ../src/");
    printf("Succeed!\nLoading Sources...");
    gROOT->ProcessLine(".L ../src/utility.cpp");
    gROOT->ProcessLine(".L ../src/var.cpp");
    gROOT->ProcessLine(".L ../src/noise.cpp");
    gROOT->ProcessLine(".L ../src/statis.cpp");
    gROOT->ProcessLine(".L ../src/polyfit.cpp");
    printf("Succeed!\nLoading Library...");
    gROOT->ProcessLine(".L /usr/lib/libfftw3.so");
    gROOT->ProcessLine(".L /usr/lib/libgslcblas.so");
    gROOT->ProcessLine(".L /usr/lib/libgsl.so");
    gROOT->ProcessLine("usleep(100000)");
    printf("Succeed!\n");
    init();
}

int init(){

    printf("Initializing ...\n");
    return 0;
}


