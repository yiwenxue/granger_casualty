#include <stdio.h>
#include <utility.h>
#include <mathematics.h>
#include <polyfit.h>
#include <noise.h>
#include <var.h>

int main(int argc, char *argv[])
{
    printf("The project of causality\n");
    auto fit = new Polyfit();
    auto fit2 = new Polyfit2();

    int times;
    double runtime;

    times= 1000;
    runtime = fit->benchmark(times);
    printf("Runtime(%8d Average): %lf\n",times,runtime);

    times = 5000;
    runtime = fit->benchmark(times);
    printf("Runtime(%8d Average): %lf\n",times,runtime);

    times = 10000;
    runtime = fit->benchmark(times);
    printf("Runtime(%8d Average): %lf\n",times,runtime);

    times= 1000;
    runtime = fit2->benchmark(times);
    printf("Runtime 2(%6d Average): %lf\n",times,runtime);

    times = 5000;
    runtime = fit2->benchmark(times);
    printf("Runtime 2(%6d Average): %lf\n",times,runtime);

    times = 10000;
    runtime = fit2->benchmark(times);
    printf("Runtime 2(%6d Average): %lf\n",times,runtime);
    return 0;
}
