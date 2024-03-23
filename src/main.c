#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "test.h"

void printUsage()
{
    printf("Usage: ./pcg [OPTIONS]\n");
    printf("Options:\n");
    printf("  --test\tTests correctness of the PCG.\n");
    printf("  --bench\tBenchmarks the PCG on conservative and aggressive parameters.\n");
}

void runBenchmarks(size_t n, size_t t, int num_trials)
{
    double time = 0;

    for (int i = 0; i < num_trials; i++)
    {
        time += bench_pcg_c4(n, t);
        printf("Done with trial %i of %i\n", i + 1, num_trials);
    }
    printf("******************************************\n");
    printf("Avg time (N=3^%zu, c=4, t=%zu): %0.4f ms\n", n, t, time / num_trials);
    printf("******************************************\n\n");
}

int main(int argc, char **argv)
{
    int num_trials = 5;

    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "--bench") == 0)
        {
            printf("******************************************\n");
            printf("Benchmarking PCG with conservative parameters (c=4, t=20)\n");
            runBenchmarks(14, 20, num_trials);
            runBenchmarks(16, 20, num_trials);
            runBenchmarks(18, 20, num_trials);

            printf("******************************************\n");
            printf("Benchmarking PCG with aggressive parameters (c=4, t=18)\n");
            runBenchmarks(14, 18, num_trials);
            runBenchmarks(16, 18, num_trials);
            runBenchmarks(18, 18, num_trials);
        }
        else if (strcmp(argv[i], "--test") == 0)
        {
            printf("******************************************\n");
            printf("Testing PCG\n");
            test_pcg();
            printf("******************************************\n");
            printf("PASS\n");
            printf("******************************************\n\n");
        }
        else
        {
            printUsage();
        }
    }

    if (argc == 1)
        printUsage();
}
