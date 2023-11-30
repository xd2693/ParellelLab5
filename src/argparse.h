#ifndef _ARGPARSE_H
#define _ARGPARSE_H

#include <getopt.h>
#include <stdlib.h>
#include <iostream>

struct options_t {
    char *in_file;          //input file path
    char *out_file;         //output file path
    int n_steps;            //steps
    double theta;           //theta
    double time_step;       //time_step
    bool visual;            //optional visualization
    bool op;                //optional optimization
};

void get_opts(int argc, char **argv, struct options_t *opts);
#endif