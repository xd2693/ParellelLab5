#ifndef _ARGPARSE_H
#define _ARGPARSE_H

#include <getopt.h>
#include <stdlib.h>
#include <iostream>

struct options_t {
    char *in_file;
    char *out_file;
    int n_steps;
    double theta;
    double time_step;
    bool visual;
};

void get_opts(int argc, char **argv, struct options_t *opts);
#endif