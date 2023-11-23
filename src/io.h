#ifndef _IO_H
#define _IO_H

#include <argparse.h>
#include <iostream>
#include <fstream>

void read_file(struct options_t* args,
               int*              n_vals,
               double***             input_vals,
               double***             output_vals);

void write_file(struct options_t*         args,
                double** output, int n_vals);

#endif