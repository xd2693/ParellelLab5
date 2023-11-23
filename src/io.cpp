#include <io.h>


void read_file(struct options_t* args,
               int*              n_vals,
               double***             input_vals,
               double***             output_vals) {

  	// Open file
	std::ifstream in;
	in.open(args->in_file);
	// Get num vals
	in >> *n_vals;
	*input_vals = (double**)malloc((*n_vals) * sizeof(double*));
	*output_vals = (double**)malloc((*n_vals) * sizeof(double*));
		
	// Alloc input and output arrays
	for (int i = 0; i < *n_vals; i++){
		(*input_vals)[i] = (double*)malloc(5 * sizeof(double));
		(*output_vals)[i] = (double*)malloc(5 * sizeof(double));
	}
	
	// Read input vals
	double value;
	for (int i = 0; i < *n_vals; i++) {
		in >> value;
		//std::cout<<value<<" ";
		for (int j = 0; j < 5; j++){
			in >> (*input_vals)[i][j];
			//in >> value;
			//std::cout<<(*input_vals)[i][j]<<" ";
		}
		//std::cout<<""<<std::endl;
	}
}

void write_file(struct options_t* args,
               	double** output, int n_vals) {
  // Open file
	std::ofstream out;
	out.open(args->out_file, std::ofstream::trunc);

	out << n_vals << std::endl;;
	// Write solution to output file
	for (int i = 0; i < n_vals; ++i) {
		out << i << " ";
		for (int j = 0; j < 5; j++){
			out << output[i][j] << " ";
		}
		out << std::endl;
		
	}

	out.flush();
	out.close();
	
}
	