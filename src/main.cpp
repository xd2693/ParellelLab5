#include <mpi.h>
#include <iostream>
#include <argparse.h>
#include <io.h>
#include <draw.h>
#include <barneshut.h>
#include <unordered_map>

using namespace std;
//using namespace std::tr1;

void run_seq(struct options_t opts, double** input, int n_vals){
    GLFWwindow* window;
    float colors[] = {1.0,0.0,0.0}; 
    if(opts.visual){       
        int result = draw_init(&window);        
        if (result<0){
            printf("fail to init window\n");
    }
    }

    for (int i = 0; i< opts.n_steps; i++){
        unordered_map<int, struct TreeNode> tree;
        tree_construct(tree, input, n_vals);
        if(opts.visual){
            glClear( GL_COLOR_BUFFER_BIT );
            for (auto const &pair : tree){
                drawOctreeBounds2D(pair.second);
            }
            for (int i = 0; i < n_vals; i++){       
                drawParticle2D(2*input[i][0]/4-1,2*input[i][1]/4-1, 0.008, colors);       
            }
            glfwSwapBuffers(window);
        }
        
        new_position(tree, input, n_vals, opts.theta, opts.time_step);
    }
       
    while(opts.visual && !glfwWindowShouldClose(window) )
    {
        glfwPollEvents();
        glfwSwapBuffers(window);
    }
}


int main(int argc, char **argv){
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    struct options_t opts;
    get_opts(argc, argv, &opts);
    
    int n_vals;
    double **input_vals, **output_vals;
    if (rank == 0){
        read_file(&opts, &n_vals, &input_vals, &output_vals);
        for(int i = 0; i< n_vals; i++){
            memcpy(output_vals[i], input_vals[i], sizeof(double) * 5);
        }
        
    }
    if (size == 1 )
        run_seq(opts, output_vals, n_vals);

    //int re = MPI_Bcast( n_vals , 1 , MPI_INT, 0 , MPI_COMM_WORLD);
        
    printf("This is process %d, n_vals = %d", rank, n_vals);
    
    //write_file(&opts, output_vals, n_vals);
    
    
    
    
    
    free(input_vals);
    free(output_vals);
    MPI_Finalize();
    return 0;
}

