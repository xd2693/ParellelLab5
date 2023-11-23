#include <mpi.h>
#include <iostream>
#include <argparse.h>
#include <io.h>
#include <draw.h>
#include <barneshut.h>
#include <unordered_map>

using namespace std;
//using namespace std::tr1;

#ifdef DEBUG
    vector<pair<double, double>> force;
    vector<pair<double, double>> debug_force;
    vector<pair<double, double>> position;
    vector<pair<double, double>> velocity;
    vector<double> mass;
#endif

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
        unordered_map<uint64_t, struct TreeNode> tree;
        tree_construct(tree, input, n_vals);
        if(opts.visual){
            glClear( GL_COLOR_BUFFER_BIT );
            for (auto const &pair : tree){
                drawOctreeBounds2D(pair.second);
            }
            for (int j = 0; j < n_vals; j++){       
                drawParticle2D(2*input[j][0]/4-1,2*input[j][1]/4-1, 0.008, colors);       
            }
            glfwSwapBuffers(window);
        }
        
        new_position(tree, input, n_vals, opts.theta, opts.time_step);
        #ifdef DEBUG
        double dfx=0, dfy=0, dp=0;
            for (int j = 0; j < n_vals; j++){
                dfx += pow((debug_force[j].first - force[j].first),2);
                dfy += pow((debug_force[j].second - force[j].second),2);
                dp += pow((position[j].first - input[j][0]),2) + pow((position[j].second - input[j][1]),2);       
                
            }
            dfx = sqrt(dfx / n_vals);
            dfy = sqrt(dfy / n_vals);
            dp = sqrt(dp / n_vals);

            printf("step %d: dfx = %f, dfy = %f, dp = %f\n", i+1, dfx, dfy, dp);
        #endif
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
            #ifdef DEBUG
                position.push_back(make_pair(output_vals[i][0], output_vals[i][1]));
                force.push_back(make_pair(0,0));
                debug_force.push_back(make_pair(0,0));
                velocity.push_back(make_pair(0,0));
                mass.push_back(output_vals[i][2]);
            #endif
        }
        
    }
    if (size == 1 )
        run_seq(opts, output_vals, n_vals);

    //int re = MPI_Bcast( n_vals , 1 , MPI_INT, 0 , MPI_COMM_WORLD);
        
    printf("This is process %d, n_vals = %d\n", rank, n_vals);
    
    //write_file(&opts, output_vals, n_vals);
#ifdef DEBUG
    printf("I am debugger %d\n", DEBUG);
#endif
    
    
    
    
    free(input_vals);
    free(output_vals);
    MPI_Finalize();
    return 0;
}

