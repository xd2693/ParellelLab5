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
    vector<pair<double, double>> new_pos;
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
        printf("Step %d tree done\n", i);
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
        double dfx=0, dfy=0, dp=0, dfxp = 0, dfyp = 0;
            for (int j = 0; j < n_vals; j++){
                double ddfx = debug_force[j].first - force[j].first;
                double debug_fx_p = fabs(ddfx/debug_force[j].first);
                dfxp += debug_fx_p;
                dfx += (ddfx * ddfx);
                double ddfy = (debug_force[j].second - force[j].second);
                double debug_fy_p = fabs(ddfy/debug_force[j].second);
                dfyp += debug_fy_p;
                dfy += (ddfy * ddfy);
                dp += pow((position[j].first - input[j][0]),2) + pow((position[j].second - input[j][1]),2);       
                // if (debug_fx_p > 0.0 || debug_fy_p > 0.0)
                    // printf("Point %d error %f%% %f%%\n", j, debug_fx_p*100, debug_fy_p*100);
            }
            dfx = sqrt(dfx / n_vals);
            dfy = sqrt(dfy / n_vals);
            dp = sqrt(dp / n_vals);

            printf("step %d: dfx = %f, dfy = %f, dp = %f dfxp = %f, dfyp = %f\n", i+1, dfx, dfy, dp, dfxp, dfyp);
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
                position.push_back(pair<double, double>(output_vals[i][0], output_vals[i][1]));
                new_pos.push_back(pair<double, double>(output_vals[i][0], output_vals[i][1]));
                //printf("DEBUG point %d %.10f %.10f %.10f %.10f\n", i, output_vals[i][0], output_vals[i][1], position.back().first, position.back().second);
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

