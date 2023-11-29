#include <mpi.h>
#include <iostream>
#include <argparse.h>
#include <io.h>
#include <draw.h>
#include <barneshut.h>
#include <unordered_map>
#include <thread>
#include <chrono>

using namespace std;
//using namespace std::tr1;

#ifdef DEBUG
    vector<pair<double, double>> force;
    vector<pair<double, double>> debug_force;
    vector<pair<double, double>> position;
    vector<pair<double, double>> new_pos;
    vector<pair<double, double>> velocity;
    vector<double> mass;
    vector<double> new_mass;
#endif

void run_seq(struct options_t opts, vector<particle>& particles, int n_vals){
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
        tree_construct(tree, particles, n_vals);
        //printf("Step %d tree done\n", i);
        if(opts.visual){
            glClear( GL_COLOR_BUFFER_BIT );
            for (auto const &pair : tree){
                drawOctreeBounds2D(pair.second);
            }
            for (int j = 0; j < n_vals; j++){       
                drawParticle2D(2*particles[j].px/4-1,2*particles[j].py/4-1, 0.008, colors);       
            }
            glfwSwapBuffers(window);
        }
        
        new_position(tree, particles, n_vals, opts.theta, opts.time_step);
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
                dp += pow((position[j].first - particles[j].px),2) + pow((position[j].second - particles[j].py),2);       
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
void parallel_parent(struct options_t opts, vector<particle>& particles, int n_vals, int size){
    GLFWwindow* window;
    float colors[] = {1.0,0.0,0.0}; 
    if(opts.visual){       
        int result = draw_init(&window);        
        if (result<0){
            printf("fail to init window\n");
        }
    }

    for (int i = 0; i< opts.n_steps; i++){
        int nCount = 0;// number of particles received from child
        for (int i = 0; i < n_vals; i++){
            MPI_Bcast( &particles[i] , sizeof(struct particle) , MPI_BYTE, 0 , MPI_COMM_WORLD);
            //printf("In process %d ", 0);
            //print_particle(particles[i]);
        }
        unordered_map<uint64_t, struct TreeNode> tree;
        //uint64_t total_weight = 0;
        tree_construct(tree, particles, n_vals);
        if(opts.visual){
            glClear( GL_COLOR_BUFFER_BIT );
            for (auto const &pair : tree){
                drawOctreeBounds2D(pair.second);
            }
            for (int j = 0; j < n_vals; j++){       
                drawParticle2D(2*particles[j].px/4-1,2*particles[j].py/4-1, 0.008, colors);       
            }
            glfwSwapBuffers(window);
        }
        while(nCount < n_vals){
            int flag;
            MPI_Status status;
            particle p;
            MPI_Iprobe( MPI_ANY_SOURCE , 1 , MPI_COMM_WORLD , &flag , &status);
            if (flag){
                MPI_Recv(&p, sizeof(struct particle), MPI_BYTE, status.MPI_SOURCE, 1, MPI_COMM_WORLD, &status);
                //p.weight = 1;
                if (p.mass < 0)
                    p.weight = 0;
                particles[p.index] = p;
                
                nCount ++;
                //printf("in process 0, %d particle received\n", nCount);
                //print_particle(p);
            }else{
                this_thread::sleep_for(chrono::milliseconds(10));
            }
        }
        
    }
    while(opts.visual && !glfwWindowShouldClose(window) )
    {
        glfwPollEvents();
        glfwSwapBuffers(window);
    }
}
void paralle_child(struct options_t opts, vector<particle>& particles, int n_vals, int size, int rank){
    

    
    
    double start_time, tree_time, end_time;
    for (int i = 0; i< opts.n_steps; i++){
        int total_weight = 0;
        int start = -1, end = 0;
        
        int workload = 0;
        int remain = 0;
        unordered_map<uint64_t, struct TreeNode> tree;
        //printf("In process %d bcast", rank);
        for(int j = 0; j< n_vals; j++){
            particle p;
            MPI_Bcast( &p , sizeof(struct particle) , MPI_BYTE, 0 , MPI_COMM_WORLD);
            total_weight += p.weight;
            //print_particle(p);
            particles[p.index] = p;
            
        }
        start_time = MPI_Wtime();
        
        tree_construct(tree, particles, n_vals);
        tree_time = MPI_Wtime();
        int temp=0;
        for (int j = 0; j < n_vals; j++){
            temp+=particles[j].weight;
            //print_particle(particles[j]);
        }
            
        //printf("total weight calculated %d, from tree %lld", temp, (long long)total_weight);
        /*for(auto const &pair : tree){
            print_node(pair.second);
        }*/
        
        int nCount = 0;
        if (opt){
            workload = total_weight / (size - 1);
            int current_weight = 0;
            vector<int> start_point(size);
            start_point[0] = 0;
            int p_id = 1;
            for(int j = 0; j < n_vals; j++){
                current_weight += particles[j].weight;
                if (current_weight >= workload ){
                    //printf("p_id %d, j %d\n, current_weight %d, weight %lld\n", p_id, j, current_weight,(long long) particles[j].weight);
                    start_point[p_id] = j + 1;

                    p_id ++;
                    if (p_id == size-1){
                        break;
                    }
                    current_weight = 0;
                    
                }                    
            }
            for (; p_id < size ; p_id++){
                start_point[p_id] = n_vals;
            }
            start = start_point[rank - 1];
        
            
            end = start_point[rank]-1;
            
        }else{
            workload = n_vals / (size - 1);
            //printf("ncount= %d\n", n_count);
            remain = n_vals % (size - 1);
            if (rank > remain){
                start = (workload + 1) * remain + workload * (rank - remain - 1);
                end = start + workload - 1;
                
            }else{
                start = (workload + 1) * (rank - 1);
                end = start + workload;
                //workload += 1;
            }
            
        }
        printf("workload=%d,process %d, start = %d, end = %d\n", workload, rank, start, end);
        MPI_Request* reqs = new MPI_Request[end - start + 1];////vector
        MPI_Status* status = new MPI_Status[end - start + 1];
        for(int j = start; j <= end; j++){
            if (particles[j].mass > 0){
                particles[j].weight = 0;
                update_position(tree, particles[j], opts.theta, opts.time_step);
            }
                           
            MPI_Isend( &particles[j], sizeof(struct particle) , MPI_BYTE , 0 , 1 , MPI_COMM_WORLD , &reqs[nCount]);
            //printf("in process %d, particle sent\n", rank);
            //print_particle(particles[j]);
            nCount++;
        }
        end_time = MPI_Wtime();
        MPI_Waitall(end - start + 1, reqs, status);
        delete[] reqs;
        delete[] status;
        printf("Child %d tree %f step time %f\n", rank, tree_time-start_time, end_time - start_time);
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
    double t1 = 0, t2  = 0;
    
    vector<particle> particles;
    if (rank == 0){
        double **input_vals, **output_vals;
        read_file(&opts, &n_vals, &input_vals, &output_vals);
        MPI_Bcast( &n_vals , 1 , MPI_INT, 0 , MPI_COMM_WORLD);

        for (int i = 0; i< n_vals; i++){
            memcpy(output_vals[i], input_vals[i], sizeof(double) * 5);
            particle p = {i, input_vals[i][0],input_vals[i][1],input_vals[i][2],input_vals[i][3],input_vals[i][4], 1};
            particles.push_back(p);            
            #ifdef DEBUG
                position.push_back(pair<double, double>(output_vals[i][0], output_vals[i][1]));
                new_pos.push_back(pair<double, double>(output_vals[i][0], output_vals[i][1]));
                //printf("DEBUG point %d %.10f %.10f %.10f %.10f\n", i, output_vals[i][0], output_vals[i][1], position.back().first, position.back().second);
                force.push_back(make_pair(0,0));
                debug_force.push_back(make_pair(0,0));
                velocity.push_back(make_pair(0,0));
                mass.push_back(output_vals[i][2]);
                new_mass.push_back(0);
            #endif
        }
        t1 = MPI_Wtime();
        if (size == 1 ){            
            run_seq(opts, particles, n_vals);
        }   
        else{
            parallel_parent(opts, particles, n_vals, size);
        }
            
        t2 = MPI_Wtime();
        printf("%f\n", t2 - t1);
        for (auto & p : particles){
            output_vals[p.index][0] = p.px;
            output_vals[p.index][1] = p.py;
            output_vals[p.index][2] = p.mass;
            output_vals[p.index][3] = p.vx;
            output_vals[p.index][4] = p.vy;
        }
        write_file(&opts, output_vals, n_vals);  
        free(input_vals);
        free(output_vals);  
    }else{
        MPI_Bcast( &n_vals , 1 , MPI_INT, 0 , MPI_COMM_WORLD);
        particles.resize(n_vals);
        paralle_child(opts, particles, n_vals, size, rank);
    }
    
    
    
    

        
    //printf("This is process %d, n_vals = %d\n", rank, n_vals);
    
#ifdef DEBUG
    printf("I am debugger %d\n", DEBUG);
#endif
    
    
    
    
    
    MPI_Finalize();
    return 0;
}

