#include <iostream>
#include <argparse.h>
#include <io.h>
#include <draw.h>
#include <barneshut.h>
#include <unordered_map>

using namespace std;
//using namespace std::tr1;

void run_seq(struct options_t opts, double** input, int n_vals){
    unordered_map<int, struct TreeNode> tree;
    
    GLFWwindow* window;
    int result = draw_init(&window);
    if (result<0){
        printf("fail to init window\n");
    }
    tree_construct(tree, input, n_vals);
       
    float colors[] = {1.0,0.0,0.0};
    
    
    glClear( GL_COLOR_BUFFER_BIT );
    for (auto const &pair : tree){
        drawOctreeBounds2D(pair.second);
    }
    for (int i = 0; i < n_vals; i++){       
        drawParticle2D(2*input[i][0]/4-1,2*input[i][1]/4-1,0.008, colors);       
    }

    
    while(!glfwWindowShouldClose(window))
    {
        glfwPollEvents();
        glfwSwapBuffers(window);
    }

}



int main(int argc, char **argv){
    struct options_t opts;
    get_opts(argc, argv, &opts);
    int n_vals;
    double **input_vals, **output_vals;
    read_file(&opts, &n_vals, &input_vals, &output_vals);
    
    
    
    
    run_seq(opts, input_vals, n_vals);

    
    free(input_vals);
    free(output_vals);
}

