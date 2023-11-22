#include<draw.h>
#include <stdio.h>

int draw_init(GLFWwindow** window){
/* OpenGL window dims */
    int width=600, height=600;
    
    if( !glfwInit() ){
        fprintf( stderr, "Failed to initialize GLFW\n" );
        return -1;
    }
    // Open a window and create its OpenGL context
    *window = glfwCreateWindow( width, height, "Simulation", NULL, NULL);
    if( *window == NULL ){
        fprintf( stderr, "Failed to open GLFW window.\n" );
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(*window); // Initialize GLEW
    if (glewInit() != GLEW_OK) {
        fprintf(stderr, "Failed to initialize GLEW\n");
        return -1;
    }
    // Ensure we can capture the escape key being pressed below
    glfwSetInputMode(*window, GLFW_STICKY_KEYS, GL_TRUE);
    return 0;
}

void drawParticle2D(double x_window, double y_window,
           double radius,
           float *colors) {
    int k = 0;
    float angle = 0.0f;
    glBegin(GL_TRIANGLE_FAN);
    glColor3f(colors[0], colors[1], colors[2]);
    glVertex2f(x_window, y_window);
    for(k=0;k<20;k++){
        angle=(float) (k)/19*2*3.141592;
        glVertex2f(x_window+radius*cos(angle),
           y_window+radius*sin(angle));
    }
    glEnd();
}
void drawOctreeBounds2D(TreeNode node) {
    
    
    printf("drawing lines! node %d\n", node.key);
    print_node(node);
    glBegin(GL_LINE_LOOP);
    // set the color of lines to be white
    glColor3f(1.0f, 1.0f, 1.0f);
    // specify the start point's coordinates
    glVertex2f(node.minX, node.minY);
    // specify the end point's coordinates
    glVertex2f(node.maxX, node.minY);
    // do the same for verticle line
    glVertex2f(node.maxX, node.maxY);
    glVertex2f(node.minX, node.maxY);
    glEnd();

}