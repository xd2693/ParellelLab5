#ifndef _DRAW_H
#define _DRAW_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <GL/glu.h>
#include <GL/glut.h>
#include <barneshut.h>

int draw_init(GLFWwindow** window);
void drawParticle2D(double x_window, double y_window,
           double radius,
           float *colors) ;
void drawOctreeBounds2D(TreeNode node);

#endif

