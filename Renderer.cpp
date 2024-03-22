#include <iostream>
#include <cmath>
#include <GL/glut.h>
#include <GL/freeglut.h>
#include <array>
#include <vector>
#include <typeinfo>

class Renderer3D{
private:
    int height=600, width=800;
public:
    static Renderer3D *A;
    std::string name = "What?";
    void init(){
        glClearColor(1.0, 1.0, 1.0, 1.0); // Grey Background
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(-5, 5, -5, 5, -5, 5); // Orthographic Projection
        glMatrixMode(GL_MODELVIEW);
        A = this;
    }
    void drawAxis(std::array<double, 3> color, std::array<double, 3> vertex1, std::array<double, 3> vertex2 ){
        glBegin(GL_LINES);
        glColor3f(color[0], color[1], color[2]);
        glVertex3f(vertex1[0], vertex1[1], vertex1[2]);
        glVertex3f(vertex2[0], vertex2[1], vertex2[2]);
        glEnd();
    }
    void drawAxes(){
        glBegin(GL_LINES);
        // X-axis (red)
        this->drawAxis({1.0, 0.0, 0.0}, {-5.0, 0.0, 0.0}, {5.0, 0.0, 0.0});
        // Y-axis (red)
        this -> drawAxis({1.0, 0.0, 0.0}, {0.0, -5.0, 0.0}, {0.0, 5.0, 0.0});
        //Z-axis (red)
        this->drawAxis({1.0, 0.0, 0.0}, {0.0, 0.0, -5.0}, {0.0, 0.0, 5.0});
    }
    Renderer3D(int *argcp, char **argv, std::string name = "OpenGL Window"){
        glutInit(argcp, argv);
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
        glutInitWindowSize(this->width, this->height);
        glutCreateWindow(name.c_str());
        this->init();
    }
    Renderer3D(){};
    void show(){
        glutMainLoop();
    }
    static void display(){
        std::cout << A->name << std::endl;
    }
};

Renderer3D this_renderer;

void display_rnd(){
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glLoadIdentity();

        // View transformation
        gluLookAt(3.0, 3.0, 3.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
        this_renderer.drawAxes();
        std::cout << this_renderer.name << std::endl;
        glutSwapBuffers();
}

int main(int argc, char **argv){
    this_renderer = Renderer3D(&argc, argv, "Custom_Name");
    //auto display = [&this_renderer](void){
    //    display_rnd(this_renderer);
    //};
    Renderer3D rnd = Renderer3D(&argc, argv);
    glutDisplayFunc(display_rnd);
    this_renderer.show();

}