#include <iostream>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

class Renderer {
private:
    GLFWwindow* window;
    float radius;
    glm::vec3 rotationAngles;

    // Callback function for mouse movement
    static void mouse_callback(GLFWwindow* window, double xpos, double ypos) {
        Renderer* renderer = static_cast<Renderer*>(glfwGetWindowUserPointer(window));
        // Adjust rotation angles based on mouse movement
        renderer->rotationAngles.x += static_cast<float>(xpos - renderer->lastX) * 0.01f;
        renderer->rotationAngles.y += static_cast<float>(ypos - renderer->lastY) * 0.01f;
        renderer->lastX = xpos;
        renderer->lastY = ypos;
    }

public:
    Renderer(int width, int height, const char* title) : radius(1.0f), rotationAngles(0.0f) {
        if (!glfwInit()) {
            std::cerr << "Failed to initialize GLFW" << std::endl;
            exit(EXIT_FAILURE);
        }

        window = glfwCreateWindow(width, height, title, NULL, NULL);
        if (!window) {
            glfwTerminate();
            std::cerr << "Failed to create window" << std::endl;
            exit(EXIT_FAILURE);
        }

        glfwMakeContextCurrent(window);
        glfwSetWindowUserPointer(window, this);
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
        glfwSetCursorPosCallback(window, mouse_callback);

        glEnable(GL_DEPTH_TEST);
    }

    ~Renderer() {
        glfwDestroyWindow(window);
        glfwTerminate();
    }

    void createSphere(float radius) {
        this->radius = radius;
    }

    void render() {
        while (!glfwWindowShouldClose(window)) {
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            int width, height;
            glfwGetFramebufferSize(window, &width, &height);
            glViewport(0, 0, width, height);
            glMatrixMode(GL_PROJECTION);
            glLoadIdentity();
            gluPerspective(45.0f, (GLfloat)width / (GLfloat)height, 0.1f, 100.0f);
            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();

            // Apply rotation based on mouse movement
            glRotatef(rotationAngles.x, 0.0f, 1.0f, 0.0f);
            glRotatef(rotationAngles.y, 1.0f, 0.0f, 0.0f);

            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // Set rendering mode to wireframe

            // Render sphere
            glColor3f(1.0f, 0.0f, 0.0f); // Red color
            glutSolidSphere(radius, 30, 30);

            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // Set rendering mode to wireframe

            glfwSwapBuffers(window);
            glfwPollEvents();
        }
    }

private:
    double lastX = 0.0, lastY = 0.0;
};

int main(int argc, char** argv) {
    const int width = 800;
    const int height = 600;

    Renderer renderer(width, height, "OpenGL Window");
    glutInit(&argc, argv);
    renderer.createSphere(1.0f);
    renderer.render();

    return 0;
}
