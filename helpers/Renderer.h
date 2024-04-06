#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <string>
#include <vector>
#include <memory>
#include"Drawables.h"


class Renderer {
private:

public:
    Shader shader;
    std::vector<Drawable*> Drawables;
    GLFWwindow* window;

    // settings
    unsigned int SCR_WIDTH = 800;
    unsigned int SCR_HEIGHT = 600;
    std::string title = "Renderer";

    // camera
    Camera camera = camera(glm::vec3(0.0f, 0.0f, 3.0f));
    float lastX = SCR_WIDTH / 2.0f;
    float lastY = SCR_HEIGHT / 2.0f;
    bool firstMouse = true;

    // timing
    float deltaTime = 0.0f;	// time between current frame and last frame
    float lastFrame = 0.0f;

    Renderer(){};
    Renderer(int inp_width, int inp_height, std::string inp_title) : {
        // glfw: initialize and configure
        // ------------------------------
        glfwInit();
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

        this->SCR_WIDTH = inp_width;
        this->SCR_HEIGHT = inp_height;
        this->title = inp_title;

        float lastX = SCR_WIDTH / 2.0f;
        float lastY = SCR_HEIGHT / 2.0f;
        bool firstMouse = true;

        glEnable(GL_DEPTH_TEST);
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    }

    // Mouse input stuff
    static void mouseCallback(GLFWwindow* window, int button, int action, int mods) {
    Renderer* renderer = static_cast<Renderer*>(glfwGetWindowUserPointer(window));
    float xpos = static_cast<float>(xposIn);
    float ypos = static_cast<float>(yposIn);

    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

    lastX = xpos;
    lastY = ypos;

    camera.ProcessMouseMovement(xoffset, yoffset);
    }

    static void cursorPositionCallback(GLFWwindow* window, double xpos, double ypos) {
    Renderer* renderer = static_cast<Renderer*>(glfwGetWindowUserPointer(window));
    if (renderer->mouseButtonDown) {
        double deltaX = xpos - renderer->lastX;
        double deltaY = ypos - renderer->lastY;
        renderer->lastX = xpos;
        renderer->lastY = ypos;

        const double sensitivity = 0.1;
        deltaX *= sensitivity;
        deltaY *= sensitivity;

        renderer->yaw += deltaX;
        renderer->pitch += deltaY;

        // Clamp pitch to avoid flipping
        if (renderer->pitch > 89.0)
            renderer->pitch = 89.0;
        if (renderer->pitch < -89.0)
            renderer->pitch = -89.0;

        // Update camera direction
        // Here, you can modify your camera transformation based on yaw and pitch
        // For simplicity, I'll just rotate the scene around the y-axis
        // For a more advanced camera system, consider using quaternions
        glLoadIdentity();
        glRotatef(static_cast<float>(renderer->yaw), 0.0f, 1.0f, 0.0f);
        glRotatef(static_cast<float>(renderer->pitch), 1.0f, 0.0f, 0.0f);
    }
    }

    ~Renderer() {
        glfwDestroyWindow(window);
        glfwTerminate();
    }


    void renderScene() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);

    for (const auto &element : this->Drawables){
        //auto element2 = element;
        element->draw();
        }

    glfwPollEvents();
    glfwSwapBuffers(window);

    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // Set rendering mode to wireframe
    // This was a drawing command previously
    //glBindVertexArray(VAO);
    //glDrawElements(GL_TRIANGLES, (GLsizei)(6 * 20 * 20), GL_UNSIGNED_INT, 0);
    //glBindVertexArray(0);

    //glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // Set rendering mode to wireframe
    }
    void addDrawable(Drawable* ToDraw){
        this->Drawables.push_back(ToDraw);
    }
    void initWindow(){
        if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        exit(EXIT_FAILURE);
    }

    this->window = glfwCreateWindow(this->width, this->height, "OpenGL Sphere", nullptr, nullptr);
    if (!window) {
        std::cerr << "Failed to create GLFW window\n";
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glfwMakeContextCurrent(this->window);
    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW\n";
        exit(EXIT_FAILURE);
    }

    glfwSetWindowUserPointer(this->window, this);
    glfwSetInputMode(this->window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

    // Set up callbacks for mouse interaction
    glfwSetMouseButtonCallback(this->window, this->mouseCallback);
    glfwSetCursorPosCallback(this->window, this->cursorPositionCallback);

    glEnable(GL_DEPTH_TEST);
    glViewport(0, 0, this->width, this->height);
    }

    void show(){
        while (!glfwWindowShouldClose(this->window)) {
        glfwPollEvents();
        this->renderScene();
        glfwSwapBuffers(this->window);
    }
    }
    std::vector<Drawable*> get_DrawableList(){
        return this->Drawables;
    }
};
