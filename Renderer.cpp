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

// Define virtual Base Class for Drawables in the Renderer
class Drawable{
public:
    virtual void draw() {std::cout << "this would be a drawing method";}
    virtual void show(){std::cout << "this would be a showing method";}
    virtual void draw2() const {std::cout << "this would be a drawing method";}

    virtual ~Drawable() {}
};

class DrawableSphere : public Drawable{
private:
    double radius;
public:
    DrawableSphere(){}
    GLuint VBO, VAO, EBO;
    DrawableSphere(double inp_radius){
        this->radius = inp_radius;
    }
    void draw() override {
    std::vector<GLfloat> vertices;
    std::vector<GLuint> indices;

    const int stacks = 20;
    const int slices = 20;
    const float radius = this->radius;

    for (int i = 0; i <= stacks; ++i) {
        float phi = static_cast<float>(i) * M_PI / stacks;
        for (int j = 0; j <= slices; ++j) {
            float theta = static_cast<float>(j) * 2.0f * M_PI / slices;
            float x = radius * std::sin(phi) * std::cos(theta);
            float y = radius * std::cos(phi);
            float z = radius * std::sin(phi) * std::sin(theta);

            vertices.push_back(x);
            vertices.push_back(y);
            vertices.push_back(z);
        }
    }

    for (int i = 0; i < stacks; ++i) {
        for (int j = 0; j < slices; ++j) {
            indices.push_back((i + 1) * (slices + 1) + j);
            indices.push_back(i * (slices + 1) + j);
            indices.push_back(i * (slices + 1) + j + 1);

            indices.push_back((i + 1) * (slices + 1) + j);
            indices.push_back(i * (slices + 1) + j + 1);
            indices.push_back((i + 1) * (slices + 1) + j + 1);
        }
    }

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(GLfloat), vertices.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), indices.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // Set rendering mode to wireframe

    glBindVertexArray(VAO);
    glDrawElements(GL_TRIANGLES, (GLsizei)(6 * 20 * 20), GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // Set rendering mode to wireframe
    }

    void draw2() const override {
        std::cout << "This is not a drawing method" << std::endl;
    }
};



class Renderer {
private:
    glm::vec3 rotationAngles;
    // Define variables for mouse interaction
    bool mouseButtonDown = false;
    double lastX = 0.0;
    double lastY = 0.0;
    double yaw = 0.0;
    double pitch = 0.0;

    int width, height;
    std::string title;

public:
    std::vector<Drawable*> Drawables;
    GLFWwindow* window;
    Renderer(){};
    Renderer(int inp_width, int inp_height, const char* inp_title) : rotationAngles(0.0f) {
        if (!glfwInit()) {
            std::cerr << "Failed to initialize GLFW" << std::endl;
            exit(EXIT_FAILURE);
        }

        this->width = inp_width;
        this->height = inp_height;
        this->title = static_cast<std::string>(inp_title);
        this->window = glfwCreateWindow(this->width, this->height, inp_title, NULL, NULL);
        if (!this->window) {
            glfwTerminate();
            std::cerr << "Failed to create window" << std::endl;
            exit(EXIT_FAILURE);
        }

        glfwMakeContextCurrent(this->window);
        glfwSetWindowUserPointer(this->window, this);
        glfwSetInputMode(this->window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);


        // Set up callbacks for mouse interaction
        glfwSetMouseButtonCallback(this->window, this->mouseCallback);
        glfwSetCursorPosCallback(this->window, this->cursorPositionCallback);

        // !!!! TODO Mouse Setup Still doesn't work correctly for some reason !!!


        glEnable(GL_DEPTH_TEST);
    }

    // Mouse input stuff
    static void mouseCallback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
    Renderer* renderer = static_cast<Renderer*>(glfwGetWindowUserPointer(window));
        if (action == GLFW_PRESS) {
            renderer->mouseButtonDown = true;
            glfwGetCursorPos(window, &renderer->lastX, &renderer->lastY);
        } else if (action == GLFW_RELEASE) {
            renderer->mouseButtonDown = false;
        }
        }
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

int main(int argc, char** argv){
    Renderer rnd = Renderer(800, 600, "Jup");
    //std::vector<Drawable*> objects;
    Drawable Base;
    DrawableSphere Derived = DrawableSphere(1.0);

    //rnd.addDrawable(&Base);
    rnd.addDrawable(&Derived);

    //objects.push_back(&Base);
    //objects.push_back(&Derived);

    //for(const auto& objPtr : objects){
    //    objPtr->draw();
    //}
    rnd.initWindow();
    for(const auto& objPtr : rnd.Drawables){
        objPtr->draw2();
    }
    rnd.renderScene();
    rnd.show();
    /*Renderer rnd = Renderer(800, 600, "Whaaaat?");
    DrawableSphere Sph = DrawableSphere(1.0);
    DrawableSphere* Sphp = &Sph;
    std::cout << typeid(Sph).name() << std::endl;
    Sph.draw();
    std::cout << typeid(*Sphp).name() << std::endl;
    std::vector<Drawable*> List;
    std::cout << typeid(List[1]).name() << std::endl;
    List.push_back(Sphp);
    //rnd.initWindow();
    //rnd.show();
    */
}