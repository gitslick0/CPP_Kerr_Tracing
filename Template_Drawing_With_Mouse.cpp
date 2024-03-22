#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <vector>
#include <cmath>

const int screenWidth = 800;
const int screenHeight = 600;

GLuint VBO, VAO, EBO;

// Define variables for mouse interaction
bool mouseButtonDown = false;
double lastX = 0.0;
double lastY = 0.0;
double yaw = 0.0;
double pitch = 0.0;

void createSphere() {
    std::vector<GLfloat> vertices;
    std::vector<GLuint> indices;

    const int stacks = 20;
    const int slices = 20;
    const float radius = 1.0f;

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
}

void renderScene() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // Set rendering mode to wireframe

    glBindVertexArray(VAO);
    glDrawElements(GL_TRIANGLES, (GLsizei)(6 * 20 * 20), GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // Set rendering mode to wireframe
}

void mouseCallback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            mouseButtonDown = true;
            glfwGetCursorPos(window, &lastX, &lastY);
        } else if (action == GLFW_RELEASE) {
            mouseButtonDown = false;
        }
    }
}

void cursorPositionCallback(GLFWwindow* window, double xpos, double ypos) {
    if (mouseButtonDown) {
        double deltaX = xpos - lastX;
        double deltaY = ypos - lastY;
        lastX = xpos;
        lastY = ypos;

        const double sensitivity = 0.1;
        deltaX *= sensitivity;
        deltaY *= sensitivity;

        yaw += deltaX;
        pitch += deltaY;

        // Clamp pitch to avoid flipping
        if (pitch > 89.0)
            pitch = 89.0;
        if (pitch < -89.0)
            pitch = -89.0;

        // Update camera direction
        // Here, you can modify your camera transformation based on yaw and pitch
        // For simplicity, I'll just rotate the scene around the y-axis
        // For a more advanced camera system, consider using quaternions
        glLoadIdentity();
        glRotatef(static_cast<float>(yaw), 0.0f, 1.0f, 0.0f);
        glRotatef(static_cast<float>(pitch), 1.0f, 0.0f, 0.0f);
    }
}

int main() {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        return -1;
    }

    GLFWwindow* window = glfwCreateWindow(screenWidth, screenHeight, "OpenGL Sphere", nullptr, nullptr);
    if (!window) {
        std::cerr << "Failed to create GLFW window\n";
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW\n";
        return -1;
    }

    glViewport(0, 0, screenWidth, screenHeight);

    createSphere();

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        renderScene();
        glfwSwapBuffers(window);
    }

    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);

    glfwTerminate();
    return 0;
}
