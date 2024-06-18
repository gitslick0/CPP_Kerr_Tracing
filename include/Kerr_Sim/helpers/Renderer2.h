#ifndef _Renderer_h_
#define _Renderer_h_

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include "ShaderClass.h"
#include "CameraClass.h"
#include "Drawables.h"
#include "./filesystem.h"
#include "Cubemap.h"

class Renderer {
public:
    Camera camera;
    Shader shader;
    GLFWwindow* window;
    float lastX, lastY;
    bool firstMouse = true;
    std::vector<Drawable*> Drawables;

    // timing
    float deltaTime = 0.1f;	// time between current frame and last frame
    float lastFrame = 0.0f;
    
    Renderer(int width, int height, const char* title) : width(width), height(height), title(title) {
        this->title = title;
        // Initialize glew to enable access to all of the OpenGL functions
        GLenum glfw_init_err = glfwInit();
        if(glfw_init_err != GLFW_TRUE){
            std::cerr << "Failed to initialize GLFW\n with error " << glfw_init_err << std::endl;
            exit(1);
        }

        // Set window Hints
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

        // Initialize window to draw in
        window = glfwCreateWindow(width, height, this->title, nullptr, nullptr);
        if (!window) {
            glfwTerminate();
            std::cerr << "Failed to create GLFW window\n";
            exit(EXIT_FAILURE);
        }

        glfwMakeContextCurrent(window);

        GLenum err = glewInit();
        if(err != GLEW_OK){
            std::cerr << "GlewInit faled with error return " << err << " and " << glewGetErrorString(err) << std::endl;
            exit(1);
        }

        // Set shader paths as defined in the compile configurations
        const char* vertexShaderPath = VERTEX_SHADER_PATH;
        const char* fragmentShaderPath = FRAGMENT_SHADER_PATH;
        std::cout << "Fragment Shader Path is " << FRAGMENT_SHADER_PATH << std::endl;

        // configure global opengl state
        // -----------------------------
        glEnable(GL_DEPTH_TEST);

        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
        glfwSetCursorPosCallback(window, cursor_position_callback);

        glfwSetWindowUserPointer(window, this);
        shader = Shader(vertexShaderPath, fragmentShaderPath);
        camera = Camera(glm::vec3(0.0f, 0.0f, 3.0f));
    }

    void addDrawable(Drawable* ToDraw){
        this->Drawables.push_back(ToDraw);
    }
    std::vector<Drawable*> getDrawables(){
        return this->Drawables;
    }

    void render() {

        Shader skyboxShader = Shader(SKYBOX_VERTEX_SHADER_PATH, SKYBOX_FRAGMENT_SHADER_PATH);
        
        float skyboxVertices[] = {
        // positions          
        -1.0f,  1.0f, -1.0f,
        -1.0f, -1.0f, -1.0f,
         1.0f, -1.0f, -1.0f,
         1.0f, -1.0f, -1.0f,
         1.0f,  1.0f, -1.0f,
        -1.0f,  1.0f, -1.0f,

        -1.0f, -1.0f,  1.0f,
        -1.0f, -1.0f, -1.0f,
        -1.0f,  1.0f, -1.0f,
        -1.0f,  1.0f, -1.0f,
        -1.0f,  1.0f,  1.0f,
        -1.0f, -1.0f,  1.0f,

         1.0f, -1.0f, -1.0f,
         1.0f, -1.0f,  1.0f,
         1.0f,  1.0f,  1.0f,
         1.0f,  1.0f,  1.0f,
         1.0f,  1.0f, -1.0f,
         1.0f, -1.0f, -1.0f,

        -1.0f, -1.0f,  1.0f,
        -1.0f,  1.0f,  1.0f,
         1.0f,  1.0f,  1.0f,
         1.0f,  1.0f,  1.0f,
         1.0f, -1.0f,  1.0f,
        -1.0f, -1.0f,  1.0f,

        -1.0f,  1.0f, -1.0f,
         1.0f,  1.0f, -1.0f,
         1.0f,  1.0f,  1.0f,
         1.0f,  1.0f,  1.0f,
        -1.0f,  1.0f,  1.0f,
        -1.0f,  1.0f, -1.0f,

        -1.0f, -1.0f, -1.0f,
        -1.0f, -1.0f,  1.0f,
         1.0f, -1.0f, -1.0f,
         1.0f, -1.0f, -1.0f,
        -1.0f, -1.0f,  1.0f,
         1.0f, -1.0f,  1.0f
    };
        // skybox VAO
        unsigned int skyboxVAO, skyboxVBO;
        glGenVertexArrays(1, &skyboxVAO);
        glGenBuffers(1, &skyboxVBO);
        glBindVertexArray(skyboxVAO);
        glBindBuffer(GL_ARRAY_BUFFER, skyboxVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(skyboxVertices), &skyboxVertices, GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);

        std::vector<std::string> faces
        {
            FileSystem::getPath("resources/textures/spacebox2/right.jpg"),
            FileSystem::getPath("resources/textures/spacebox2/left.jpg"),
            FileSystem::getPath("resources/textures/spacebox2/top.jpg"),
            FileSystem::getPath("resources/textures/spacebox2/bottom.jpg"),
            FileSystem::getPath("resources/textures/spacebox2/front.jpg"),
            FileSystem::getPath("resources/textures/spacebox2/back.jpg")
        };

        unsigned int cubemapTexture = loadCubemap(faces);

        skyboxShader.use();
        skyboxShader.setInt("skybox", 0);

        while (!glfwWindowShouldClose(window)) {
            // per-frame time logic
            // --------------------
            float currentFrame = static_cast<float>(glfwGetTime());
            float time = glfwGetTime();
            deltaTime = currentFrame - lastFrame;
            lastFrame = currentFrame;

            glClearColor(0.0f, 0.5f, 0.9f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            processInput(window);

            shader.use();
            // pass projection matrix to shader (note that in this case it could change every frame 
            glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)width / (float)height, 0.1f, 100.0f);
            shader.setMat4("projection", projection);

            // camera/view transformation
            glm::mat4 view = camera.GetViewMatrix();
            shader.setMat4("view", view);

            int jj = 0;
            std::cout << "length of drawable list (within shader) " << this->Drawables.size() << std::endl;

            for (const auto &element : this->Drawables){
            std::cout << "currently at element " << jj+1 << std::endl;
            if( auto p = dynamic_cast<DrawableDisk*>( element ) ){
                std::cout << "about to draw the disk";
            }
            if( auto p = dynamic_cast<DrawableLine*>( element ) ){
                std::cout << "about to draw the line";
            }
            std::cout << element << std::endl;
            element->update(time);
            std::cout << "update function called " << std::endl;
            shader.setMat4("model", element->getModelMatrix());

            if( auto p = dynamic_cast<DrawableBLSphere*>( element )){
                std::cout << "about to draw BH ";
            }

            //auto element2 = element;
            element->draw();
            jj++;
            std::cout << "drew element " << jj << std::endl;
            }


            // draw skybox as last
            glDepthFunc(GL_LEQUAL);  // change depth function so depth test passes when values are equal to depth buffer's content
            skyboxShader.use();
            view = glm::mat4(glm::mat3(view)); // remove translation from the view matrix
            skyboxShader.setMat4("view", view);
            skyboxShader.setMat4("projection", projection);
            // skybox cube
            glBindVertexArray(skyboxVAO);
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_CUBE_MAP, cubemapTexture);
            glDrawArrays(GL_TRIANGLES, 0, 36);
            glBindVertexArray(0);
            glDepthFunc(GL_LESS); // set depth function back to default
            
            // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
            // -------------------------------------------------------------------------------

            glfwSwapBuffers(window);
            glfwPollEvents();
        }

        // optional: de-allocate all resources once they've outlived their purpose:
        // ------------------------------------------------------------------------
        glDeleteVertexArrays(1, &skyboxVAO);
        glDeleteBuffers(1, &skyboxVBO);

        glfwTerminate();
    }

    static void cursor_position_callback(GLFWwindow* window, double xposIn, double yposIn) {
        Renderer* renderer = static_cast<Renderer*>(glfwGetWindowUserPointer(window));
        if (renderer) {
            float xpos = static_cast<float>(xposIn);
            float ypos = static_cast<float>(yposIn);

            if (renderer->firstMouse)
            {
                renderer->lastX = xpos;
                renderer->lastY = ypos;
                renderer->firstMouse = false;
            }

            float xoffset = xpos - renderer->lastX;
            float yoffset = renderer->lastY - ypos; // reversed since y-coordinates go from bottom to top

            renderer->lastX = xpos;
            renderer->lastY = ypos;

            renderer->camera.ProcessMouseMovement(xoffset, yoffset);
            // std::cout << "cursor position callback called" << std::endl;
            // renderer->camera.Yaw = 0.3;
            // std::cout << renderer->camera.Yaw << std::endl;
            // Implement method to rotate view based on mouse movement
            // This is beyond the scope of this example
        }
    }

    void processInput(GLFWwindow* win) {
        if (glfwGetKey(win, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
            glfwSetWindowShouldClose(win, true);
        }

        if (glfwGetKey(win, GLFW_KEY_W) == GLFW_PRESS){
            camera.ProcessKeyboard(FORWARD, deltaTime);
        }
        if (glfwGetKey(win, GLFW_KEY_S) == GLFW_PRESS)
            camera.ProcessKeyboard(BACKWARD, deltaTime);
        if (glfwGetKey(win, GLFW_KEY_A) == GLFW_PRESS)
            camera.ProcessKeyboard(LEFT, deltaTime);
        if (glfwGetKey(win, GLFW_KEY_D) == GLFW_PRESS)
            camera.ProcessKeyboard(RIGHT, deltaTime);

        // Implement camera movement using WASD keys
        // This is beyond the scope of this example
    }


    ~Renderer() {
        glfwDestroyWindow(window);
        glfwTerminate();
    }

private:
    int width;
    int height;
    const char* title;
};

#endif