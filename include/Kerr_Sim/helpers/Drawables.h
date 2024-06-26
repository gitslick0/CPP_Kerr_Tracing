#ifndef _Drawables_h_
#define _Drawables_h_

#include<vector>
#include<iostream>
#include<GL/glew.h>
#include<GL/glu.h>
#include<glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include"coordinates.h"

class Drawable{
private:
    glm::vec3 global_position = glm::vec3(0.0f, 0.0f, 0.0f);
    glm::mat4 model = glm::mat4(1.0f);
public:
    virtual void draw(){std::cout << "this is the abstract drawing method of Drawables" << std::endl;}
    virtual void show(){std::cout << "this is the abstract showing method of Drawables" << std::endl;}
    virtual void draw2(){std::cout << "this is the virtual testing function draw2 of Drawables" << std::endl;}

    virtual ~Drawable(){}
    void setGlobalPosition(glm::vec3 inp_global_position){
        this->model = glm::mat4(1.0f); // reset position (if was set before)
        this->model = glm::translate(this->model, inp_global_position);
        this->global_position = inp_global_position;
    }
    glm::vec3 getGlobalPosition(){
        return global_position;
    }
    glm::mat4 getModelMatrix(){
        return this->model;
    }

    void setModelMatrix(glm::mat4 inp_model){
        this->model = inp_model;
    }
    
    // implement update of the model matrix or other here
    //std::cout << "this is the abstract update method for Drawables" << std::endl;
    virtual void update(float time) = 0;
};

class DrawableBLSphere : public Drawable{
    // Implements a Drawable Object that has constant radius in Boyer Lindquist Coordinates
private:
    double radius, a_spin;
    glm::vec3 color = glm::vec3(0.0f, 0.0f, 0.0f);
public:
    DrawableBLSphere(){}
    GLuint VBO, VAO, EBO;
    const int stacks = 100, slices = 100;
    DrawableBLSphere(double inp_radius, double inp_a_spin ){
        this->radius = inp_radius;
        this->a_spin = inp_a_spin;
    }
    std::vector<GLfloat> createVertices(){
        std::vector<GLfloat> vertices;
        const int stacks = this->stacks;
        const int slices = this->slices;

        for (int i = 0; i <= stacks; ++i) {
            float phi = static_cast<float>(i) * M_PI / stacks;
            for (int j = 0; j <= slices; ++j) {
                float theta = static_cast<float>(j) * 2.0f * M_PI / slices;
                CoordVec3 BL_Coords = CoordVec3(radius, theta, phi, false, true);
                CoordVec3 Cart_Coords = BL_Coords.to_cart(a_spin);

                vertices.push_back(Cart_Coords.x);
                vertices.push_back(Cart_Coords.z);
                vertices.push_back(Cart_Coords.y);

                // push back color
                vertices.push_back(this->color.x); vertices.push_back(this->color.y); vertices.push_back(this->color.z);
            }
        }
        return vertices;
    }

    std::vector<GLuint> createIndices(){
        std::vector<GLuint> indices;
        const int stacks = this->stacks;
        const int slices = this->slices;

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

        return indices;
    }

    void draw() override {
    
    std::vector<GLfloat> vertices = this->createVertices();
    std::vector<GLuint> indices = this->createIndices();

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(GLfloat), vertices.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), indices.data(), GL_STATIC_DRAW);

    // Set Position Data
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(0);

    // Set Color Data
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(GLfloat), (void*)(3*sizeof(GLfloat)));
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // Set rendering mode to wireframe

    glBindVertexArray(VAO);
    glDrawElements(GL_TRIANGLES, (GLsizei)(6 * stacks * slices), GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // Set rendering mode to wireframe
    }

    void draw2() override {
        std::cout << "This is not a drawing method" << std::endl;
    }

    void update(float time) override {
        glm::mat4 updated_model_matrix = glm::mat4(1.0f);
        updated_model_matrix = glm::rotate(updated_model_matrix, 3.0f*time, glm::vec3(0.0, 1.0, 0.0));
        updated_model_matrix = glm::translate(updated_model_matrix, this->getGlobalPosition());
        this->setModelMatrix(updated_model_matrix);
    }

    void setColor(glm::vec3 inp_color){
        this->color = inp_color;
    }
    glm::vec3 getColor(){
        return this->color;
    }
};

class DrawableSphere : public Drawable{
private:
    double radius;
    glm::vec3 color = glm::vec3(1.0f, 1.0f, 1.0f);
public:
    DrawableSphere(){}
    GLuint VBO, VAO, EBO;
    const int stacks = 100, slices = 100;
    DrawableSphere(double inp_radius){
        this->radius = inp_radius;
    }
    std::vector<GLfloat> createVertices(){
        std::vector<GLfloat> vertices;
        const int stacks = this->stacks;
        const int slices = this->slices;
        const float radius = this->radius;

        for (int i = 0; i <= stacks; ++i) {
            float phi = static_cast<float>(i) * M_PI / stacks;
            for (int j = 0; j <= slices; ++j) {
                float theta = static_cast<float>(j) * 2.0f * M_PI / slices;
                float x = radius * std::sin(phi) * std::cos(theta);
                float y = radius * std::sin(phi) * std::sin(theta);
                float z = radius * std::cos(phi);

                vertices.push_back(x);
                vertices.push_back(z);
                vertices.push_back(y);

                // push back color
                vertices.push_back(this->color.x); vertices.push_back(this->color.y); vertices.push_back(this->color.z);
            }
        }
        return vertices;
    }

    std::vector<GLuint> createIndices(){
        std::vector<GLuint> indices;
        const int stacks = this->stacks;
        const int slices = this->slices;

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

        return indices;
    }

    void draw() override {
    
    std::vector<GLfloat> vertices = this->createVertices();
    std::vector<GLuint> indices = this->createIndices();

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(GLfloat), vertices.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), indices.data(), GL_STATIC_DRAW);

    // Set Position Data
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(0);

    // Set Color Data
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(GLfloat), (void*)(3*sizeof(GLfloat)));
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // Set rendering mode to wireframe

    glBindVertexArray(VAO);
    glDrawElements(GL_TRIANGLES, (GLsizei)(6 * stacks * slices), GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // Set rendering mode to wireframe
    }

    void draw2() override {
        std::cout << "This is not a drawing method" << std::endl;
    }

    void update(float time) override {
        glm::mat4 updated_model_matrix = glm::mat4(1.0f);
        updated_model_matrix = glm::rotate(updated_model_matrix, 3.0f*time, glm::vec3(0.0, 1.0, 0.0));
        updated_model_matrix = glm::translate(updated_model_matrix, this->getGlobalPosition());
        this->setModelMatrix(updated_model_matrix);
    }

    void setColor(glm::vec3 inp_color){
        this->color = inp_color;
    }
    glm::vec3 getColor(){
        return this->color;
    }
};

class DrawableDisk : public Drawable{
private:
    double r_in, r_out;
    GLuint VBO, VAO, EBO;
    const int slices = 120;
    double height = 0.0;
public:
    bool is_disk = true;
    DrawableDisk(){};
    DrawableDisk(double inp_r_in, double inp_r_out, double inp_height){
        this->r_in = inp_r_in;
        this->r_out = inp_r_out;
        this->height = inp_height;
    }
   std::vector<GLfloat> createVertices(){
        std::vector<GLfloat> vertices;
        const int slices = this->slices;
        const float radius_in = this->r_in;
        const float radius_out = this->r_out;

        for (int i = 0; i <= slices; ++i) {
            float phi = static_cast<float>(i) * 2*M_PI / slices;
                float x_in = radius_in*std::cos(phi);
                float y_in = radius_in*std::sin(phi);
                float z_in = this->height;

                float x_out = radius_out*std::cos(phi);
                float y_out = radius_out*std::sin(phi);
                float z_out = this->height;


                vertices.push_back(x_in); vertices.push_back(z_in); vertices.push_back(y_in);
                vertices.push_back(std::abs(std::cos(phi/2+0.6))); vertices.push_back(std::abs(std::cos(phi/2+0.6))); vertices.push_back(std::abs(0.7f*std::cos(phi/2+0.6)));
                vertices.push_back(x_out); vertices.push_back(z_out); vertices.push_back(y_out);
                vertices.push_back(0.0); vertices.push_back(0.0); vertices.push_back(0.0);
            }
        return vertices;
    }

    std::vector<GLuint> createIndices(){
        std::vector<GLuint> indices;
        const int slices = this->slices;

        for (int i = 0; i < slices; ++i) {
            for (int j = 0; j < slices; ++j) {
                indices.push_back(i*2);
                indices.push_back((i*2)+1);
                indices.push_back((i+1)*2 + 1);

                indices.push_back(i*2);
                indices.push_back((i+1)*2+1);
                indices.push_back((i+1)*2);
            }
        }

        return indices;
    }

    void draw() override {

    std::vector<GLfloat> vertices = this->createVertices();
    std::vector<GLuint> indices = this->createIndices();

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(GLfloat), vertices.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), indices.data(), GL_STATIC_DRAW);

    // Set Position Data
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(0);

    // Set Color Data
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(GLfloat), (void*)(3*sizeof(GLfloat)));
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // Set rendering mode to wireframe

    glBindVertexArray(VAO);
    glDrawElements(GL_TRIANGLES, (GLsizei)(6 * slices * slices), GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // Set rendering mode to wireframe
    }

    void draw2() override {
        std::cout << "This is not a drawing method" << std::endl;
    }

    void update(float time) override {
        glm::mat4 updated_model_matrix = glm::mat4(1.0f);
        updated_model_matrix = glm::rotate(updated_model_matrix, 3.0f*time, glm::vec3(0.0, 1.0, 0.0));
        updated_model_matrix = glm::translate(updated_model_matrix, this->getGlobalPosition());
        this->setModelMatrix(updated_model_matrix);
        /*glm::mat4 updated_model_matrix = glm::mat4(1.0f);
        updated_model_matrix = this->getModelMatrix();
        updated_model_matrix = glm::translate(updated_model_matrix, -this->getGlobalPosition());
        updated_model_matrix = glm::rotate(updated_model_matrix, 0.01f*time, glm::vec3(0.0, 1.0, 0.0));
        updated_model_matrix = glm::translate(updated_model_matrix, -this->getGlobalPosition());
        this->setModelMatrix(updated_model_matrix);*/
    }
};

class DrawableLine : public Drawable{
private:
    GLuint VBO, VAO, EBO;
    //glm::vec3 color = glm::vec3(1.0f, 1.0f, 1.0f);
    std::vector<CoordVec3> coordinates;
    std::vector<glm::vec3> colors;
    int n_draw;
    std::vector<float> time_ray;
public:
    DrawableLine(){};
    DrawableLine(std::vector<CoordVec3> inp_vertices){
        // No specific colors specified, use white everywhere
        coordinates = inp_vertices;
        for(int j=0; j<int(inp_vertices.size()); j++){
            time_ray.push_back(static_cast<float>(j));
            glm::vec3 color = glm::vec3(1.0f, 1.0f, 1.0f);
            colors.push_back(color);
        }
    }
    DrawableLine(std::vector<CoordVec3> inp_vertices, std::vector<float> inp_time_ray){
        coordinates = inp_vertices;
        time_ray = inp_time_ray;
        colors = std::vector<glm::vec3>(inp_vertices.size(), glm::vec3(1.0f, 1.0f, 1.0f));
    }
    DrawableLine(std::vector<CoordVec3> inp_vertices, std::vector<float> inp_time_ray, std::vector<glm::vec3> inp_colors){
        coordinates = inp_vertices;
        time_ray = inp_time_ray;
        colors=inp_colors;
    }
    std::vector<GLfloat> createVertices(){
        std::vector<GLfloat> Vertices;
        for(int i = 0; i < int(coordinates.size()); i++){
            Vertices.push_back(coordinates[i].x); Vertices.push_back(coordinates[i].z); Vertices.push_back(coordinates[i].y);
            // Set Color
            Vertices.push_back(this->colors[i].x); Vertices.push_back(this->colors[i].y); Vertices.push_back(this->colors[i].z);
        }
        return Vertices;
    }

    std::vector<GLuint> createIndices(){
        std::vector<GLuint> indices;
        for(int i=0; i < int(coordinates.size()) - 1; i++){
            indices.push_back(i); indices.push_back(i+1);
        }
        return indices;
    }
    
    void draw() override {
        std::vector<GLfloat> vertices = this->createVertices();
        std::vector<GLuint> indices = this->createIndices();

        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        glGenBuffers(1, &EBO);

        glBindVertexArray(VAO);

        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(GLfloat), vertices.data(), GL_STATIC_DRAW);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), indices.data(), GL_STATIC_DRAW);

        // Set Position Data
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (GLvoid*)0);
        glEnableVertexAttribArray(0);

        // Set Color Data
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(GLfloat), (void*)(3*sizeof(GLfloat)));
        glEnableVertexAttribArray(1);

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);

        //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // Set rendering mode to wireframe

        glBindVertexArray(VAO);
        //glDrawElements(GL_LINES, (GLsizei)(3*vertices.size()), GL_UNSIGNED_INT, 0);
        glDrawElements(GL_LINES, (GLsizei)(3*this->n_draw), GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // Set rendering mode to wireframe
    }
    void update(float time) override {
        /* This function would take "Photon Time" for the rendering of the photon lines
        auto lower = std::lower_bound(time_ray.begin(), time_ray.end(), (time-5)*2);
        auto index = std::distance(time_ray.begin(), lower);
        this->n_draw = index;
        */
        //double size = static_cast<double>(coordinates.size());
        this->n_draw = std::max(static_cast<int>(coordinates.size()/20.0 * fmod((time - 5), 20.0)), 1);
        //std::cout << fmod((time - 5), 20) << std::endl;
        //std::cout << this->n_draw << std::endl;
    }
};

#endif