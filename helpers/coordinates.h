#ifndef _coordinates_h_
#define _coordinates_h_


#include <iostream>
#include <vector>

class CoordVec4{
private:
    int length = 4;
    bool cartesian = true;
    bool bl = false;
public:
    float x,y,z,t;
    CoordVec4(float xIn, float yIn, float zIn, float tIn):x(xIn), y(yIn), z(zIn), t(tIn){};
    CoordVec4(){};
    void display(){
        std::cout << "[" << x << ", " << y << ", " << z << "]" << std::endl;
    }

    CoordVec4 operator+(const CoordVec4& other) const {
        return CoordVec4(x+other.x, y+other.y, z+other.z, t+other.t);
    }

    CoordVec4& operator+=(const CoordVec4& other){
        this->x+= other.x;
        this->y+= other.y;
        this->z+= other.z;
        this->t+= other.t;
        return *this;
    }

    CoordVec4 operator*(const float scalar) const {
        return CoordVec4(x*scalar, y*scalar, z*scalar, t*scalar);
    }
    CoordVec4 operator*(const double scalar) const {
            return CoordVec4(x*scalar, y*scalar, z*scalar, t*scalar);
    }
    bool is_bl(){
        return bl;
    }
    bool is_cart(){
        return cartesian;
    }

    CoordVec4 to_cart(float a_spin){
        CoordVec4 coords;
        float x_c, y_c, z_c, t_c;
        if (false == this->cartesian && true == this->bl){
            // Converts a blcoordinate vector to a cartesian coordinate vector (r, theta, phi, t) -> (x, y, z, t);
            x_c = std::sqrt(this->x * this->x + a_spin*a_spin) * std::sin(this->y)*std::cos(this->z);
            y_c = std::sqrt(this->x * this->x + a_spin*a_spin) * std::sin(this->y)*std::sin(this->z);
            z_c = this->x * std::cos(this->y);
            t_c = this->t;
            // update the cordinate kind
            this->cartesian = true; this->bl = false;
        }
        else if (true == this->cartesian && false == this->bl){
            x_c = this->x; y_c = this->y; z_c = this->z; t_c = this->t;
            this->cartesian = true; this->bl=false;
        }
        else {
            std::cerr << "Coordinate Vector is either both cartesian and bl or neither cartesian nor bl! Cannot convert" << std::endl;
            exit(1);
        }
        coords = CoordVec4(x_c,y_c,z_c,t_c);
        return coords;
    }

    CoordVec4 to_bl(float a_spin){
        CoordVec4 coords;
        float radius, theta, phi, t_b;
        if (false == this->bl && true == this->cartesian){
            // Converts a cartesian vector to a bl coordinate vector (x, y, z, t) -> (r, theta, phi, t);
            float x = this->x, y = this->y, z = this->z;
            float R1 = std::sqrt(pow(z,4) + 2*pow(y,2)*pow(z,2) + 2*pow(x,2)*pow(z,2) + 2*pow(a_spin,2)*pow(z,2) + pow(y,4) + 2*pow(x,2)*pow(y,2) - 2*pow(a_spin,2)*pow(y,2) + 
                        pow(x,4) - 2*pow(a_spin,2)*pow(x,2) + pow(a_spin,4));
            radius = std::sqrt(R1 + pow(x,2)+pow(y,2)+pow(z,2)-pow(a_spin,2))/std::sqrt(2);
            theta = std::acos(z/radius);
            phi = std::atan(y/x);
            t_b = this->t;
            // update the cordinate kind
            this->cartesian = true; this->bl = false;
        }
        else if (true == this->bl && false == this->cartesian){
            radius = this->x; theta = this->y; phi = this->z; t_b = this->t;
            this->cartesian = true; this->bl=false;
        }
        else {
            std::cerr << "Coordinate Vector is either both cartesian and bl or neither cartesian nor bl! Cannot convert" << std::endl;
            exit(1);
        }
        coords = CoordVec4(radius, theta, phi, t_b);
        return coords;
    }
};

class CoordVec3{
private:
    int length = 3;
public:
    bool cartesian = true;
    bool bl = false;
    float x,y,z;
    CoordVec3(float xIn, float yIn, float zIn):x(xIn), y(yIn), z(zIn){};
    CoordVec3(){};
    void display(){
        std::cout << "[" << x << ", " << y << ", " << z << "]" << std::endl;
    }

    CoordVec3 operator+(const CoordVec3& other) const{
        return CoordVec3(x+other.x, y+other.y, z+other.z);
    }

    CoordVec3& operator+=(const CoordVec3& other){
        this->x+= other.x;
        this->y+= other.y;
        this->z+= other.z;
        return *this;
    }

    CoordVec3 operator*(const float scalar) const{
        return CoordVec3(x*scalar, y*scalar, z*scalar);
    }
    CoordVec3 operator*(const double scalar) const{
            return CoordVec3(x*scalar, y*scalar, z*scalar);
    }
    CoordVec3 to_cart(float a_spin){
        CoordVec3 coords;
        float x_c, y_c, z_c;
        if (false == this->cartesian && true == this->bl){
            // Converts a blcoordinate vector to a cartesian coordinate vector (r, theta, phi, t) -> (x, y, z, t);
            x_c = std::sqrt(this->x * this->x + a_spin*a_spin) * std::sin(this->y)*std::cos(this->z);
            y_c = std::sqrt(this->x * this->x + a_spin*a_spin) * std::sin(this->y)*std::sin(this->z);
            z_c = this->x * std::cos(this->y);
            // update the cordinate kind
            this->cartesian = true; this->bl = false;
        }
        else if (true == this->cartesian && false == this->bl){
            x_c = this->x; y_c = this->y; z_c = this->z;
            this->cartesian = true; this->bl=false;
        }
        else {
            std::cerr << "Coordinate Vector is either both cartesian and bl or neither cartesian nor bl! Cannot convert" << std::endl;
            exit(1);
        }
        coords = CoordVec3(x_c,y_c,z_c);
        return coords;
    }

    CoordVec3 to_bl(float a_spin){
        CoordVec3 coords;
        float radius, theta, phi;
        if (false == this->bl && true == this->cartesian){
            // Converts a cartesian vector to a bl coordinate vector (x, y, z, t) -> (r, theta, phi, t);
            float x = this->x, y = this->y, z = this->z;
            float R1 = std::sqrt(pow(z,4) + 2*pow(y,2)*pow(z,2) + 2*pow(x,2)*pow(z,2) + 2*pow(a_spin,2)*pow(z,2) + pow(y,4) + 2*pow(x,2)*pow(y,2) - 2*pow(a_spin,2)*pow(y,2) + 
                        pow(x,4) - 2*pow(a_spin,2)*pow(x,2) + pow(a_spin,4));
            radius = std::sqrt(R1 + pow(x,2)+pow(y,2)+pow(z,2)-pow(a_spin,2))/std::sqrt(2);
            theta = std::acos(z/radius);
            phi = std::atan(y/x);
            // update the cordinate kind
            this->cartesian = true; this->bl = false;
        }
        else if (true == this->bl && false == this->cartesian){
            radius = this->x; theta = this->y; phi = this->z;
            this->cartesian = true; this->bl=false;
        }
        else {
            std::cerr << "Coordinate Vector is either both cartesian and bl or neither cartesian nor bl! Cannot convert" << std::endl;
            exit(1);
        }
        coords = CoordVec3(radius, theta, phi);
        return coords;
    }
};

#endif