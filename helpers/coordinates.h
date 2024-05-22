#ifndef _coordinates_h_
#define _coordinates_h_


#include <iostream>
#include <vector>
#include <cmath>

class CoordVec4{
private:
    int length = 4;
public:
    bool cartesian = true;
    bool bl = false;
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
    CoordVec3(float xIn, float yIn, float zIn, bool cartIn, bool blIn):x(xIn), y(yIn), z(zIn), cartesian(cartIn), bl(blIn){};
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
        }
        else if (true == this->cartesian && false == this->bl){
            x_c = this->x; y_c = this->y; z_c = this->z;
        }
        else {
            std::cerr << "Coordinate Vector is either both cartesian and bl or neither cartesian nor bl! Cannot convert" << std::endl;
            exit(1);
        }
        coords = CoordVec3(x_c,y_c,z_c, true, false);
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
        }
        else if (true == this->bl && false == this->cartesian){
            radius = this->x; theta = this->y; phi = this->z;
            this->cartesian = true; this->bl=false;
        }
        else {
            std::cerr << "Coordinate Vector is either both cartesian and bl or neither cartesian nor bl! Cannot convert" << std::endl;
            exit(1);
        }
        coords = CoordVec3(radius, theta, phi, false, true);
        return coords;
    }
    float cart_norm(float a_spin = 0.0){
        float norm = -1;
        if (this->cartesian){
            norm = std::sqrt(this->x*this->x + this->y*this->y + this->z*this->z);
        }
        else if (!this->cartesian && this->bl)
        {
            CoordVec3 cartesian = this->to_cart(a_spin);
            norm = std::sqrt(cartesian.x*cartesian.x + cartesian.y*cartesian.y + cartesian.z*cartesian.z);
        }
        else
        {
            std::cerr << "Coordinate was neither bl nor cartesian or both at once" << std::endl;
            norm = -1;
        }
        return norm;
    }
};

class MomentumVector : public CoordVec4{
private:
    float pr = this->x;
    float pth = this->y;
    float pph = this->z;
    float pt = this->t;
    
public:
    bool _check_momentum(CoordVec4 inp_momentum){
        float length = pow(inp_momentum.x,2)+pow(inp_momentum.y,2)+pow(inp_momentum.z,2) - pow(inp_momentum.t,2);
        return (length == 0.0);
    }
    MomentumVector(float inp_pr, float inp_pth, float inp_pph, float inp_pt){
        bool valid = this->_check_momentum(CoordVec4(inp_pr, inp_pth, inp_pph, inp_pt));
        if (!valid) {
            std::cerr << "The input momenta do not belong to a photon" << std::endl;
            std::cerr << "Adjusting Energy to meet photon requirements" << std::endl;
        }
        this->x = inp_pr; this->pr = inp_pr;
        this->y = inp_pth; this->pth = inp_pth;
        this->z = inp_pph; this->pph = inp_pph;
        this->t = std::sqrt(pow(inp_pr,2) + pow(inp_pth, 2) + pow(inp_pph,2));
    }
    MomentumVector(float inp_pr, float inp_pth, float inp_pph){
        this->x = inp_pr; this->pr = inp_pr;
        this->y = inp_pth; this->pth = inp_pth;
        this->z = inp_pph; this->pph = inp_pph;
        this->t = std::sqrt(pow(inp_pr,2) + pow(inp_pth, 2) + pow(inp_pph,2));
    }
};

class SourcePosition : public CoordVec3{
    // Another "Alias" Class. This encapsulates the relevant data of the source position for many of the YNOGK functions.
    // As Kerr spacetime is time invariant and polar rotation invariant, t and phi coordinates are not necessary for many calculations
    // Contains: radius, sinobs, muobs | (radius, sin(theta), cos(theta)) where radius and theta are taken as the BL coordinates of the source
private:

public:
    float radius;
    float sinobs; float muobs;
    SourcePosition(CoordVec3 inp){
        if(inp.bl && !inp.cartesian){
        this->x = inp.x; this->y = inp.y; this->z = inp.z;
        this->radius = inp.x; this->sinobs = std::sin(inp.y); this->muobs = std::cos(inp.z);
        }
        else if (inp.cartesian && !inp.bl){
            std::cerr << "CoordVec is cartesian! Need BL or additional a_spin information!" << std::endl;
            this->x = inp.x; this->y = inp.y; this->z = inp.z; 
        }
        else{
            std::cerr << "CoordVec neither cartesian nor BL or both simultaneously!" << std::endl;
        }
    }
    SourcePosition(CoordVec4 inp) : SourcePosition(CoordVec3(inp.x, inp.y, inp.z)) {}
};

#endif