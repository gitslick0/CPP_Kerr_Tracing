#include"./helpers/Renderer2.h"

int main() {
    const int WIDTH = 800;
    const int HEIGHT = 600;
    const char* TITLE = "OpenGL Window";

    Renderer renderer(WIDTH, HEIGHT, TITLE);

    DrawableSphere Sphere = DrawableSphere(3.0);
    DrawableDisk Disk = DrawableDisk(5.0, 20.0, 0.0);
    DrawableSphere Source = DrawableSphere(0.2);

    std::vector<CoordVec3> points = {CoordVec3(5.0, 3.0, 3.0), CoordVec3(4.8, 3.8, 3.8), CoordVec3(4.6, 4.3, 4.3)};
    DrawableLine Line = DrawableLine(points);
    points[0].display();

    renderer.addDrawable(&Sphere);
    renderer.addDrawable(&Disk);
    renderer.addDrawable(&Source);
    renderer.addDrawable(&Line);

    Disk.setGlobalPosition(glm::vec3(0.0, 0.0, 0.0));
    Source.setColor(glm::vec3(1.0, 1.0, 0.6));
    Source.setGlobalPosition(glm::vec3(4.0, 3.0, 0.0));
    Sphere.setColor(glm::vec3(0.0, 0.0, 0.0));

    // Example: Render the scene
    renderer.render();

    return 0;
}