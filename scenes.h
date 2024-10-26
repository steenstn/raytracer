#ifndef SCENES
#define SCENES_h
#include <vector>

#include "cube.h"
#include "material.h"
#include "matrix4.h"
#include "mesh.h"
#include "plane.h"
#include "sphere.h"
#include "triangle.h"
#include "vector.h"
#include "vertex.h"

float makeRandom(void) { return (float)rand() / (float)RAND_MAX; }

Vector randomMixed(Vector mix) {
  float r = makeRandom();
  float g = makeRandom();
  float b = makeRandom();

  r = (r + mix.x) / 2;
  g = (g + mix.y) / 2;
  b = (b + mix.z) / 2;

  return Vector(r, g, b);
}

Vector mixWhite(Vector v) {
  return Vector((v.x + 1) / 2, (v.y + 1) / 2, (v.z + 1) / 2);
}

struct Camera {
  Vector position;
  float focusLength;
  float depthOfField;
};
struct Scene {
  std::vector<Mesh *> meshes;
  Camera camera;
  Vector skyColor;
};

Scene lonelyBalls() {
  Scene scene;
  Material mirror = Material(Vector(0.91, 0.91, 0.91));
  mirror.reflective = true;
  scene.meshes.push_back(new Sphere(-8, -5, 0, 8, mirror));
  scene.meshes.push_back(
      new Sphere(4, 0, 3, 3, Material(Vector(0.12, 0.24, 0.86))));

  scene.meshes.push_back(
      new Sphere(6.2, 1, 10, 2, Material(Vector(0.2, 0.64, 0.16))));

  scene.meshes.push_back(
      new Plane(0, 1, 0, Vector(0, -40, 0),
                Material(Vector(1, 1, 1), Vector(0.8, 0.8, 0.8))));
  scene.meshes.push_back(
      new Plane(0, -1, 0, Vector(0, 3, 0),
                Material(Vector(0, 0, 0), Vector(0.8, 0.8, 0.8))));

  scene.camera.position = Vector(0, -2.0, 35); // Starting point
  scene.camera.focusLength = 35;
  scene.camera.depthOfField = 0.9;
  scene.skyColor = Vector(1.0, 1.0, 1.0);
  return scene;
}

Scene fiveBallsOfColor() {
  Scene scene;
  Material mirror = Material(mixWhite(Vector(1, 0.47, 0.47)));
  mirror.reflective = true;

  Material greenMirror = Material(mixWhite(Vector(0.6, 1.47, 0.47)));
  mirror.reflective = true;

  scene.meshes.push_back(new Sphere(
      30, -22, -40, 25, Material(mixWhite(Vector(0.64, 0.94, 0.94)))));
  scene.meshes.push_back(new Sphere(0, -7, 0, 10, mirror));

  scene.meshes.push_back(
      new Sphere(-5, 1, 20, 2, Material(mixWhite(Vector(1.5, 1.5, 0.2)))));
  scene.meshes.push_back(
      new Sphere(3, 2.5, 20, 1, Material(mixWhite(Vector(0.23, 0.6, 2)))));

  scene.meshes.push_back(new Sphere(-18, -4, 2, 7, greenMirror));

  scene.meshes.push_back(new Sphere(
      -20, -30, 6, 4, Material(Vector(30, 30, 30), Vector(1, 1, 1))));

  scene.meshes.push_back(
      new Plane(0, 1, 0, Vector(0, -40, 0),
                Material(Vector(0.9, 0.9, 1), Vector(1, 1, 1))));
  scene.meshes.push_back(
      new Plane(0, -1, 0, Vector(0, 3, 0),
                Material(Vector(0, 0, 0), Vector(0.8, 0.8, 0.8))));

  scene.camera.position = Vector(0, -3, 40);
  scene.camera.depthOfField = 0.9;
  scene.camera.focusLength = 35;
  scene.skyColor = Vector(1.0, 1.0, 1.0);
  return scene;
}

Scene threeBallsWithGlass() {
  Scene scene;
  scene.meshes.push_back(
      new Sphere(-3.5, -5.0, 2.0, 2.0,
                 Material(Vector(5.0, 5.0, 5.0), Vector(1.0, 1.0, 1.0))));
  Material yellow = Material(Vector(1.0, 0.6, 0.1));
  yellow.reflective = true;
  scene.meshes.push_back(new Sphere(-2.0, -1.0, -3.0, 2.0, yellow));
  Material glass = Material(mixWhite(Vector(0.2, 0.5, 1.0)));
  glass.transparent = true;
  scene.meshes.push_back(new Sphere(1.0, 0.0, 0.8, 1.0, glass));
  scene.meshes.push_back(new Sphere(3.0, -2.0, -3.0, 3.0,
                                    Material(mixWhite(Vector(0.8, 0.2, 0.2)))));
  scene.meshes.push_back(
      new Plane(0.0, 1.0, 0.0, Vector(0.0, -1000.0, 0.0),
                Material(Vector(0.9, 0.9, 1.0), Vector(0.9, 0.9, 1.0))));
  scene.meshes.push_back(new Plane(0.0, -1.0, 0.0, Vector(0.0, 1.0, 0.0),
                                   Material(mixWhite(Vector(0.2, 0.3, 0.2)))));
  scene.camera.depthOfField = 0.1;
  scene.camera.focusLength = 7;
  scene.camera.position = Vector(0.0, -0.5, 7.0);
  return scene;
}

Scene glassBalls() {
  Scene scene;
  scene.meshes.push_back(
      new Plane(0, -1, 0, Vector(0, 1, 0), Material(Vector(0.8, 0.8, 0.8))));
  scene.meshes.push_back(
      new Plane(0.0, 1.0, 0.0, Vector(0.0, -8.0, 0.0),
                Material(Vector(2, 2, 2), Vector(0.9, 0.9, 1.0))));
  scene.meshes.push_back(new Plane(1.0, 0.0, 0.0, Vector(-6.0, 0.0, 0.0),
                                   Material(Vector(0.2, 0.8, 1.0))));
  scene.meshes.push_back(new Plane(-1.0, 0.0, 0.0, Vector(6.0, 0.0, 0.0),
                                   Material(Vector(0.8, 0.2, 1.0))));
  scene.meshes.push_back(new Plane(0.0, 0.0, 1.0, Vector(0.0, 0.0, -8.0),
                                   Material(Vector(0.8, 0.8, 0.8))));
  scene.meshes.push_back(new Plane(0.0, 0.0, -1.0, Vector(0.0, 0.0, 15.0),
                                   Material(Vector(0.9, 0.9, 0.9))));

  // scene.meshes.push_back(new Sphere(0, 0, 0, 1, Material(Vector(0.8, 0.2,
  // 0.25))));

  // Glass balls
  Material glass = Material(Vector(0.7, 0.7, 0.7));
  glass.transparent = true;
  glass.reflective = true;
  Material mirror = Material(Vector(1, 1, 1));
  mirror.reflective = true;

  scene.meshes.push_back(new Sphere(0, 0, 0, 1, glass));
  scene.meshes.push_back(new Sphere(-2.2, 0, 0, 1, glass));
  scene.meshes.push_back(new Sphere(2.2, 0, 0, 1, glass));

  // scene.meshes.push_back(new Sphere(-1, -3, 0, 1, Material(Vector(30, 30,
  // 30), Vector(1, 1, 1))));
  scene.camera.depthOfField = 0.0;
  scene.camera.focusLength = 7;
  scene.camera.position = Vector(0.0, -2, 12.0);
  return scene;
}

Scene woo() {
  Scene scene;
  scene.meshes.push_back(
      new Sphere(4.0, 0.0, 0.0, 4.0, Material(Vector(0.1, 0.31, 0.91))));
  scene.meshes.push_back(
      new Sphere(-4.0, 0.0, 0.0, 4.0, Material(Vector(0.1, 0.51, 0.21))));

  scene.camera.position = Vector(0.0, 0.0, 30.0);

  scene.camera.focusLength = 30;
  scene.camera.depthOfField = 0.3;
  return scene;
}

Scene molecule() {
  Scene scene;
  scene.meshes.push_back(
      new Sphere(0.0, 0.0, 0.0, 1.0, Material(Vector(0.2, 0.3, 1.0))));

  Vector position = Vector(0.0, 0.0, 0.0);

  for (int i = 0; i < 400; i++) {
    Vector new_direction =
        Vector(makeRandom() - makeRandom(), makeRandom() - makeRandom(),
               makeRandom() - makeRandom());
    new_direction.normalize();
    position = position + new_direction;
    scene.meshes.push_back(new Sphere(position.x, position.y, position.z, 1.0,
                                      makeRandom() > 0.2
                                          ? Material(Vector(0.2, 0.3, 1.0))
                                          : Material(Vector(0.5, 0.14, 0.43))));
  }
  scene.camera.position = Vector(0.0, -0.5, 30.0);
  scene.camera.focusLength = 30;
  scene.camera.depthOfField = 0.3;
  scene.skyColor = Vector(1.0, 1.0, 1.0);
  return scene;
}

#endif
