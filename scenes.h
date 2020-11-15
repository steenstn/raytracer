#ifndef SCENES
#define SCENES_h
#include <vector>

#include "vector.h"
#include "sphere.h"
#include "mesh.h"
#include "material.h"
#include "plane.h"
#include "vertex.h"
#include "triangle.h"
#include "cube.h"
#include "matrix4.h"

float random(void) {
	return (float)rand() / (float)RAND_MAX;
}

Vector randomMixed(Vector mix) {
	float r = random();
	float g = random();
	float b = random();

	r = (r + mix.x) / 2;
	g = (g + mix.y) / 2;
	b = (b + mix.z) / 2;

	return Vector(r, g, b);
}

Vector mixWhite(Vector v) {
	return Vector((v.x + 1) / 2, (v.y + 1) / 2, (v.z + 1) / 2);
}

struct Scene {
	std::vector<Mesh*> meshes;
	Vector camera;
	float focusLength;
	float depthOfField;
};

Scene lonelyBalls() {
	Scene scene;
	Material mirror = Material(Vector(0.91, 0.91, 0.91));
	mirror.reflective = true;
	scene.meshes.push_back(new Sphere(-8, -5, 0, 8, mirror));
	scene.meshes.push_back(new Sphere(4, 0, 3, 3, Material(Vector(0.12, 0.24, 0.86))));

	scene.meshes.push_back(new Sphere(6.2, 1, 10, 2, Material(Vector(0.2, 0.64, 0.16))));

	scene.meshes.push_back(new Plane(0, 1, 0, Vector(0, -40, 0), Material(Vector(1, 1, 1), Vector(0.8, 0.8, 0.8))));
	scene.meshes.push_back(new Plane(0, -1, 0, Vector(0, 3, 0), Material(Vector(0, 0, 0), Vector(0.8, 0.8, 0.8))));

	scene.camera = Vector(0, -2.0, 35); // Starting point
	scene.focusLength = 35;
	scene.depthOfField = 0.9;
	return scene;
}

Scene fiveBallsOfColor() {
	Scene scene;
	Material mirror = Material(mixWhite(Vector(1, 0.47, 0.47)));
	mirror.reflective = true;

	Material greenMirror = Material(mixWhite(Vector(0.6, 1.47, 0.47)));
	mirror.reflective = true;

	scene.meshes.push_back(new Sphere(30, -22, -40, 25, Material(mixWhite(Vector(0.64, 0.94, 0.94)))));
	scene.meshes.push_back(new Sphere(0, -7, 0, 10, mirror));

	scene.meshes.push_back(new Sphere(-5, 1, 20, 2, Material(mixWhite(Vector(1.5, 1.5, 0.2)))));
	scene.meshes.push_back(new Sphere(3, 2.5, 20, 1, Material(mixWhite(Vector(0.23, 0.6, 2)))));

	scene.meshes.push_back(new Sphere(-18, -4, 2, 7, greenMirror));

	scene.meshes.push_back(new Sphere(-20, -30, 6, 4, Material(Vector(30, 30, 30), Vector(1, 1, 1))));

	scene.meshes.push_back(new Plane(0, 1, 0, Vector(0, -40, 0), Material(Vector(0.9, 0.9, 1), Vector(1, 1, 1))));
	scene.meshes.push_back(new Plane(0, -1, 0, Vector(0, 3, 0), Material(Vector(0, 0, 0), Vector(0.8, 0.8, 0.8))));

	scene.camera = Vector(0, -3, 40);
	scene.depthOfField = 0.9;
	scene.focusLength = 35;
	return scene;
}

#endif