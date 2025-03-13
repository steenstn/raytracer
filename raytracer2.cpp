#include <ctime>
#include <time.h>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
// #include <omp.h>

#include "material.h"
#include "mesh.h"
#include "scenes.h"
#include "sphere.h"
#include "vector.h"

//const int SCREENWIDTH = 800;
//const int SCREENHEIGHT = 600;
const int SCREENWIDTH = 800;
const int SCREENHEIGHT = 600;

const double aaFactor = 1; // Antialias factor

int numFrames = 1;
int numRays = 10; // Number of rays/iteration

const int WIDTH = SCREENWIDTH * aaFactor;
const int HEIGHT = SCREENHEIGHT * aaFactor;
unsigned char img[WIDTH * HEIGHT * 3]; // Slutbilden sparad som en lång jävla sträng

short picture[SCREENWIDTH][SCREENHEIGHT][3];
double pictureAA[WIDTH][HEIGHT][3];

int bounces = 0; // Number of reflected rays
Vector BLACK = Vector();
Vector WHITE = Vector(1.0,1.0,1.0);
Vector ambientColor = WHITE;

int maxBounces = 50; // Maximum number of bounces allowed

std::vector<Mesh *> theMeshes; // All of the meshes

void savebmp(const char *filename, int w, int h);
void superSample(float numPasses);
Vector shootRay(Vector s, Vector d, int index);
Vector shootRefractedRay(Vector s, Vector d, int index, float n1);

Vector shootRefractedRay2(Vector s, Vector d, int index, float n1);
float step(float a, float x) { return (float)(x >= a); }
float clamp(float x, float a, float b) { // Clamp x between a and b(THE CLAMPS!)
  return (x < a ? a : (x > b ? b : x));
}

struct oSphere {
  Vector position;
  float radius;
};

struct oColor {
  float r;
  float g;
  float b;
};

struct ObjectHit {
  int index;
  Vector position;
};

const int num_objects = 2;
oSphere *all_objects = new oSphere[num_objects];
oColor *all_colors = new oColor[num_objects];

struct oSpheres {
  float radius[num_objects];
  Vector position[num_objects];
};

struct oColors {
  float r[num_objects];
  float g[num_objects];
  float b[num_objects];
};

struct ObjectHit ray_sphere_intersection2(oSphere* all_spheres, int num_spheres, Vector start, Vector direction);
Vector shoot_ray2(oSphere *spheres, int num_spheres, Vector start, Vector direction);

Vector get_normal(Vector position, struct oSphere *sphere);

int main(void) {

  //  srand((unsigned)time(0));
  Scene scene;
  scene.meshes.push_back(
    new Sphere(0.0, 0.0, 0.0, 1.0, Material(Vector(0.2, 0.3, 1.0))));

  all_objects[0].position = Vector(0,5,-10);
  all_objects[0].radius = 10;

  all_colors[0].r = 0.2;
  all_colors[0].g = 0.3;
  all_colors[0].b = 1.0;

  all_objects[1].position = Vector(2,0,0);
  all_objects[1].radius = 1;

  all_colors[1].r = 0.6;
  all_colors[1].g = 0.5;
  all_colors[1].b = 0.1;

  Vector da_position = Vector(0.0, 0.0, 0.0);

  scene.camera.position = Vector(0.0, -0.5, 10.0);
  scene.camera.focusLength = 30;
  scene.camera.depthOfField = 0.0;
  scene.skyColor = Vector(1.0, 1.0, 1.0);

  theMeshes = scene.meshes;
  Vector s = scene.camera.position;
  float DoF = scene.camera.depthOfField;

  std::cout << "Number of meshes: " << theMeshes.size() << std::endl << std::endl;
  std::cout << "Number of rays/pixel: " << numRays << std::endl;
  // Vector s(0,0,15); // Starting point
  Vector s2;
  float focusLength = scene.camera.focusLength;
  int totalRays = 0;

  int numberPasses = 0;
  for (;;) {
    // cout << "Samples: " << totalRays << endl;
    // s=s+Vector(0,0,-0.1);

    float xmax = 5, ymax = 5;
    clock_t t;
    t = clock();
    // #pragma omp parallel for
    for (int screenY = 0; screenY < HEIGHT; screenY++) {
      for (int screenX = 0; screenX < WIDTH; screenX++) {
        Vector endColor, dir;
        float x, y;
        x = (float)(screenX * 6) / (float)WIDTH - 3.0;
        y = (float)(screenY * 6) * (float)HEIGHT / (float)WIDTH /
          (float)HEIGHT -
          3.0 * (float)HEIGHT / (float)WIDTH;

        dir = Vector(x / xmax, y / ymax, -1); // Direction
        dir.normalize();
        for (int i = 0; i < numRays; i++) {
          bounces = 0;
          s2 = s + Vector(2.0 * makeRandom() - 1, 2.0 * makeRandom() - 1,
                          2.0 * makeRandom() - 1) * DoF;
          Vector dir2;
          Vector position2 = s + dir * focusLength;

          dir2 = position2 - s2;
          dir2.normalize();

          endColor = endColor + shoot_ray2(all_objects, num_objects, s2, dir2); // Fire it up
        }
        endColor = endColor / numRays;
        /*if(endColor.x>1)
            endColor.x=1;
        if(endColor.y>1)
            endColor.y=1;
        if(endColor.z>1)
            endColor.z=1;*/

        pictureAA[screenX][screenY][0] += endColor.x;
        pictureAA[screenX][screenY][1] += endColor.y;
        pictureAA[screenX][screenY][2] += endColor.z;
      }
      if ((int)screenY % 20 == 0)
        std::cout << screenY << " / " << HEIGHT << std::endl;
    }
    t = clock() -t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC;
    std::cout << "Time: " << time_taken;
    totalRays += numRays;
    std::cout << "Total rays: " << totalRays << std::endl;
    numberPasses++;
    superSample(numberPasses);
    // stringstream ss;
    // ss<<totalRays;
    // string out="bilder/aa"+ss.str()+".bmp";
    std::string out = "result.bmp";
    char outchar[20];
    for (int a = 0; a < out.length(); a++) {
      outchar[a] = out[a];
    }
    outchar[out.length()] = '\0';

    savebmp(outchar, SCREENWIDTH, SCREENHEIGHT);
  }
}

Vector shoot_ray2(oSphere* all_spheres, int num_spheres, Vector start, Vector direction) {
  int max_bounces = 50; 

  Vector resulting_color = Vector();
  
  float fraction = 1.0f;
  for(int num_bounces = 0; num_bounces < max_bounces; num_bounces++) {

    struct ObjectHit hit =
      ray_sphere_intersection2(all_spheres, num_spheres, start, direction);

    if (hit.index == -1) {
      break;
    }

    Vector new_direction = Vector(2 * makeRandom() - 1, 2 * makeRandom() -1, 2*makeRandom() -1);

    struct oSphere hit_sphere = all_objects[hit.index];
   
    Vector hit_normal =  hit.position - hit_sphere.position;
    hit_normal.normalize();
    new_direction = new_direction.cross(hit_normal);
    new_direction.normalize();

    float eps1 = makeRandom() * 6.28318; // 2*PI
    float eps2 = sqrtf(makeRandom());

    float x = cosf(eps1) * eps2;
    float y = sinf(eps1) * eps2;
    float z = sqrtf(1.0f - eps2 * eps2);
    Vector ssx = new_direction * x + hit_normal.cross(new_direction) * y + hit_normal * z;
    ssx.normalize();
    direction = ssx;
    start = hit.position;

    Vector this_color = Vector(all_colors[hit.index].r, all_colors[hit.index].g, all_colors[hit.index].b);

    Vector emittance = hit.index == 1 ? Vector(4,4,4) : Vector(0,0,0);
    resulting_color = emittance + resulting_color * this_color; 
  }
  return resulting_color;
}


Vector shootRay(Vector s, Vector d, int index) {
  if (bounces > maxBounces) {
    return ambientColor;
  }
  bounces++;
  d.normalize();
  Vector n, c, v, endMovement, pos;
  float wee, spec, dot = 0;
  float t[2];
  Vector endColor;
  Vector refractionDir = d;
  Vector refractionPos;
  float distances = 9999999;
  bool hit = false;

  for (int i = 0; i < theMeshes.size(); i++) // Find closest intersection
  {
    if (i != index)
      if (theMeshes.at(i)->checkIntersection(s, d, distances, pos) == true) {
        index = i;
        hit = true;
      }
  }

  if (hit) {
    Material m = theMeshes.at(index)->material;
    if (m.emittance.x > 0 || m.emittance.y > 0 || m.emittance.z > 0)
      return m.emittance;

    n = theMeshes.at(index)->getNormal(pos);
    if (index == 0) {
      //       m.reflectance=Vector(0.75,0.75,0.75)*step(0,sin(4*pos.x)*cos(4*pos.z))+Vector(0.05,0.05,0.05);
      //       // Checkers floor
    }
    if (m.reflectance.sum() > ((float)rand() / (float)RAND_MAX) * 3.0 ||
      m.transparent == true || m.reflective == true) {
      if (!m.reflective && !m.transparent)
        m.reflectance = m.reflectance * (3.0 / m.reflectance.sum());

      Vector emittance = m.emittance;
      Vector reflected, refracted, explicitLight, empty;

      //   // pick a makeRandom direction from here and keep going
      Vector newDir = Vector(2 * makeRandom() - 1, 2 * makeRandom() - 1,
                             2 * makeRandom() - 1);

      newDir = newDir.cross(theMeshes.at(index)->getNormal(pos));
      newDir.normalize();

      float eps1 = makeRandom() * 3.14159 * 2.0f;
      float eps2 = sqrtf(makeRandom());

      float x = cosf(eps1) * eps2;
      float y = sinf(eps1) * eps2;
      float z = sqrtf(1.0f - eps2 * eps2);
      Vector tempnormal = theMeshes.at(index)->getNormal(pos);
      Vector ssx = newDir * x + tempnormal.cross(newDir) * y + tempnormal * z;
      ssx.normalize();

      if (m.transparent == true) {
        refracted = refracted + shootRefractedRay(pos, d, index, 1);

      } else if (m.reflective == true && makeRandom() > 0.6) {
        reflected =
          reflected + shootRay(pos, d - n * 2 * (n.dot(d)), -1) * m.spec;
      } else {
        reflected = reflected + shootRay(pos, ssx, -1);
      }
      return (reflected + refracted) * m.reflectance;
    } else {
      return ambientColor;
    }
  }
  //   float nohitdot=d.dot(Vector(0,-1,0));
  //     nohitdot=sqrt(clamp(nohitdot,0,1));
  return ambientColor;
}

ObjectHit ray_sphere_intersection2(oSphere* spheres, int num_spheres, Vector start, Vector direction) {

  float shortest_distance = 999999;
  bool hit = false;
  int hit_index = -1;

  for (int i = 0; i < num_spheres; i++) {
    Vector c = spheres[i].position;
    float radius = spheres[i].radius;
    Vector v = (start - c);
    float v_dot_direction = v.dot(direction);

    float wee = v_dot_direction * v_dot_direction -
      (v.x * v.x + v.y * v.y + v.z * v.z - radius * radius);
    if(wee <= 0.0) {
      continue;
    }

    float dot_product = v_dot_direction * -1.0;
    float wee_sqrt = sqrt(wee);
    float intersection1 = dot_product + wee_sqrt;
    float intersection2 = dot_product - wee_sqrt;

    Vector end_position;

    if (intersection1 < intersection2 && intersection1 > 0.000001) {
      end_position = direction * intersection1;
      hit = true;
    } else if (intersection2 < intersection1 && intersection2 > 0.000001) {
      end_position = direction * intersection2;
      hit = true;
    } else {
      continue;
    }

    float length = end_position.length();
    if (hit == true && length < shortest_distance) {
      hit_index = i;
      shortest_distance = length;
    }
  }

  if (hit) {
    return {hit_index, start + direction * shortest_distance};
  } else {
    return {-1, Vector()};
  }
}


Vector get_normal(Vector position, oSphere *sphere) {
  Vector normal = (position-sphere->position);// / (position - sphere->position).length();
  normal.normalize();
  return normal;
}

Vector shootRefractedRay(Vector s, Vector d, int index, float n1) {
  float n2 = 1.5;

  float n = n1 / n2;
  Vector normal = theMeshes.at(index)->getNormal(s);
  float cosI = normal.dot(d);

  float sinT2 = n * n * (1.0 - cosI * cosI);

  Vector refracted = d * n + normal * (n * cosI - sqrt(1 - sinT2));

  refracted.normalize();
  Vector newStart = s + refracted * 0.0001;

  float distances = 9999999;
  bool hit = false;
  int index2 = 0;
  Vector pos = newStart;
  for (int i = 0; i < theMeshes.size(); i++) // Find closest intersection
  {
    if (theMeshes.at(i)->checkIntersection(newStart, refracted, distances,
                                           pos) == true) {
      index2 = i;
      hit = true;
    }
  }

  if (!hit) {
    return Vector();
  }

  Vector normalOut = theMeshes.at(index2)->getNormal(pos) * -1;

  float cosOut = normalOut.dot(refracted);
  float sinT2Out = n * n * (1.0 - cosOut * cosOut);
  if (sinT2Out > 1) {
    return Vector();
  }
  Vector refractedOut =
    refracted * n + normalOut * (n * cosOut - sqrt(1.0 - sinT2Out));
  refractedOut.normalize();

  return shootRay(pos + refracted * 0.0001, refractedOut, -1);
}

void superSample(float numPasses) {
  double endR, endG, endB;
  int endX = 0, endY = 0;
  for (int y = 0; y < HEIGHT; y += aaFactor) {

    for (int x = 0; x < WIDTH; x += aaFactor) {
      endR = endG = endB = 0;

      for (int a = 0; a < aaFactor; a++) // x-movement
      {
        for (int b = 0; b < aaFactor; b++) // y-movement
        {
          double red = pictureAA[x + a][y + b][0];
          double green = pictureAA[x + a][y + b][1];
          double blue = pictureAA[x + a][y + b][2];

          endR += red * 255 / numPasses;
          endG += green * 255 / numPasses;
          endB += blue * 255 / numPasses;
        }
      }

      endR /= (aaFactor * aaFactor);
      endG /= (aaFactor * aaFactor);
      endB /= (aaFactor * aaFactor);
      endR = clamp(endR, 0, 255);
      endG = clamp(endG, 0, 255);
      endB = clamp(endB, 0, 255);
      picture[endX][endY][0] = endR;
      picture[endX][endY][1] = endG;
      picture[endX][endY][2] = endB;

      endX++;
    }
    endX = 0;
    endY++;
  }
}

void savebmp(const char *filename, int w, int h) {
  /*std::ofstream theFile;
  theFile.open(filename);
  // Write header
  theFile << "P6 " << w << " " << h << " " << 255 << "\n";

  for(int i=0;i<SCREENHEIGHT;i++)
{
  for(int j=0;j<SCREENWIDTH;j++)
  {
      theFile << (unsigned char)picture[j][i][2] << " ";
      theFile << (unsigned char)picture[j][i][1] << " ";
      theFile << (unsigned char)picture[j][i][0] << " ";

  }
          theFile << "\n";
}
  theFile.close();*/

  int i;
  FILE *f;
  int filesize = 54 + 3 * w * h;

  unsigned char bmpfileheader[14] = {'B', 'M', 0, 0,  0, 0, 0,
    0,   0,   0, 54, 0, 0, 0};
  unsigned char bmpinfoheader[40] = {40, 0, 0, 0, 0, 0, 0,  0,
    0,  0, 0, 0, 1, 0, 24, 0};
  unsigned char bmppad[3] = {0, 0, 0};

  bmpfileheader[2] = (unsigned char)(filesize);
  bmpfileheader[3] = (unsigned char)(filesize >> 8);
  bmpfileheader[4] = (unsigned char)(filesize >> 16);
  bmpfileheader[5] = (unsigned char)(filesize >> 24);

  bmpinfoheader[4] = (unsigned char)(w);
  bmpinfoheader[5] = (unsigned char)(w >> 8);
  bmpinfoheader[6] = (unsigned char)(w >> 16);
  bmpinfoheader[7] = (unsigned char)(w >> 24);
  bmpinfoheader[8] = (unsigned char)(h);
  bmpinfoheader[9] = (unsigned char)(h >> 8);
  bmpinfoheader[10] = (unsigned char)(h >> 16);
  bmpinfoheader[11] = (unsigned char)(h >> 24);

  // f = fopen("img.raw","wb");
  // fwrite(img,3,w*h,f);
  // fclose(f);

  // Make the int-array into a char-array for the bmp
  int counter = 0;
  for (int i = 0; i < SCREENHEIGHT; i++) {
    for (int j = 0; j < SCREENWIDTH; j++) {
      img[counter] = (unsigned char)picture[j][i][2];
      img[counter + 1] = (unsigned char)picture[j][i][1];
      img[counter + 2] = (unsigned char)picture[j][i][0];
      counter = counter + 3;
    }
  }

  f = fopen(filename, "wb+");
  fwrite(bmpfileheader, 1, 14, f);
  fwrite(bmpinfoheader, 1, 40, f);

  for (i = 0; i < h; i++) {
    fwrite(img + (w * (h - i - 1) * 3), 3, w, f);
    fwrite(bmppad, 1, (4 - (w * 3) % 4) % 4, f);
  }
  fclose(f);
}
