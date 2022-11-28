#include <iostream>
#include "../tga/tgaimage.h"
#include "../model/geometry.h"
#include "../model/model.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
const int width  = 800;
const int height = 800;

/*Lesson 1: Bresenham’s Line Drawing Algorithm*/

// Attempt 1
/*Simplest attempt: define two points in the coord. plane and incrementally 
    set pixel value at small intervals moving from one point to the next
    */
void line1(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color)
{
    for (float t = 0.; t < 1.; t+=0.01)
    {
        int x = x0 + (x1 - x0)*t;
        int y = y0 + (y1 - y0)*t;
        image.set(x, y, color);
    }
}
/*Problem: choice of constant `t` is arbitrarily chosen*/

// Attempt 2
/*We can choose t based on the number of pixels to be drawn. */
void line2(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color)
{
    for (int x = x0; x <= x1; x++)
    {
        float t = (x-x0) / (float)(x1-x0);
        int y = y0*(1.-t) + y1*t;
        image.set(x, y, color);
    }
}
/*Problem: one line has holes since its height is greater than its width.
    also, the result of a line depends on the order of points provided, when it shouldn't*/

// Final Attempt + Optimized
/*The error variable gives us the distance to the best straight line from our current (x, y) pixel. 
    Each time error is greater than one pixel, we increase (or decrease) y by one,
     and decrease the error by one as well.*/
/*We can also eliminate floating points division by modifying our "error" variable to be 
    equivalent to error * dx * 2.*/

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color)
{
    bool steep = false;
    if (std::abs(x0 - x1) < std::abs(y0 - y1)) //if line is steep, transpose image
    {
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    if (x0 > x1) // if first x is larger, flip image to make left-to-right
    {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }
    int dx = x1-x0;
    int dy = y1-y0;
    float derror2 = std::abs(dy)*2;
    float error2 = 0;
    
    int y = y0;
    for (int x = x0; x <= x1; x++)
    {
        if (steep) {
            image.set(y, x, color); //untranspose if line was steep
        } else {
            image.set(x, y, color);
        }
        error2 += derror2;
        if (error2 > dx) {
            y += (y1 > y0 ? 1 : -1); //increment or decrement depending on y values
            error2 -= dx*2;
        }
    }
}


/* Basic Wireframe Rendering */
/* We can use a basic file format to define vertices and render a wireframe model. 
    The .OBJ file is a geometry definition file format.
    The "Model" class reads a Wavefront OBJ file and stores vertices and faces.
    The "geometry.h" header contains basic vector types and operations for storing and 
    manipulating vectors.
    */
void wireframe(Model &model, TGAImage &image, TGAColor color)
{
    for (int i=0; i < model.nfaces(); i++)
    {
        std::vector<int> face = model.face(i);
        
        for (int j=0; j < 3; j++) //x,y,z vertices
        {
            // Get vertices at the face index for face j 
            Vec3f v0 = model.vert( face[j] );
            Vec3f v1 = model.vert( face[(j+1) % 3] );
            // Scale coordinates///
            /*Add 1 to each x,y coordinate (since )*/
            int x0 = (v0.x + 1) * width / 2.;
            int y0 = (v0.y + 1) * height / 2.;
            int x1 = (v1.x + 1) * width / 2.;
            int y1 = (v1.y + 1) * height / 2.;
            line(x0, y0, x1, y1, image, color);
        }
    }
}

////////////////////////////////////////////////////////

void line_ex()
{
    TGAImage image(100, 100, TGAImage::RGB);
    line(13, 20, 80, 40, image, white); 
    line(20, 13, 40, 80, image, red); 
    line(80, 40, 13, 20, image, red);

	image.write_tga_file("line_example.tga");
}


void wireframe_ex(char *filename, char* outname)
{
    Model model(filename);
    TGAImage image(width, height, TGAImage::RGB);
    wireframe(model, image, green);

    // image.flip_vertically(); //make origin the left-bottom corner
    image.write_tga_file(outname);
}


int main(int argc, char** argv)
{
    char* output;
    if (argc==3)
        output = argv[2];
    else 
        output = "out.tga";

    wireframe_ex(argv[1], output);
}
