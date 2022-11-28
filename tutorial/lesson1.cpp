#include <iostream>
#include "../tga/tgaimage.h"
#include "../model/geometry.h"
#include "../model/model.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
const int width  = 800;
const int height = 800;

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
/*                      Lesson 1: Bresenham’s Line Drawing Algorithm                          */
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////



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




/* Final Attempt + Optimized */
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
    if(steep) { //individual loops is an optimization tactic: https://github.com/ssloy/tinyrenderer/issues/28
        for(int x = x0; x <= x1; ++x) {
            image.set(y, x, color);
            error2 += derror2;
            if (error2 > dx) {
                y += (y1>y0 ? 1 : -1);
                error2 -= dx*2;
            }
        }
    } else {
        for(int x = x0; x <= x1; ++x) {
            image.set(x, y, color);
            error2 += derror2;
            if (error2 > dx) {
                y += (y1>y0 ? 1 : -1);
                error2 -= dx*2;
            }
        }
    }}


/* Basic Wireframe Rendering */

/* We can use a basic file format to define vertices and render a wireframe model. 
    The .OBJ file is a geometry definition file format.
    The "Model" class reads a Wavefront OBJ file and stores vertices and faces.
    The "geometry.h" header contains basic vector types and operations for storing and 
    manipulating vectors.

    How does this function work?
    -We cycle through the list of "faces" in our model. For each face, we iterate (j) over the
     x, y, and z vertices/coordinates in order to draw lines.
     We will draw lines from vertex j -> vertex (j+1)%3, which means:
        vertex 0 -> vertex 1, vertex 1 -> vertex 2, vertex 2 -> vertex 0
    -These vertices are pulled from the model by indexing with the face index - remember, a face 
     is a list of vertex indexes (at least, the first number after each space) so at each face, we 
     can access its relevant vertices by indexing the 0th, 1st, and 2nd vertex within the face.
    -Now that we have the vertices stored in our vec structs (geometry.h) we can scale the x and y
     coordinates. The tutorial adds 1 to each (I'm not sure why this is, OBJ files are unitless so
     pit is difficult to say for the "african_head" example. It may be just to keep the coordinates
     positive, since they appear to be between -1 and 1. Removing the +1 results in an image that is
     not correctly positioned.). (Also, translating all the vertices by a constant scalar will not
     affect the output shape, just its position in the image).
     Then we scale the coords. by width or height /2 so that the the output is centered within our image.
     Thus, after scaling the coordinates they are equivalent to a pixel within our image.
    -Also note that, for now, we are doing nothing with the z dimension.
    */
void wireframe(Model &model, TGAImage &image, TGAColor color)
{
    for (int i=0; i < model.nfaces(); i++)
    {
        std::vector<int> face = model.face(i);
        
        for (int j=0; j < 3; j++) //x,y,z vertices
        {
            // Get list of vertices at the face index 
            Vec3f v0 = model.vert( face[j] );
            Vec3f v1 = model.vert( face[(j+1) % 3] );
            // Scale coordinates
            int x0 = (v0.x + 1.) * width / 2.;
            int y0 = (v0.y + 1.) * height / 2.;
            int x1 = (v1.x + 1.) * width / 2.;
            int y1 = (v1.y + 1.) * height / 2.;
            line(x0, y0, x1, y1, image, color);
        }
    }
}

////////////////////////////////////////////////////////

void line_ex()
{
    //Test line function
    TGAImage image(100, 100, TGAImage::RGB);
    line(13, 20, 80, 40, image, white); 
    line(20, 13, 40, 80, image, red); 
    line(80, 40, 13, 20, image, red);

	image.write_tga_file("line_example.tga");
}



void wireframe_ex(const char *filename, const char* outname)
{
    //Make wireframe model 
    Model model(filename);
    TGAImage image(width, height, TGAImage::RGB);
    wireframe(model, image, green);

    // image.flip_vertically(); //make origin the left-bottom corner
    image.write_tga_file(outname);
}


int main(int argc, char** argv)
{
    const char* output;
    if (argc==3)
        output = argv[2];
    else 
        output = "out.tga";

    wireframe_ex(argv[1], output);
}
