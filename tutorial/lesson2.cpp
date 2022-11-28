#include <iostream>
#include "../tga/tgaimage.h"
#include "../model/geometry.h"
#include "../model/model.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
const int width  = 200;
const int height = 200;

// Reusing from lesson 1:
void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color);


////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
/*                Lesson 2: Triangle Rasterization and back face culling                       */
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////


/*Part 1: Moving from Line Wireframe to Filled-in Triangles
    OpenGL triangulates almost any polygon so it is a very useful exercise.
    We simply draw lines connecting 3 points, making sure to sort them first so
    that the image does not depend on the order that the vertices were passed in. 
 */

void basic_triangle2d(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color)
{
    // sort
    if (t0.y > t1.y) std::swap(t0, t1);
    if (t0.y > t2.y) std::swap(t0, t2);
    if (t1.y > t2.y) std::swap(t1, t2);
    // draw
    line(t0.x, t0.y, t1.x, t1.y, image, color);
    line(t1.x, t1.y, t2.x, t2.y, image, color);
    line(t2.x, t2.y, t0.x, t0.y, image, color);
}

// Old-School: Line Sweeping
/*Draw 2d triangles and fill them in by drawing horizontal line-segments between left and
  right triangle boundaries.
  First consider that one side of a triangle is boundary A and the other 2 sides make up boundary B.
  Unless the triangle's base is completely flat, boundary B will actually have two parts so we need to
  think of how to fill each section.
  We can split the triangle into two parts by cutting it horizontally at the vertex that breaks up boundary
  B and fill each section separetly.
  */
void fill_triangle2d(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color)
{
    // sort
    if (t0.y > t1.y) std::swap(t0, t1);
    if (t0.y > t2.y) std::swap(t0, t2);
    if (t1.y > t2.y) std::swap(t1, t2);
    int total_height = t2.y - t0.y;
    // Triangle Section 1
    for (int y = t0.y; y <= t1.y; y++)
    {
        int segment_height = t1.y - t0.y+1;
        // Determine coefficients for line drawing
        float alpha = (float)(y - t0.y) / total_height;
        float beta;
        if (segment_height > 0) 
            beta = (float)(y - t0.y) / segment_height;
        else beta = 0.0f;
        // At the given y position, determine the x positions
        // which correspond to boundary A and boundary B. Then
        // sweep horizontally between the xs and color the pixel
        Vec2i A = t0 + (t2-t0)*alpha; //x pos. on bound. A
        Vec2i B = t0 + (t1-t0)*beta;  //x pos. on bound. B sec. 1
        if (A.x > B.x) std::swap(A, B);
        for (int j=A.x; j <= B.x; j++)
            image.set(j, y, color);
    }

    // Triangle Section 2
    for (int y = t1.y; y <= t2.y; y++)
    {
        int segment_height = t2.y - t1.y+1;
        float alpha = (float)(y - t0.y) / total_height;
        float beta;
        if (segment_height > 0) 
            beta = (float)(y - t1.y) / segment_height;
        else beta = 0.0f;
        Vec2i A = t0 + (t2-t0)*alpha; //x pos. on bound. A
        Vec2i B = t1 + (t2-t1)*beta;  //x pos. on bound. B sec. 2
        if (A.x > B.x) std::swap(A, B);
        for (int j = A.x; j <= B.x; j++)
            image.set(j, y, color);
    }
}





void triangle_ex(TGAImage &image)
{   
    // Triangles
    Vec2i t0[3] = {Vec2i(10, 70),   Vec2i(50, 160),  Vec2i(70, 70)};  //flat base
    Vec2i t1[3] = {Vec2i(180, 50),  Vec2i(150, 1),   Vec2i(70, 180)}; 
    Vec2i t2[3] = {Vec2i(180, 150), Vec2i(120, 160), Vec2i(130, 180)}; 

    fill_triangle2d(t0[0], t0[1], t0[2], image, green);
    fill_triangle2d(t1[0], t1[1], t1[2], image, white);
    fill_triangle2d(t2[0], t2[1], t2[2], image, red);
}

int main()
{
    TGAImage image(width, height, TGAImage::RGB);
    triangle_ex(image);


    image.write_tga_file("out.tga");
}




////////////////////////////
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
