#include <iostream>
#include <vector>
#include "../tga/tgaimage.h"
#include "../model/geometry.h"
#include "../model/model.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
const int width  = 800;
const int height = 800;

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

// Bounding Box and Barycentric Coordinates

/*Insted of sweeping between the bondaries of a triangle, we can compute a triangle's
  bounding *box* and compute the barycentric coordinates of each pixel in that box. If
  any component of the coords. is negative, the pixel is outside the triangle.
  -Barycentric coordinates: a coordinate system in which the location of a point is specified
   in reference to a "simplex" - which is equal to a triangle for points in a 2d plane.
   *Wikipedia: the barycentric coords. of one point can be interpreted as masses at the vertices of
   the simplex (a triangle here, so 3 vertices) such that the point is the center of the mass. 
   The masses can be zero or negative, and they are all positive IF AND ONLY IF the point is
   inside the simplex (triange).
   -We can find the barycentric coords for point P by placing weights (1-u-v, u, v) at the triangle's
    3 vertices (A,B,C) -> P = (1-u-v)A + uB + vC -> P = A + uAB + vAC
    We have vectors AB (connecting vertex A to B), AC (connecting A to C), and AP (connecting A to P,
    which may or may not be within the triangle). 
    We need to find 2 numbers, u and v, which satisfy
     {uAB_x + vAC_x + PA_x = 0   &   uAB_y + vAC_y + PA_y = 0} , a linear system of equations.
     Written in matrix notation, we need vector [u, v, 1]'s product with col vector [AB_x, AC_x, PA_x] to = 0
     while simultaneously it's product with col vector [AB_y, AC_y, PA_y] must also = 0. 
     In other words, we are looking for the intersection of these two lines, which means we can just compute
     their cross-product.

  New rasterization routine:
  -Iterate through the pixels of a bounding box, calculating the barycentric coordinates for that point.
  -If it has at least one negative component, the pixel is outside the triangle.
*/

Vec3f barycentric(Vec2i *pts, Vec2i P)
{
    /*Compute cross-product to calculate barycentric coords*/
    Vec3f u = Vec3f(pts[2].x-pts[0].x, pts[1].x-pts[0].x, pts[0].x-P.x //Cx-Ax, Bx-Ax, Ax-Px
              )^ Vec3f(pts[2].y - pts[0].y, pts[1].y - pts[0].y, pts[0].y - P.y); //Cy-Ay, By-Ay, Ay-Py 
    
    /* 'pts' and 'P' are integer vectors, so if abs(u.z) < 0, then u.z == 0, and the current point
        'P' is not in the triangle. Thus, we should return a vector with some negative coords.*/
    if (std::abs(u.z) < 1) return Vec3f(-1, -1, -1);

    /*Otherwise, return the vector of "weights" (equivalent to 1-u-v, u, v)*/
    return Vec3f( 1.0F - (u.x+u.y) / u.z,
                               u.y / u.z,
                               u.x / u.z);
}


/*Now to draw a triangle, we first compute a bounding box.
    The bounding box is described by two points: the bottom left and upper right vertices.
    To find these corners, we iterate through the vertices of the triangle and choose min & max
    coords. 
*/
void triangle(Vec2i *pts, TGAImage &image, TGAColor color)
{
    Vec2i bboxmin(image.width() - 1, image.height() - 1);
    Vec2i bboxmax(0, 0);
    Vec2i clamp(image.width() - 1, image.height() - 1); //clip bb to image dims if triangle extends outside image
    // Compute bounding box
    for (int i = 0; i < 3; i++)
    {
        bboxmin.x = std::max(0, std::min(bboxmin.x, pts[i].x));
        bboxmin.y = std::max(0, std::min(bboxmin.y, pts[i].y));

        bboxmax.x = std::min(clamp.x, std::max(bboxmax.x, pts[i].x));
        bboxmax.y = std::min(clamp.y, std::max(bboxmax.y, pts[i].y));
    }
    Vec2i P;
    // For all points within the bounding box, computer the barycentric
    // coordinates. If any are less than 0, the current point is not in
    // the triangle, and we should skip.
    for (P.x = bboxmin.x; P.x < bboxmax.x; P.x++)
        for (P.y = bboxmin.y; P.y < bboxmax.y; P.y++)
        {
            Vec3f bc = barycentric(pts, P);
            if (bc.x < 0 || bc.y < 0 || bc.z < 0) continue;
            image.set(P.x, P.y, color);
        }
}

// Flat Shading Rendering
/*To test our triangle drawing/filling, let's render the african_head.obj model and randomly
    shade each triangle.
    We first load in the model and iterate through all its faces, recalling that for each face is really
    an list of 3 indexes to 3 corresponding vertices (3 pairs of x,y,z coords) which make up a triangle. 
    Then we scale the x and y coords to fit the screen and draw them with the triangle function. 
*/
// In main just call: flat_color_shading(image)
void flat_color_shading(Model &model, TGAImage &image)
{
    for (int i = 0; i < model.nfaces(); i++)
    {
        Vec2i screen_coords[3];
        for (int j = 0; j < 3; j++)
        {
            Vec3f world_coords = model.vert(i, j); //The actual vertices of the current face
            screen_coords[j] = Vec2i(                 //The coords scaled to be shown on screen
                (world_coords.x+1.)*width/2., (world_coords.y+1.)*height/2.
            );
        }
        triangle(screen_coords, image, TGAColor(rand()%255, rand()%255, rand()%255, 255));
    }
}


// Polygon Shading & back-face culling
/*
    Instead of randomly coloring the triangles, we can shade them based on a predent light source.
    Key to this is the simple fact that a" polygon is illuminated most brightly when it is orthogonal
    to the light direction". In practice, the intensity of illumination can be calculated simply as the
    scalar product of the light vector and the *normal* to the given triangle, where the normal is simply
    the cross-product of the triangle's 2 sides (a surface normal is just a unit vector that is perpendicular
    to a surface at a specific spot).

    Back-face culling:
    Note that when we compute the itensity of light upon the polygon as the scalar product of the normal and
    the light vector, it can be negative, which indicates that light is actually coming from behind the polygon
    (it's hitting the back of the face). In this case, we do not draw the triangle.

    Gamma correction:
    To display different levels of illumination, we modify the color value of the triangle with the intensity level.
    In reality, (128, 128, 128) is not actually half as bright as (255, 255, 255), since humans perceive light and 
    color in a non-linear manner. Gamma encoding of images corrects for this, but we ignore that for now.
*/
void illuminate_model(Model &model, TGAImage &image)
{
    Vec3f light_dir(0, 0, -1);

    for (int i = 0; i < model.nfaces(); i++)
    {
        Vec2i screen_coords[3]; 
        Vec3f world_coords[3]; 
        for (int j = 0; j < 3; j++)  
        {
            Vec3f v = model.vert(i, j);
            screen_coords[j] = Vec2i( (v.x+1.)*width/2., (v.y+1.)*height/2. );
            world_coords[j] = v;
        }
        // Compute normal
        Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1]-world_coords[0]); //side 1 cross-product side 2
        n.normalize();
        float intensity = n*light_dir;
        if (intensity > 0)
        {
            triangle(screen_coords, image, TGAColor(intensity*255, intensity*255, intensity*255, 255));
        }
    }
}



//  (Basic triangle drawing testers)
void linesweep_triangle_ex(TGAImage &image)
{   
    // Triangles
    Vec2i t0[3] = {Vec2i(10, 70),   Vec2i(50, 160),  Vec2i(70, 70)};  //flat base
    Vec2i t1[3] = {Vec2i(180, 50),  Vec2i(150, 1),   Vec2i(70, 180)}; 
    Vec2i t2[3] = {Vec2i(180, 150), Vec2i(120, 160), Vec2i(130, 180)}; 

    fill_triangle2d(t0[0], t0[1], t0[2], image, green);
    fill_triangle2d(t1[0], t1[1], t1[2], image, white);
    fill_triangle2d(t2[0], t2[1], t2[2], image, red);
}


void bb_triangle_ex(TGAImage &image)
{
    // Triangle vertices
    Vec2i pts[3] = {Vec2i(10, 10), Vec2i(100, 30), Vec2i(190, 160)};
    // Draw triangle via bounding box method
    triangle(pts, image, red);
}





int main()
{
    TGAImage image(width, height, TGAImage::RGB);
    Model model("./obj/african_head.obj");

    illuminate_model(model, image);
    image.write_tga_file("out.tga");
}





////////////////////////////////////////////////////////////////////////////////////
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
    }
}
