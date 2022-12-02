#include <iostream>
#include <vector>
#include "../tga/tgaimage.h"
#include "../model/geometry.h"
#include "../model/model.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
const int width  = 1024;
const int height = 1024;

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
/*                Lesson 3: Removal of Hidden Faces (Z Buffering)                             */
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

/*If we just draw all triangles in our model from back to front, the front facets will
erase the back ones, and this ay technically work (called the painter's algorithm). In
practice, it is highly computational and we may not be sure what order to render the triangles.

Instead, we can "lose a dimension" - we will lose z, or the depth dimension. When a triangle
is projected on the screen, the depth (z-value) of a generated pixel is in the projected screen
is compared to the value already stored in the z-buffer, and it replaces that pixel only if the
new value is closer to the "camera". 
*/

/* First, let's implement a helper function, which maps the model's (or "world's") coordinate
 scale to the image's scale (just as we've done repeatedly in prior lessons)*/ 
Vec3f world2screen(Vec3f v)
{
    return Vec3f(
        int((v.x + 1.) * width/2. + .5),
        int((v.y + 1.) * height/2. + .5),
        v.z
    );
}


// Implementation
/*
    We will update our "triangle" rasterization function to incorporate the z-buffer.
    Again, we start by defining a bounding box 
*/
Vec3f barycentric(Vec3f *pts, Vec3f P)
{
    /*Compute cross-product to calculate barycentric coords*/
    Vec3f u = Vec3f(pts[2].x-pts[0].x, pts[1].x-pts[0].x, pts[0].x-P.x //Cx-Ax, Bx-Ax, Ax-Px
              )^ Vec3f(pts[2].y - pts[0].y, pts[1].y - pts[0].y, pts[0].y - P.y); //Cy-Ay, By-Ay, Ay-Py 
    
    /* u.z is actually an integer, and thus if it is zero then the point is not within
    the triangle defined by 'pts', and so we should return the vector of weights only
    if it u.z > 0 */
    if (std::abs(u.z) > 1e-2) 
        return Vec3f( 1.0f - (u.x+u.y) / u.z,
                                   u.y / u.z,
                                   u.x / u.z);
    // Otherwise return negative coords so triangle is not drawn
    return Vec3f(-1, -1, -1);
}

/* Triangle Rasterizer
    Now we have some updates to our triangle function. First we compute the bounding box as 
    before, initializing it's coords with maximum positive and negative float values and then
    replacing those values where the triangles 'pts', iterating over all 3 of the triangle's 
    vertices to ensure we use the bottom-left corner and top-right corner for the box's points.

    Then, we iterate through all points in the bounding box and for each:
        We calculate the barycentric coordinates, and if any are < 0 the pixel is not
        rendered since it point falls outside of the triangle.
        If it does fall within the triangle:
            We position the current point P's depth at the weighted sum of the triangle's
            vertices' z-dimensions, where the weights are the barycentric coordinates.**
            Then, if the current points depth is closer to the camera (ie, has a larger value
            than) the z-buffer, we render the pixel and we update the zbuffer so that the new
            closest depth at that pixel is that of current point P. 
            If the depth is NOT closer than the z-buffer, we do not render this point since there
            are other pixels which cover it up.

            **more on the computation of the point's depth: 
                We know the triangle's vertices' Z-coordinates, so to compute the actual depth of 
                the current point P, we need to interpolate between all the vertices' z-coordinates.
                We use the barycentric coordinates to do so:
                "If the barycentric coordinates are used to compute the position of a point located
                 on the triangle using the triangle vertices, we can interpolate any other data defined
                 at the triangle's vertices (like for example the color) in the exact same way.
                 In other words, barycentric coordinates are used to interpolate vertex data across the
                 triangle's surface (the technique can be applied to any data type, float, color, etc.). "
                 - https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/barycentric-coordinates
*/

void triangle(Vec3f *pts, float *zbuffer, TGAImage &image, TGAColor color)
{
    Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max() );
    Vec2f bboxmax( -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max() );
    Vec2f clamp(image.width() - 1, image.height() - 1); //clip bb to image dims if triangle extends outside image
    // Compute bounding box
    for (int i = 0; i < 3; i++) //for triangle vertex i
        for (int j = 0; j < 2; j++) //for coord (of x,y) j
        {
            bboxmin[j] = std::max(0.0f,     std::min(bboxmin[j], pts[i][j]));
            bboxmax[j] = std::max(clamp[j], std::max(bboxmax[j], pts[i][j]));
        }
    Vec3f P;
    // For all points within the bounding box, computer the barycentric
    // coordinates. If any are less than 0, the current point is not in
    // the triangle, and we should skip. Otherwise, determine if the
    // point's z dimension (depth) is closer to camera then the z buffer,
    // and only render the point if it is.
    for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++)
        for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++)
        {
            Vec3f bc = barycentric(pts, P);
            if (bc.x < 0 || bc.y < 0 || bc.z < 0) continue;
            P.z = 0;
            for (int i = 0; i < 3; i++)
                P.z += pts[i].z * bc[i];

            int idx = P.x + P.y*width;
            if ( zbuffer[idx] < P.z )
            {
                zbuffer[idx] = P.z;
                image.set(P.x, P.y, color);
            }
        }
}


// Function to illuminate model while incorporating z-buffering
// Good info on minimal lighting: https://math.hws.edu/graphicsbook/c7/s2.html
void illuminate_model(Model &model, TGAImage &image)
{
    /*Initialize zbuffer with largest possible negative depth (so we
        start rendering infront of the zbuffer*/
    float *zbuffer = new float[width*height];
    for (int i=width*height; i--; zbuffer[i] = -std::numeric_limits<float>::max());

    /*Specify light vector for illumniation*/
    Vec3f light_dir(0, 0, -1);

    /*Iterate through the model's triangles ("faces") and render via the
     bounding-box + barycentrix-coords method, while also checking the z-buffer to
     remove "hidden" faces.
      */
    for (int i = 0; i < model.nfaces(); i++)
    {
        Vec3f screen_coords[3]; 
        Vec3f world_coords[3]; 
        for (int j = 0; j < 3; j++)  
        {
            Vec3f v = model.vert( i, j );
            screen_coords[j] = world2screen(v); //scale "world" coords to fit screen
            world_coords[j] = v;
        }
        // Compute triangle's normal and use to calculate light intensity fr the current face
        Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1]-world_coords[0]); //side 1 cross-product side 2
        n.normalize();
        float intensity = n*light_dir;
        if (intensity > 0)
        {
            triangle(screen_coords, zbuffer, image,
                     TGAColor(intensity*225, intensity*193, intensity*110, 255)
                    );
        }
    }
}


///////////////// PART 2 ////////////////////////////////////////////////
/*Attempting to add Texture
    -Mapping in the "diffuse map" using the vertex-texture ("vt") coords.
    The tutorial does not do it with barycentric coord method (yet), nonetheless I tried (unsuccessfully).
    -As before, we create a bounding box and calculate the barycentric coordinates for each point to determine
    if the point is in the triangle.
    -Now we also use the barycentric points to to calculate the uv points which index into our texture map.
     The barycentric coordinates define the position of our current point P in relation to a triangle - 
        specifically the triangle we have provided and wish to draw. As above, barycentric coordinates allow us 
        to interpolate across that triangle, this time for the purpose of pulling in a color value from the texture map.**

    **More on mapping our square texture onto a triangle:
        https://computergraphics.stackexchange.com/questions/1866/how-to-map-square-texture-to-triangle:
        "First, we find the barycentric coordinates of P. Barycentric coordinates represent how much weight each vertex
         contributes to the point, and can be used to interpolate any value which is known at the vertices across the face
         of a triangle...
         Once you have [barycentric] coordinates, you can determine the texture coordinates of P by interpolating the values 
         at the vertices using the barycentric coordinates as weights"

*/
Vec3f barycentric(Vec3i *pts, Vec3i P)
{
    Vec3f u = Vec3f(pts[2].x-pts[0].x, pts[1].x-pts[0].x, pts[0].x-P.x //Cx-Ax, Bx-Ax, Ax-Px
              )^ Vec3f(pts[2].y - pts[0].y, pts[1].y - pts[0].y, pts[0].y - P.y); //Cy-Ay, By-Ay, Ay-Py 
    if (std::abs(u.z) > 1e-2) 
        return Vec3f( 1.0f - (u.x+u.y) / u.z,
                                   u.y / u.z,
                                   u.x / u.z);
    return Vec3f(-1, -1, -1);
}

void texture_triangle(Vec3i *pts, int *zbuffer,  Vec2i *uvPts, float intensity, TGAImage &diffusemap, TGAImage &image)
{
    Vec2i bboxmin( std::numeric_limits<int>::max(),  std::numeric_limits<int>::max() );
    Vec2i bboxmax( -std::numeric_limits<int>::max(), -std::numeric_limits<int>::max() );
    Vec2i clamp(image.width() - 1, image.height() - 1); //clip bb to image dims if triangle extends outside image
    // Compute bounding box
    for (int i = 0; i < 3; i++) //for each triangle vertex, i
        for (int j = 0; j < 2; j++) //for each coord (x & y), j
        {
            bboxmin[j] = std::max(0,     std::min(bboxmin[j], pts[i][j]));
            bboxmax[j] = std::max(clamp[j], std::max(bboxmax[j], pts[i][j]));
        }
    Vec3i P;
    // For all points within the bounding box, computer the barycentric
    // coordinates. If any are less than 0, the current point is not in
    // the triangle, and we should skip. Otherwise, determine if the
    // point's z dimension (depth) is closer to camera then the z buffer,
    // and only render the point if it is.
    // Use the barycentric coords. to calculate depth and to interpolate
    // the uv points.
    for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++)
        for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++)
        {
            Vec3f bc = barycentric(pts, P);
            if (bc.x < 0 || bc.y < 0 || bc.z < 0) continue;
            P.z = 0;
            for (int i = 0; i < 3; i++)
                P.z += pts[i].z * bc[i];

            int idx = P.x + P.y*width;
            if ( zbuffer[idx] < P.z )
            {
                zbuffer[idx] = P.z;
                Vec2f uv(0,0);
                for (int i = 0; i < 3; i++)
                {
                    uv.x += uvPts[i].x * bc[i];
                    uv.y += uvPts[i].y * bc[i];
                }

                TGAColor color = diffusemap.get(uv.x, uv.y);

                image.set(
                    P.x, P.y, 
                    TGAColor( *(color.r) * intensity, *(color.g) * intensity,
                              *(color.b) * intensity, 255)
                );
            }
        }
}


// Adding Texture
/*Here we get our diffuse map (the texture map) which was loaded by the model,
    then initialize the z-buffer and light direction.
    Then we began iterating through the faces in our model. For each face,
    we get the triangle coordinates and scale them to screen size. We also get
    the coordinates for our texture map.
    Then we compute the triangle's normal to calculate the light intensity,
    and only draw the triangle if the intensity is positive (back-face culling).
    The triangle rasterization method needs the triangle scren coordinates, the
    texture coordinates for mapping each pixel to the correct color value, the intensity
    level for modifying the color based on light intensity, the zbuffer to determine depth
    (ie, should each pixel be drawn or is it hidden by other faces), the diffuse texture map to
    sample color/texture from, and the image to draw onto. 
*/
void render_texture_model(Model &model, TGAImage &image)
{
    /*Load in model's diffuse map for texture mapping*/
    TGAImage diffuse_map = model.diffuse();

    /*Initialize zbuffer with smallest possible depth (so we
        start rendering infront of the zbuffer*/
    int *zbuffer = new int[width*height];
    for (int i=0; i < width*height; i++)
        zbuffer[i] = std::numeric_limits<int>::min();

    /*Specify light vector for illumniation*/
    Vec3f light_dir(0, 0, -1);

    /*Iterate through the model's triangles ("faces") and render via the
     bounding-box + barycentrix-coords method, while also checking the z-buffer to
     remove "hidden" faces.
      */
    for (int i = 0; i < model.nfaces(); i++)
    {
        Vec3i screen_coords[3]; 
        Vec3f world_coords[3]; 
        Vec2i tex_coords[3];
        for (int j = 0; j < 3; j++)  
        {
            Vec3f v = model.vert( i, j );
            Vec2f uv = model.uv( i, j );
            screen_coords[j] = world2screen(v); //scale "world" coords to fit screen
            world_coords[j] = v;
            tex_coords[j] = {(int)(uv.x * diffuse_map.width()), (int)(uv.y * diffuse_map.height())};
        }
        // Compute triangle's normal and use to calculate light intensity fr the current face
        Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1]-world_coords[0]); //side 1 cross-product side 2
        n.normalize();
        float intensity = n*light_dir;

        if (intensity > 0)
        {
            texture_triangle(screen_coords, zbuffer, tex_coords, intensity, diffuse_map, image);
        }
    }
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{ 
    const char* mod;
    if (argc==1)
        mod = "../obj/african_head.obj";     
    else
        mod = argv[1];

    Model model(mod);
    TGAImage image(width, height, TGAImage::RGB);

    render_texture_model(model, image);
    image.write_tga_file("out.tga");
}







// How to implement textured triangle with the line sweeping method

void linesweep_texture_triangle(Vec3i* pts, Vec2i* uv, float intensity, int *zbuffer, TGAImage &diffusemap, TGAImage &image) {
    if (pts[0].y==pts[1].y && pts[0].y==pts[2].y) return; // degenerate triangle
    // Sort points so point 0 is bottom left and point 2 is top right
    if (pts[0].y>pts[1].y) { std::swap(pts[0], pts[1]); std::swap(uv[0], uv[1]); }
    if (pts[0].y>pts[2].y) { std::swap(pts[0], pts[2]); std::swap(uv[0], uv[2]); }
    if (pts[1].y>pts[2].y) { std::swap(pts[1], pts[2]); std::swap(uv[1], uv[2]); }

    // Line-sweeping method in one loop
    int total_height = pts[2].y-pts[0].y;
    for (int i=0; i<total_height; i++) 
    {
        bool second_half = i>pts[1].y-pts[0].y || pts[1].y==pts[0].y; //if current y position is in the second segment
        int segment_height = second_half ? pts[2].y-pts[1].y : pts[1].y-pts[0].y;
        float alpha = (float)i/total_height;
        float beta  = (float)(i-(second_half ? pts[1].y-pts[0].y : 0))/segment_height; 
        Vec3i A   =               pts[0]  + Vec3f(pts[2]-pts[0])*alpha;
        Vec3i B   = second_half ? pts[1]  + Vec3f(pts[2]-pts[1])*beta   : pts[0]  + Vec3f(pts[1]-pts[0]  )*beta;
        // Get the current A-side and B-side coords in the uv-texture map (diffuse map)
        Vec2i uvA =               uv[0] +      (uv[2]-uv[0])*alpha;
        Vec2i uvB = second_half ? uv[1] +      (uv[2]-uv[1])*beta       : uv[0] +      (uv[1]-uv[0])*beta;
        if (A.x>B.x) { std::swap(A, B); std::swap(uvA, uvB); }
        // For the given boundary positions, sweep through the pixels between them,
        // calculating a pixel poisition for the image (P) and the texture map (uvP).
        for (int j=A.x; j<=B.x; j++) 
        {
            float phi = B.x==A.x ? 1.0f : (float)(j-A.x)/(float)(B.x-A.x);
            Vec3i   P = Vec3f(A) + Vec3f(B-A)*phi;
            Vec2i uvP =     uvA +   (uvB-uvA)*phi;
            int idx = P.x+P.y*width;
            if (zbuffer[idx]<P.z) 
            {
                zbuffer[idx] = P.z;
                TGAColor color = diffusemap.get(uvP.x, uvP.y); //use the texture map coords to get the color value
                image.set(P.x, P.y, TGAColor(
                    *(color.r)*intensity, *(color.g)*intensity, *(color.b)*intensity, 255
                    ));
            }
        }
    }
}


