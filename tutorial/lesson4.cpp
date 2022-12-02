#include <vector>
#include <cmath>
#include <limits>
#include "../tga/tgaimage.h"
#include "../model/geometry.h"
#include "../model/model.h"


////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
/*                Lesson 4: Rendering in Perspective (Projection)                             */
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

/*
    In previous lessons we rendered our model in orthographic projection by simply forgetting 
    about the z-coordinate (see: http://www.cim.mcgill.ca/~langer/557/4-slides-9pp.pdf)

    In general we can take a vector in N-dimensional space and embed it into a higher N+1
    dimensional space by converting it to have homogeneous coordinates - which means appending a new
    coordinate which, in this case, can always equals 1.  
    In this N+1 dimensional space, we can represent linear transformations as matrices, which means
    we can use them to easily transform our coordinate vectors, AND this has the added but important benefit
    that we can pre-compute a series of linear transformations as a single composite transformation.
    To project our N+1 coordinate vector back into N-dimensional space, we simply divide all of our coordinates
    by the last dimension.
    So for example, if our final 3D coordinates were [x, y, z], we can project them back to 2D as [x/z, y/z].

    Now to apply to full 3D, as in computer graphics:
    We have a point (x,y,z). We convert it to homogeneous coordinates, (x,y,z,1), apply linear transformations
    (for example to modify the projection), and then divide x, y, and z by whatever the last dimension has become.

    A basic example:
    take column vector [x,y,z] and convert to [x,y,z,1].
    apply linear transformation: 
    [1 0  0   0     [x          [x
     0 1  0   0      y      =    y
     0 0  1   0      z           z
     0 0 -1/c 0]     1]         1 - z/c]
     
     Then we convert back to 3D by doing: [x / (1-z/c), y / (1-z/c), z / (1-z/c)] 


*/

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
const int width  = 800;
const int height = 800;
const int depth  = 255;

Vec3f light_dir(0,0,-1);
Vec3f camera(0,0,3);



// First some functions to convert between vectors and matrices, which we've added to geometry
/*In m2v we map from 4 dimensional space to 3 dimensional space by dividing each xyz coord. by
    the value of the final dimension.
*/
Vec3f m2v(Matrix m) 
{
    return Vec3f(m[0][0]/m[3][0], m[1][0]/m[3][0], m[2][0]/m[3][0]);
}

/*Convert a 3D vector into 4D with 1 in the 4th dim- ie, homogeneous coordinates*/
Matrix v2m(Vec3f v)
{
    Matrix m(4, 1);
    m[0][0] = v.x;
    m[1][0] = v.y;
    m[2][0] = v.z;
    m[3][0] = 1.f;
    return m;

}

/*
    Wikipedia: "the viewport is an area (typically rectangular) expressed in rendering-device-specific coordinates,
     e.g. pixels for screen coordinates, in which the objects of interest are going to be rendered. Clipping to the 
     world-coordinates window is usually applied to the objects before they are passed through the window-to-viewport
     transformation. For a 2D object, the latter transformation is simply a combination of translation and scaling, 
     the latter not necessarily uniform.[1] An analogy of this transformation process based on traditional photography 
     notions is to equate the world-clipping window with the camera settings and the variously sized prints that can be 
     obtained from the resulting film image as possible viewports."

     Good info: https://phoenix.goucher.edu/~kelliher/s2000/cs320/mar27.html
*/
Matrix viewport(int x, int y, int w, int h)
{
    Matrix m = Matrix::identity(4);
    m[0][3] = x+w/2.f;
    m[1][3] = y+h/2.f;
    m[2][3] = depth/2.f;

    m[0][0] = w/2.f;
    m[1][1] = h/2.f;
    m[2][2] = depth/2.f;
    return m;
}

Matrix projection()
{
    Matrix mat = Matrix::identity(4);
    mat[3][2] = -1.f / camera.z;
    return mat;
}


/*Texture mapped triangle, from lesson 3*/
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


/*
    Take model coordinates, map to world coordinates, use linear transformations to change the perspective (thus modifying the coordinates),
    and then map those coordinates to the screen.
    This is summed up by the line below: `m2v(ViewPort*Projection*v2m(v))`,
        reading right to left:
         we convert triangle vertex v to matrix form which also converts to homogeneous coords (adds another dimension),
         then multiply by the projection matrix to change the projection,
         then multiply by the viewport to convert from the projected "world" coords. to the current viewport,
            which is a transformation in itself,
         then we convert back to vector format with `m2v` which also divides each xyz coord by the added 4th dimension,
            thus projection from 4D homogeneous coordinates back to 3D 
*/
void project_and_draw(Model &model, int *zbuffer, TGAImage& image)
{
    /*Get diffuse texture map from model*/
    TGAImage diffusemap = model.diffuse();


    Matrix ViewPort = viewport(width/8, height/8, width*3/4, height*3/4);
    Matrix Projection = projection();

    for (int i = 0; i < model.nfaces(); i++)
    {
        Vec3i screen_coords[3]; 
        Vec3f world_coords[3]; 
        Vec2i uv_coords[3];
        for (int j = 0; j < 3; j++)  
        {
            Vec3f v = model.vert( i, j );
            screen_coords[j] = m2v(ViewPort*Projection*v2m(v));
            world_coords[j] = v;
        }
        // Compute triangle's normal and use to calculate light intensity fr the current face
        Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1]-world_coords[0]); //side 1 cross-product side 2
        n.normalize();
        float intensity = n*light_dir;

        // Draw triangle if not back-of-face
        if (intensity > 0)
        {
            for (int k=0; k<3; k++)
            {
                Vec2f uv = model.uv( i, k );
                uv_coords[k] = {(int)(uv.x * diffusemap.width()), (int)(uv.y * diffusemap.height())};
            }

            texture_triangle(screen_coords, zbuffer, uv_coords, intensity, diffusemap, image);
        }
    }
}



int main(int argc, char** argv)
{
    const char* mod;
    if (argc==1)
        mod = "../obj/african_head.obj";     
    else
        mod = argv[1];

    Model model(mod);
    TGAImage image(width, height, TGAImage::RGB);

    /*Initialize zbuffer with smallest possible depth (so we
    start rendering infront of the zbuffer*/
    int *zbuffer = new int[width*height];
    for (int i=0; i < width*height; i++)
        zbuffer[i] = std::numeric_limits<int>::min();

    project_and_draw(model, zbuffer, image);
    image.write_tga_file("out.tga");
}