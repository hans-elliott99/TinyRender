#include <vector>
#include <cmath>
#include <limits>
#include "../tga/tgaimage.h"
#include "../model/geometry.h"
#include "../model/model.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
const int width  = 800;
const int height = 800;
const int depth  = 255;

Vec3f light_dir(0,0,-1);
Vec3f eye(0,0,3);
Vec3f center(0,0,0);


////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
/*                Lesson 4: Rendering in Perspective (Projection)                             */
/*                Lesson 5: Moving the Camera                                                 */
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

/*
Summary: Our models (characters, for example) are created in their own local frame (object coordinates). 
         They are inserted into a scene expressed in world coordinates. The transformation from one to another 
         is made with matrix Model. Then, we want to express it in the camera frame (eye coordinates), the
         transformation is called ModelView. Then, we deform the scene to create a perspective deformation with
         Projection matrix (lesson 4), this matrix transforms the scene to so-called clip coordinates.
         Finally, we draw the scene, and the matrix transforming clip coordinates to the screen coordinates is called Viewport.
*/



/*
    LESSON 4
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
    
    Again, we are using the "perpsecitve projection" here. Before we used an orthogonal projection for projecting our 3D points onto a 2D
    plane in such a way that the projection preserves the orthogonality (perpendicularity) of the axes on the 2D plane. This means
    that the projected coordinates maintain the same angles and relationships as the original 3D coordinates, allowing for accurate
    representation of 3D shapes in a 2D space.
    In a perspective projection, the 3D coordinates of a point are mapped onto a 2D plane in such a way that objects that are farther away
    from the viewer appear smaller on the 2D plane. This creates the illusion of depth and distance, and allows for a more realistic representation
    of 3D scenes. 



*/



// First some functions to convert between vectors and matrices - which we've added to geometry
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
    First, the VIEWPORT:
    In prior lessons we would convert from the model's coordinates (or "world" coords.) to screen coordiantes by
    scaling them like so: screen_coords[j] = Vec2i((v.x+1.)*width/2., (v.y+1.)*height/2.);
    But we can do this using matrix multiplication . The code below creates the matrix:
        [ w/2   0   0   x+(w/2)
        0     h/2 0   y+(h/2)
        0     0   d/2     d/2
        0     0   0         1]
    which maps the "bi-unit" cube [-1,1]*[-1,1]*[-1,1] (which is what our model gets from OBJ) to the screen 
    cube [x, x+w]*[y, y+h]*[0, d] (d is the resolution of z buffer)

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
    // mat[3][2] = -1.f / eye.z;
    mat[3][2] = -1.f/(eye-center).norm();
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
            Vec2i uv_coords[3];
            for (int k=0; k<3; k++)
            {
                Vec2f uv = model.uv( i, k );
                uv_coords[k] = {(int)(uv.x * diffusemap.width()), (int)(uv.y * diffusemap.height())};
            }

            texture_triangle(screen_coords, zbuffer, uv_coords, intensity, diffusemap, image);
        }
    }
}



/*
    LESSON 5
    The "lookat" function - or, changing where the camera is looking in the scene.   
    (see: https://www.scratchapixel.com/lessons/mathematics-physics-for-computer-graphics/lookat-function)
    To change what the camera is looking at in the scene, we need to define a point in space
    where the camera exists (from) and a point in space where it looks (to).  
    Now, we can adjust what the camera "sees" by applying appropriate linear transformations to
    the coordinates of our model. As always, we will use a 4x4 matrix to perform these linear
    transformations on our (cartesian) coordinate system. The matrix will take the form:
        [Right_x   Right_y   Right_z    0
        Up_x      Up_y      Up_z        0
        Forward_x Forward_y Forward_z   0
        T_x       T_y       T_z         1]
    where right is the x-axis, up is the y-axis, forward is the z-axis, and T is the added dimension for 
    homogenous coordinates.
    To compute each row of the above matrix, we perform the following steps:
        1 - Compute the forward (z) axis: the forward axis is aligned between points "from" and "to", and thus
            we just need to calculate (From-To) and normalize the resulting vector.
        2 - Compute the right (x) axis:  Cartesian coordinates are defined by three unit vectors that are perpendicular to each other.
            If we take two vectors A and B, they can be seen as lying in a plane. And the cross product of 
            these two vectors creates a third vector C perpendicular to that plane and thus perpendicular to both A and B. We can use 
            this property to create the right vector. The idea being to use some arbitrary vector and calculate the cross vector between
            the forward vector and this arbitrary vector. 
            In this case, we know the arbitrary vector should be pointing "up", since that is perpendicular to forward.
        3 - Compute the Up (y) axis: Since we have two orthogonal vectors, the forward and right vector, computing the cross product between 
            these two vectors will just give us the missing third vector.
        4 - Row 4: replace the first three coefficients of the row with the coordinates of the from point. (according to scratchapixels)**
        
        Put simply: we want to draw a scene with a camera situated in point e (eye), the camera should be pointed to the point c (center) in such way that
         a given vector u (up) is to be vertical in the final render.
         
         ** Note on step 4: in the implementation below, we do things slightly differently. This step is a translation of the origin to 
                the point of viewer e (eye) - which is, we can recall, all we really wanted. We change the basis of our cartesian coordinate
                system such that the new origin is at the current "eye" or "viewer", which allows us to effectively move whatt we're looking
                at in the scene.

    Implementation:
    First we will also experiment with an alternative form of shading: GOURARD SHADING
    Our OBJ files also contain normal vectors ("vn ") which specify the normal vector at each vertex. So, rather than calculating a single normal vector for
    each triangle as we have done before, we can shade the model by interpolating between the normal vectors provided at each vertex.
*/

Matrix lookat(Vec3f eye, Vec3f center, Vec3f up)
{
    Vec3f z = (eye-center).normalize();
    Vec3f x = cross(up, z).normalize();
    Vec3f y = cross(z, x).normalize();

    Matrix res = Matrix::identity(4);
    for (int i=0; i<3; i++)
    {
        res[0][i] = x[i];
        res[1][i] = y[i];
        res[2][i] = z[i];
        res[i][3] = -center[i];
    }
    return res;
}




/*Similair to the texture triangle, but now we interpolate between the triangle's vertice's
  normal vectors using the barycentric coordinates, which creates a smoothed shading look.
  We can still map in the color values from the diffuse texture map, and the result looks really
  good.
*/
void gourard_triangle(Vec3i *pts, int *zbuffer, Vec2i *uvPts,  float *intensity, TGAImage &diffusemap, TGAImage &image)
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
            if (zbuffer[idx] < P.z)
            {
                zbuffer[idx] = P.z;
                
                Vec2f uv(0,0);
                float bc_intens = 0.f;
                for (int i = 0; i < 3; i++)
                {
                    uv.x += uvPts[i].x * bc[i];
                    uv.y += uvPts[i].y * bc[i];
                    bc_intens += intensity[i] * bc[i];
                }
                bc_intens = std::abs(bc_intens); // final intensity floats were < 0, this -> correct render

                TGAColor color = diffusemap.get(uv.x, uv.y);

                image.set(
                    P.x, P.y, 
                    TGAColor( *(color.r) * bc_intens, *(color.g) * bc_intens,
                              *(color.b) * bc_intens, 255)
                );
            }
        }
}


/*To establish our projection and view, we compute the transformation matrices and convert them
    into a composite transformation. Then we iterate through the faces in the model, grabbing the
    triangle coordinates and using the composite transformation to adjust them - ie, change the camera
    position (ModelView), adjust the Projection, and convert the coordinates from world to screen.
    Then we do our usual light intensity and texture mapping routine and rasterize the triangle.
*/
void gourard_draw(Model &model, int *zbuffer, TGAImage& image)
{
    TGAImage diffusemap = model.diffuse();

    Matrix ModelView = lookat(eye, center, Vec3f(0,1,0)); //eye, center, up
    Matrix ViewPort = viewport(width/8, height/8, width*3/4, height*3/4);
    Matrix Projection = projection();

    //Precompute the transformations (composite transformaton)
    Matrix z = (ViewPort*Projection*ModelView);

    for (int i = 0; i < model.nfaces(); i++)
    {
        Vec3i screen_coords[3]; 
        Vec3f world_coords[3]; 
        float intensity[3];
        for (int j = 0; j < 3; j++)  
        {
            Vec3f v = model.vert( i, j );
            screen_coords[j] = m2v(z*v2m(v));
            world_coords[j] = v;
            intensity[j] = model.normal(i, j)*light_dir;
        }
        // Compute triangle's normal and use to calculate light intensity for the current face (just for back face culling)
        Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1]-world_coords[0]); //side 1 cross-product side 2
        n.normalize();
        float bf_intensity = n*light_dir;

        // Draw triangle if not back-of-face (this is not necessary for the output in this case, but it does look slightly cleaner)
        if (bf_intensity > 0)
        {
            Vec2i uv_coords[3];
            for (int k=0; k<3; k++)
            {
                Vec2f uv = model.uv( i, k );
                uv_coords[k] = {(int)(uv.x * diffusemap.width()), (int)(uv.y * diffusemap.height())};
            }


            gourard_triangle(screen_coords, zbuffer, uv_coords, intensity, diffusemap, image);
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

    // project_and_draw(model, zbuffer, image); //perspective project
    gourard_draw(model, zbuffer, image); //perspective project + gourard shading + texture mapping

    std::cout << "Saving output.";
    image.write_tga_file("out.tga");
}