#include <iostream>
#include <vector>
#include "../tga/tgaimage.h"
#include "../model/geometry.h"
#include "../model/model.h"

/*
    Lesson 4-5 Reworked:
    Refactored/simplified geometry code (vectors simply store floats).
    Fixed some bugs.  
    Overall pipeline:
        Read in model's faces, vertex coordinates ("world" coords), vertex normals (for shading), and texture (uv) coordinates.
        For each face (triangle), apply transformations to the vertex coords to: 
            A) move the camera (ie, change the ModelView, ie, lookat)
            B) project the 3d coordinates to 2d (ie, perspective projection)
            C) convert from model's "world" coordinates to screen coordinates (ie, viewport).
        Pass all this to the triangle rasterizer which:
            A) computes a bounding box given the screen coordinates
            B) determines if each point in the BB is in the triangle based on barycentric coords
            C) if in the triangle, the pixel can be drawn if its z coord is greater than its z-buffer value:
                - calculate the z/depth coordinate using the barycentric coords & compare to current z-buff value
            D) if drawing, use the uv coordinates to pull in the color from the texture map (diffusemap)
                - we have texture coordinates which correspond to the triangle's vertices, so
                - interpolate across the triangle using the barycentric coords.
            C) finally, use Gourard shading to modify the texture-map color for shading
                - we have normal coordinates which correspond to the triangle's vertices, so
                - interpolate across the triangle using the barycentric coordinates.
*/

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
const int width  = 800;
const int height = 800;
const int depth = 255;

vec3 light_dir = vec3(0,0,-1).normalize();
vec3 eye(1,1,3);
vec3 center(0,0,0);


vec3 world2screen(vec3 v)
{
    return vec3(
        int((v.x + 1.) * width/2. + .5),
        int((v.y + 1.) * height/2. + .5),
        v.z
    );
}

// Convert vector to matrix, with homogeneous coords
mat<4,1> v2m(vec3 v)
{
    mat<4,1> m;
    m[0][0] = v.x;
    m[1][0] = v.y;
    m[2][0] = v.z;
    m[3][0] = 1.f;
    return m;
}

// Convert from matrix, with homog coords, back to 3d vec
vec3 m2v(mat<4,1> m)
{
    return vec3(m[0][0]/m[3][0], m[1][0]/m[3][0], m[2][0]/m[3][0]);
}


// Viewport: map model's "world coords." to the screen coords
mat<4,4> viewport(int x, int y, int w, int h)
{
    mat<4,4> m = mat<4,4>().identity();
    m[0][3] = x+w/2.f;
    m[1][3] = y+h/2.f;
    m[2][3] = depth/2.f;

    m[0][0] = w/2.f;
    m[1][1] = h/2.f;
    m[2][2] = depth/2.f;
    return m;
}

// Lookat matrix : change where the camera is looking
mat<4,4> lookat(vec3 eye, vec3 center, vec3 up)
{
    vec3 z = (eye-center).normalize();
    vec3 x = cross(up, z).normalize();
    vec3 y = cross(z, x).normalize();

    mat<4,4> res = mat<4,4>().identity();
    for (int i=0; i<3; i++)
    {
        res[0][i] = x[i];
        res[1][i] = y[i];
        res[2][i] = z[i];
        res[i][3] = -center[i];
    }
    return res;
}



// Projection matrix - convert to perspective projection 
mat<4,4> projection(vec3 eye, vec3 center)
{
    mat<4,4> m = mat<4,4>().identity(); 
    // m[3][2] = -1.f / eye.z;
    m[3][2] = -1.f / (eye-center).norm();
    return m;
}



vec3 barycentric(vec3 *pts, vec2 P)
{
    /*Compute cross-product to calculate barycentric coords*/
    vec3 u = cross(vec3(pts[2].x-pts[0].x, pts[1].x-pts[0].x, pts[0].x-P.x), //Cx-Ax, Bx-Ax, Ax-Px
                   vec3(pts[2].y - pts[0].y, pts[1].y - pts[0].y, pts[0].y - P.y)); //Cy-Ay, By-Ay, Ay-Py 
    
    /* 'pts' and 'P' are integer vectors, so if abs(u.z) < 0, then u.z == 0, and the current point
        'P' is not in the triangle. Thus, we should return a vector with some negative coords.*/
    if (std::abs(u.z) < 1) return vec3(-1, -1, -1);

    /*Otherwise, return the vector of "weights" (equivalent to 1-u-v, u, v)*/
    return vec3( 1.0F - (u.x+u.y) / u.z,
                               u.y / u.z,
                               u.x / u.z);
}


/*Now to draw a triangle, we first compute a bounding box.
    The bounding box is described by two points: the bottom left and upper right vertices.
    To find these corners, we iterate through the vertices of the triangle and choose min & max
    coords. 
*/

void triangle(vec3 *pts, vec2 *uvPts, int *zbuffer, TGAImage &diffusemap, TGAImage &image, float *intensity)
{
    vec2 bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max() );
    vec2 bboxmax( -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max() );
    vec2 clamp(image.width() - 1, image.height() - 1); //clip bb to image dims if triangle extends outside image
    // Compute bounding box
    for (int i = 0; i < 3; i++) //for triangle vertex i
        for (int j = 0; j < 2; j++) //for coord (of x,y) j
        {
            bboxmin[j] = std::max(0.0f,     std::min(bboxmin[j], pts[i][j]));
            bboxmax[j] = std::max(clamp[j], std::max(bboxmax[j], pts[i][j]));
        }
    vec3 P;
    // For all points within the bounding box, computer the barycentric
    // coordinates. If any are less than 0, the current point is not in
    // the triangle, and we should skip. Otherwise, determine if the
    // point's z dimension (depth) is closer to camera then the z buffer,
    // and only render the point if it is.
    for (int x = (int)bboxmin.x; x <= (int)bboxmax.x; x++)     /*cast coords to int or else we will miss faces/have holes*/
        for (int y = (int)bboxmin.y; y <= (int)bboxmax.y; y++)
        {
            vec3 bc = barycentric(pts, vec2(x,y));
            if (bc.x < 0 || bc.y < 0 || bc.z < 0) continue;
            int z = 0;
            for (int i = 0; i < 3; i++)
                z += pts[i].z * bc[i];

            int idx = x + y*width;
            if ( zbuffer[idx] < z )
            {
                zbuffer[idx] = z;
                
                vec2 uv;
                float bc_intens = 0.f;
                for (int i = 0; i < 3; i++)
                {
                    uv.x += uvPts[i].x * bc[i];
                    uv.y += uvPts[i].y * bc[i];
                    bc_intens += intensity[i] * bc[i];
                }
                bc_intens = std::abs(bc_intens);
                TGAColor color = diffusemap.get(uv.x, uv.y) * bc_intens;

                image.set(x, y, color);
            }
        }
}


/*Gourard Shading, Texture Mapping*/
void render_model(Model &model, TGAImage &image)
{
    TGAImage diffusemap = model.diffuse();
    /*Initialize zbuffer with largest possible negative depth (so we
        start rendering infront of the zbuffer*/
    int *zbuffer = new int[width*height];
    for (int i=width*height; i--; zbuffer[i] = -std::numeric_limits<int>::max());

    /*Transformation*/
    mat<4,4> ModelView = lookat(eye, center, /*up=*/vec3(0,1,0));
    mat<4,4> ViewPort = viewport(width/8, height/8, width*3/4, height*3/4);
    mat<4,4> Projection = projection(eye, center);

    mat<4,4> z = (ViewPort*Projection*ModelView); 
    std::cout << z; //transformation matrix

    for (int i = 0; i < model.nfaces(); i++)
    {
        vec3 screen_coords[3]; 
        vec3 world_coords[3]; 
        float intensity[3];
        for (int j = 0; j < 3; j++)  
        {
            vec3 v = model.vert( i, j );
            screen_coords[j] = m2v(z*v2m(v)); //scale "world" coords to fit screen
            world_coords[j] = v;
            intensity[j] = model.normal(i, j) * light_dir;
        }
        // We don't want to use this BF culling method anymore, doesn't work as expected
        // vec3 n = cross(world_coords[2] - world_coords[0], world_coords[1]-world_coords[0]); //side 1 cross-product side 2
        // n.normalize();
        // float bf_intensity = n*light_dir;
        // if (bf_intensity > 0)
        {
            vec2 uv_coords[3];
            for (int k=0; k<3; k++)
            {
                vec2 uv = model.uv( i, k );
                uv_coords[k] = {uv.x * diffusemap.width(), uv.y * diffusemap.height()};
            }

            triangle(screen_coords, uv_coords, zbuffer, diffusemap, image, intensity);
        }
    }
}



int main(int argc, char** argv)
{
    TGAImage image(width, height, TGAImage::RGB);

    Model model = Model("../obj/african_head.obj");
    if (argc > 1)
        Model model = Model(argv[1]);

    render_model(model, image);
    std::cout <<"saving\n";
    image.write_tga_file("out.tga");
}


/*Notes: */
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
/*                Lesson 4: Rendering in Perspective (Projection)                             */
/*                Lesson 5: Moving the Camera                                                 */
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

/*
Summary: Our models are created in their own local frame (object coordinates). 
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

/*  VIEWPORT:
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




/*  TEXTURE MAPPING:
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
