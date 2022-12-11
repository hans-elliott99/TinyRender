////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
/*                Lesson 6: Shaders for the software renderer                             */
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <iostream>

#include "../src/tgaimage.h"
#include "../src/model.h"
#include "../src/geometry.h"
#include "../src/graphics.h"

const int width  = 800;
const int height = 800;
const int depth = 255;

vec3 light_dir(1,1,1);
// vec3       eye(0,1,3);
vec3    center(0,0,0);
vec3        up(0,1,0);


/*
    Shaders:
    First off, let's note that we have created some new source code in the "graphics" files. This
    is our equivalent to the OpenGL lib, and it contains functions and objects which we use for 
    rendering our graphics (lookat, projection, viewport, triangle, etc...)
    The "shaders" are to be used for manipulation by the programmer and they interact with hard-coded graphics methods.
    In our simplified graphics pipeling, the programmer is responsible for reading in the model's data and for deciding how 
    to transform it (via the Shaders). The actual rasterization of the geometry happens with our triangle rasterization function.

        - the *vertex shader* transforms the model's coordinates and prepares the data for the fragment shader:
            First it embeds our 3d coordinates into 4d space (homoegeneous coordinates) and applies the linear transformations.
            Just as in previous lessons, we use the ModelView matrix to conver the model's coordinates to world coordinates,
            use the (perspective) projection matrix to alter the projection, and the ViewPort matrix to 'move' the camera.
            It also reads in the vertex-intensities from the vertex normal data to be used by the fragment shader - they are
            stored in a variable called varying_itensity - "varying" is a keyword in the OpenGL Shading Language which signals that
            the variable stores data which should be interpolated inside the triangle (ie, with barycentric coords)
            - the vertex shader is called for each face before the rasterizer

        - the *fragment shader* determines the color of each pixel, and determines if we should draw or discard the pixel.
            It calculates the intensity for the current pixel by inerpolating the vertex normals with the barycentric coords,
            and adjusts the colors as such.
            Additionally, it could interpolate the diffuse texture map to determine the color of each pixel, or interpolate other
            kinds of data (see further below).
            - the fragment shader is called for each pixel by the rasterizer
*/

struct GourardShader_1 : public IShader
{
    vec3     varying_intensity; //normal coords
    mat<2,3> varying_uv;        //texture coords

    Model* model_ptr;
    mat<4,4> ViewPort, Projection, ModelView, Z;

    GourardShader_1(mat<4,4> viewport, mat<4,4> projection, mat<4,4> modelview) : 
        ViewPort(viewport),
        Projection(projection),
        ModelView(modelview),
        Z(viewport*projection*modelview)
        {std::cout << Z;}

    /*Vertex Shader*/
    vec4 vertex(Model &model, int iface, int nthvert) {
        model_ptr = &model;
        vec4 vertex                = embed<4>(model.vert(iface, nthvert)); //read in 3d vertex from model and convert to homogeneous coords
        varying_intensity[nthvert] = std::max(0.f, model.normal(iface, nthvert)*light_dir); // get diffuse lighting intensity
        varying_uv.set_col(nthvert, model.uv(iface, nthvert));

        return  Z * vertex; //transformations
    }
    /*Fragment Shader*/
    bool fragment(vec3 baryc, TGAColor &color) {
        float intensity = varying_intensity*baryc;       /*dot product interpolates the intensity*/
        vec2 uv = varying_uv*baryc;                      /*interpolate uv coordinates for current pixel*/
        color = model_ptr->diffuse(uv) * intensity;
        
        return false;
    }
};

/*Using intensity 'buckets' with a constant color leads to a cool effect*/
struct ToonShader : public IShader
{
    vec3 varying_intensity;
    mat<4,4> ViewPort, Projection, ModelView, Z;

    ToonShader(mat<4,4> viewport, mat<4,4> projection, mat<4,4> modelview) : 
        ViewPort(viewport),
        Projection(projection),
        ModelView(modelview),
        Z(viewport*projection*modelview)
        {std::cout << Z;}

    /*Vertex Shader*/
    vec4 vertex(Model &model, int iface, int nthvert) {
        vec4 vertex = embed<4>(model.vert(iface, nthvert)); //read in 3d vertex from model and convert to homogeneous coords
        varying_intensity[nthvert] = std::max(0.f, model.normal(iface, nthvert)*light_dir); // get diffuse lighting intensity

        return  Z * vertex; //transformations
    }
    /*Fragment Shader*/
    bool fragment(vec3 baryc, TGAColor &color) {
        float intensity = varying_intensity*baryc;       /*dot product interpolates the intensity*/
        if      (intensity>.85) intensity =   1;
        else if (intensity>.60) intensity = .80;
        else if (intensity>.45) intensity = .60;
        else if (intensity>.30) intensity = .45;
        else if (intensity>.15) intensity = .30;
        else intensity = 0;

        color = TGAColor(15, 252, 3) * intensity;        
        return false;
    }
};


/*Phong Shading
    In Phong shading, we use the vertex normals associated with the current face (from the OBJ) and interpolate them
    across the triangle. It differs from Gourard shading in that the light calculation is computer per fragment rather
    than per vertex, as in Gourard.

*/
struct PhongShader : public IShader
{
    mat<3,3> varying_norm; 
    mat<2,3> varying_uv;        //texture coords

    Model* model_ptr;
    mat<4,4> ViewPort, Projection, ModelView, Z;

    PhongShader(mat<4,4> viewport, mat<4,4> projection, mat<4,4> modelview) : 
        ViewPort(viewport),
        Projection(projection),
        ModelView(modelview),
        Z(viewport*projection*modelview)
        {}

    /*Vertex Shader*/
    vec4 vertex(Model &model, int iface, int nthvert) {
        model_ptr = &model;
        vec4 vertex                = embed<4>(model.vert(iface, nthvert)); //read in 3d vertex from model and convert to homogeneous coords
        varying_uv.set_col(nthvert, model.uv(iface, nthvert));
        varying_norm.set_col(nthvert,
         proj<3>( (Projection*ModelView).invert_transpose() * embed<4>(model.normal(iface, nthvert)) ) 
        );
        return  Z * vertex; //transformations
    }
    /*Fragment Shader*/
    bool fragment(vec3 baryc, TGAColor &color) {
        vec3 bn = (varying_norm*baryc).normalize();
        vec2 uv = varying_uv*baryc; 

        float diff = std::max(0.f, bn*light_dir);       /*dot product interpolates the intensity*/
        color = model_ptr->diffuse(uv) * diff;
        
        return false;
    }
};



/*SIMPLE NORMAL MAPPING
    (more on normal mapping below)
    
    We can use our varying texture coordinates to pull in other information too. 
    For example, if we have a normal mapping texture ('african_head_nm.tga' or 'african_head_nm_tangent'.tga) we can compute light
    intensity by sampling from this map rather than interpolating across the triange (before we had to interpolate since we only had vector normals).

    First, note that we store uniform_M (Projection matrix X ModelView matrix) and uniform_Mit (the inverse transpose of that).
    In OpenGL Shading Language, 'uniform' is a keyword which denotes that the variable stores a constant value.
    We use these to transform the sample normal vector and light direction so they match our model's new coord system  
    
    *Why use the inverse transpose for the texture normals?* - that was explained at the end of lesson 5:
        If a model's coordinates are transformed by an affine mapping M, then its normal vectors must be transformed with a
        mapping equal to the inverse tranpose of M.
    Essentially, since the 'normals' we use for light intensity are based on the surface normals of the faces (triangles) that
     make up our model, when we rescale our model's coordinates, the new faces (triangles) need to have their surface normals 
     re-computed.
*/


struct GourardShader : public IShader
{
    mat<2,3> varying_uv;        //texture coords
    mat<4,4> uniform_M;         //Projection*ModelView)
    mat<4,4> uniform_Mit;       //(Projection*ModelView).invert_transpose()
    mat<4,4> transform;         //ViewPort*Projection*ModelView

    Model*   model_ptr;


    GourardShader(mat<4,4> viewport, mat<4,4> projection, mat<4,4> modelview) : 
        uniform_M  ( projection*modelview ),
        uniform_Mit( (projection*modelview).invert_transpose() ),
        transform  ( viewport*projection*modelview )
        {}

    /*Vertex Shader*/
    vec4 vertex(Model &model, int iface, int nthvert) {
        model_ptr = &model;
        vec4 vertex                = embed<4>(model.vert(iface, nthvert)); //read in 3d vertex from model and convert to homogeneous coords
        varying_uv.set_col(nthvert, model.uv(iface, nthvert));

        return  transform * vertex; 
    }

    /*Fragment Shader*/
    bool fragment(vec3 baryc, TGAColor &color) {
        vec2 uv = varying_uv*baryc;                      /*interpolate uv coordinates for current pixel*/

        /*Use the uv coords (texture coords) to sample from the normal texture map at position uv.x, uv.y */
        /*We also transform the uv coordinates and light direction based on the ModelView & Projection matrices*/
        vec3 norm  = proj<3>( uniform_Mit*embed<4>(model_ptr->normal(uv)) ).normalize();
        vec3 light = proj<3>( uniform_M*embed<4>(light_dir) ).normalize();
        float intensity = std::max(0.f, norm*light);

        /*Finally, sample the color from the diffuse texture map and alter by intensity level*/    
        color = model_ptr->diffuse(uv) * intensity;
        return false;
    }
};

/*SPECULAR MAPPING

    We can use the Phong reflection model to make the reflection of light look more realistic. Bui Tuong Phong realized that shiny surfaces
    have small, intense specular highlights while dull surfaces have large highlights which fade out more gradually. We can approximate this
    effect by considering final lighting as a weighted sum of ambient lighting (which is constant per scene), diffuse lighting, and specular
    lighting.
    For diffuse lighting, we've been computing light intensity as the cosine (dot product) of the light vector (l) and the surface normal vector (n).
    For specular lighting, we're interested in the cosine of the angle between vectors r (reflected light direction) and v (view direction).
    If this relationship reveals that we *can* see the light reflected by this pixel, the pixel is illuminated.
    We can calculate the 'reflected light direction', r, for each pixel as r = 2n<n,l>-l (where all vars are normalized).
    But we also want to model the effect where glossy surfaces tend to reflect light in one direction much more than in others. Thus, we can
    raise 'r' to a power. That power is actually extracted from a specular texture map, the data of which effectively tells us how glossy the
    given point is.

    So our vertex shader hasnt changed, but for each pixel our fragment shader now:
    1. interpolates the uv coordinates
    2. uses the uv coordinates to sample from the normal texture map (and transforms the normal and light vectors)
    3. calculates the reflection direction based on the normal and light vectors
    4. computes specular lighting by raising the reflection direction to a power sampled from the 
        specular texture map (using the uv coords)
    5. compute diffuse lighting by multiplying the normal vector by the light direction
    6. compute the final RGB pixel color by using a weighted sum of ambient light, diffuse light, and specular light, applied
        to the color sampled from the diffuse texture map.

    Finally, a note on the coefficients used in the Phong Reflection Model:
        Below, a coefficient of 0.5 is used for ambient light, and 1.0 is applied to the diffuse lighting and 0.9 is appled to the specular.
        Technically, these coefficients are meant to sum to 1.
*/
struct GourardSpecular : public IShader
{
    mat<2,3> varying_uv;        //texture coords
    mat<4,4> uniform_M;         //Projection*ModelView)
    mat<4,4> uniform_Mit;       //(Projection*ModelView).invert_transpose()
    mat<4,4> transform;         //ViewPort*Projection*ModelView

    Model*   model_ptr;


    GourardSpecular(mat<4,4> viewport, mat<4,4> projection, mat<4,4> modelview) : 
        uniform_M  ( projection*modelview ),
        uniform_Mit( (projection*modelview).invert_transpose() ),
        transform  ( viewport*projection*modelview )
        {}

    /*Vertex Shader*/
    vec4 vertex(Model &model, int iface, int nthvert) {
        model_ptr = &model;
        vec4 vertex                = embed<4>(model.vert(iface, nthvert));
        varying_uv.set_col(nthvert, model.uv(iface, nthvert));

        return  transform * vertex; 
    }

    /*Fragment Shader*/
    bool fragment(vec3 baryc, TGAColor &color) {
        vec2 uv = varying_uv*baryc;                     

        vec3 norm  = proj<3>( uniform_Mit*embed<4>(model_ptr->normal(uv)) ).normalize();
        vec3 light = proj<3>( uniform_M*embed<4>(light_dir) ).normalize();
        vec3 refl  = (2.f*norm*(norm*light) - light).normalize();
        
        float specular = pow( std::max(refl.z, 0.f), model_ptr->specular(uv) );
        float diffuse  = std::max(0.f, norm*light);

        TGAColor c = model_ptr->diffuse(uv);
        for (int i: {0,1,2})
            color[i] = std::min<float>( 0.5 + c[i]*(diffuse + .9*specular), 255 );

        return false;
    }
};


// 6b:



/*Can swap in and shader*/
void render_pipeline(Model &model, TGAImage &image, TGAImage &zbuffer)
{
    float cam[3];
    std::cout << "Camera position (try: 0 1 3): ";
    for (int i : {0,1,2})
        std::cin >> cam[i];
    vec3 eye(cam[0], cam[1], cam[2]);

    mat<4,4> modelview = lookat_mat(eye, center, up);
    mat<4,4> viewport = viewport_mat(width/8, height/8, width*3/4, height*3/4);
    mat<4,4> proj = projection_mat( -1.f / (eye-center).norm() ); 
    light_dir.normalize();

    GourardSpecular shader(viewport, proj, modelview);

    for (int i=0; i<model.nfaces(); i++) {
        vec4 screen_coords[3];
        for (int j : {0, 1, 2}) {
            screen_coords[j] = shader.vertex(model, i, j);
        }
        triangle(screen_coords, shader, image, zbuffer);
    }
}






int main(int argc, char** argv)
{
    TGAImage image(width, height, TGAImage::RGB);
    TGAImage zbuffer(width, height, TGAImage::GRAYSCALE);

    Model model = Model("../obj/african_head.obj");
    if (argc > 1)
        Model model = Model(argv[1]);

    render_pipeline(model, image, zbuffer);
    std::cout <<"saving\n";
    image.write_tga_file("out.tga");
    zbuffer.write_tga_file("zbuf.tga");
}
