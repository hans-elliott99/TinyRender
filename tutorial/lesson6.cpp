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

const int width  = 1200;
const int height = 1200;
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
    
    To clarify, we sample a color from our normal texture map, but the RGB components correspond to the X, Y, and Z coordinates,
    respectively, of the surface normal for the given pixel.
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


/*TANGENT SPACE NORMAL MAPPING
    We can use a tangent space (aka, Darboux or surface space) normal map instead of an RGB one. A major advantage of using
    a tangent space map is that it would deform correctly as our model mutated (for example,
    if we animated it), where as we would need a a new RGB normal map for every frame.  
    Wikipedia (normal mapping):
    - The orientation of coordinate axes differs depending on the space in which the normal map was encoded. 
      A straightforward implementation encodes normals in object-space, so that red, green, and blue components 
      correspond directly with X, Y, and Z coordinates. In object-space the coordinate system is constant.
    - However object-space normal maps cannot be easily reused on multiple models, as the orientation of the surfaces differ. 
      Since color texture maps can be reused freely, and normal maps tend to correspond with a particular texture map, 
      it is desirable for artists that normal maps have the same property.
    - Normal map reuse is made possible by encoding maps in tangent space. 
      The tangent space is a vector space which is tangent to the model's surface. 
      The coordinate system varies smoothly (based on the derivatives of position with respect to texture coordinates)
      across the surface.

    It also boost performance by separating changes relative to the object from changes relative to the World Space:
    - ex: rotating an object in tangent space leaves the relationship between vertices the same, while rotating the
        object in World space means all pixels need to be recalculated to take in changes for lighting, shading, etc.
    Video lecture: https://www.youtube.com/watch?v=19dey0OYPYo&ab_channel=MichaelLoeser
    How to use it:
        1. Calculate normal vector (per vertex) : which is just the "up" vector, perpendicular to the surface (always positive)
        2. Calculate tangent vector : can be chosen as desired, but must be parallel to the surface and perpendicular to the normal vector
        3. Calculate bitangent (or, bi-normal) vector : perpendicular to both the normal and tangent vector.
        4. Using the Tangent (T), Bitangent (B), and Normal (N) vectors, construct the 3x3 translation matrix:
            M_tan = [[Tx, Bx, Nx], [Ty, By, Ny], [Tz, Bz, Nz]]
            - so to move a vertex to tangent space: Position_tan = Position_world X M_tan

        5. To return to world space, we use the inverse Tangent Space matrix, but since each vector is perpendicular (orthogonal) to the other,
           this is equivalent to the transpose of the matrix (easier to calculate)
           M_world = [[Tx, Ty, Tz], [Bx, By, Bz], [Nx, Ny, Nz]]
    The normal maps are generally blue because their RGB values at each pixel make up the normal vector (x,y,z). 
    Normal vectors have a positive z value, which is the B color value. 

        
    To use the tangent space normal map we sample a normal vector from the texture map and
    transform its coordinates from Darboux to the global system coordinates (ie, world coords).  

    In the shader below, we compute a Darboux basis as a triplet of vectors (i,j,n) where n is the normal vector
    and i and j are computed as below. (See tutorial for matching equations/matrices: 
        https://github.com/ssloy/tinyrenderer/wiki/Lesson-6bis:-tangent-space-normal-mapping)   
    Once we have the Darboux basis packaged into the change-of-basis matrix B, we read in the normal
    from the texture map and multiply with B to change the basis from Darboux space to the world/global coordinates.
*/

struct TangentNormalShader : public IShader
{
    mat<2,3> varying_uv;        //texture coords
    mat<4,3> varying_tri;       //triangle coords
    mat<3,3> varying_norm;      //vertex normals to be inerpolated
    mat<3,3> ndc_tri;           //triangle in normalized device coordinates

    mat<4,4> uniform_M;         //Projection*ModelView)
    mat<4,4> uniform_Mit;       //(Projection*ModelView).invert_transpose()
    mat<4,4> ViewPort;          //ViewPort*Projection*ModelView

    Model*   model_ptr;

    TangentNormalShader(mat<4,4> viewport, mat<4,4> projection, mat<4,4> modelview) : 
        uniform_M  ( projection*modelview ),
        uniform_Mit( (projection*modelview).invert_transpose() ),
        ViewPort   ( viewport )
        { }

    /*Vertex Shader*/
    vec4 vertex(Model &model, int iface, int nthvert) {
        model_ptr = &model;
        varying_uv    .set_col(nthvert, model.uv(iface, nthvert));
        varying_norm  .set_col( nthvert, proj<3>(uniform_Mit * embed<4>(model.normal(iface, nthvert), 0.f)) );
        
        vec4 m_vertex  = uniform_M * embed<4>(model.vert(iface, nthvert));
        varying_tri   .set_col(nthvert, m_vertex);
        ndc_tri       .set_col(nthvert, proj<3>(m_vertex/m_vertex[3]));

        return ViewPort*m_vertex;
    }

    /*Fragment Shader*/
    bool fragment(vec3 baryc, TGAColor &color) {
        vec2 uv = varying_uv*baryc;                      /*interpolate uv coordinates for current pixel*/
        vec3 bn = (varying_norm*baryc).normalize();

        mat<3,3> A;
        A[0] = ndc_tri.col(1) - ndc_tri.col(0);
        A[1] = ndc_tri.col(2) - ndc_tri.col(0);
        A[2] = bn;

        mat<3,3> Ai = A.invert();
        vec3 i = Ai * vec3(varying_uv[0][1] - varying_uv[0][0], varying_uv[0][2] - varying_uv[0][0], 0);
        vec3 j = Ai * vec3(varying_uv[1][1] - varying_uv[1][0], varying_uv[1][2] - varying_uv[1][0], 0);

        mat<3,3> B;
        B.set_col(0, i.normalize());
        B.set_col(1, j.normalize());
        B.set_col(2, bn);

        vec3 norm  = (B * model_ptr->normal(uv)).normalize();
        // vec3 light = proj<3>(uniform_M * embed<4>(light_dir, 0.f)).normalize();  //what tutorial does
        vec3 light = (B * light_dir).normalize(); //tutorial doesnt translate light with B, but my result looks better w it

        float diffuse  = std::max(0.f, norm*light);

        // (For specular lighting, set to true - can test with and without)
        bool spec = false;
        if (spec) {
            vec3 refl  = (2.f*norm*(norm*light) - light).normalize(); //for specular lighting
            float specular = pow( std::max(refl.z, 0.f), model_ptr->specular(uv) );
            TGAColor c = model_ptr->diffuse(uv);
            for (int i: {0,1,2})
                color[i] = std::min<float>( 0.5 + c[i]*(diffuse + .3*specular), 255 );
        } else 
            color = model_ptr->diffuse(uv)*diffuse;
        return false;
    }
};




/*Can swap in any shader*/
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

    TangentNormalShader shader(viewport, proj, modelview);

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
