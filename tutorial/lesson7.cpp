#include <iostream>

#include "../src/tgaimage.h"
#include "../src/model.h"
#include "../src/geometry.h"
#include "../src/graphics.h"

const int width  = 1200;
const int height = 1200;

vec3 light_dir(1,1,1);
// vec3       eye(0,1,3);
vec3    center(0,0,0);
vec3        up(0,1,0);


std::vector<float> shadowbuffer(width*height, -std::numeric_limits<float>::max());

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
/*                                  Lesson 7: Shadow Mapping                                  */
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////


/*
Summary:
    *Rendering hard shadows.*
    Up until now shading has been local - we either shaded a particular face or a particular pixel
    based on light direction and the normal vector at the vertices of the face (or at the pixel).
    Now, we will do 2 passes of rendering: 
    - first, we do the render while placing the camera (eye) at the position of the light source,
      so we can determine what parts are lit vs hidden by the light.  
    - then we do our typical render while taking into account the visibility informaition from the first pass.
      
    In the first pass (DepthShader below) we simply retrieve the model's coordinates, transform them (with the lookat matrix
    adjusted so that the camera, or eye, is the light direction), and then use the fragment shader to determine the degree to 
    which each pixel can be "seen" by the light source.  
    Most importantly, the "zbuffer" for this pass, which we call the shadowbuffer, is saved for use by the next pass of rendering.
    Since this now encodes data about which pixels are visible from the light source, it will allow us to determine whether a pixel
    should be in shadow or not.    

    **Also note (in render_pipeline) we save the transformation matrix which consists of the ViewPort*Projection*ModelView where
    the model-view is based from the light source. This is saved in the MDepth variable which is passed into the Shader and used in
    the creation of a uniform variable, MShadow.

    In the second pass (Shader below) we reuse a basic shader from the last lesson (with diffuse and specular lighting) and 
    add just a bit of code to incorporate the shading.  
    First, the vertext shader remains unchaged - receiving the current vertex and using transformation matrices to convert
    the model to screen coordinates.
    In the fragment shader, we calculate the "shadowbuffer points" by converting our screen coordinates (varying_tri)
    to shadowbuffer coordinates. The specific matrix transformation used to do this is MDepth (as described above)
    times the inverse of the standard transformations (M) for the actual render - (Viewport*Projection*ModelView).invert() or M.invert().  
    Why? 
    Well matrix Mdepth is a transformation from object space to the shadowbuffer space, and matrix M
    is a transformation from from object coordinates to screen coordinates.  
    Thus, the M.inverse() converts from screen coordinates to object coordinates, and MDepth * M.inverse() converts from
    screen coordinates to shadowbuffer coordinates.

    Then, we can use the x and y shadowbuffer coordinates to index into the shadowbuffer array and encode the 'shadow' coefficient
    based on whether or not the shadowbuffer value at that point is less than the z (or depth) dimension of the current pixel (in SB coordinates).
    In other words, "to determine whether the current pixel is lit or no it suffices to compare its z-coordinate with the value we stored in the shadow buffer."

    There is one issue: "z-fighting", which wikipedia explains well:
    "[Z-fighting] occurs when two or more primitives have very similar distances to the camera. This would cause them to have near-similar or identical values
    in the z-buffer, which keeps track of depth. This then means that when a specific pixel is being rendered, it is ambiguous which one of the two primitives are
    drawn in that pixel because the z-buffer cannot distinguish precisely which one is farther from the other. If one pixel was unambiguously closer, the less close
    one could be discarded."
    It is easy to tell when this is happening because the render has a terrible zerbra-like pattern across the model's faces.
    One way this can be fixed is by moving the primitives (triangles in our case) apart, which we achieve by adding a constant value to the z-coordinate of our
    shadowbuffer point.  
    

*/


struct DepthShader : public IShader
{
    mat<3,3> varying_tri;

    mat<4,4> ViewPort, Projection, ModelView, Z;

    DepthShader(mat<4,4> viewport, mat<4,4> projection, mat<4,4> modelview) : 
        ViewPort(viewport),
        Projection(projection),
        ModelView(modelview),
        Z(viewport*projection*modelview),
        varying_tri()
        {}

    vec4 vertex(Model& model, int iface, int nthvert) {
        vec4 vertex = embed<4>( model.vert(iface, nthvert) );
        vertex = Z * vertex; 
        varying_tri.set_col(nthvert, proj<3>(vertex/vertex[3]) );

        return vertex;
    }

    bool fragment(vec3 baryc, TGAColor &color) {
        vec3 p = varying_tri * baryc;
        color = TGAColor(255, 255, 255) * (p.z/g_DEPTH); /*g_DEPTH=255 (graphics.h)*/
        return false;
    }
};



/*Reuse our final shader from lesson 6, but add in shadows*/
struct Shader : public IShader
{
    mat<2,3> varying_uv;        //texture coords
    mat<3,3> varying_tri;       //triangle coords

    mat<4,4> uniform_M;         //Projection*ModelView
    mat<4,4> uniform_Mit;       //(Projection*ModelView).invert_transpose()
    mat<4,4> uniform_MShadow;   //MShadow*(Viewport*Projection*ModelView).invert()
    mat<4,4> ViewPort;          //ViewPort*Projection*ModelView

    Model*   model_ptr;

    Shader(mat<4,4> viewport, mat<4,4> projection, mat<4,4> modelview, mat<4,4> shadow_transform) : 
        uniform_M      ( projection*modelview ),
        uniform_Mit    ( (projection*modelview).invert_transpose() ),
        uniform_MShadow( shadow_transform * (viewport*projection*modelview).invert() ),
        ViewPort       ( viewport )
        { }

    /*Vertex Shader*/
    vec4 vertex(Model &model, int iface, int nthvert) {
        model_ptr = &model;
        varying_uv.set_col(nthvert, model.uv(iface, nthvert));
        vec4 m_vertex  = ViewPort * uniform_M * embed<4>(model.vert(iface, nthvert));
        varying_tri.set_col(nthvert, proj<3>(m_vertex / m_vertex[3]));

        return m_vertex;
    }

    /*Fragment Shader*/
    bool fragment(vec3 baryc, TGAColor &color) {
        /* NEW: Shadows */
        vec4 sb_p = uniform_MShadow * embed<4>(varying_tri*baryc);     //get shadowbuffer points
        sb_p = sb_p / sb_p[3];                                         //convert back to 3d (perspective project)
        int idx = int(sb_p[0]) + int(sb_p[1]) * width;                 //idx into the shadowbuffer array
        float shadow = 0.3 + 0.7 *(shadowbuffer[idx] < sb_p[2]+40); //use coefs to avoid z-fighting
        
        /* old */
        vec2 uv    = varying_uv*baryc;                      /*interpolate uv coordinates for current pixel*/
        vec3 norm  = proj<3>(uniform_Mit * embed<4>(model_ptr->normal(uv))).normalize();
        vec3 light = proj<3>(uniform_M * embed<4>(light_dir)).normalize();
        
        float diffuse  = std::max(0.f, norm*light);

        // (For specular lighting, set to true - can test with and without)
        bool spec = true;
        if (spec) {
            vec3 refl  = (2.f*norm*(norm*light) - light).normalize(); //for specular lighting
            float specular = pow( std::max(refl.z, 0.f), model_ptr->specular(uv) );
            TGAColor c = model_ptr->diffuse(uv);

        /* NEW : add shadow to final color calculation. */
            for (int i: {0,1,2})
                color[i] = std::min<float>( 0.5 + c[i]*shadow*(diffuse + .9*specular), 255 ); //include shadow value, determines if pixel is in light or not.

        } else 
            color = model_ptr->diffuse(uv) * diffuse;
        return false;
    }
};




/*Can swap in any shader*/
void render_pipeline(Model &model)
{
    light_dir.normalize();
    mat<4,4> viewport = viewport_mat(width/8, height/8, width*3/4, height*3/4);

    float cam[3];
    std::cout << "Camera position (try: 0 1 3): ";
    for (int i : {0,1,2})
        std::cin >> cam[i];
    vec3 eye(cam[0], cam[1], cam[2]);


    /*FIRST PASS: render shadows*/
    TGAImage depth(width, height, TGAImage::RGB);
    // shadowbuffer initialized above    
    mat<4,4> depth_lookat = lookat_mat(light_dir, center, up);
    mat<4,4> depth_proj   = projection_mat(0);


    DepthShader depthshader(viewport, depth_proj, depth_lookat);

    for (int i=0; i<model.nfaces(); i++) {
        vec4 screen_coords[3];
        for (int j : {0, 1, 2}) {
            screen_coords[j] = depthshader.vertex(model, i, j);
        }
        triangle(screen_coords, depthshader, depth, shadowbuffer);
    }

    std::cout << "saving depth\n";
    depth.write_tga_file("depth.tga");

    /*compute matrix to use in conversion of screen coordinates to shadowbuffer index*/
    mat<4,4> MDepth = viewport*depth_proj*depth_lookat;


    /*SECOND PASS: render the rest of the model*/
    TGAImage image(width, height, TGAImage::RGB);
    std::vector<float> zbuffer(width*height, -std::numeric_limits<float>::max());
    mat<4,4> modelview    = lookat_mat(eye, center, up);
    mat<4,4> proj         = projection_mat( -1.f / (eye-center).norm() ); 

    /*pass in the Mshadow matrix*/
    Shader shader(viewport, proj, modelview, MDepth);

    for (int i=0; i<model.nfaces(); i++) {
        vec4 screen_coords[3];
        for (int j : {0, 1, 2}) {
            screen_coords[j] = shader.vertex(model, i, j);
        }
        triangle(screen_coords, shader, image, zbuffer);
    }

    std::cout << "saving render\n";
    image.write_tga_file("out.tga");

}



int main(int argc, char** argv)
{

    Model model = Model("../obj/african_head.obj");
    if (argc > 1)
        Model model = Model(argv[1]);

    render_pipeline(model);
}
