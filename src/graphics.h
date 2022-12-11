#pragma once

#include "tgaimage.h"
#include "geometry.h"

const int g_DEPTH = 255;

struct IShader
{
    virtual ~IShader() {};
    virtual vec4 vertex(Model &model, int iface, int nthvert) = 0;

    /*Return bool so the fragment shader can determine if we draw a pixel*/
    virtual bool fragment(vec3 baryc, TGAColor &color) = 0;

};


mat<4,4> viewport_mat(int x, int y, int w, int h)
{
    mat<4,4> m = mat<4,4>().identity();
    m[0][3] = x+w/2.f;
    m[1][3] = y+h/2.f;
    m[2][3] = g_DEPTH/2.f;

    m[0][0] = w/2.f;
    m[1][1] = h/2.f;
    m[2][2] = g_DEPTH/2.f;
    return m;
};

mat<4,4> projection_mat(float coeff /*coeff = -1/c*/)
{
    mat<4,4> m = mat<4,4>().identity(); 
    m[3][2] = coeff; //-1.0f / (eye-center).norm();
    return m;
};


mat<4,4> lookat_mat(vec3 eye, vec3 center, vec3 up)
{
    vec3 z = (eye-center).normalize();
    vec3 x = cross(up, z).normalize();
    vec3 y = cross(z, x).normalize();

    mat<4,4> ModelView = mat<4,4>().identity();
    for (int i=0; i<3; i++)
    {
        ModelView[0][i] = x[i];
        ModelView[1][i] = y[i];
        ModelView[2][i] = z[i];
        ModelView[i][3] = -center[i];
    }
    return ModelView;
};


vec3 barycentric(vec2 vA, vec2 vB, vec2 vC, vec2 P)
{
    /*Compute cross-product to calculate barycentric coords*/
    // vec3 u = cross(vec3(vC.x - vA.x, vB.x - vA.x, vA.x - P.x),
    //                vec3(vC.y - vA.y, vB.y - vA.y, vA.y - P.y));
    vec3 s[2];
    for (int i : {0, 1}) {
        s[i][0] = vC[i]-vA[i];
        s[i][1] = vB[i]-vA[i];
        s[i][2] = vA[i]-P[i];
    }
    vec3 u = cross(s[0], s[1]);
    
    // u.z is stored as a float but is truly an integer, so if it is < 1 it is
    // equivalent to 0, and thus the current point P is not in triangle ABC 
    /*Otherwise, return the vector of "weights" (equivalent to 1-u-v, u, v)*/
    if (std::abs(u[2]) < 1) return vec3(-1,-1,-1);

    return vec3(  1.0f - (u.x+u.y) / u.z,
                                u.y / u.z,
                                u.x / u.z);
};


void triangle(vec4 *pts, IShader &shader, TGAImage &image, TGAImage &zbuffer)
{
    vec3 pts3d[3];
    for (int i: {0,1,2})
        pts3d[i] = proj<3>(pts[i]/pts[i][3]);

    vec2 bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max() );
    vec2 bboxmax( -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max() );
    vec2 clamp(image.width() - 1, image.height() - 1); //clip bb to image dims if triangle extends outside image
    // Compute bounding box
    for (int i = 0; i<3; i++) //for triangle vertex i
        for (int j = 0; j<2; j++) //for coord (of x,y) j
        {
            bboxmin[j] = std::min(bboxmin[j], pts3d[i][j]);
            bboxmax[j] = std::max(bboxmax[j], pts3d[i][j]);
        }
    // For all points within the bounding box, computer the barycentric
    // coordinates. If any are less than 0, the current point is not in
    // the triangle, and we should skip. Otherwise, determine if the
    // point's z dimension (depth) is closer to camera then the z buffer,
    // and only render the point if it is.
    TGAColor color;
    for (int x = (int)bboxmin.x; x <= (int)bboxmax.x; x++)     /*cast coords to int or else we will miss faces/have holes*/
        for (int y = (int)bboxmin.y; y <= (int)bboxmax.y; y++)
        {
            vec3 bc = barycentric(proj<2>(pts3d[0]), proj<2>(pts3d[1]), proj<2>(pts3d[2]),
                                  vec2(x,y));
            float z = pts[0][2] * bc.x + pts[1][2] * bc.y + pts[2][2] * bc.z;
            float w = pts[0][3] * bc.x + pts[1][3] * bc.y + pts[2][3] * bc.z;
            int frag_depth = std::max( 0, std::min(g_DEPTH, int(z / w+0.5)) );
            
            if (bc.x < 0 || bc.y < 0 || bc.z < 0 || zbuffer.get(x, y)[0]>frag_depth) continue;
            bool discard = shader.fragment(bc, color);
            if (!discard) {
                zbuffer.set(x,y, TGAColor(frag_depth));
                image.set(x,y,color);
            }
        }
};
