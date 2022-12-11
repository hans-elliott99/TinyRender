#include <vector>
#include "geometry.h"
#include "tgaimage.h"

#pragma once

class Model
{
private:
    std::vector<vec3> verts {};
    std::vector<vec2> tex_coord {};
    std::vector<vec3> norms {};

    std::vector<int> facet_v {};
    std::vector<int> facet_t {};
    std::vector<int> facet_n {};
    TGAImage diffusemap{};             // diffuse color texture
    TGAImage normalmap{};              // normal map texture
    TGAImage specularmap{};            // specular map texture

    bool diffuse_success{false};

    void load_texture(const std::string filename, const std::string suffix, TGAImage &img);

public:
    Model(const char *filename);
    ~Model();
    int nverts() const;
    int nfaces() const;

    vec3 normal(const int iface, const int nthvert) const; // per triangle corner normal vertex
    vec3 normal(const vec2 &uv) const;                     // fetch the normal vector from the normal map texture
    vec3 vert(const int i) const;
    vec3 vert(const int iface, const int nthvert) const;
    vec2 uv(const int iface, const int nthvert) const;
    
    TGAColor diffuse(vec2 uv)  const;
    TGAImage get_diffuse() {return diffusemap; };
    float specular(vec2 uv) const;

    std::vector<int> face(int idx);
};