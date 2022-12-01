#include <vector>
#include "geometry.h"
#include "../tga/tgaimage.h"

#pragma once

class Model
{
private:
    std::vector<Vec3f> verts {};
    std::vector<Vec2f> tex_coord {};
    std::vector<Vec3f> norms {};

    std::vector<int> facet_v {};
    std::vector<int> facet_t {};
    std::vector<int> facet_n {};
    TGAImage diffusemap{};         // diffuse color texture
    // TGAImage normalmap{};          // normal map texture
    // TGAImage specularmap{};        // specular map texture

    void load_texture(const std::string filename, const std::string suffix, TGAImage &img);

public:
    Model(const char *filename);
    ~Model();
    int nverts() const;
    int nfaces() const;

    Vec3f normal(const int iface, const int nthvert) const; // per triangle corner normal vertex
    Vec3f normal(const Vec2f &uv) const;                     // fetch the normal vector from the normal map texture
    Vec3f vert(const int i) const;
    Vec3f vert(const int iface, const int nthvert) const;
    Vec2f uv(const int iface, const int nthvert) const;
    
    const TGAImage& diffuse()  const { return diffusemap;  }
    TGAColor diffuse_get(Vec2i uv) {return diffusemap.get(uv.x, uv.y); }
    // const TGAImage& specular() const { return specularmap; }

    std::vector<int> face(int idx);
};