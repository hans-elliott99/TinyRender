#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "model.h"

/* Wavefront OBJ Model
 * https://en.wikipedia.org/wiki/Wavefront_.obj_file
 */

// If the line's first 2 characters are "v" we have a list of geometric vertices.
/* We will assume that each vertex list contains just 3 elements - x, y, and z.
    * (They could also contain a 4th, "[w]", which is optional and defaults to 1.0)
*/
//Else if they're "f " we have a polygonal face element...
/* First number after each space - the number of the vertex in the array that we have read before.
   Second number - the texture coordinates of this vertex of this triangle.
   Third number - Normals...
    * For example, f 1193/1240/1193 1180/1227/1180 1179/1226/1179 says that 1193, 1180 and 1179 vertices form a triangle.
    * Note that in obj files indexes start from 1, meaning you should look for 1192, 1179 and 1178 vertices given C++ zero-indexing. 
    * For vector 1193, the texture coordinates are at index 1240 of the "vt" vectors.
*/


Model::Model(const char *filename)
{
    std::ifstream in;
    in.open( filename, std::ifstream::in );
    if (in.fail()) return;

    std::string line;
    while (!in.eof())
    {
        // Get a line of text from the file and write to the input-string-stream
        // using std::string::c_str(), which returns a pointer to an array that contains a null-terminated
        // seq of chars representing the current value of the std::string object.
        std::getline(in, line);
        std::istringstream iss( line.c_str() );

        char trash;                          //temporary var which we can write data to.
        if (!line.compare(0, 2, "v "))
        {
            iss >> trash;
            vec3 v;
            for (int i=0; i < 3; i++) iss >> v[i];
            verts.push_back(v);
        } else if (!line.compare(0, 3, "vn ")) 
        {
            iss >> trash >> trash;
            vec3 n;
            for (int i=0;i<3;i++) iss >> n[i];
            norms.push_back( n.normalize() );
        } else if (!line.compare(0, 3, "vt "))
        {
            iss >> trash >> trash;
            vec2 uv;
            for (int i=0;i<2;i++) iss >> uv[i];
            tex_coord.push_back( {uv.x, 1-uv.y} );
        } else if (!line.compare(0, 2, "f "))
        {
            int f, t, n;
            iss >> trash;
            int cnt = 0;
            while (iss >> f >> trash >> t >> trash >> n)
            {
                facet_v.push_back(--f); //decrement since obj files are 1 indexed
                facet_t.push_back(--t);
                facet_n.push_back(--n);
                cnt++;
            }
            if (cnt != 3)
            {
                std::cerr << "Error: the obj file must be triangulated." << '\n';
                in.close();
                return;
            }
        }
    }
    in.close();
    std::cerr << " #v: " << nverts() << " #f: " << nfaces() 
              << " #vt: " << tex_coord.size() << " #vn: " << norms.size()
             << '\n';
    
    // Load texture
    load_texture(filename, "_diffuse.tga",    diffusemap);
    load_texture(filename, "_nm_tangent.tga", normalmap  ); //_nm_tangent
    load_texture(filename, "_spec.tga",       specularmap);

}
// https://stackoverflow.com/questions/9158894/differences-between-c-string-and-compare


Model::~Model() {}

int Model::nverts() const
{
    return verts.size();
}

int Model::nfaces() const
{
    return facet_v.size()/3;
}

vec3 Model::vert(const int i) const
{
    return verts[i];
}

vec3 Model::vert(const int iface, const int nthvert) const
{
    return verts[ facet_v[iface*3 + nthvert] ];
}


vec2 Model::uv(const int iface, const int nthvert) const
{
    return tex_coord[ facet_t[iface*3 + nthvert] ];
}

vec3 Model::normal(const int iface, const int nthvert) const {
    return norms[facet_n[iface*3 + nthvert]];
}

vec3 Model::normal(const vec2 &uv) const {
    // Sample the color from the normalmap and convert to normal values
    TGAColor c = normalmap.get(uv[0]*normalmap.width(), uv[1]*normalmap.height());
    return vec3(c[2], c[1], c[0]) * 2.f/255. - vec3(1,1,1);
}



TGAColor Model::diffuse(vec2 uv) const {
    vec2 _uv(uv[0] * diffusemap.width(), uv[1] * diffusemap.height());
    return diffusemap.get(_uv[0], _uv[1]);
}

float Model::specular(vec2 uv) const {
    vec2 _uv(uv[0]*specularmap.width(), uv[1]*specularmap.height());
    return specularmap.get(uv[0], uv[1])[0] / 1.f;
}



void Model::load_texture(std::string filename, const std::string suffix, TGAImage &img)
{
    size_t dot = filename.find_last_of(".");
    if (dot == std::string::npos) return; //no file extension
    std::string texfile = filename.substr(0, dot) + suffix;
    diffuse_success = img.read_tga_file(texfile.c_str());

    std::cerr << "texture file " << texfile << " | loading " //<< texfile.c_str()
              << (diffuse_success ? " -success" : " -fail") << std::endl;
    
}



/*For backwards compatibility, get a basic triangle "face"*/
std::vector<int> Model::face(int idx)
{
    // the 3 face inds that index the vertices are stored sequentially
    std::vector<int> f = {
        facet_v[idx*3],
        facet_v[idx*3 + 1],
        facet_v[idx*3 + 2]
    };
    return f;
}