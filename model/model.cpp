#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "model.h"

/* Wavefront OBJ Model
 * https://en.wikipedia.org/wiki/Wavefront_.obj_file
 */
Model::Model(const char *filename) : verts_(), faces_()
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
        // If the line's first 2 characters are "v " we have a list of geometric vertices.
        /* We will assume that each vertex list contains just 3 elements - x, y, and z.
         * (They could also contain a 4th, "[w]", which is optional and defaults to 1.0)
        */
        if (!line.compare(0, 2, "v "))
        {
            iss >> trash;
            Vec3f v;
            for (int i=0; i < 3; i++) iss >> v.raw[i];
            verts_.push_back(v);

        //Else if they're "f " we have a polygonal face element...
        /* We care about the first number after each space - the number of the vertex in the array that we have read before.
         * For example, f 1193/1240/1193 1180/1227/1180 1179/1226/1179 says that 1193, 1180 and 1179 vertices form a triangle.
         * Note that in obj files indexes start from 1, meaning you should look for 1192, 1179 and 1178 vertices respectively. 
        */
        } else if (!line.compare(0, 2, "f "))
        {
            std::vector<int> f;
            int itrash, idx;
            iss >> trash;
            while (iss >> idx >> trash >> itrash >> trash >> itrash)
            {
                idx--;            //wavefront obj is 1-indexed so decrement
                f.push_back(idx);
            }
            faces_.push_back(f);
        }
    }
    std::cerr << "# No. of Verts: " << verts_.size() << " No. of Faces: " << faces_.size() << '\n';
}
// https://stackoverflow.com/questions/9158894/differences-between-c-string-and-compare


Model::~Model() {}

int Model::nverts()
{
    return (int)verts_.size();
}

int Model::nfaces()
{
    return (int)faces_.size();
}

std::vector<int> Model::face(int idx)
// Extract faces at index
{
    return faces_[idx];
}

Vec3f Model::vert(int idx)
// Extract vertices at index
{
    return verts_[idx];
}


