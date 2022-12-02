/* implement some basic vector geometry */
#include <iostream>
#include <cmath>
#include <assert.h>
#include <vector>

#pragma once

template <class T> struct Vec2
{
    union
    {
        struct {T u, v;};
        struct {T x, y;};
        T raw[2];
    };
    Vec2() : u(0), v(0) {}
    Vec2(T _u, T _v) : u(_u), v(_v) {}
	inline T&      operator[](const int i)       {assert(i>=0 && i < 2); return raw[i]; } 		
	inline T& 	   operator[](const int i) const {assert(i>=0 && i < 2); return raw[i]; } 				
    inline Vec2<T> operator+(const Vec2<T> &V) const { return Vec2<T>(u + V.u, v + V.v); }
    inline Vec2<T> operator-(const Vec2<T> &V) const { return Vec2<T>(u - V.u, v - V.v); }
    inline Vec2<T> operator *(float f)          const { return Vec2<T>(u*f, v*f); }
	template <class > friend std::ostream& operator<<(std::ostream& s, Vec2<T>& v);
};


template <class T> struct Vec3 
{
	union {
		struct {T x, y, z;};
		struct {T ivert, iuv, inorm; };
		T raw[3];
	};
	Vec3<T>() : x(0), y(0), z(0) {}
	Vec3<T>(T _x, T _y, T _z) : x(_x),y(_y),z(_z) {}
	template <class Y> Vec3<T>(const Vec3<Y> &v);

	Vec3<T> & operator=(const Vec3<T> &v) {
        if (this != &v) {
            x = v.x;
            y = v.y;
            z = v.z;
        }
        return *this;
    }
	inline T&      operator[](const int i)       {assert(i>=0 && i < 3); return raw[i]; } 		
	inline T& 	   operator[](const int i) const {assert(i>=0 && i < 3); return raw[i]; } 				
	inline Vec3<T> operator ^(const Vec3<T> &v) const { return Vec3<T>(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x); } //cross-product with external vec
	inline Vec3<T> operator +(const Vec3<T> &v) const { return Vec3<T>(x+v.x, y+v.y, z+v.z); }
	inline Vec3<T> operator -(const Vec3<T> &v) const { return Vec3<T>(x-v.x, y-v.y, z-v.z); }
	inline Vec3<T> operator *(float f)          const { return Vec3<T>(x*f, y*f, z*f); }
	inline T       operator *(const Vec3<T> &v) const { return x*v.x + y*v.y + z*v.z; }
	float norm() const { return std::sqrt(x*x+y*y+z*z); }
	Vec3<T> &normalize(T l=1) { *this = (*this)*(l/norm()); return *this; }
	template <class > friend std::ostream& operator<<(std::ostream& s, Vec3<T>& v);
};


typedef Vec2<float> Vec2f;
typedef Vec2<int>   Vec2i;
typedef Vec3<float> Vec3f;
typedef Vec3<int>   Vec3i;

template <class t> std::ostream& operator<<(std::ostream& s, Vec2<t>& v) {
	s << "(" << v.x << ", " << v.y << ")\n";
	return s;
}

template <class t> std::ostream& operator<<(std::ostream& s, Vec3<t>& v) {
	s << "(" << v.x << ", " << v.y << ", " << v.z << ")\n";
	return s;
}



const int DEFAULT_ALLOC=4;

class Matrix {
    std::vector<std::vector<float> > m;
    int rows, cols;
public:
    Matrix(int r=DEFAULT_ALLOC, int c=DEFAULT_ALLOC);
    inline int nrows();
    inline int ncols();

    static Matrix identity(int dimensions);
    std::vector<float>& operator[](const int i);
    Matrix operator*(const Matrix& a);
    Matrix transpose();
    Matrix inverse();

    friend std::ostream& operator<<(std::ostream& s, Matrix& m);
};
