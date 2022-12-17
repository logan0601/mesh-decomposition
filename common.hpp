#ifndef VEC_H
#define VEC_H

#include <iostream>
#include <cmath>
#include <vector>

using std::vector;


class Utils
{
public:
    static int MAX_FACES;   // Max Number Faces

    static double EPS;      // Math Calculation eps(Global EPS)

    static double ETA;      // For Edge Weight Computing. Default: 0.2
    static double DELTA;    // Balance Between Angle & Geodesic Distance

    static double UPPER;    // Max Double

    static int NUM_SEED;    // Try Initial Seed Number

    static double ANG_EPS;  // Stopping Angle eps
    static double DIS_EPS;  // Stopping Distance eps
    static double AVG_EPS;  // Stopping Average Distance Ratio eps;

    static double PROB_EPS; // Probability eps

};


class Vec
{
public:
    double x, y, z;

    // Init
    Vec(double x_=0.0, double y_=0.0, double z_=0.0): x(x_), y(y_), z(z_) { }
    Vec(const Vec& v) { x=v.x; y=v.y; z=v.z; }

    // Vector Operator
    Vec operator+(const Vec& v) const { return Vec(x+v.x, y+v.y, z+v.z); }
    Vec operator-(const Vec& v) const { return Vec(x-v.x, y-v.y, z-v.z); }
    Vec operator*(const Vec& v) const { return Vec(x*v.x, y*v.y, z*v.z); }

    // Scalar Operator
    Vec operator+(double s) const { return Vec(x+s, y+s, z+s); }
    Vec operator-(double s) const { return Vec(x-s, y-s, z-s); }
    Vec operator*(double s) const { return Vec(x*s, y*s, z*s); }
    Vec operator/(double s) const { return Vec(x/s, y/s, z/s); }

    // Assign Operator
    Vec& operator+=(const Vec& v) { return *this = *this + v; }
    Vec& operator-=(const Vec& v) { return *this = *this - v; }
    Vec& operator+=(double s) { return *this = *this + s; }
    Vec& operator-=(double s) { return *this = *this - s; }
    Vec& operator*=(double s) { return *this = *this * s; }
    Vec& operator/=(double s) { return *this = *this / s; }

    // Vector Operator Product
    double dot(const Vec& b) const { return x*b.x + y*b.y + z*b.z; }
    Vec cross(const Vec& b) const { return Vec(y*b.z-z*b.y, z*b.x-x*b.z, x*b.y-y*b.x); }

    // Compare Operator
    double max() { return x>y&&x>z ? x : (y>z ? y : z); }
    double min() { return x<y&&x<z ? x : (y<z ? y : z); }

    // Length Operator
    double len() { return sqrt(x*x + y*y + z*z); }
    double len2() { return x*x + y*y + z*z; }
    Vec& norm() { return *this = *this / this->len(); }

    // Index Operator
    double operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }

    static double dot(const Vec& a, const Vec& b);
    static Vec cross(const Vec& a, const Vec& b);
};

inline std::ostream& operator<<(std::ostream& os, const Vec& v)
{
    os << "<" << v.x << ", " << v.y << ", " << v.z << ">";
    return os;
}


class Edge
{
public:
    int vid[2];
    int fid;
    double angle;
    double ang_dis;
    double geo_dis;
    double tot_dis;

    Edge(int v1_, int v2_, int f_, double ang_, double ang_dis_, double geo_dis_);
};


class Face
{
public:
    int vid[3];
    int label;
    Vec center;
    Vec normal;
    std::vector<Edge> nbrs;

    Face(int v1_, int v2_, int v3_, const Vec& x_, const Vec& y_, const Vec& z_);
};

#endif // VEC_H
