#ifndef COMMON_H
#define COMMON_H

#include <cmath>
#include <vector>
#include <string>

using std::vector;
using std::string;


class Vec
{
public:
    double x, y, z;

    // Init
    Vec(double x_=0, double y_=0, double z_=0);
    Vec(const Vec& v);

    // Vector Op
    Vec operator+(const Vec& v) const;
    Vec operator-(const Vec& v) const;
    Vec operator*(const Vec& v) const;

    // Scalar Op
    Vec operator+(double s) const;
    Vec operator-(double s) const;
    Vec operator*(double s) const;
    Vec operator/(double s) const;

    // Assign Op
    Vec& operator+=(const Vec& v);
    Vec& operator-=(const Vec& v);
    Vec& operator+=(double s);
    Vec& operator-=(double s);
    Vec& operator*=(double s);
    Vec& operator/=(double s);

    // Length Op
    double len() const;
    double len2() const;
    void norm();

    // Index Op
    double operator[](int i) const;

    // Product Op
    static double dot(const Vec& a, const Vec& b);
    static Vec cross(const Vec& a, const Vec& b);
};


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
    vector<Edge> edges;

    Face(int v1_, int v2_, int v3_, const Vec& x_, const Vec& y_, const Vec& z_);
};


#endif // COMMON_H
