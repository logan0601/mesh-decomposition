#include "common.h"


Vec::Vec(double x_, double y_, double z_) { x=x_; y=y_; z=z_; }
Vec::Vec(const Vec& v) { x=v.x; y=v.y; z=v.z; }

Vec Vec::operator+(const Vec& v) const { return Vec(x+v.x, y+v.y, z+v.z); }
Vec Vec::operator-(const Vec& v) const { return Vec(x-v.x, y-v.y, z-v.z); }
Vec Vec::operator*(const Vec& v) const { return Vec(x*v.x, y*v.y, z*v.z); }

Vec Vec::operator+(double s) const { return Vec(x+s, y+s, z+s); }
Vec Vec::operator-(double s) const { return Vec(x-s, y-s, z-s); }
Vec Vec::operator*(double s) const { return Vec(x*s, y*s, z*s); }
Vec Vec::operator/(double s) const { return Vec(x/s, y/s, z/s); }

Vec& Vec::operator+=(const Vec& v) { return *this = *this + v; }
Vec& Vec::operator-=(const Vec& v) { return *this = *this - v; }
Vec& Vec::operator+=(double s) { return *this = *this + s; }
Vec& Vec::operator-=(double s) { return *this = *this - s; }
Vec& Vec::operator*=(double s) { return *this = *this * s; }
Vec& Vec::operator/=(double s) { return *this = *this / s; }

double Vec::len() const { return sqrt(x*x + y*y + z*z); }
double Vec::len2() const { return x*x + y*y + z*z; }
void Vec::norm() { *this = *this / this->len(); }

double Vec::dot(const Vec& a, const Vec& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
Vec Vec::cross(const Vec& a, const Vec& b) { return Vec(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x); }

Edge::Edge(int v1_, int v2_, int f_, double ang_, double ang_dis_, double geo_dis_)
{
    vid[0] = v1_;
    vid[1] = v2_;
    fid = f_;
    angle = ang_;
    ang_dis = ang_dis_;
    geo_dis = geo_dis_;
}

Face::Face(int v1_, int v2_, int v3_, const Vec& x_, const Vec& y_, const Vec& z_)
{
    vid[0] = v1_;
    vid[1] = v2_;
    vid[2] = v3_;
    label = 0;
    center = (x_ + y_ + z_) / 3.0;
    normal = Vec::cross(z_ - x_, y_ - x_);
    if(normal.len() >= 1e-10) normal.norm();
}
