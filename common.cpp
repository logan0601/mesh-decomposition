#include "common.hpp"

int Utils::MAX_FACES = 5000;

double Utils::EPS = 1e-10;

double Utils::ETA = 0.2;
double Utils::DELTA = 0.8;

double Utils::UPPER = 1e300;

int Utils::NUM_SEED = 20;

double Utils::ANG_EPS = 0.3;
double Utils::DIS_EPS = 10;
double Utils::AVG_EPS = 0.3;

double Utils::PROB_EPS = 0.2;


double Vec::dot(const Vec& a, const Vec& b)
{
    return a.x*b.x + a.y*b.y + a.z*b.z;
}


Vec Vec::cross(const Vec& a, const Vec& b)
{
    return Vec(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
}


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
    if(normal.len() >= Utils::EPS) normal = normal.norm();
}
