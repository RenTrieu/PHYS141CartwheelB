/*
 * Header File: OctTree.h
 * Author: Darren Trieu Nguyen
 * Last Modified: 3-4-20
 */

#include "nBodyLeapint.h"

/* Constants */
#define heurThreshold 1.0

/* Particle Declaration */
class Particle {
    public:
    double x;
    double y;
    double z;
    double mass;
    double ax = 0.0;
    double ay = 0.0;
    double az = 0.0;

    Particle(double xi, double yi, double zi, double m)
    : x(xi), y(yi), z(zi), mass(m) { }
};

/* OctTreeNode Declaration */          
// Based off of http://mathandcode.com/2016/02/16/quadtree.html
// Just ported to C++ and extended to 3 dimensions 
class OctTreeNode {
    public:
    /* Six bounds on the octtree node rectangular prism */
    double xs; double xl;
    double ys; double yl;
    double zs; double zl;

    /* The center */
    double cx;
    double cy;
    double cz;

    /* Center of Masses */
    double cmx = 0;
    double cmy = 0;
    double cmz = 0;
    double totalMass = 0;

    /* Pointers to 8 subnodes */
    OctTreeNode* nsss = nullptr;
    OctTreeNode* nssl = nullptr;
    OctTreeNode* nsls = nullptr;
    OctTreeNode* nsll = nullptr;
    OctTreeNode* nlss = nullptr;
    OctTreeNode* nlsl = nullptr;
    OctTreeNode* nlls = nullptr;
    OctTreeNode* nlll = nullptr;

    /* Array of pointers of the 8 subnodes */
    OctTreeNode** subNodeList;
    int subListSize;

    /* Pointer to a particle */
    Particle* particle = nullptr;

    /* OctTreeNode Constructor Prototype */
    OctTreeNode(double smallx, double bigx, double smally, double bigy, 
                double smallz, double bigz);

    /* Add particle function */
    void addParticle(Particle* p);

    /* Calculate acceleration function */
    int calculateAccel(Particle* p);
};
