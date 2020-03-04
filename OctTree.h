/*
 * Header File: OctTree.h
 * Author: Darren Trieu Nguyen
 * Last Modified: 3-3-20
 */

/* Particle Declaration */
class Particle {
    public:
    double x;
    double y;
    double z;

    Particle(double xi, double yi, double zi) : x(xi), y(yi), z(zi) { }
};

/* OctTreeNode Declaration */          
// Based off of http://mathandcode.com/2016/02/16/quadtree.html
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

    /* Pointers to 8 subnodes */
    OctTreeNode* nsss = nullptr;
    OctTreeNode* nssl = nullptr;
    OctTreeNode* nsls = nullptr;
    OctTreeNode* nsll = nullptr;
    OctTreeNode* nlss = nullptr;
    OctTreeNode* nlsl = nullptr;
    OctTreeNode* nlls = nullptr;
    OctTreeNode* nlll = nullptr;

    /* Pointer to a particle */
    Particle* particle = nullptr;

    /* OctTreeNode Constructor Prototype */
    OctTreeNode(double smallx, double bigx, double smally, double bigy, 
                double smallz, double bigz);

    /* Add particle function */
    void addParticle(Particle* p);
};
