/*
 * Program: OctTree.cpp
 * Author: Darren Trieu Nguyen
 * Last Modified: 3-3-20
 */
#include <algorithm>
#include <cstring>
#include <cmath>
#include "OctTree.h"

using namespace std;

/* OctTreeNode Constructor */
OctTreeNode::OctTreeNode(double smallx, double bigx,
            double smally, double bigy,
            double smallz, double bigz) : xs(smallx), xl(bigx),
                                          ys(smally), yl(bigy),
                                          zs(smallz), zl(bigz) {
    this->cx = (xs+xl) / 2.0;
    this->cy = (ys+yl) / 2.0;
    this->cz = (zs+zl) / 2.0;
}

/* OctTreeNode addParticle Function */
void OctTreeNode::addParticle(Particle* p) {
    if (this->nsss == nullptr) {
        if (this->particle == nullptr) {
            /* No particle or subnodes detected, adding particle to node */
            this->particle = p;
        }
        else {
            /* No subnodes detected, but particle detected */
            Particle* tmp = this->particle;

            /* If the particles are in the same space, exclude new particle */
            if (abs(tmp->x - p->x) < 0.001 && abs(tmp->y - p->y) < 0.001) {
                return;
            }

            /* Creating subnodes */
            this->nsss = new OctTreeNode(xs, ys, zs, cx, cy, cz);
            this->nssl = new OctTreeNode(xs, ys, cz, cx, cy, zl);
            this->nsls = new OctTreeNode(xs, cy, zs, cx, yl, cz);
            this->nsll = new OctTreeNode(xs, cy, cz, cx, yl, zl);
            this->nlss = new OctTreeNode(cx, ys, zs, xl, cy, cz);
            this->nlsl = new OctTreeNode(cx, ys, cz, xl, cy, zl);
            this->nlls = new OctTreeNode(cx, cy, zs, xl, yl, cz);
            this->nlll = new OctTreeNode(cx, cy, cz, xl, yl, zl);
            this->particle = nullptr;

            /* Recursively call addParticle to sort the particles 
               into the corresponding subnodes */
            this->addParticle(tmp);
            this->addParticle(p);
        }
    }
    else {
        /* Using Math And Code's binary trick */
        // TODO: Check binary arithmetic
        int a = (p->x >= cx) ? 1:0;
        int b = (p->y >= cy) ? 1:0;
        int c = (p->z >= cz) ? 1:0;
        int d=a|(b<<1)|(c<<2);
        if (c == 0) {
            this->nsss->addParticle(p);
        }
        else if (c == 1) {
            this->nlss->addParticle(p);
        }
        else if (c == 2) {
            this->nsls->addParticle(p);
        }
        else if (c == 3) {
            this->nlls->addParticle(p);
        }
        else if (c == 4) {
            this->nssl->addParticle(p);
        }
        else if (c == 5) {
            this->nlsl->addParticle(p);
        }
        else if (c == 6) {
            this->nsll->addParticle(p);
        }
        else if (c == 7) {
            this->nlll->addParticle(p);
        }
        return;
    }
}

