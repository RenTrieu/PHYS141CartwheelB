/*
 * Program: OctTree.cpp
 * Author: Darren Trieu Nguyen
 * Last Modified: 3-4-20
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
            cmx = p->x * p->mass;
            cmy = p->y * p->mass;
            cmz = p->z * p->mass;
            totalMass = p->mass;

            this->particle = p;
        }
        else {
            /* No subnodes detected, but particle detected */
            Particle* tmp = this->particle;

            /* If the particles are in the same space, exclude new particle */
            if ((abs(tmp->x - p->x) < 0.001) && (abs(tmp->y - p->y) < 0.001)) {
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

            cmx = 0;
            cmy = 0;
            cmz = 0;
            totalMass = 0;

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

        /* Center of mass calculation */
        cmx += p->x * p->mass;
        cmy += p->y * p->mass;
        cmz += p->z * p->mass;
        totalMass = p->mass;

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

/* 
 * OctTreeNode calculateAccel Function 
 * Parameters: Particle * p - pointer to particle for which to calculate 
 *                            acceleration
 * Returns: 0 if acceleration is calculated, 
 *          1 if acceleration is not calculated
 */
int OctTreeNode::calculateAccel(Particle* p) {
    if (this->nsss == nullptr) {
        /* Empty subnode and empty particle case */
        if (this->particle == nullptr) {
            return 1;
        }
        /* Empty subnode, but existing particle case */
        else {
            double distVal = (this->particle->x - p->x)
                             *(this->particle->x - p->x)
                             +(this->particle->y - p->y)
                             *(this->particle->y - p->y)
                             +(this->particle->z - p->z)
                             *(this->particle->z - p->z) + 0.00001;
            distVal = distVal * distVal * distVal;
            distVal = 1.0f / sqrtf(distVal);

            /* Summing up acceleration */
            p->ax += (this->particle->x - p->x)*GCONST
                     *this->particle->mass*distVal;
            p->ay += (this->particle->y - p->y)*GCONST
                     *this->particle->mass*distVal;
            p->az += (this->particle->z - p->z)*GCONST
                     *this->particle->mass*distVal;
            return 0;
        }
    }
    /* Subnodes exist case */
    else {
        /* Heuristic Analysis */
        double avgSize = (this->xl-this->xs 
                          + this->yl-this->ys
                          + this->zl-this->zs) / 3.0;
        double distVal = sqrtf((this->particle->x - p->x)
                               *(this->particle->x - p->x)
                               +(this->particle->y - p->y)
                               *(this->particle->y - p->y)
                               +(this->particle->z - p->z)
                               *(this->particle->z - p->z));
        /* Recursive Case */
        if ((avgSize / distVal) >= heurThreshold) {
            this->nsss->calculateAccel(p);
            this->nssl->calculateAccel(p);
            this->nsls->calculateAccel(p);
            this->nsll->calculateAccel(p);
            this->nlss->calculateAccel(p);
            this->nlsl->calculateAccel(p);
            this->nlls->calculateAccel(p);
            this->nlll->calculateAccel(p);
            return 0;
        }
        /* Non-Recursive Center of Mass Case */
        else {
            this->nsss->totalMass*(this->nss->cmx - this->p)
            return 0;
        }
    }
    return 1;
}

