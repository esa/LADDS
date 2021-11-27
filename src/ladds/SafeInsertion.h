//
// Created by albert on 27.11.21.
//
#pragma once

#include <autopas/AutoPasDecl.h>

#include "Particle.h"

#include <iostream>


namespace SafeInsertion {



std::vector<Particle> insert(autopas::AutoPas<Particle>& autopas,std::vector<Particle>& newSatellites,double cutoff){
    std::vector<Particle> delayedInsertion = {};

    for (auto &nSat : newSatellites) {
        //only insert satellites, if they have a reasonable distance to other satellites
        bool collisionFree = true;
        double collisionRadius = 2*cutoff;
        std::array<double,3> boxSpan = {collisionRadius,collisionRadius,collisionRadius};
        std::array<double,3> lowCorner = autopas::utils::ArrayMath::sub(nSat.getPosition(),boxSpan);
        std::array<double,3> highCorner = autopas::utils::ArrayMath::add(nSat.getPosition(),boxSpan);
        for(auto iter = autopas.getRegionIterator(lowCorner, highCorner);iter.isValid();++iter) {
            std::array<double,3> diff = autopas::utils::ArrayMath::sub(nSat.getPosition(),iter->getPosition());
            collisionFree = collisionFree && (autopas::utils::ArrayMath::dot(diff,diff) > collisionRadius*collisionRadius);
        }
        if(collisionFree) {
            autopas.addParticle(nSat);
        } else {
            std::cout << "delayed!" << std::endl;
            delayedInsertion.push_back(nSat);
        }
    }
    return delayedInsertion;
}


}