/**
 * @file TypeDefinitions.h
 * @author F. Gratl
 * @date 06.12.21
 */

#pragma once

#include <autopas/AutoPasDecl.h>

#include "ladds/particle/Particle.h"

using AutoPas_t = autopas::AutoPas<LADDS::Particle>;

const std::string LADDS_SPD_LOGGER_NAME = "laddsLog";