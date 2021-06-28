#include <iostream>
#include "Debris.h"
#include "autopas/AutoPasDecl.h"
#include "Logger.h"

// Declare the main AutoPas class as extern template instantiation. It is instantiated in AutoPasClass.cpp.
extern template class autopas::AutoPas<Debris>;

int main() {

  Logger logger;

  autopas::AutoPas<Debris> autopas;

  constexpr size_t numDebris = 2;

  autopas.setBoxMin({0.,0.,0.});
  autopas.setBoxMax({10.,10.,10.});
  autopas.setCutoff(1);
  autopas.init();

  for (size_t i = 0; i < numDebris; ++i) {
    autopas.addParticle(Debris{{static_cast<double>(i),0,0}, {0.,0.,0.}, i});
  }

  for (const auto &d : autopas) {
    logger.log(Logger::Level::info, d.toString());
  }
  return 0;
}
