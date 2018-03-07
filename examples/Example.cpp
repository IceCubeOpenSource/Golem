#include "DataLoader.h"
#include "EventConverter.h"
#include "PhysicsModel.h"

#include <Golem/Golem.h>

int main() {
  analysis::DataLoader data_loader("resources/monte_carlo/simulation.hdf5", "resources/data/observation.hdf5");
  analysis::EventConverter event_converter("resources/monte_carlo/generation.lic","resources/cross_sections/");
  analysis::PhysicsModel physics_model;

  golem::Golem<analysis::PhysicsModel,phys_tools::likelihood::SAYLikelihood> golem(physics_model);

  phys_tools::likelihood::likelihoodPoint p = golem.MinLLH();

  return 0;
}
