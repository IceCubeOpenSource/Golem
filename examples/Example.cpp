#include <deque>
#include "DataLoader.h"
#include "EventConverter.h"
#include "PhysicsModel.h"

#include <Golem/Golem.h>

int main() {
    analysis::DataLoader data_loader("resources/monte_carlo/simulation.hdf5", "resources/data/observation.hdf5");
    analysis::EventConverter event_converter("resources/monte_carlo/generation.lic","resources/cross_sections/");
    analysis::PhysicsModel physics_model;

    std::deque<analysis::DataLoader::Event> raw_simulation_events = data_loader.GetSimulationEvents();
    std::deque<analysis::DataLoader::Event> raw_observation_events = data_loader.GetDataEvents();

    golem::Golem<analysis::PhysicsModel,phys_tools::likelihood::SAYLikelihood> golem(physics_model);

    std::deque<analysis::PhysicsModel::Event> simulation_events;
    for(const analysis::DataLoader::Event e : raw_simulation_events) {
        simulation_events.push_back(event_converter.ConvertEvent(e));
    }
    std::deque<analysis::PhysicsModel::Event> observation_events;
    for(const analysis::DataLoader::Event e : raw_observation_events) {
        observation_events.push_back(event_converter.ConvertEvent(e));
    }

    golem.SetObservationEvents(observation_events);
    golem.SetSimulationEvents(simulation_events);

    phys_tools::likelihood::likelihoodPoint p = golem.MinLLH();

    const std::vector<double>& parameter_values = p.params;
    for(double val : parameter_values) {
        std::cout << val << std::endl;
    }

    return 0;
}
