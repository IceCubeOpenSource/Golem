#ifndef __GOLEM_H__
#define __GOLEM_H__

#include <vector>
#include <deque>

#include <PhysTools/histogram.h>
#include <PhysTools/bin_type.h>
#include <PhysTools/likelihood/likelihood.h>

namespace golem {

template<typename PhysicsModel, typedef Likelihood>
class Golem {
    typedef PhysicsModel::Event Event;
    typedef PhysicsModel::WeighterMaker WeighterMaker;
    typedef PhysicsModel::UncertaintyWeighter UncertaintyWeighter;
    typedef PhysicsModel::HistogramSet HistogramSet;
    typedef PhysicsModel::Prior Prior;
    static constexpr unsigned int NParameters = PhysicsModel::NParameters;

    typedef phys_tools::likelihood::detail::SwitchableWeighter<phys_tools::likelihood::SimpleDataWeighter, decltype(WeighterMaker(std::vector<double>))> DataWeighter;

    typedef phys_tools::likelihood::LikelihoodProblem<std::reference_wrapper<const Event>, HistogramSet, DataWeighter, phys_tools::likelihood::detail::WeighterCollection<WeighterMaker, UncertaintyWeighter>, Prior, Likelihood, NParameters> LType;

    static constexpr dataWeighterIndex = 0;
    static constexpr asimovWeighterIndex = 1;

    PhysicsModel model;
    HistogramSet observationHistogram;
    HistogramSet simulationHistogram;
    phys_tools::ParameterSet params;
    WeighterMaker WM;
    UncertaintyWeighter UW;
    std::unique_ptr<LType> likelihoodProblem;

    std::deque<Event> observation;
    std::deque<Event> simulation;

    std::mt19937 rng;

public:
    Golem(PhysicsModel model):
            model(model),
            observationHistogram(model.MakeHistogramSet()),
            simulationHistogram(model.MakeHistogramSet()),
            params(model.MakeParameterSet()),
            WM(model.MakeWeighterMaker()),
            UW(model.MakeUncertaintyWeighter()) {
    }

    // Parameter functions
    phys_tools::ParameterSet & GetParameters();

    // LLH functions
    phys_tools::likelihoodPoint MinLLH() const;
    double EvalLLH(std::vector<double> params) const;

    // Observation / Expectation functions
    void SetObservationEvents(std::deque<Event> events) {
        observation = std::move(events);
        for(const Event& e : observation)
            model.AddEventToHistogram(observationHistogram, e);
        if(!observation.empty() && !simulation.empty())
            MakeLikelihoodProblem();
    }

    void SetSimulationEvents(std::deque<Event> events) {
        simulation = std::move(events);
        for(const Event& e : simulation)
            model.AddEventToHistogram(simulationHistogram, e);
        if(!observation.empty() && !simulation.empty())
            MakeLikelihoodProblem();
    }

    const std::deque<Event> & GetObservationEvents() const {
        return observation;
    }

    const std::deque<Event> & GetSimulationEvents() const {
        return simulation;
    }

    // Distribution functions
    HistogramSet GetExpectationDistribution(std::vector<double>) const;
    HistogramSet GetSimulationDistribution(std::vector<double>) const;

    // Fun stuff
    //MakeHistogramOfArbitraryProperty(SomeHistogramType histogram, SomeParameterPointerType parameter_pointer)

    // Test functions
    void SetUpAsimov(std::vector<double> parameters) {
        if(simulation.empty())
            throw std::runtime_error("Simulation cannot be empty for Asimov test.");
        // make the weighter for parameters
        auto weighter = WM(parameters);

        // Set the asimov weighter
        likelihoodProblem.dataWeighter.setWeighter(asimovWeigherIndex);
        std::get<asimovWeighterIndex>(likelihoodProblem.dataWeighter.implementations) = weighter;

        // overwrite the observation in the Golem
        observation = simulation;
        observationHistogram = simulationHistogram;

        // overwrite the observation in the LP
        likelihoodProblem.observation = simulationHistogram;
    }

    std::deque<Event> GetRealization(std::vector<double> parameters) const {
        std::vector<double> weights;
        auto weighter = WM(parameters);
        double expectation=0;
        double weight;
        for(const Event& e : simulation) {
            weight = weighter(e);
            weights.push_back(weight);
            expectation += weight;
        }
        return phys_tools::likelihood::generateSample(weights, simulation, expectation, rng);
    }

    // RNG functions
    void SetSeed(unsigned int seed) {
        rng.seed(seed);
    }

private:
    void MakeLikelihoodProblem() {
        likelihoodProblem.reset(); // delete the LP
        //TODO make the likelihood problem again...

    }
}

} // namespace golem

#endif
