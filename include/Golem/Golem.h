#ifndef __GOLEM_H__
#define __GOLEM_H__

#include <vector>
#include <deque>

#include <PhysTools/histogram.h>
#include <PhysTools/bin_type.h>
#include <PhysTools/likelihood/likelihood.h>

namespace golem {

template<typename PhysicsModel>
class Golem {
    PhysicsModel model;
    phys_tools::ParameterSet params;

    typedef PhysicsModel::Event Event;
    typedef PhysicsModel::WeighterMaker WeighterMaker;
    typedef PhysicsModel::HistogramSet HistogramSet;

    Golem(PhysicsModel model);

    // Parameter functions
    phys_tools::ParameterSet & GetParameters();
    
    //LLH functions
    phys_tools::likelihoodPoint MinLLH() const;
    double EvalLLH(std::vector<double> params) const;

    //Observation / Expectation functions
    void SetObservationEvents(std::deque<Event> events);
    void SetSimulationEvents(std::deque<Event> events);
    const std::deque<Event> & GetObservationEvents() const;
    const std::deque<Event> & GetSimulationEvents() const;

    //Distribution functions
    HistogramSet GetExpectationDistribution(std::vector<double>) const;
    HistogramSet GetSimulationDistribution(std::vector<double>) const;

    //Fun stuff
    //MakeHistogramOfArbitraryProperty(SomeHistogramType histogram, SomeParameterPointerType parameter_pointer)

    //Test functions
    void SetUpAsimov();
    std::deque<Event> GetRealization(std::vector<double>) const;
}

} // namespace golem

#endif
