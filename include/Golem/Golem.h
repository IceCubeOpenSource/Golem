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

    typedef phys_tools::likelihood::SwitchableWeighter<phys_tools::likelihood::SimpleDataWeighter, decltype(WeighterMaker(std::vector<double>))> DataWeighter;

    typedef phys_tools::likelihood::LikelihoodProblem<std::reference_wrapper<const Event>, HistogramSet, DataWeighter, phys_tools::likelihood::detail::WeighterCollection<WeighterMaker, UncertaintyWeighter>, Prior, Likelihood, NParameters> LType;

    PhysicsModel model;
    phys_tools::ParameterSet params;
    WeighterMaker WM;
    std::unique_ptr<LType> likelihoodProblem;

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
