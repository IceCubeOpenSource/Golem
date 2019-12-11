#ifndef GOLEM_H
#define GOLEM_H

#include <vector>
#include <deque>
#include <memory>

#include <PhysTools/histogram.h>
#include <PhysTools/bin_types.h>
#include <PhysTools/likelihood/likelihood.h>

namespace golem {

template<typename PhysicsModel_,typename Likelihood_>
class Golem {
    typedef PhysicsModel_ PhysicsModel;
    typedef Likelihood_ Likelihood;
    typedef typename PhysicsModel::Event Event;
    typedef typename PhysicsModel::WeighterMaker WeighterMaker;
    typedef typename PhysicsModel::UncertaintyWeighter UncertaintyWeighter;
    typedef typename PhysicsModel::HistogramSet HistogramSet;
    typedef typename PhysicsModel::Prior Prior;

    static constexpr unsigned int NParameters = PhysicsModel::NParameters;

    typedef phys_tools::likelihood::detail::SwitchableWeighter<phys_tools::likelihood::simpleDataWeighter, decltype(std::declval<WeighterMaker>()(std::declval<std::vector<double>>()))> DataWeighter;

    typedef phys_tools::likelihood::LikelihoodProblem<std::reference_wrapper<const Event>, HistogramSet, DataWeighter, phys_tools::likelihood::detail::WeighterCollection<WeighterMaker, UncertaintyWeighter>, Prior, Likelihood, NParameters> LType;

    static constexpr unsigned int dataWeighterIndex = 0;
    static constexpr unsigned int asimovWeighterIndex = 1;

    PhysicsModel model;
    HistogramSet observationHistogram;
    HistogramSet simulationHistogram;
    phys_tools::ParameterSet params;
    WeighterMaker WM;
    UncertaintyWeighter UW;
    phys_tools::likelihood::detail::WeighterCollection<WeighterMaker, UncertaintyWeighter> WC;
    Prior prior;
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
            WM(model.MakeWeighterMaker(params)),
            UW(model.MakeUncertaintyWeighter()),
            WC(WM,UW),
            prior(model.MakePrior(params)) {
        MakeLikelihoodProblem();
    }

    // Parameter functions
    phys_tools::ParameterSet & GetParameters() {
        return params;
    }

    // LLH functions
    phys_tools::likelihood::likelihoodPoint MinLLH(double changeTolerance=1e-6, unsigned int historySize=5) const {
        if(!likelihoodProblem)
            throw std::runtime_error("LikelihoodProblem needs to exist before the minimizer can be called.");
        phys_tools::lbfgsb::LBFGSB_Driver minimizer(params);
        minimizer.setChangeTolerance(changeTolerance);
        minimizer.setHistorySize(historySize);
        minimizer.minimize(phys_tools::likelihood::BFGS_Function<LType>(*likelihoodProblem));
        phys_tools::likelihood::likelihoodPoint result;
        result.likelihood = minimizer.minimumValue();
        result.params = minimizer.minimumPosition();
        return result;
    }

    double EvalLLH(std::vector<double> params) const {
        if(!likelihoodProblem)
            throw std::runtime_error("LikelihoodProblem needs to exist before thread count can be set.");
        return -likelihoodProblem->evaluateLikelihood(params);
    }

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

    //TODO Distribution functions
    //HistogramSet GetObservationDistribution() const;
    //
    template<typename T>
    static std::tuple<T> groupHistograms(T t) {
        return std::tuple<T>(t);
    }

    template<typename T, typename... Args>
    static decltype(std::tuple_cat(std::declval<T>(), groupHistograms(std::declval<Args>()...))) groupHistograms(T t, Args... args) {
        return std::tuple_cat(std::tuple<T>(t), groupHistograms(args...));
    }

    template<typename DataType, typename T>
    static std::tuple<T> getNumericHistogram(T t) {
        return std::tuple<phys_tools::histograms::histogram<T::dimensions>, DataType>(phys_tools::histograms::histogram<T::dimensions, DataType>());
    }

    template<typename DataType, typename T, typename... Args>
    static decltype(std::tuple_cat(std::declval<T>(), getNumericHistogram(std::declval<Args>()...))) getNumericHistogram(T t, Args... args) {
        return std::tuple_cat(std::tuple<phys_tools::histograms::histogram<T::dimensions>, DataType>(phys_tools::histograms::histogram<T::dimensions, DataType>()),
                getNumericHistogram<DataType>(args...));
    }

    //decltype(otherhorriblething<HistogramSet>) GetSimulationDistribution(std::vector<double> parameters) const {
    //phys_tools::histograms::histogram<std::tuple_element<0,HistogramSet>::type::dimensions,double> GetSimulationDistribution(std::vector<double> parameters) const {

    //}

    /*
    decltype(std::declval<HistogramSet>()) GetSimulationDistribution(std::vector<double> parameters) const {
        auto weighter = WM(parameters);
        //decltype(simulationHistogram.getAxis(0));
    }
    */

    //TODO Fun stuff
    //MakeHistogramOfArbitraryProperty(SomeHistogramType histogram, SomeParameterPointerType parameter_pointer)

    // Test functions
    void SetUpAsimov(std::vector<double> parameters) {
        if(simulation.empty())
            throw std::runtime_error("Simulation cannot be empty for Asimov test.");
        // make the weighter for parameters
        auto weighter = WM(parameters);

        // Set the asimov weighter
        likelihoodProblem.dataWeighter.setWeighter(asimovWeighterIndex);
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

    void SetThreadCount(unsigned int nThreads) {
        if(likelihoodProblem)
            likelihoodProblem->setEvaluationThreadCount(nThreads);
        else
            throw std::runtime_error("LikelihoodProblem needs to exist before thread count can be set.");
    }

private:
    void MakeLikelihoodProblem() {
        likelihoodProblem.reset( new LType(
                phys_tools::likelihood::makeLikelihoodProblem<std::reference_wrapper<const Event>,NParameters>(
                    observationHistogram,
                    {simulationHistogram},
                    prior,
                    {0.0},
                    DataWeighter(phys_tools::likelihood::detail::SwitchableWeighter<phys_tools::likelihood::simpleDataWeighter, decltype(WM(std::declval<std::vector<double>>()))>(phys_tools::likelihood::simpleDataWeighter(), WM(params.getParameterValues()))),
                    WC,
                    Likelihood(),
                    params.getParameterValues()
                    )
                )
            );
    }
};

} // namespace golem

#endif
