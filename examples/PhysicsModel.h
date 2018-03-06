#ifndef GOLEM_EXAMPLE_PHYSICS_MODEL_H
#define GOLEM_EXAMPLE_PHYSICS_MODEL_H

class PhysicsModel {
private:
public:
    ///The PhysicsModel constructor.
    ///This might take properties that you want to build into your model, but 
    ///this simple example doesn't have any. 
    PhysicsModel(){}

    ///This is the number of parameters the model has.
    static constexpr unsigned int NParameters = 1;

    ///This is the representaion of one observed or simulated event used for
    ///fitting/evaluating this model. It should contain all of the observables
    ///which are relevant to the model, and any information needed to compute
    ///the weights needed to find the expectation for the model.
    struct Event {
        double primaryEnergy;
        double primaryZenith;
        double energy;
        double zenith;
        double weight;
    }

    ///This is a template alias for the general type of histograms we will use.
    ///The only remaining template parameter is the number of dimensions 
    ///(observables) the histogram will have.
    template<unsigned int Dimension>
    using histogram=phys_tools::histograms::histogram<Dimension, phys_tools::likelihood::entryStoringBin<std::reference_wrapper<const Event>>>;
    
    ///This type is the full collection of data histograms that will be used to
    ///evaluate the model. There will be one HistogramSet for the observation
    ///and one for the expectation. 
    using HistogramSet = std::tuple<
        histogram<2>
        >;

    struct WeighterMaker {
        template<typename DataType>
            std::function<DataType(const Event&)> operator()(const std::vector<DataType>& params) const {
                assert(parameters.size() == NParameters);
                cachedValueWeighter<Event,DataType,double> w(Event::weight);
                DataType normalization = p.extractParameter("normalization",parameters)
                    return w*normalization;
            }
    };

    struct UncertaintyWeighter {
        template<typename DataType>
            DataType operator()(const Event& e, const DataType& w) const {
                return w*w;
            }
    };

    struct Prior {
        const ParameterSet& p;
        phys_tools::likelihood::GaussianPrior normalization;
        Prior(const ParameterSet& p): p(p), normalization(1.,0.1){}
        template<typename DataType>
            DataType operator()(const std::vector<DataType>& parameters) const {
                assert(parameters.size() == NParameters);
                return normalization(p.extractParameter("normalization",parameters));
            }
    };

    void AddEventToHistogram(HistogramSet& h, const Event& e) const {
        h.add(e.energy, cos(e.zenith), phys_tools::histograms::amount(std::cref(e)));
    }

    HistogramSet MakeHistogramSet() const {
        using phys_tools::likelihood::entryStoringBin;
        using HistType = phys_tools::histograms::histogram<2, entryStoringBin<std::reference_wrapper<const Event>>>;
        LogarithmicAxis energy_axis(0,0.1);
        LinearAxis cos_zenith_axis(0,0.1);
        HistogramSet h = std::make_tuple(HistType(energy_axis,cos_zenith_axis));
        return h;
    }

    phys_tools::ParameterSet MakeParameterSet() const {
        phys_tools::ParameterSet parameters;
        parameters.addParameter("normalization");
        parameters.setParameterLowerLimit("normalization", 0);
        parameters.setParameterValue("normalization", 1);
        return parameters;
    }

    WeighterMaker MakeWeigherMaker() {
        return WeighterMaker();
    }

    UncertaintyWeighter MakeUncertaintyWeighter() {
        return UncertaintyWeighter();
    }
};

#endif //GOLEM_EXAMPLE_PHYSICS_MODEL_H
