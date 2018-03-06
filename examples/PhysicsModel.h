#ifndef GOLEM_EXAMPLE_PHYSICS_MODEL_H
#define GOLEM_EXAMPLE_PHYSICS_MODEL_H

class PhysicsModel {
private:
public:
    PhysicsModel(){}

    static constexpr unsigned int NParameters = 1;

    struct Event {
        double primaryEnergy;
        double primaryZenith;
        double energy;
        double zenith;
        double weight;
    }

    typedef std::tuple<
        phys_tools::histograms::histogram<2, phys_tools::likelihood::entryStoringBin<std::reference_wrapper<const Event>>>
        > HistogramSet;

    struct WeighterMaker {
        WeighterMaker() {}
        template<typename DataType>
            std::function<DataType(const Event&)> operator()(const std::vector<DataType>& params) const {
                assert(parameters.size() == NParameters);
                cachedValueWeighter<Event,DataType,double> w(Event::weight);
                DataType normalization = p.extractParameter("normalization",parameters)
                    return w*normalization;
            }
    };

    struct UncertaintyWeighter {
        UncertaintyWeighter() {}
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
