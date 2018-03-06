#ifndef __PHYSICS_MODEL_H_
#define __PHYSICS_MODEL_H_

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
      template<typename DataType>
        std::function<DataType(const Event&)> operator()(const std::vector<DataType>& params) const {
          assert(parameters.size() == NParameters);
          cachedValueWeighter<Event,DataType,double> w(Event::weight);
          return w*params[0];
        }
    };

    struct UncertaintyWeighter {
      template<typename DataType>
        DataType operator()(const Event& e, const DataType& w) const {
          return w*w;
        }
    };

    struct Prior {
      phys_tools::likelihood::GaussianPrior normalization;
      Prior(): normalization(1.,0.1) {}
      template<typename DataType>
        DataType operator()(const std::vector<DataType>& parameters) const {
          assert(parameters.size() == NParameters);
          return normalization(parameters[0]);;
        }
    };

    void AddEventToHistogram(HistogramSet& h, const Event& e) const {
      h.add(e.energy, cos(e.zenith), phys_tools::histograms::amount(std::cref(e)));
    }
};

#endif

