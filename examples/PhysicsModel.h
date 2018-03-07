#ifndef GOLEM_EXAMPLE_PHYSICS_MODEL_H
#define GOLEM_EXAMPLE_PHYSICS_MODEL_H

#include <tuple>
#include <deque>
#include <vector>
#include <functional>
#include <PhysTools/histogram.h>
#include <PhysTools/likelihood/likelihood.h>
#include <PhysTools/likelihood/physics_weighters.h>
#include <PhysTools/optimization/ParameterSet.h>

namespace analysis {

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
    };

    ///This is the object which, given a particular set of model parameter values,
    ///makes an object which will weight the simulated events to that particular
    ///model.
    struct WeighterMaker {
        const phys_tools::ParameterSet& p;

        WeighterMaker(const phys_tools::ParameterSet& p): p(p) {}

        ///This is the function that will be called create a weighting object
        ///for a particular model. DataType is the type of the parameter values
        ///(which is often but not always `double`) and also the type that the
        ///calculated weights should have.
        template<typename DataType>
        std::function<DataType(const Event&)> operator()(const std::vector<DataType>& parameters) const {
            //Check for mistakes; make sure we got the right number of parameters
            assert(parameters.size() == NParameters);
            //In this example we suppose that all per-event weight information
            //is in the `weight` member of `Event`, so we just make a helper
            //object to fetch that.
            phys_tools::cachedValueWeighter<Event,DataType,double> w(&Event::weight);
            //The free parameter in the model is an overall normalization.
            DataType normalization = p.extractParameter("normalization",parameters);
            //So, we just scale all per-event weights by that factor overall.
            return normalization*w;
        }
    };

    ///This type is used to estimate the uncertainty on the weight associated
    ///with a simualted event.
    struct UncertaintyWeighter {
        UncertaintyWeighter() {}

        ///In the typical case the uncertainty (the variance) is just the square
        ///of the weight. However, in case something more complex is required
        ///this function is also given the event itself, from which it could get
        ///additional information.
        template<typename DataType>
        DataType operator()(const Event& e, const DataType& w) const {
            return w*w;
        }
    };

    ///This is the type which is used to represent any priors on the model
    ///parameters.
    struct Prior {
        const phys_tools::ParameterSet& p;
        phys_tools::likelihood::GaussianPrior normalization;

        ///Make gaussian prior for the normalization with a preferred (mean)
        ///value of 1.0 and a width (standard deviation) of 0.1.
        Prior(const phys_tools::ParameterSet& p): p(p), normalization(1.,0.01){}

        ///For a given set of model parameters evaluate the log of the prior
        template<typename DataType>
        DataType operator()(const std::vector<DataType>& parameters) const {
            //Check for mistakes; make sure we got the right number of parameters
            assert(parameters.size() == NParameters);
            //We already made the gaussian prior subobject for the normalization,
            //so we just call it with the normalization value that was given to us.
            return normalization(p.extractParameter("normalization",parameters));
        }
    };

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

    ///Having defined the type of our set of histograms (which fixes the number
    ///of them and their dimensions), we also need this function which will make
    ///an instance of that type, implementing the particular binning we want to
    ///use in each histogram.
    HistogramSet MakeHistogramSet() const {
        //Here we will just hard-code our binning choices. We could also have
        //member variables and setter functions to allow these to be changed at
        //runtime.
        phys_tools::histograms::LogarithmicAxis energy_axis(0,0.1); //ten bins per decade, based at 10^0
        phys_tools::histograms::LinearAxis cos_zenith_axis(0,0.1); //bin with a width of 0.1, based at 0
        //Create our one histogram with our preferred axes.
        HistogramSet h = std::make_tuple(histogram<2>(energy_axis,cos_zenith_axis));
        return h;
    }

    ///This function acts as an interface between our event type and our set of
    ///histograms: It defines how the propoerties stored in `Event` are used to
    ///select where that event belongs in the histogram(s).
    void AddEventToHistogram(HistogramSet& h, const Event& e) const {
        using phys_tools::histograms::amount;
        //We only have one histogram in our set, so all events go into it.
        //We chose to order our observables as reconstructed energy then cosine
        //of reconstructed zenith, so we pass those values to the histogram. We
        //pass a reference to the event itself as the 'amount' to put the event
        //itself into the histogram.
        std::get<0>(h).add(e.energy, cos(e.zenith), amount(std::cref(e)));
    }

    ///This function creates and provides an object describing all of the
    ///parameters in our model: Their names, any limits they have, default
    ///values, etc. The order in which we define the parameters here is the order
    ///which will be used for all vectors of parameter values, but we can use
    ///the convenience functions of the ParameterSet to refer to parameters by
    ///their names in almost all circumstances.
    phys_tools::ParameterSet MakeParameterSet() const {
        phys_tools::ParameterSet parameters;
        parameters.addParameter("normalization");
        //We can use limits to enforce simple physical constarints, like that a
        //flux normalization can never be negative.
        parameters.setParameterLowerLimit("normalization", 0);
        //This is the default or seed value for this parameter.
        parameters.setParameterValue("normalization", 1);
        //We could keep adding parameters here, but this example only needs one.
        return parameters;
    }

    ///This function constructs an instance of our WeighterMaker type. The Golem
    ///object will pass it a reference to the ParameterSet being used, but we
    ///can pass extra informtion of our own to the WeighterMaker constructor if
    ///we want to.
    WeighterMaker MakeWeighterMaker(const phys_tools::ParameterSet& params) {
        return WeighterMaker(params);
    }

    ///This function constructs an instance of our MakeUncertaintyWeighter type.
    ///It is assumed that the UncertaintyWeighter doesn't know or care about
    ///model parameters.
    UncertaintyWeighter MakeUncertaintyWeighter() {
        return UncertaintyWeighter();
    }

    ///This function makes a Prior instance. Like the WeighterMaker it probably
    ///wants to know about the model parameters, so this information will be
    ///passed back from the Golem, and we could pass extra information of or own
    ///to its constructor if we had any.
    Prior MakePrior(const phys_tools::ParameterSet& params) {
        return Prior(params);
    }
};

} // namespace analysis

#endif //GOLEM_EXAMPLE_PHYSICS_MODEL_H
