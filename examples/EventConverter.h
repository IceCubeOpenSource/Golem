#ifndef GOLEM_EXAMPLE_EVENT_CONVERTER_H
#define GOLEM_EXAMPLE_EVENT_CONVERTER_H

#include <deque>
#include <vector>
#include <LeptonWeighter/Weighter.h>

#include "DataLoader.h"
#include "PhysicsModel.h"

namespace analysis {

/// \brief Event Converter class.
class EventConverter {
private:
    std::unique_ptr<LW::Weighter> weighter;
protected:
    /// \brief Gets a data loader event and returns its weight.
    /// @param e_input Event from the data loader.
    /// \details Transfer properties from the DataLoader::Event to a LeptonWeighter::Event and apply the weighter.
    double GetEventWeight(DataLoader::Event& e){
        LW::Event lw_e {
          e.primary_type,
          e.final_state_particle_0,
          e.final_state_particle_1,
          e.interaction_x,
          e.interaction_y,
          e.primaryEnergy,
          e.primaryAzimuth,
          e.primaryZenith,
          e.x,
          e.y,
          e.z,
          e.radius,
          e.total_column_depth
        };
        return (*weighter)(lw_e);
    }
public:
    /// \brief Event converter constructor.
    /// @param path_to_lic_file Path to Lepton Injector Configuration (LIC) file.
    /// @param path_to_cross_sections Path to the folder where cross section splines are located.
    /// \details Generated the generator and cross section weighters, which then are assambled
    /// into the LeptonWeighter::Weighter.
    EventConverter(std::string path_to_lic_file, std::string path_to_cross_sections){
        std::vector<std::shared_ptr<LW::Generator>> generators = LW::MakeGeneratorsFromLICFile(path_to_lic_file);
        std::shared_ptr<LW::CrossSectionFromSpline> xs = std::make_shared<LW::CrossSectionFromSpline>(path_to_cross_sections + "/dsdxdy_nu_CC_iso.fits",
                path_to_cross_sections + "/dsdxdy_nubar_CC_iso.fits",
                path_to_cross_sections + "/dsdxdy_nu_NC_iso.fits",
                path_to_cross_sections + "/dsdxdy_nubar_NC_iso.fits");
        // Constructs and sets the weighter
        weighter.reset(new LW::Weighter(xs,generators));
    }

    /// \brief Converts one data loader event structure to a physics model event structure.
    /// @param e_input Event from the data loader.
    /// \details This function takes one event from the data loader and extracts
    /// the properties needed for the analysis. To do so it dumps the content into
    /// PhysicsModel::Event structure. In this particular example the analizer is only
    /// needs energy and zenith information. Additionally the weight is calculated by
    /// LeptonWeighter and store in the object.
    PhysicsModel::Event ConvertEvent(DataLoader::Event& e_input){
        PhysicsModel::Event e_output;

        e_output.energy = e_input.energy;
        e_output.zenith = e_input.zenith;
        e_output.primaryEnergy = e_input.primaryEnergy;
        e_output.primaryZenith = e_input.primaryZenith;
        e_output.weight = GetEventWeight(e_input);
        return e_output;
    }
};

} // namespace analysis

#endif // GOLEM_EXAMPLE_EVENT_CONVERTER_H
