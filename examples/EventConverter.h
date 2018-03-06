#ifndef GOLEM_EXAMPLE_EVENT_CONVERTER_H
#define GOLEM_EXAMPLE_EVENT_CONVERTER_H

#include <LeptonWeighter/Weighter.h>

#include "DataLoader.h"
#include "PhysicsModel.h"

class EventConverter {
  private:
    std::unique_ptr<LW::Weighter> weighter;
  protected:
    double GetEventWeight(DataLoader::Event& e_input){
      LW::Event lw_e {
        e.primaryType,
        e.final_state_particle_0,
        e.final_state_particle_1,
        e.intX,
        e.intY,
        e.primaryEnergy,
        e.primaryAzimuth,
        e.primaryZenith,
        e.x,
        e.y,
        e.z,
        e.r,
        e.totalColumnDepth
      };
      return weighter(lw_e);
    }
  public:
    EventConverter(std::string path_to_lic_file, std::string path_to_cross_sections){
      std::vector<std::shared_ptr<LW::Generator>> generators = LW::MakeGeneratorsFromLICFile(path_to_lic_file);
      std::shared_ptr<LW::CrossSectionFromSpline> xs = std::make_shared<LW::CrossSectionFromSpline>(path_to_cross_sections + "/dsdxdy_nu_CC_iso.fits",
                                                                                                    path_to_cross_sections + "/dsdxdy_nubar_CC_iso.fits",
                                                                                                    path_to_cross_sections + "/dsdxdy_nu_NC_iso.fits",
                                                                                                    path_to_cross_sections + "/dsdxdy_nubar_NC_iso.fits");
      weighter.reset(new LW::Weighter(xs,generators));
    }

    PhysicsModel::Event ConvertEvent(DataLoader::Event& e_input){
      PhysicsModel::Event e_output;

      e_output.energy = e_input.energy;
      e_output.zenith = e_input.zenith;
      e_output.primaryEnergy = e_input.primaryEnergy;
      e_output.primaryZenith = e_input.primaryZenith;
      e_output.weight = GetEventWeight(e_input);
      return e_output;
    }
}

#endif
