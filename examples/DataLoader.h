#ifndef GOLEM_EXAMPLE_DATA_LOADER_H
#define GOLEM_EXAMPLE_DATA_LOADER_H

#include <deque>
#include <PhysTools/tableio.h>
#include <LeptonWeighter/ParticleType.h>

namespace analysis {

class DataLoader {
private:
    std::string simulation_data_path;
    std::string observation_data_path;
public:
    DataLoader(std::string simulation_data_path, std::string observation_data_path):
	simulation_data_path(simulation_data_path), observation_data_path(observation_data_path) {}

    struct Event {
        double primaryEnergy;
        double primaryZenith;
        double primaryAzimuth;
        double energy;
        double zenith;
        double azimuth;
        LW::ParticleType final_state_particle_0;
        LW::ParticleType final_state_particle_1;
        LW::ParticleType primary_type;
        double interaction_x;
        double interaction_y;
        double total_column_depth;
        double radius;
        double x;
        double y;
        double z;
    };
protected:
    void readFile(const std::string& filePath,
                  std::function<void(phys_tools::tableio::RecordID,Event&)> action) const;
public:
    std::deque<Event> GetSimulationEvents() const;
    std::deque<Event> GetDataEvents() const;
};

} // namespace analysis

#endif // GOLEM_EXAMPLE_DATA_LOADER_H

