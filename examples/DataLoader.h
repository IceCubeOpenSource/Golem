#ifndef GOLEM_EXAMPLE_DATA_LOADER_H
#define GOLEM_EXAMPLE_DATA_LOADER_H

#include <deque>
#include <vector>
#include <PhysTools/tableio.h>
#include <LeptonWeighter/ParticleType.h>
#include <hdf5.h>
#include <hdf5_hl.h>

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

    template<typename CallbackType>
    void readFile(const std::string& filePath, CallbackType action) const{
        using namespace phys_tools::cts;
        phys_tools::tableio::H5File h5file(filePath);
        if(!h5file)
            throw std::runtime_error("Unable to open "+filePath);
        std::set<std::string> tables=phys_tools::tableio::getTables(h5file,"/");
        if(tables.empty())
            throw std::runtime_error(filePath+" contains no tables");
        std::map<phys_tools::tableio::RecordID,Event> intermediateData;

        using particle = phys_tools::tableio::TableRow<
            phys_tools::tableio::field<double,CTS("totalEnergy")>,
            phys_tools::tableio::field<double,CTS("zenith")>,
            phys_tools::tableio::field<double,CTS("azimuth")>,
            phys_tools::tableio::field<double,CTS("finalStateX")>,
            phys_tools::tableio::field<double,CTS("finalStateY")>,
            phys_tools::tableio::field<int,CTS("finalType1")>,
            phys_tools::tableio::field<int,CTS("finalType2")>,
            phys_tools::tableio::field<int,CTS("initialType")>,
            phys_tools::tableio::field<double,CTS("totalColumnDepth")>,
            phys_tools::tableio::field<double,CTS("radius")>,
            phys_tools::tableio::field<double,CTS("z")>
                >;

        if(tables.count("EventProperties")){
          phys_tools::tableio::readTable<particle>(h5file, "EventProperties", intermediateData,
                    [](const particle& p, Event& e){
                    e.primaryEnergy=p.get<CTS("totalEnergy")>();
                    e.primaryZenith=p.get<CTS("zenith")>();
                    e.primaryAzimuth=p.get<CTS("azimuth")>();
                    e.energy=p.get<CTS("totalEnergy")>();
                    e.zenith=p.get<CTS("zenith")>();
                    e.azimuth=p.get<CTS("azimuth")>();
                    e.interaction_x=p.get<CTS("finalStateX")>();
                    e.interaction_y=p.get<CTS("finalStateY")>();
                    e.final_state_particle_0=static_cast<LW::ParticleType>(p.get<CTS("finalType1")>());
                    e.final_state_particle_1=static_cast<LW::ParticleType>(p.get<CTS("finalType2")>());
                    e.primary_type=static_cast<LW::ParticleType>(p.get<CTS("initialType")>());
                    e.total_column_depth=p.get<CTS("totalColumnDepth")>();
                    e.radius=p.get<CTS("radius")>();
                    e.z=p.get<CTS("z")>();
                    });
        }

        for(std::map<phys_tools::tableio::RecordID,Event>::value_type& item : intermediateData)
            action(item.first,item.second);
    }
public:
    std::deque<Event> GetSimulationEvents() const;
	std::deque<Event> GetDataEvents() const;
};

} // namespace analysis

#endif // GOLEM_EXAMPLE_DATA_LOADER_H

