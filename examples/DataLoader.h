#ifndef GOLEM_EXAMPLE_DATA_LOADER_H
#define GOLEM_EXAMPLE_DATA_LOADER_H

#include <PhysTools/tableio.h>
#include <LeptonWeighter/ParticleType.h>

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
        LW::ParticleType initialType;
        double finalStateX;
        double finalStateY;
        double totalColumnDepth;
        double radius;
        double z;
    }

    template<typename CallbackType>
    void readFile(const std::string& filePath, CallbackType action){
        H5File h5file(filePath);
        if(!h5file)
            throw std::runtime_error("Unable to open "+filePath);
        std::set<std::string> tables;
        H5Giterate(h5file,"/",NULL,&collectTableNames,&tables);
        if(tables.empty())
            throw std::runtime_error(filePath+" contains no tables");
        std::map<RecordID,Event> intermediateData;

        using particle = TableRow<
            field<double,CTS("totalEnergy")>,
            field<double,CTS("zenith")>,
            field<double,CTS("azimuth")>,
            field<double,CTS("finalStateX")>,
            field<double,CTS("finalStateY")>,
            field<int,CTS("finalType1")>,
            field<int,CTS("finalType2")>,
            field<int,CTS("initialType")>,
            field<double,CTS("totalColumnDepth")>,
            field<double,CTS("radius")>,
            field<double,CTS("z")>
                >;

        if(tables.count("EventProperties")){
            readTable<particle>(h5file, "EventProperties", intermediateData,
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

        for(std::map<RecordID,Event>::value_type& item : intermediateData)
            action(item.first,item.second);
    }

    std::deque<Event> GetSimulationEvents(){
        std::deque<Event> simulation_events;
        try {
            readFile(simulation_data_path,
                    [&](RecordID id, Event& e){
                    simulation_events.push_back(e);
                    }
                    );
        } catch ( std::exception & ex){
            std::cerr << ex.what() << std::endl;
        }
        return simulation_events;
    }

    std::deque<Event> GetDataEvents(){
        std::deque<Event> observation_events;
        try {
            readFile(observation_data_path,
                    [&](RecordID id, Event& e){
                    observation_events.push_back(e);
                    }
                    );
        } catch ( std::exception & ex){
            std::cerr << ex.what() << std::endl;
        }
        return observation_events;
    }
};
