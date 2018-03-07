#include "DataLoader.h"

namespace analysis {
	
void DataLoader::readFile(const std::string& filePath,
                          std::function<void(phys_tools::tableio::RecordID,Event&)> action) const{
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

std::deque<DataLoader::Event> DataLoader::GetSimulationEvents() const{
	std::deque<Event> simulation_events;
	try {
		readFile(simulation_data_path,
		        [&](phys_tools::tableio::RecordID id, Event& e){
		        simulation_events.push_back(e);
		        });
	} catch (std::exception & ex){
		std::cerr << ex.what() << std::endl;
	}
	return simulation_events;
}

std::deque<DataLoader::Event> DataLoader::GetDataEvents() const{
	std::deque<Event> observation_events;
	try {
		readFile(observation_data_path,
		        [&](phys_tools::tableio::RecordID id, Event& e){
		        observation_events.push_back(e);
		        });
	} catch (std::exception & ex){
		std::cerr << ex.what() << std::endl;
	}
	return observation_events;
}
	
} // namespace analysis
