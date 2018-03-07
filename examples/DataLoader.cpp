#include "DataLoader.h"

namespace analysis {

std::deque<DataLoader::Event> DataLoader::GetSimulationEvents() const{
	std::deque<Event> simulation_events;
	try {
		readFile(simulation_data_path,
				[&](phys_tools::tableio::RecordID id, Event& e){
				simulation_events.push_back(e);
				}
				);
	} catch ( std::exception & ex){
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
				}
				);
	} catch ( std::exception & ex){
		std::cerr << ex.what() << std::endl;
	}
	return observation_events;
}
	
} // namespace analysis
