#!/bin/sh

if [ "$#" -lt 2 ]; then
	echo "Usage: create_project name path"
	exit 1
fi

PROJ_NAME=$1
PROJ_PATH=$2

PROJ_DIR="${PROJ_PATH}/${PROJ_NAME}"

echo "Generating $PROJ_DIR"

mkdir -p "$PROJ_DIR"
if [ "$?" -ne 0 ]; then
	echo "Failed to create project directory $PROJ_DIR" 1>&2
	exit 1
fi

mkdir -p "${PROJ_DIR}/include/${PROJ_NAME}"
if [ "$?" -ne 0 ]; then
	echo "Failed to create project include directory ${PROJ_DIR}/include/${PROJ_NAME}" 1>&2
	exit 1
fi

mkdir -p "${PROJ_DIR}/src"
if [ "$?" -ne 0 ]; then
	echo "Failed to create project source directory ${PROJ_DIR}/src" 1>&2
	exit 1
fi

#===============================================================================
# *PhysicsModel.h
#===============================================================================

PHYSICSMODEL_HEADER="${PROJ_DIR}/include/${PROJ_NAME}/${PROJ_NAME}PhysicsModel.h"

echo "#ifndef ${PROJ_NAME}_PHYSICSMODEL_H
#define ${PROJ_NAME}_PHYSICSMODEL_H

#include <deque>
#include <functional>
#include <tuple>
#include <vector>
#include <PhysTools/histogram.h>
#include <PhysTools/likelihood/likelihood.h>
#include <PhysTools/likelihood/physics_weighters.h>
#include <PhysTools/optimization/ParameterSet.h>

class ${PROJ_NAME}PhysicsModel{
public:
	//Change the constructor to take any additional arguments you want
	${PROJ_NAME}PhysicsModel();

	//!!! Change this to the number of parameters in your model
	static constexpr unsigned int NParameters = 0;

	struct Event{
		//!!! Put variables here for your observables and 
		//    anything you need to compute your weights
	};

	struct WeighterMaker{
		const ParameterSet& p;

		WeighterMaker(const ParameterSet& p): p(p) {}

		template<typename DataType>
		std::function<DataType(const Event&)> operator()(const std::vector<DataType>& params) const {
			//!!! Put the code to implement your weight calculation here
		}
	};

	struct UncertaintyWeighter{
		UncertaintyWeighter() {}
		
		template<typename DataType>
		DataType operator()(const Event& e, const DataType& w) const {
			return w*w;
		}
	};

	struct Prior{
		const ParameterSet& p;
		//Add any other members you need here

		Prior(const ParameterSet& p):p(p)
		//Initialize any members you added here
		{}

		template<typename DataType>
		DataType operator()(const std::vector<DataType>& parameters) const {
			//Compute and return any prior you want here.
			return 1.;
		}
	};

	template<unsigned int Dimension>
	using histogram=phys_tools::histograms::histogram<Dimension, phys_tools::likelihood::entryStoringBin<std::reference_wrapper<const Event>>>;

	using HistogramSet = std::tuple<
		//!!! Put the types of histograms you want to fit here
	>;

	HistogramSet MakeHistogramSet() const;

	void AddEventToHistogram(HistogramSet& h, const Event& e) const;
	
	phys_tools::ParameterSet MakeParameterSet() const;

	WeighterMaker MakeWeigherMaker(const ParameterSet& params) const;

	UncertaintyWeighter MakeUncertaintyWeighter() const;

	Prior MakePrior(const ParameterSet& params) const;

private:
	//Add any member variables here
};

#endif //${PROJ_NAME}_PHYSICSMODEL_H
" > $PHYSICSMODEL_HEADER
if [ "$?" -ne 0 ]; then
	echo "Failed to write template header file $PHYSICSMODEL_HEADER" 1>&2
	exit 1
fi

#===============================================================================
# *PhysicsModel.cpp
#===============================================================================
PHYSICSMODEL_IMPL="${PROJ_DIR}/src/${PROJ_NAME}PhysicsModel.cpp"
echo "#include <${PROJ_NAME}PhysicsModel.h>

//Change the constructor to initialize any member variables that you add
${PROJ_NAME}PhysicsModel::${PROJ_NAME}PhysicsModel(){}

${PROJ_NAME}PhysicsModel::HistogramSet 
${PROJ_NAME}PhysicsModel::MakeHistogramSet() const{
	//!!! Construct your histograms here
	
	HistogramSet h = std::make_tuple( /* insert histograms */ );
	return h;
}

void ${PROJ_NAME}PhysicsModel::AddEventToHistogram(HistogramSet& h, const Event& e) const{
	//!!! Write code here to put e into h
}

phys_tools::ParameterSet ${PROJ_NAME}PhysicsModel::MakeParameterSet() const{
	phys_tools::ParameterSet parameters;

	//!!! Define your parameters here
	//The number of parameters you add should match NParameters

	return parameters;
}

${PROJ_NAME}PhysicsModel::WeighterMaker
${PROJ_NAME}PhysicsModel::MakeWeigherMaker(const ParameterSet& params) const{
	//Change this code if your WeighterMaker takes any other arguments
	return WeighterMaker(params);
}

${PROJ_NAME}PhysicsModel::UncertaintyWeighter 
${PROJ_NAME}PhysicsModel::MakeUncertaintyWeighter() const{
	return UncertaintyWeighter();
}

${PROJ_NAME}PhysicsModel::Prior 
${PROJ_NAME}PhysicsModel::MakePrior(const ParameterSet& params) const{
	//Change this code if your Prior takes any other arguments
	return Prior(params);
}
" > $PHYSICSMODEL_IMPL

#===============================================================================
# *DataLoader.h
#===============================================================================

DATALOADER_HEADER="${PROJ_DIR}/include/${PROJ_NAME}/${PROJ_NAME}DataLoader.h"

echo "#ifndef ${PROJ_NAME}_DATALOADER_H
#define ${PROJ_NAME}_DATALOADER_H

#include <deque>
#include <PhysTools/tableio.h>
#include <LeptonWeighter/ParticleType.h>

class ${PROJ_NAME}DataLoader{
public:
	struct Event{
		//Define the properties you want for each event you will read in
	};
	${PROJ_NAME}DataLoader();
protected:
	void readFile(const std::string& filePath,
                  std::function<void(phys_tools::tableio::RecordID,Event&)> action) const;
public:
	std::deque<Event> GetSimulationEvents() const;
	std::deque<Event> GetDataEvents() const;
};

#endif //${PROJ_NAME}_DATALOADER_H
" > $DATALOADER_HEADER

#===============================================================================
# *DataLoader.cpp
#===============================================================================
DATALOADER_IMPL="${PROJ_DIR}/src/${PROJ_NAME}DataLoader.cpp"
echo "#include <${PROJ_NAME}DataLoader.h>

${PROJ_NAME}DataLoader::${PROJ_NAME}DataLoader(){}

void 
${PROJ_NAME}DataLoader::readFile(const std::string& filePath,
  std::function<void(phys_tools::tableio::RecordID,Event&)> action) const{
	using namespace phys_tools::cts;
	phys_tools::tableio::H5File h5file(filePath);
	if(!h5file)
		throw std::runtime_error(\"Unable to open \"+filePath);
	std::set<std::string> tables=phys_tools::tableio::getTables(h5file,\"/\");
	if(tables.empty())
		throw std::runtime_error(filePath+\" contains no tables\");
	std::map<phys_tools::tableio::RecordID,Event> intermediateData;

	//Put code here to read all tables of interest from the file

	for(std::map<phys_tools::tableio::RecordID,Event>::value_type& item : intermediateData)
		action(item.first,item.second);
}

std::deque<Event> ${PROJ_NAME}DataLoader::GetSimulationEvents() const{
	//Write code here to read whatever files contain the simulated data
}

std::deque<Event> ${PROJ_NAME}DataLoader::GetDataEvents() const{
	//Write code here to read whatever files contain the observed data
}

" > $DATALOADER_IMPL
