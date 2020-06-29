#ifndef GECODESYSTEM_HPP_
#define GECODESYSTEM_HPP_

#include "CellwiseOdeSystemInformation.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include <cmath>
#include <iostream>
#include "AbstractOdeSystem.hpp"

class GecOdeSystem : public AbstractOdeSystem

{
	public:

//	friend class boost::serialization::access;

//	template<class Archive>
//	void serialize(Archive & archive, const unsigned int version)
//	{
//		archive & boost::serialization::base_object<AbstractOdeSystem>(*this);
//	}

	GecOdeSystem(std::vector<double> stateVariables=std::vector<double>());

	~GecOdeSystem();

	void Init();

	void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);
};


#include "SerializationExportWrapper.hpp"

CHASTE_CLASS_EXPORT(GecOdeSystem)

namespace boost
{
	namespace serialization
	{

		/**

		* Serialize information required to construct the system.

		*/

		template<class Archive>

		inline void save_construct_data(

		Archive & ar, const GecOdeSystem * t, const BOOST_PFTO unsigned int file_version)

		{
			const std::vector<double> state_variables = t->rGetConstStateVariables();

			ar & state_variables;
		}

		/**

		* De-serialize constructor parameters and initialise a system.

		*/

		template<class Archive>

		inline void load_construct_data(

		Archive & ar, GecOdeSystem * t, const unsigned int file_version)
		{

			std::vector<double> state_variables;
			ar & state_variables;

			// Invoke inplace constructor to initialise instance

			::new(t)GecOdeSystem(state_variables);
		}
	}
}

GecOdeSystem::GecOdeSystem(std::vector<double> stateVariables)
                : AbstractOdeSystem(3)
{

	mpSystemInfo.reset(new CellwiseOdeSystemInformation<GecOdeSystem>);

	Init();

	if (stateVariables != std::vector<double>())
	{
		SetStateVariables(stateVariables);
	}
}

GecOdeSystem::~GecOdeSystem()
{

}

void GecOdeSystem::Init()
{
	//initialize member variables
}

void GecOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
	// Species names map to
	// Create named variables for current concentration/molecule number

	double X = rY[0];

	double Y1 = rY[1];

	double Y2 = rY[2];

	// Assign derivatives

	rDY[0] =0;
	rDY[1] = 10*X*Y1 - 0.01*Y1*Y2;
	rDY[2] = 0.01*Y1*Y2 - 10*Y2;

}

/*bool GecOdeSystem::CalculateStoppingEvent(double time, const std::vector<double>& rY)

{

}*/

/*double GecOdeSystem::CalculateRootFunction(double time, const std::vector<double>& rY)

{

}*/

// Assign initial conditions

template<>

void CellwiseOdeSystemInformation<GecOdeSystem>::Initialise()

{
	// Names, units & initial conditions for state variables.

	this->mVariableNames.push_back("X");
	this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0); // will be filled in later

	this->mVariableNames.push_back("Y1");
	this->mVariableUnits.push_back("dimensionless");
	this->mInitialConditions.push_back(0.0); // will be filled in later

	this->mVariableNames.push_back("Y2");
	this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0); // will be filled in later

	this->mInitialised = true;
}

// Serialization for Boost >= 1.36

#include "SerializationExportWrapperForCpp.hpp"

CHASTE_CLASS_EXPORT(GecOdeSystem)

#endif /*GECODESYSTEM_HPP_*/
