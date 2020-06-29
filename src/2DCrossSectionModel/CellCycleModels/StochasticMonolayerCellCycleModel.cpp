#include "PetscTools.hpp"
#include "StochasticMonolayerCellCycleModel.hpp"

AbstractCellCycleModel* StochasticMonolayerCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell cycle model
    StochasticMonolayerCellCycleModel* p_model = new StochasticMonolayerCellCycleModel();

    // Set the values of the new cell cycle model's member variables
    p_model->SetGeneration(mGeneration);
    p_model->SetMaxTransitGenerations(mMaxTransitGenerations);
    p_model->SetCellProliferativeType(mCellProliferativeType);

    return p_model;
}

void StochasticMonolayerCellCycleModel::SetG1Duration()
{
    assert(mpCell!=NULL);

    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    switch (mpCell->GetCellCycleModel()->GetCellProliferativeType())
    {
		case STEM:
			mG1Duration = 1 + 4*p_gen->ranf(); // U[1,5] according to Meineke
			break;
		case TRANSIT:
			mG1Duration = 1 + 2*p_gen->ranf(); // U[1,3] according to Meineke
			break;
		case DIFFERENTIATED:
			mG1Duration = DBL_MAX;
			break;
		default:
			NEVER_REACHED;
    }
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(StochasticMonolayerCellCycleModel)
