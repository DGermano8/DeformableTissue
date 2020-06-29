#include "UniformDistributionSimpleWntCellCycleModel.hpp"
#include "PetscTools.hpp"
#include "Debug.hpp"
#include "SimulationTime.hpp"

UniformDistributionSimpleWntCellCycleModel::UniformDistributionSimpleWntCellCycleModel()
    : mUseCellProliferativeTypeDependentG1Duration(false),
      mWntStemThreshold(0.8),
      mWntTransitThreshold(0.65),
      mWntLabelledThreshold(0.65)
{
}


AbstractCellCycleModel* UniformDistributionSimpleWntCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell cycle model
    UniformDistributionSimpleWntCellCycleModel* p_model = new UniformDistributionSimpleWntCellCycleModel();

    /*
     * Set each member variable of the new cell-cycle model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new cell-cycle model's member variables (namely
     * mBirthTime, mCurrentCellCyclePhase, mReadyToDivide) will already have been
     * correctly initialized in its constructor.
     *
     * Note 2: one or more of the new cell-cycle model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new cell-cycle model.
     */
    p_model->SetBirthTime(mBirthTime);
    p_model->SetDimension(mDimension);
    p_model->SetMinimumGapDuration(mMinimumGapDuration);
    p_model->SetStemCellG1Duration(mStemCellG1Duration);
    p_model->SetTransitCellG1Duration(mTransitCellG1Duration);
    p_model->SetSDuration(mSDuration);
    p_model->SetG2Duration(mG2Duration);
    p_model->SetMDuration(mMDuration);
    p_model->SetUseCellProliferativeTypeDependentG1Duration(mUseCellProliferativeTypeDependentG1Duration);
    p_model->SetWntStemThreshold(mWntStemThreshold);
    p_model->SetWntTransitThreshold(mWntTransitThreshold);
    p_model->SetWntLabelledThreshold(mWntLabelledThreshold);

    return p_model;
}


void UniformDistributionSimpleWntCellCycleModel::SetUseCellProliferativeTypeDependentG1Duration(bool useCellProliferativeTypeDependentG1Duration)
{
    mUseCellProliferativeTypeDependentG1Duration = useCellProliferativeTypeDependentG1Duration;
}


void UniformDistributionSimpleWntCellCycleModel::SetG1Duration()
{
    assert(mpCell != NULL);

    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
    {
        if (mUseCellProliferativeTypeDependentG1Duration)
        {
			mG1Duration = 1 + 4*p_gen->ranf(); // U[1,5] according to Meineke
		}
		else
		{
			// Normally stem cells should behave just like transit cells in a Wnt simulation
			mG1Duration = 1 + 2*p_gen->ranf(); // U[1,3] according to Meineke
		}
    }
    else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
    {
		mG1Duration = 1 + 2*p_gen->ranf(); // U[1,3] according to Meineke
    }
    else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
		mG1Duration = DBL_MAX;
    }
    else
    {
        NEVER_REACHED;
    }

    // Check that the normal random deviate has not returned a small or negative G1 duration
    if (mG1Duration < mMinimumGapDuration)
    {
        mG1Duration = mMinimumGapDuration;
    }
}


double UniformDistributionSimpleWntCellCycleModel::GetWntLevel()
{
    assert(mpCell != NULL);
    double level = 0;

    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);            
            if (level > 1.0)
            {
            	PRINT_VARIABLE(level);
            }
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);
            break;
        }
        default:
            NEVER_REACHED;
    }
    return level;
}


WntConcentrationType UniformDistributionSimpleWntCellCycleModel::GetWntType()
{
    WntConcentrationType wnt_type;
    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            wnt_type = WntConcentration<DIM>::Instance()->GetType();
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            wnt_type = WntConcentration<DIM>::Instance()->GetType();
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            wnt_type = WntConcentration<DIM>::Instance()->GetType();
            break;
        }
        default:
            NEVER_REACHED;
    }
    return wnt_type;
}


void UniformDistributionSimpleWntCellCycleModel::UpdateCellCyclePhase()
{
    // The cell can divide if the Wnt concentration >= wnt_division_threshold
    double wnt_division_threshold = DBL_MAX;
    
    // Set up under what level of Wnt stimulus a cell will divide
    double healthy_threshold = mWntTransitThreshold;
    double labelled_threshold = mWntLabelledThreshold;

    if (mpCell->GetMutationState()->IsType<WildTypeCellMutationState>())
    {
        wnt_division_threshold = healthy_threshold;
    }
    else if (mpCell->GetMutationState()->IsType<ApcOneHitCellMutationState>())
    {
        // should be less than healthy values
        wnt_division_threshold = 0.77*healthy_threshold;
    }
    else if (mpCell->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())
    {
        // less than above value
        wnt_division_threshold = 0.155*healthy_threshold;
    }
    else if (mpCell->GetMutationState()->IsType<ApcTwoHitCellMutationState>())
    {
        // should be zero (no Wnt-dependence)
        wnt_division_threshold = 0.0;
    }
    else if (mpCell->GetMutationState()->IsType<StromalCellMutationState>())
    {
        // should be zero (no Wnt-dependence)
        wnt_division_threshold = 0.0;
    }
    else
    {
        NEVER_REACHED;
    }

    if (mpCell->HasCellProperty<CellLabel>())
    {
        wnt_division_threshold = labelled_threshold;
    }

    double wnt_level = GetWntLevel();
    WntConcentrationType wnt_type = GetWntType();

    // Set the cell type to TRANSIT if the Wnt stimulus exceeds wnt_division_threshold - don't want 
    // any stromal cells to suddenly become transit though
    if ( (wnt_level >= wnt_division_threshold) && (!mpCell->GetMutationState()->IsType<StromalCellMutationState>()) )
    {
        // For a RADIAL Wnt type, override the cell type to StemCellProliferativeType if the Wnt stimulus exceeds a higher threshold
        if ((wnt_type == RADIAL) && (wnt_level > mWntStemThreshold))
        {
            /*
             * This method is usually called within a CellBasedSimulation, after the CellPopulation
             * has called CellPropertyRegistry::TakeOwnership(). This means that were we to call
             * CellPropertyRegistry::Instance() here when setting the CellProliferativeType, we
             * would be creating a new CellPropertyRegistry. In this case the cell proliferative
             * type counts, as returned by AbstractCellPopulation::GetCellProliferativeTypeCount(),
             * would be incorrect. We must therefore access the CellProliferativeType via the cell's
             * CellPropertyCollection.
             */
            boost::shared_ptr<AbstractCellProperty> p_stem_type =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<StemCellProliferativeType>();
            mpCell->SetCellProliferativeType(p_stem_type);
        }
		else
		{
			boost::shared_ptr<AbstractCellProperty> p_transit_type =
				mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
			mpCell->SetCellProliferativeType(p_transit_type);
		}
    }
	else
	{
		// The cell is set to have DifferentiatedCellProliferativeType and so in G0 phase
		boost::shared_ptr<AbstractCellProperty> p_diff_type =
			mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
		mpCell->SetCellProliferativeType(p_diff_type);
	}

    AbstractSimpleCellCycleModel::UpdateCellCyclePhase();
}

void UniformDistributionSimpleWntCellCycleModel::InitialiseDaughterCell()
{
    WntConcentrationType wnt_type = GetWntType();

    if (wnt_type == RADIAL)
    {
        boost::shared_ptr<AbstractCellProperty> p_transit_type =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_transit_type);
    }

    AbstractSimpleCellCycleModel::InitialiseDaughterCell();
}

bool UniformDistributionSimpleWntCellCycleModel::CanCellTerminallyDifferentiate()
{
    return false;
}

double UniformDistributionSimpleWntCellCycleModel::GetWntStemThreshold()
{
    return mWntStemThreshold;
}

void UniformDistributionSimpleWntCellCycleModel::SetWntStemThreshold(double wntStemThreshold)
{
    assert(wntStemThreshold <= 1.0);
    assert(wntStemThreshold >= 0.0);
    mWntStemThreshold = wntStemThreshold;
}

double UniformDistributionSimpleWntCellCycleModel::GetWntTransitThreshold()
{
    return mWntTransitThreshold;
}

void UniformDistributionSimpleWntCellCycleModel::SetWntTransitThreshold(double wntTransitThreshold)
{
//    assert(wntTransitThreshold <= 1.0);
//    assert(wntTransitThreshold >= 0.0);
    mWntTransitThreshold = wntTransitThreshold;
}

double UniformDistributionSimpleWntCellCycleModel::GetWntLabelledThreshold()
{
    return mWntLabelledThreshold;
}

void UniformDistributionSimpleWntCellCycleModel::SetWntLabelledThreshold(double wntLabelledThreshold)
{
    assert(wntLabelledThreshold <= 1.0);
    assert(wntLabelledThreshold >= 0.0);
    mWntLabelledThreshold = wntLabelledThreshold;
}

void UniformDistributionSimpleWntCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile <<  "\t\t\t<WntStemThreshold>"<< mWntStemThreshold << "</WntStemThreshold>\n";
    *rParamsFile <<  "\t\t\t<WntTransitThreshold>"<< mWntTransitThreshold << "</WntTransitThreshold>\n";
    *rParamsFile <<  "\t\t\t<WntLabelledThreshold>"<< mWntLabelledThreshold << "</WntLabelledThreshold>\n";

    // Call direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(UniformDistributionSimpleWntCellCycleModel)
