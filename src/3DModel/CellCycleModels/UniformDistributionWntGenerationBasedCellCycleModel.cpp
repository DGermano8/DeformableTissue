#include "UniformDistributionWntGenerationBasedCellCycleModel.hpp"
#include "PetscTools.hpp"

UniformDistributionWntGenerationBasedCellCycleModel::UniformDistributionWntGenerationBasedCellCycleModel()
    : mUseCellProliferativeTypeDependentG1Duration(false),
      mWntStemThreshold(0.8),
      mWntTransitThreshold(0.65),
      mWntLabelledThreshold(0.65),
      mGeneration(0),
      mMaxTransitGenerations(2)
{
}

AbstractCellCycleModel* UniformDistributionWntGenerationBasedCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell cycle model
    UniformDistributionWntGenerationBasedCellCycleModel* p_model = new UniformDistributionWntGenerationBasedCellCycleModel();

    // Set the values of the new cell cycle model's member variables
    p_model->SetDimension(mDimension);
    p_model->SetCellProliferativeType(mCellProliferativeType);
    p_model->SetMaxTransitGenerations(mMaxTransitGenerations);

    return p_model;
}


void UniformDistributionWntGenerationBasedCellCycleModel::SetGeneration(unsigned generation)
{
    mGeneration = generation;
}

unsigned UniformDistributionWntGenerationBasedCellCycleModel::GetGeneration() const
{
    return mGeneration;
}

void UniformDistributionWntGenerationBasedCellCycleModel::SetMaxTransitGenerations(unsigned maxTransitGenerations)
{
    mMaxTransitGenerations = maxTransitGenerations;
}

unsigned UniformDistributionWntGenerationBasedCellCycleModel::GetMaxTransitGenerations() const
{
    return mMaxTransitGenerations;
}

void UniformDistributionWntGenerationBasedCellCycleModel::SetUseCellProliferativeTypeDependentG1Duration(bool useCellProliferativeTypeDependentG1Duration)
{
    mUseCellProliferativeTypeDependentG1Duration = useCellProliferativeTypeDependentG1Duration;
}


void UniformDistributionWntGenerationBasedCellCycleModel::SetG1Duration()
{
    assert(mpCell != NULL);

    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    switch (mCellProliferativeType)
    {
        case STEM:
            if (mUseCellProliferativeTypeDependentG1Duration)
            {
                mG1Duration = 1 + 4*p_gen->ranf(); // U[1,5] according to Meineke
            }
            else
            {
                // Normally stem cells should behave just like transit cells in a Wnt simulation
                mG1Duration = 1 + 2*p_gen->ranf(); // U[1,3] according to Meineke
            }
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

    // Check that the normal random deviate has not returned a small or negative G1 duration
    if (mG1Duration < mMinimumGapDuration)
    {
        mG1Duration = mMinimumGapDuration;
    }
}


double UniformDistributionWntGenerationBasedCellCycleModel::GetWntLevel()
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


WntConcentrationType UniformDistributionWntGenerationBasedCellCycleModel::GetWntType()
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


void UniformDistributionWntGenerationBasedCellCycleModel::UpdateCellCyclePhase()
{
    // The cell can divide if the Wnt concentration >= wnt_division_threshold and if
	// its generation is less than the max
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

    unsigned generation = GetGeneration();

    // Set the cell type to TRANSIT if the Wnt stimulus exceeds wnt_division_threshold
    // and it hasn't reached the maximum number of generations
    if ( (wnt_level >= wnt_division_threshold) && (generation < mMaxTransitGenerations) )
    {
        CellProliferativeType cell_type = TRANSIT;

        // For a RADIAL Wnt type, override the cell type to STEM if the Wnt stimulus exceeds a higher threshold
        if ( (wnt_type == RADIAL) && (wnt_level > mWntStemThreshold) )
        {
            cell_type = STEM;
        }

        mCellProliferativeType = cell_type;
    }
    else
    {
        // The cell is DIFFERENTIATED and so in G0 phase
        mCellProliferativeType = DIFFERENTIATED;
    }
    AbstractSimpleCellCycleModel::UpdateCellCyclePhase();
}

void UniformDistributionWntGenerationBasedCellCycleModel::ResetForDivision()
{
    mGeneration++;
    if (mGeneration > mMaxTransitGenerations)
    {
        mCellProliferativeType = DIFFERENTIATED;
    }
    if (mCellProliferativeType == STEM)
    {
        mGeneration = 0;
    }
    AbstractSimpleCellCycleModel::ResetForDivision();
}

void UniformDistributionWntGenerationBasedCellCycleModel::InitialiseDaughterCell()
{
    WntConcentrationType wnt_type = GetWntType();

    SetGeneration(0);

    if (wnt_type == RADIAL)
    {
        mCellProliferativeType = TRANSIT;
    }

    /*
     * If the parent cell is a stem cell then its generation was reset
     * to zero when ResetForDivision() was called. The daughter cell's
     * generation must therefore be incremented here.
     */
    if (mGeneration == 0)
    {
        mGeneration = 1;
    }
    /*
     * In generation-based cell-cycle models, the daughter cell
     * is always of type transit or differentiated.
     */
    mCellProliferativeType = TRANSIT;
    if (mGeneration > mMaxTransitGenerations)
    {
        mCellProliferativeType = DIFFERENTIATED;
    }

    AbstractSimpleCellCycleModel::InitialiseDaughterCell();
}

bool UniformDistributionWntGenerationBasedCellCycleModel::CanCellTerminallyDifferentiate()
{
    return false;
}

double UniformDistributionWntGenerationBasedCellCycleModel::GetWntStemThreshold()
{
    return mWntStemThreshold;
}

void UniformDistributionWntGenerationBasedCellCycleModel::SetWntStemThreshold(double wntStemThreshold)
{
    assert(wntStemThreshold <= 1.0);
    assert(wntStemThreshold >= 0.0);
    mWntStemThreshold = wntStemThreshold;
}

double UniformDistributionWntGenerationBasedCellCycleModel::GetWntTransitThreshold()
{
    return mWntTransitThreshold;
}

void UniformDistributionWntGenerationBasedCellCycleModel::SetWntTransitThreshold(double wntTransitThreshold)
{
    assert(wntTransitThreshold <= 1.0);
    assert(wntTransitThreshold >= 0.0);
    mWntTransitThreshold = wntTransitThreshold;
}

double UniformDistributionWntGenerationBasedCellCycleModel::GetWntLabelledThreshold()
{
    return mWntLabelledThreshold;
}

void UniformDistributionWntGenerationBasedCellCycleModel::SetWntLabelledThreshold(double wntLabelledThreshold)
{
    assert(wntLabelledThreshold <= 1.0);
    assert(wntLabelledThreshold >= 0.0);
    mWntLabelledThreshold = wntLabelledThreshold;
}

void UniformDistributionWntGenerationBasedCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile <<  "\t\t\t<WntStemThreshold>"<< mWntStemThreshold << "</WntStemThreshold>\n";
    *rParamsFile <<  "\t\t\t<WntTransitThreshold>"<< mWntTransitThreshold << "</WntTransitThreshold>\n";
    *rParamsFile <<  "\t\t\t<WntLabelledThreshold>"<< mWntLabelledThreshold << "</WntLabelledThreshold>\n";

    // Call direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(UniformDistributionWntGenerationBasedCellCycleModel)
