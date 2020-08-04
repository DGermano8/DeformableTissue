#include "SimpleWntWithContactInhibitionCellCycleModel.hpp"
#include "CellwiseData.hpp"
#include "PetscTools.hpp"

SimpleWntWithContactInhibitionCellCycleModel::SimpleWntWithContactInhibitionCellCycleModel()
    : mUseCellProliferativeTypeDependentG1Duration(false),
      mWntStemThreshold(0.8),
      mWntTransitThreshold(0.65),//(0.65),
      mWntLabelledThreshold(0.65)
{
}


AbstractCellCycleModel* SimpleWntWithContactInhibitionCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell cycle model
	SimpleWntWithContactInhibitionCellCycleModel* p_model = new SimpleWntWithContactInhibitionCellCycleModel();

    // Set the values of the new cell cycle model's member variables
    p_model->SetDimension(mDimension);
    p_model->SetCellProliferativeType(mCellProliferativeType);
    p_model->SetDimension(2);

    return p_model;
}


void SimpleWntWithContactInhibitionCellCycleModel::SetUseCellProliferativeTypeDependentG1Duration(bool useCellProliferativeTypeDependentG1Duration)
{
    mUseCellProliferativeTypeDependentG1Duration = useCellProliferativeTypeDependentG1Duration;
}


void SimpleWntWithContactInhibitionCellCycleModel::SetG1Duration()
{
    assert(mpCell != NULL);

    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    switch (mCellProliferativeType)
    {
        case STEM:
            if (mUseCellProliferativeTypeDependentG1Duration)
            {
                mG1Duration = 1 + 4*p_gen->ranf();
            }
            else
            {
                // Normally stem cells should behave just like transit cells in a Wnt simulation
                mG1Duration = 1 + 2*p_gen->ranf();
            }
            break;
        case TRANSIT:
            mG1Duration = 1 + 2*p_gen->ranf();
            break;
        case DIFFERENTIATED:
            mG1Duration = DBL_MAX;
            break;
        default:
            NEVER_REACHED;
    }

    // Check that the uniform random deviate has not returned a small G1 duration
    if (mG1Duration < mMinimumGapDuration)
    {
        mG1Duration = mMinimumGapDuration;
    }
}


double SimpleWntWithContactInhibitionCellCycleModel::GetWntLevel()
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


WntConcentrationType SimpleWntWithContactInhibitionCellCycleModel::GetWntType()
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


void SimpleWntWithContactInhibitionCellCycleModel::UpdateCellCyclePhase()
{
    // The cell can divide if the Wnt concentration >= wnt_division_threshold and also if the cell
	// area is sufficiently large

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

    // Set the cell type to TRANSIT if the Wnt stimulus exceeds wnt_division_threshold
    if (wnt_level >= wnt_division_threshold)
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

    // Get cell area
    double cell_area = CellwiseData<2>::Instance()->GetValue(mpCell, 0);;
//    switch (mDimension)
//    {
//        case 1:
//        {
//            const unsigned DIM = 1;
//            cell_area = CellwiseData<2>::Instance()->GetValue(mpCell, 0);
//            break;
//        }
//        case 2:
//        {
//            const unsigned DIM = 2;
//            cell_area = CellwiseData<2>::Instance()->GetValue(mpCell, 0);
//            break;
//        }
//        case 3:
//        {
//            const unsigned DIM = 3;
//            cell_area = CellwiseData<2>::Instance()->GetValue(mpCell, 0);
//            break;
//        }
//        default:
//            NEVER_REACHED;
//    }

    ///\todo Fix this usage of cell mutation state (see #1145)
    boost::shared_ptr<AbstractCellMutationState> p_wild_type(new WildTypeCellMutationState);
    mpCell->SetMutationState(p_wild_type);

    // Use the cell area to determine whether the cell should move out of G1
    if (mCurrentCellCyclePhase == G_ONE_PHASE)
    {
        // Update G1 duration based on cell area
        double dt = SimulationTime::Instance()->GetTimeStep();

        double equilibrium_area = sqrt(3.0/4.0);
//        double fraction_of_equilibrium_area_at_which_cell_enters_quiescence = 0.75;
//        double quiescent_area = equilibrium_area * fraction_of_equilibrium_area_at_which_cell_enters_quiescence;
        double quiescent_area = equilibrium_area * mQuiescentAreaFraction;

        if (cell_area < quiescent_area)
        {
            mG1Duration += dt;

            ///\todo Fix this usage of cell mutation state (see #1145)
//            mpCell->AddCellProperty(CellPropertyRegistry::Instance()->Get<CellLabel>());
        }
    }

    double time_since_birth = GetAge();
    assert(time_since_birth >= 0);

    if (mpCell->GetCellCycleModel()->GetCellProliferativeType()==DIFFERENTIATED)
    {
        mCurrentCellCyclePhase = G_ZERO_PHASE;
    }
    else if ( time_since_birth < GetMDuration() )
    {
        mCurrentCellCyclePhase = M_PHASE;
    }
    else if ( time_since_birth < GetMDuration() + mG1Duration)
    {
        mCurrentCellCyclePhase = G_ONE_PHASE;
    }
    else if ( time_since_birth < GetMDuration() + mG1Duration + GetSDuration())
    {
        mCurrentCellCyclePhase = S_PHASE;
    }
    else if ( time_since_birth < GetMDuration() + mG1Duration + GetSDuration() + GetG2Duration())
    {
        mCurrentCellCyclePhase = G_TWO_PHASE;
    }
    AbstractSimpleCellCycleModel::UpdateCellCyclePhase();
}

void SimpleWntWithContactInhibitionCellCycleModel::InitialiseDaughterCell()
{
    WntConcentrationType wnt_type = GetWntType();

    if (wnt_type == RADIAL)
    {
        mCellProliferativeType = TRANSIT;
    }

    AbstractSimpleCellCycleModel::InitialiseDaughterCell();
}

bool SimpleWntWithContactInhibitionCellCycleModel::CanCellTerminallyDifferentiate()
{
    return false;
}

double SimpleWntWithContactInhibitionCellCycleModel::GetWntStemThreshold()
{
    return mWntStemThreshold;
}

void SimpleWntWithContactInhibitionCellCycleModel::SetWntStemThreshold(double wntStemThreshold)
{
    assert(wntStemThreshold <= 1.0);
    assert(wntStemThreshold >= 0.0);
    mWntStemThreshold = wntStemThreshold;
}

double SimpleWntWithContactInhibitionCellCycleModel::GetWntTransitThreshold()
{
    return mWntTransitThreshold;
}

void SimpleWntWithContactInhibitionCellCycleModel::SetWntTransitThreshold(double wntTransitThreshold)
{
    assert(wntTransitThreshold <= 1.0);
    assert(wntTransitThreshold >= 0.0);
    mWntTransitThreshold = wntTransitThreshold;
}

double SimpleWntWithContactInhibitionCellCycleModel::GetWntLabelledThreshold()
{
    return mWntLabelledThreshold;
}

void SimpleWntWithContactInhibitionCellCycleModel::SetWntLabelledThreshold(double wntLabelledThreshold)
{
    assert(wntLabelledThreshold <= 1.0);
    assert(wntLabelledThreshold >= 0.0);
    mWntLabelledThreshold = wntLabelledThreshold;
}

void SimpleWntWithContactInhibitionCellCycleModel::SetQuiescentAreaFraction(double quiescentAreaFraction)
{
    mQuiescentAreaFraction = quiescentAreaFraction;
}

void SimpleWntWithContactInhibitionCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile <<  "\t\t\t<WntStemThreshold>"<< mWntStemThreshold << "</WntStemThreshold>\n";
    *rParamsFile <<  "\t\t\t<WntTransitThreshold>"<< mWntTransitThreshold << "</WntTransitThreshold>\n";
    *rParamsFile <<  "\t\t\t<WntLabelledThreshold>"<< mWntLabelledThreshold << "</WntLabelledThreshold>\n";
    *rParamsFile <<  "\t\t\t<QuiescentAreaFraction>"<< mQuiescentAreaFraction << "</QuiescentAreaFraction>\n";

    // Call direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(SimpleWntWithContactInhibitionCellCycleModel)
