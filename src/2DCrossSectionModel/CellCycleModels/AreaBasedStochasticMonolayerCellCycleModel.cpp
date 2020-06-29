#include "PetscTools.hpp"
#include "AreaBasedStochasticMonolayerCellCycleModel.hpp"
#include "Debug.hpp"

AreaBasedStochasticMonolayerCellCycleModel::AreaBasedStochasticMonolayerCellCycleModel()
    : StochasticMonolayerCellCycleModel()
{
}

void AreaBasedStochasticMonolayerCellCycleModel::UpdateCellCyclePhase()
{
    // Get cell area
    double cell_area;
    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            cell_area = CellwiseData<DIM>::Instance()->GetValue(mpCell, 0);
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            cell_area = CellwiseData<DIM>::Instance()->GetValue(mpCell, 0);
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            cell_area = CellwiseData<DIM>::Instance()->GetValue(mpCell, 0);
            break;
        }
        default:
            NEVER_REACHED;
    }

//    PRINT_VARIABLE(SimulationTime::Instance()->GetTime());
    PRINT_VARIABLE(cell_area);
//    PRINT_VARIABLE(mpCell->GetCellCycleModel()->GetCellProliferativeType());

    //G_ZERO_PHASE = 0, G_ONE_PHASE = 1, S_PHASE = 2, G_TWO_PHASE = 3, M_PHASE = 4;

//    mCurrentCellCyclePhase  = GetCurrentCellCyclePhase();
//    PRINT_VARIABLE(mCurrentCellCyclePhase);

    if (mCurrentCellCyclePhase == G_ONE_PHASE)
    {
        // Update G1 duration based on cell area
        double dt = SimulationTime::Instance()->GetTimeStep();

        double equilibrium_area = sqrt(3.0)*0.5;
//        double fraction_of_equilibrium_area_at_which_cell_enters_quiescence = 0.75;
//        double quiescent_area = equilibrium_area * fraction_of_equilibrium_area_at_which_cell_enters_quiescence;
        double quiescent_area = equilibrium_area * mQuiescentAreaFraction;

//        PRINT_VARIABLE(quiescent_area);
//        PRINT_VARIABLE((cell_area < quiescent_area));

        if (cell_area < quiescent_area)
        {
            mG1Duration += dt;
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
}

AbstractCellCycleModel* AreaBasedStochasticMonolayerCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell cycle model
    AreaBasedStochasticMonolayerCellCycleModel* p_model = new AreaBasedStochasticMonolayerCellCycleModel();

    // Set the values of the new cell cycle model's member variables
    p_model->SetGeneration(mGeneration);
    p_model->SetMaxTransitGenerations(mMaxTransitGenerations);
    p_model->SetCellProliferativeType(mCellProliferativeType);
    p_model->SetQuiescentAreaFraction(mQuiescentAreaFraction);
    p_model->SetDimension(2);

    return p_model;
}

void AreaBasedStochasticMonolayerCellCycleModel::SetQuiescentAreaFraction(double quiescentAreaFraction)
{
    mQuiescentAreaFraction = quiescentAreaFraction;
}

//bool AreaBasedStochasticMonolayerCellCycleModel::ReadyToDivide()
//{
//    assert(mpCell != NULL);
//
//    // Get cell area
//    double cell_area;
//    switch (mDimension)
//    {
//        case 1:
//        {
//            const unsigned DIM = 1;
//            cell_area = CellwiseData<DIM>::Instance()->GetValue(mpCell, 0);
//            break;
//        }
//        case 2:
//        {
//            const unsigned DIM = 2;
//            cell_area = CellwiseData<DIM>::Instance()->GetValue(mpCell, 0);
//            break;
//        }
//        case 3:
//        {
//            const unsigned DIM = 3;
//            cell_area = CellwiseData<DIM>::Instance()->GetValue(mpCell, 0);
//            break;
//        }
//        default:
//            NEVER_REACHED;
//    }
//
//    ///\todo Fix this usage of cell mutation state (see #1145)
//    boost::shared_ptr<AbstractCellMutationState> p_wild_type(new WildTypeCellMutationState);
//    mpCell->SetMutationState(p_wild_type);
//
//    if (mCurrentCellCyclePhase == G_ONE_PHASE)
//    {
//        // Update G1 duration based on cell area
//        double dt = SimulationTime::Instance()->GetTimeStep();
//
//        double equilibrium_area = sqrt(3.0)*0.5;
////        double fraction_of_equilibrium_area_at_which_cell_enters_quiescence = 0.75;
////        double quiescent_area = equilibrium_area * fraction_of_equilibrium_area_at_which_cell_enters_quiescence;
//        double quiescent_area = equilibrium_area * mQuiescentAreaFraction;
//
//    if (!mReadyToDivide)
//    {
//    	UpdateCellCyclePhase();
//
//        if ( (mCurrentCellCyclePhase != G_ZERO_PHASE) &&
//             (GetAge() >= GetMDuration() + GetG1Duration() + GetSDuration() + GetG2Duration()) &&
//             (cell_area > quiescent_area))
//        {
//            mReadyToDivide = true;
//        }
//    }
//    return mReadyToDivide;
//}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(AreaBasedStochasticMonolayerCellCycleModel)
