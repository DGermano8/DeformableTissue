/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "PetscTools.hpp"

StochasticDurationGenerationBasedCellCycleModel::StochasticDurationGenerationBasedCellCycleModel()
{
}

StochasticDurationGenerationBasedCellCycleModel::StochasticDurationGenerationBasedCellCycleModel(const StochasticDurationGenerationBasedCellCycleModel& rModel)
   : AbstractSimpleGenerationalCellCycleModel(rModel)
{
    /*
     * The member variables mGeneration and mMaxTransitGeneration are
     * initialized in the AbstractSimpleGenerationalCellCycleModel
     * constructor.
     *
     * The member variables mCurrentCellCyclePhase, mG1Duration,
     * mMinimumGapDuration, mStemCellG1Duration, mTransitCellG1Duration,
     * mSDuration, mG2Duration and mMDuration are initialized in the
     * AbstractPhaseBasedCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mG1Duration is (re)set as soon as InitialiseDaughterCell()
     * is called on the new cell-cycle model.
     */
}

AbstractCellCycleModel* StochasticDurationGenerationBasedCellCycleModel::CreateCellCycleModel()
{
    return new StochasticDurationGenerationBasedCellCycleModel(*this);
}

//AbstractCellCycleModel* StochasticDurationGenerationBasedCellCycleModel::CreateCellCycleModel()
//{
    // Create a new cell-cycle model
//    StochasticDurationGenerationBasedCellCycleModel* p_model = new StochasticDurationGenerationBasedCellCycleModel();

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
     *
     * Note 3: the member variable mDimension remains unset, since this cell-cycle
     * model does not need to know the spatial dimension, so if we were to call
     * SetDimension() on the new cell-cycle model an exception would be triggered;
     * hence we do not set this member variable.
     */
/*  
    p_model->SetBirthTime(mBirthTime);
    p_model->SetCellProliferativeType(mCellProliferativeType);
    p_model->SetMinimumGapDuration(mMinimumGapDuration);
    p_model->SetStemCellG1Duration(mStemCellG1Duration);
    p_model->SetTransitCellG1Duration(mTransitCellG1Duration);
    p_model->SetSDuration(mSDuration);
    p_model->SetG2Duration(mG2Duration);
    p_model->SetMDuration(mMDuration);
    p_model->SetGeneration(mGeneration);
    p_model->SetMaxTransitGenerations(mMaxTransitGenerations);

    return p_model;
*/
//}

void StochasticDurationGenerationBasedCellCycleModel::SetG1Duration()
{
    assert(mpCell != NULL);

    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    switch (mCellProliferativeType)
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

void StochasticDurationGenerationBasedCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output

    // Call method on direct parent class
    AbstractSimpleGenerationBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(StochasticDurationGenerationBasedCellCycleModel)
