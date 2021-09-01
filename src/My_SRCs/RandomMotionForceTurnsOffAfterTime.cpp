#include "RandomMotionForceTurnsOffAfterTime.hpp"

template<unsigned DIM>
RandomMotionForceTurnsOffAfterTime<DIM>::RandomMotionForceTurnsOffAfterTime()
    : AbstractForce<DIM>(),
	  mMovementParameter(0.01),
      mTurnOffAfterTime(DOUBLE_UNSET)
{
}

template<unsigned DIM>
RandomMotionForceTurnsOffAfterTime<DIM>::~RandomMotionForceTurnsOffAfterTime()
{
}

template<unsigned DIM>
void RandomMotionForceTurnsOffAfterTime<DIM>::SetMovementParameter(double movementParameter)
{
    assert(movementParameter > 0.0);
    mMovementParameter = movementParameter;
}

template<unsigned DIM>
double RandomMotionForceTurnsOffAfterTime<DIM>::GetMovementParameter()
{
    return mMovementParameter;
}

template<unsigned DIM>
void RandomMotionForceTurnsOffAfterTime<DIM>::SetTurnOffAfterTime(double turnOffAfterTime)
{
    assert(turnOffAfterTime > 0.0);
    mTurnOffAfterTime = turnOffAfterTime;
}

template<unsigned DIM>
double RandomMotionForceTurnsOffAfterTime<DIM>::GetTurnOffAfterTime()
{
    return mTurnOffAfterTime;
}

template<unsigned DIM>
void RandomMotionForceTurnsOffAfterTime<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    double dt = SimulationTime::Instance()->GetTimeStep();

    if( SimulationTime::Instance()->GetTime() <= mTurnOffAfterTime)
    {
        double moving_mag = 0.0;
        if( mTurnOffAfterTime - SimulationTime::Instance()->GetTime() <= 1)
        {
            moving_mag = mMovementParameter*(mTurnOffAfterTime - SimulationTime::Instance()->GetTime());
        }
        else
        {
            moving_mag = mMovementParameter;
        }
        // Iterate over the nodes
        for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
            node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
            ++node_iter)
        {
                    c_vector<double, DIM> force_contribution;
            for (unsigned i=0; i<DIM; i++)
            {
                /*
                * The force on this cell is scaled with the timestep such that when it is
                * used in the discretised equation of motion for the cell, we obtain the
                * correct formula
                *
                * x_new = x_old + sqrt(2*D*dt)*W
                *
                * where W is a standard normal random variable.
                */
                double xi = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();

                force_contribution[i] = (sqrt(2.0*moving_mag*dt)/dt)*xi;
            }
            node_iter->AddAppliedForceContribution(force_contribution);
        }
    }
}

template<unsigned DIM>
void RandomMotionForceTurnsOffAfterTime<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MovementParameter>" << mMovementParameter << "</MovementParameter> \n";
    *rParamsFile << "\t\t\t<TurnOffAfterTime>" << mTurnOffAfterTime << "</TurnOffAfterTime> \n";

    // Call direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class RandomMotionForceTurnsOffAfterTime<1>;
template class RandomMotionForceTurnsOffAfterTime<2>;
template class RandomMotionForceTurnsOffAfterTime<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RandomMotionForceTurnsOffAfterTime)
