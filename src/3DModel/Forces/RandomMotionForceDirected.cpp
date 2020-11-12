#include "RandomMotionForceDirected.hpp"

template<unsigned DIM>
RandomMotionForceDirected<DIM>::RandomMotionForceDirected()
    : AbstractForce<DIM>(),
	  mMovementParameter(0.01)
{
}

template<unsigned DIM>
RandomMotionForceDirected<DIM>::~RandomMotionForceDirected()
{
}

template<unsigned DIM>
void RandomMotionForceDirected<DIM>::SetMovementParameter(double movementParameter)
{
    assert(movementParameter > 0.0);
    mMovementParameter = movementParameter;
}

template<unsigned DIM>
double RandomMotionForceDirected<DIM>::GetMovementParameter()
{
    return mMovementParameter;
}

template<unsigned DIM>
void RandomMotionForceDirected<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    double dt = SimulationTime::Instance()->GetTimeStep();

    if(SimulationTime::Instance()->GetTime() < 10.0)
    {
        // Iterate over the nodes
        for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
            node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
            ++node_iter)
        {
            c_vector<double, DIM> force_contribution = zero_vector<double>(DIM);;
            for (unsigned i=0; i<DIM-1; i++)
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

                force_contribution[i] = (sqrt(2.0*mMovementParameter*dt)/dt)*xi;
            }
            node_iter->AddAppliedForceContribution(force_contribution);
        }
    }
}

template<unsigned DIM>
void RandomMotionForceDirected<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MovementParameter>" << mMovementParameter << "</MovementParameter> \n";

    // Call direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class RandomMotionForceDirected<1>;
template class RandomMotionForceDirected<2>;
template class RandomMotionForceDirected<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RandomMotionForceDirected)
