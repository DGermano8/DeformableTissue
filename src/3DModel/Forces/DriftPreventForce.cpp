#include "DriftPreventForce.hpp"

template<unsigned DIM>
DriftPreventForce<DIM>::DriftPreventForce()
    : AbstractForce<DIM>(),
	  mMovementParameter(0.01),
      mTissueMiddle(0.0)
{
}

template<unsigned DIM>
DriftPreventForce<DIM>::~DriftPreventForce()
{
}

template<unsigned DIM>
void DriftPreventForce<DIM>::SetMovementParameter(double movementParameter)
{
    assert(movementParameter > 0.0);
    mMovementParameter = movementParameter;
}

template<unsigned DIM>
double DriftPreventForce<DIM>::GetMovementParameter()
{
    return mMovementParameter;
}

template<unsigned DIM>
void DriftPreventForce<DIM>::SetTissueMiddle(double tissue_middle)
{
	mTissueMiddle = tissue_middle;
}

template<unsigned DIM>
double DriftPreventForce<DIM>::GetTissueMiddle()
{
	return mTissueMiddle;
}

template<unsigned DIM>
void DriftPreventForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    double dt = SimulationTime::Instance()->GetTimeStep();

   	double tissue_max_height = -10.0;
	double tissue_min_height = 10.0;
	double average_height = 0.0;

    double tissue_middle = GetTissueMiddle();

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // First, create and store a copy of this real node and cell
        unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        c_vector<double, DIM> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

		if (real_node_location[2] > tissue_max_height)
		{
			tissue_max_height = real_node_location[2];
		}
		if (real_node_location[2] < tissue_min_height)
		{
			tissue_min_height = real_node_location[2];
		}
		average_height = average_height + real_node_location[2];

	}
   	average_height = average_height/(rCellPopulation.GetNumRealCells());


    // Iterate over the nodes
    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {

        c_vector<double, 3> force_due_to_drift;
		force_due_to_drift[0] = 0.0; force_due_to_drift[1] = 0.0;

        force_due_to_drift[2] = -1.0*(average_height - tissue_middle)/(SimulationTime::Instance()->GetTimeStep());

        node_iter->AddAppliedForceContribution(force_due_to_drift);

    }
}

template<unsigned DIM>
void DriftPreventForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MovementParameter>" << mMovementParameter << "</MovementParameter> \n";

    // Call direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class DriftPreventForce<1>;
template class DriftPreventForce<2>;
template class DriftPreventForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DriftPreventForce)
