#include "TissueSlabRigidBoundaryCondition.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "StromalCellMutationState.hpp"
#include "Debug.hpp"

TissueSlabRigidBoundaryCondition::TissueSlabRigidBoundaryCondition(AbstractCellPopulation<2>* pCellPopulation)
    : AbstractCellPopulationBoundaryCondition<2>(pCellPopulation),
    mMaxHeightForPinnedCells(DOUBLE_UNSET)
{
}

TissueSlabRigidBoundaryCondition::~TissueSlabRigidBoundaryCondition()
{
}

void TissueSlabRigidBoundaryCondition::ImposeBoundaryCondition(const std::vector< c_vector<double, 2> >& rOldLocations)
{
    // Iterate over all nodes associated with real cells to update their positions
    // according to any tissue boundary conditions - don't want the top row of cells to be affected by bc
    for (AbstractCellPopulation<2>::Iterator cell_iter = this->mpCellPopulation->Begin();
    	 cell_iter != this->mpCellPopulation->End();
    	 ++cell_iter)
    {
		assert(this->mpCellPopulation->IsCellAssociatedWithADeletedLocation(*cell_iter) == false);

		// Get index of node associated with cell
        unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);

        // Get pointer to this node
        Node<2>* p_node = this->mpCellPopulation->GetNode(node_index);

        // Pin the bottom cells
        if (p_node->rGetLocation()[1] < mMaxHeightForPinnedCells)
        {
        	// Get old node location
        	c_vector<double, 2> old_node_location = rOldLocations[node_index];

        	// Return node to old location
        	p_node->rGetModifiableLocation()[0] = old_node_location[0];
        	p_node->rGetModifiableLocation()[1] = old_node_location[1];
        }

        assert(p_node->rGetLocation()[1] >= 0.0);

        // Boundary conditions on the vertical edges:
        if ((cell_iter->GetCellCycleModel()->GetCellProliferativeType()==DIFFERENTIATED))
        {
            // Any cell that has moved beyond the LH boundary of the box should be moved back
            if (p_node->rGetLocation()[0] < 0.0)
            {
                p_node->rGetModifiableLocation()[0] = 0.0;
            }

            double width = mCellPopulationWidth;

            // Any cell that has moved beyond the RH boundary of the box should be moved back
            if (p_node->rGetLocation()[0] > width)
            {
                p_node->rGetModifiableLocation()[0] = width;
            }

            assert(p_node->rGetLocation()[0] >= 0.0);
            assert(p_node->rGetLocation()[0] <= width);
        }
    }
}

bool TissueSlabRigidBoundaryCondition::VerifyBoundaryCondition()
{
	bool boundary_condition_satisfied = true;

	/*
	 * Here we verify that the boundary condition is still satisfied by simply
	 * checking that no cells lies below the y=0 boundary.
	 */
    for (AbstractCellPopulation<2>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        // Get index of node associated with cell
        unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);

        // Get pointer to this node
        Node<2>* p_node = this->mpCellPopulation->GetNode(node_index);

        // If this (real) node lies below the y=0 boundary, or outside of the periodic boundaries, break and return false
        if ( (cell_iter->GetCellCycleModel()->GetCellProliferativeType()==DIFFERENTIATED) && ( (p_node->rGetLocation()[1] < 0.0) || (p_node->rGetLocation()[0] < 0.0)
        		|| (p_node->rGetLocation()[0] > mCellPopulationWidth) ) )
        {
        	boundary_condition_satisfied = false;
        	break;
        }
    }

    return boundary_condition_satisfied;
}

void TissueSlabRigidBoundaryCondition::SetMaxHeightForPinnedCells(double maxHeightForPinnedCells)
{
	mMaxHeightForPinnedCells = maxHeightForPinnedCells;
}

double TissueSlabRigidBoundaryCondition::GetMaxHeightForPinnedCells()
{
	return mMaxHeightForPinnedCells;
}

void TissueSlabRigidBoundaryCondition::SetCellPopulationWidth(double cellPopulationWidth)
{
	mCellPopulationWidth = cellPopulationWidth;
}

double TissueSlabRigidBoundaryCondition::GetCellPopulationWidth()
{
	return mCellPopulationWidth;
}

void TissueSlabRigidBoundaryCondition::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<MaxHeightForPinnedCells>" << mMaxHeightForPinnedCells << "</MaxHeightForPinnedCells>\n";
    *rParamsFile << "\t\t<CellPopulationWidth>" << mCellPopulationWidth << "</CellPopulationWidth>\n";

    // Call method on direct parent class
	AbstractCellPopulationBoundaryCondition<2>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(TissueSlabRigidBoundaryCondition)
