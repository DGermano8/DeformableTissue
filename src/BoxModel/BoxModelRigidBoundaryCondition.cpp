#include "BoxModelRigidBoundaryCondition.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "StromalCellMutationState.hpp"
#include "RandomNumberGenerator.hpp"

BoxModelRigidBoundaryCondition::BoxModelRigidBoundaryCondition(AbstractCellPopulation<2>* pCellPopulation)
    : AbstractCellPopulationBoundaryCondition<2>(pCellPopulation),
    mMaxHeightForPinnedCells(DOUBLE_UNSET),
    mUseJiggledBottomCells(false),
    mCellPopulationWidth(DOUBLE_UNSET)
{
}

BoxModelRigidBoundaryCondition::~BoxModelRigidBoundaryCondition()
{
}

void BoxModelRigidBoundaryCondition::ImposeBoundaryCondition(const std::vector< c_vector<double, 2> >& rOldLocations)
{
    // Iterate over all nodes associated with real cells to update their positions
    // according to any tissue boundary conditions - don't want the top row of cells to be affected by bc
    for (AbstractCellPopulation<2>::Iterator cell_iter = this->mpCellPopulation->Begin();
    	 cell_iter != this->mpCellPopulation->End();
    	 ++cell_iter)
    {
		assert(this->mpCellPopulation->IsCellAssociatedWithADeletedLocation(*cell_iter) == false);

		// We only apply these boundary conditions to the stromal cells
        if (cell_iter->GetMutationState()->IsType<StromalCellMutationState>())
        {
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

			// Boundary conditions on the vertical edges of the stromal box

			// Any cell that has moved beyond the LH boundary of the box should be moved back
			if (p_node->rGetLocation()[0] < 0.0)
			{
				p_node->rGetModifiableLocation()[0] = 0.0;

				if (mUseJiggledBottomCells)
				{
					// Here we give the cell a push to the right
					p_node->rGetModifiableLocation()[0] = 0.05*RandomNumberGenerator::Instance()->ranf();
				}

			}

			// Any cell that has moved beyond the RH boundary of the box should be moved back
			if (p_node->rGetLocation()[0] > mCellPopulationWidth)
			{
				p_node->rGetModifiableLocation()[0] = mCellPopulationWidth;

				if (mUseJiggledBottomCells)
				{
					// Here we give the cell a push to the left
					p_node->rGetModifiableLocation()[0] = mCellPopulationWidth - 0.05*RandomNumberGenerator::Instance()->ranf();
				}
			}

			assert(p_node->rGetLocation()[0] >= 0.0);
			assert(p_node->rGetLocation()[0] <= mCellPopulationWidth);
        }
    }
}

bool BoxModelRigidBoundaryCondition::VerifyBoundaryCondition()
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
        if ( (cell_iter->GetMutationState()->IsType<StromalCellMutationState>())
        		&& ( (p_node->rGetLocation()[1] < 0.0) || (p_node->rGetLocation()[0] < 0.0)
        		|| (p_node->rGetLocation()[0] > mCellPopulationWidth) ) )
        {
        	boundary_condition_satisfied = false;
        	break;
        }
    }

    return boundary_condition_satisfied;
}

void BoxModelRigidBoundaryCondition::SetMaxHeightForPinnedCells(double maxHeightForPinnedCells)
{
	mMaxHeightForPinnedCells = maxHeightForPinnedCells;
}

double BoxModelRigidBoundaryCondition::GetMaxHeightForPinnedCells()
{
	return mMaxHeightForPinnedCells;
}

void BoxModelRigidBoundaryCondition::SetCellPopulationWidth(double cellPopulationWidth)
{
	mCellPopulationWidth = cellPopulationWidth;
}

double BoxModelRigidBoundaryCondition::GetCellPopulationWidth()
{
	return mCellPopulationWidth;
}

void BoxModelRigidBoundaryCondition::SetUseJiggledBottomCells(bool useJiggledBottomCells)
{
	mUseJiggledBottomCells = useJiggledBottomCells;
}

bool BoxModelRigidBoundaryCondition::GetUseJiggledBottomCells()
{
	return mUseJiggledBottomCells;
}

void BoxModelRigidBoundaryCondition::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<MaxHeightForPinnedCells>" << mMaxHeightForPinnedCells << "</MaxHeightForPinnedCells>\n";
    *rParamsFile << "\t\t<UseJiggledBottomCells>" << mUseJiggledBottomCells << "</UseJiggledBottomCells>\n";
    *rParamsFile << "\t\t<CellPopulationWidth>" << mCellPopulationWidth << "</CellPopulationWidth>\n";

    // Call method on direct parent class
	AbstractCellPopulationBoundaryCondition<2>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(BoxModelRigidBoundaryCondition)
