#include "FixedEpithelialBoundary3d.hpp"

    FixedEpithelialBoundary3d::FixedEpithelialBoundary3d(AbstractCellPopulation<3>* pCellPopulation, double maxHeightForPinnedCells, double cellPopulationWidth, double cellPopulationDepth)
        : AbstractCellPopulationBoundaryCondition<3>(pCellPopulation),
          mMaxHeightForPinnedCells(maxHeightForPinnedCells),
          mCellPopulationWidth(cellPopulationWidth),
          mCellPopulationDepth(cellPopulationDepth)
    {
//        mpStaticCastCellPopulation = static_cast<MeshBasedCellPopulation<3>*>(&mpCellPopulation);
    	assert((this->mpCellPopulation == NULL));
    }

    double FixedEpithelialBoundary3d::GetMaxHeightForPinnedCells()
    {
    	return mMaxHeightForPinnedCells;
    }

    double FixedEpithelialBoundary3d::GetCellPopulationWidth()
    {
    	return mCellPopulationWidth;
    }

    double FixedEpithelialBoundary3d::GetCellPopulationDepth()
    {
    	return mCellPopulationDepth;
    }

    void ImposeBoundaryCondition()
    {
        // Iterate over all nodes associated with real cells to update their positions
        // according to any tissue boundary conditions

        for (AbstractCellPopulation<3>::Iterator cell_iter = this->mpCellPopulation->Begin();
             cell_iter != this->mpCellPopulation->End();
             ++cell_iter)
        {
			assert(this->mpCellPopulation->IsCellAssociatedWithADeletedLocation(*cell_iter) == false);

			// Get index of node associated with cell
			unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
//			assert(!(this->mpCellPopulation->IsGhostNode(node_index)));

			// Get pointer to this node
            Node<3>* p_node = this->mpCellPopulation->GetNode(node_index);

//	        if ( (mpStaticCastCellPopulation->GetCellUsingLocationIndex(node_index)->GetMutationState()->IsType<StromalCellMutationState>())
//	        		&& ( (p_node->rGetLocation()[2] < mMaxHeightForPinnedCells)
//	        				|| (p_node->rGetLocation()[0] < 0.5)
//	        				|| (p_node->rGetLocation()[0] > mCellPopulationWidth-0.5)
//	        				|| (p_node->rGetLocation()[1] < 0.5)
//	        				|| (p_node->rGetLocation()[1] > mCellPopulationDepth-0.5) ) )

			// Pin the cells which are at the edges or on the bottom layer (both stromal and wildtype)
			if ( (p_node->rGetLocation()[2] < mMaxHeightForPinnedCells)
							|| (p_node->rGetLocation()[0] < 0.5)
							|| (p_node->rGetLocation()[0] > mCellPopulationWidth-0.5)
							|| (p_node->rGetLocation()[1] < 0.5)
							|| (p_node->rGetLocation()[1] > mCellPopulationDepth-0.5)  )
			{
				// Get old node location
				c_vector<double, 3> old_node_location = rOldLocations[node_index];

				// Return node to old location
				p_node->rGetModifiableLocation()[0] = old_node_location[0];
				p_node->rGetModifiableLocation()[1] = old_node_location[1];
				p_node->rGetModifiableLocation()[2] = old_node_location[2];
			}
		}
    }

    bool VerifyBoundaryCondition()
    {
        bool condition_satisfied = true;

        for (AbstractCellPopulation<3>::Iterator cell_iter = this->mpCellPopulation->Begin();
             cell_iter != this->mpCellPopulation->End();
             ++cell_iter)
        {
            c_vector<double, 3> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
            double x_coordinate = cell_location(0);
            double y_coordinate = cell_location(1);
            double z_coordinate = cell_location(2);

            if ( (x_coordinate < 0.5) || (x_coordinate > mCellPopulationWidth-0.5)
            		|| (y_coordinate < 0.5) || (y_coordinate > mCellPopulationDepth-0.5) ) // Note, don't think we can check the bottom layer are pinned
            {
                condition_satisfied = false;
                break;
            }
        }
        return condition_satisfied;
    }

    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t<MaxHeightForPinnedCells>" << mMaxHeightForPinnedCells << "</MaxHeightForPinnedCells>\n";
    	*rParamsFile << "\t\t<CellPopulationWidth>" << mCellPopulationWidth << "</CellPopulationWidth>\n";
        *rParamsFile << "\t\t<CellPopulationDepth>" << mCellPopulationDepth << "</CellPopulationDepth>\n";

        AbstractCellPopulationBoundaryCondition<3>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
    }
};

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(FixedEpithelialBoundary3d)
