#include "PeriodicStromalBoxBoundaryCondition3dIncreasingDomain.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "StromalCellMutationState.hpp"

PeriodicStromalBoxBoundaryCondition3dIncreasingDomain::PeriodicStromalBoxBoundaryCondition3dIncreasingDomain(AbstractCellPopulation<3>* pCellPopulation)
    : AbstractCellPopulationBoundaryCondition<3>(pCellPopulation),
    mCellPopulationWidth(DOUBLE_UNSET),
    mCellPopulationDepth(DOUBLE_UNSET),
    mMaxHeightForPinnedCells(DOUBLE_UNSET),
	mIncreaseDomainTime(DBL_MAX),
	mEndIncreaseDomainTime(0.0),
	mMultiplyerIncreaseDomainTime(1.0)
{
}

PeriodicStromalBoxBoundaryCondition3dIncreasingDomain::~PeriodicStromalBoxBoundaryCondition3dIncreasingDomain()
{
}
 //void PeriodicStromalBoxBoundaryCondition3dIncreasingDomain::ImposeBoundaryCondition(const std::vector< c_vector<double, 3> >& rOldLocations)
void PeriodicStromalBoxBoundaryCondition3dIncreasingDomain::ImposeBoundaryCondition(const std::map<Node<3>*, c_vector<double, 3> >& rOldLocations)
{
	double sim_time = SimulationTime::Instance()->GetTime();
	double gotCellPopulationWidth = GetCellPopulationWidth();
	double gotCellPopulationDepth = GetCellPopulationDepth();
	if( sim_time < mIncreaseDomainTime )
	{
		gotCellPopulationWidth = GetCellPopulationWidth();
		gotCellPopulationDepth = GetCellPopulationDepth();
	}
	else if( sim_time >= mIncreaseDomainTime && sim_time < mEndIncreaseDomainTime)
	{
		gotCellPopulationWidth = GetCellPopulationWidth() +  mMultiplyerIncreaseDomainTime*(sim_time - mIncreaseDomainTime)/mEndIncreaseDomainTime ;
		gotCellPopulationDepth = GetCellPopulationDepth() +  sqrt(0.75)*mMultiplyerIncreaseDomainTime*(sim_time - mIncreaseDomainTime)/mEndIncreaseDomainTime ;
	}
	else if(sim_time >= mEndIncreaseDomainTime && mEndIncreaseDomainTime > 0.0)
	{
		gotCellPopulationWidth = GetCellPopulationWidth() +  mMultiplyerIncreaseDomainTime*(mEndIncreaseDomainTime - mIncreaseDomainTime)/mEndIncreaseDomainTime;
		gotCellPopulationDepth = GetCellPopulationDepth() +  sqrt(0.75)*mMultiplyerIncreaseDomainTime*(mEndIncreaseDomainTime - mIncreaseDomainTime)/mEndIncreaseDomainTime ;
	}

	for (AbstractCellPopulation<3>::Iterator cell_iter = this->mpCellPopulation->Begin();
		 cell_iter != this->mpCellPopulation->End();
		 ++cell_iter)
	{
		assert(this->mpCellPopulation->IsCellAssociatedWithADeletedLocation(*cell_iter) == false);

		// Get index of node associated with cell
        unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);

        // Get pointer to this node
        Node<3>* p_node = this->mpCellPopulation->GetNode(node_index);
		
		if (this->mpCellPopulation->GetCellUsingLocationIndex(node_index)->GetMutationState()->IsType<StromalCellMutationState>())
		{
			// Pin the bottom cells 
			if (p_node->rGetLocation()[2] < mMaxHeightForPinnedCells) 
			{
				// Get old node location
				//c_vector<double, 3> old_node_location = rOldLocations[node_index];
				typename std::map<Node<3>*, c_vector<double, 3> >::const_iterator it = rOldLocations.find(p_node);
    			c_vector<double, 3> old_node_location = it->second;
	
				// Return node to old location
				p_node->rGetModifiableLocation()[0] = old_node_location[0];
				p_node->rGetModifiableLocation()[1] = old_node_location[1];
				// p_node->rGetModifiableLocation()[2] = 0.0;
				p_node->rGetModifiableLocation()[2] = mMaxHeightForPinnedCells;
			}
	
			// Want to enforce a periodic boundary box on only the stromal cells, such that cells which move outside the 
			// boundaries of the box get moved to the opposite side (as in the cylindrical mesh)
			if (p_node->rGetLocation()[0] < 0.0)
			{
				/* Move to the opposite edge */
				 p_node->rGetModifiableLocation()[0] = p_node->rGetLocation()[0] + gotCellPopulationWidth;
			}
			else if (p_node->rGetLocation()[0] > gotCellPopulationWidth)
			{
				/* Move to the opposite edge */
				p_node->rGetModifiableLocation()[0] = p_node->rGetLocation()[0] - gotCellPopulationWidth;
			}
			if (p_node->rGetLocation()[1] < 0.0)
			{
				/* Move to the opposite edge */
				 p_node->rGetModifiableLocation()[1] = p_node->rGetLocation()[1] + gotCellPopulationDepth;
			}
			else if (p_node->rGetLocation()[1] > gotCellPopulationDepth)
			{
				/* Move to the opposite edge */
				p_node->rGetModifiableLocation()[1] = p_node->rGetLocation()[1] - gotCellPopulationDepth;
			}
		}
	}
}

bool PeriodicStromalBoxBoundaryCondition3dIncreasingDomain::VerifyBoundaryCondition()
{
	double sim_time = SimulationTime::Instance()->GetTime();
	double gotCellPopulationWidth = GetCellPopulationWidth();
	double gotCellPopulationDepth = GetCellPopulationDepth();

	if( sim_time < mIncreaseDomainTime )
	{
		gotCellPopulationWidth = GetCellPopulationWidth();
		gotCellPopulationDepth = GetCellPopulationDepth();
	}
	else if( sim_time >= mIncreaseDomainTime && sim_time < mEndIncreaseDomainTime)
	{
		gotCellPopulationWidth = GetCellPopulationWidth() + mMultiplyerIncreaseDomainTime*(sim_time - mIncreaseDomainTime)/mEndIncreaseDomainTime ;
		gotCellPopulationDepth = GetCellPopulationDepth() + sqrt(0.75)*mMultiplyerIncreaseDomainTime*(sim_time - mIncreaseDomainTime)/mEndIncreaseDomainTime ;
	}
	else if(sim_time >= mEndIncreaseDomainTime  && mEndIncreaseDomainTime > 0.0)
	{
		gotCellPopulationWidth = GetCellPopulationWidth() + mMultiplyerIncreaseDomainTime*(mEndIncreaseDomainTime - mIncreaseDomainTime)/mEndIncreaseDomainTime;
		gotCellPopulationDepth = GetCellPopulationDepth() + sqrt(0.75)*mMultiplyerIncreaseDomainTime*(mEndIncreaseDomainTime - mIncreaseDomainTime)/mEndIncreaseDomainTime ;
	}

	bool boundary_condition_satisfied = true;

	/*
	 * Here we verify that the boundary condition is still satisfied by simply
	 * checking that no cells lies below the Z=0 boundary.
	 */
    for (AbstractCellPopulation<3>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        // Get index of node associated with cell
        unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);

        // Get pointer to this node
        Node<3>* p_node = this->mpCellPopulation->GetNode(node_index);

        // If this node lies below the z=0 boundary(NO! - Dom), or outside of the periodic boundaries, break and return false
		// if ( (cell_iter->GetMutationState()->IsType<StromalCellMutationState>()) && 
        // 		((p_node->rGetLocation()[2] < 0.0) || (p_node->rGetLocation()[0] < 0.0) || (p_node->rGetLocation()[0] > mCellPopulationWidth)
        // 		|| (p_node->rGetLocation()[1] < 0.0) || (p_node->rGetLocation()[1] > mCellPopulationDepth)) )

        // If this node is outside of the periodic boundaries, break and return false
        if ( (cell_iter->GetMutationState()->IsType<StromalCellMutationState>()) && 
        		  ((p_node->rGetLocation()[0] < 0.0) || (p_node->rGetLocation()[0] > gotCellPopulationWidth)
        		|| (p_node->rGetLocation()[1] < 0.0) || (p_node->rGetLocation()[1] > gotCellPopulationDepth)) )
        {
        	boundary_condition_satisfied = false;
        	break;
        }
    }

    return boundary_condition_satisfied;
}

void PeriodicStromalBoxBoundaryCondition3dIncreasingDomain::SetMaxHeightForPinnedCells(double maxHeightForPinnedCells)
{
	mMaxHeightForPinnedCells = maxHeightForPinnedCells;
}

double PeriodicStromalBoxBoundaryCondition3dIncreasingDomain::GetMaxHeightForPinnedCells()
{
	return mMaxHeightForPinnedCells;
}

void PeriodicStromalBoxBoundaryCondition3dIncreasingDomain::SetCellPopulationWidth(double cellPopulationWidth)
{
	mCellPopulationWidth = cellPopulationWidth;
}

double PeriodicStromalBoxBoundaryCondition3dIncreasingDomain::GetCellPopulationWidth()
{
	return mCellPopulationWidth;
}

void PeriodicStromalBoxBoundaryCondition3dIncreasingDomain::SetCellPopulationDepth(double cellPopulationDepth)
{
	mCellPopulationDepth = cellPopulationDepth;
}

double PeriodicStromalBoxBoundaryCondition3dIncreasingDomain::GetCellPopulationDepth()
{
	return mCellPopulationDepth;
}

void PeriodicStromalBoxBoundaryCondition3dIncreasingDomain::SetIncreaseDomainTime(double increaseDomainTime)
{
	mIncreaseDomainTime = increaseDomainTime;
}

double PeriodicStromalBoxBoundaryCondition3dIncreasingDomain::GetIncreaseDomainTime()
{
	return mIncreaseDomainTime;
}

void PeriodicStromalBoxBoundaryCondition3dIncreasingDomain::SetEndIncreaseDomainTime(double endIncreaseDomainTime)
{
	mEndIncreaseDomainTime = endIncreaseDomainTime;
}

double PeriodicStromalBoxBoundaryCondition3dIncreasingDomain::GetEndIncreaseDomainTime()
{
	return mEndIncreaseDomainTime;
}

void PeriodicStromalBoxBoundaryCondition3dIncreasingDomain::SetMultiplyerIncreaseDomainTime(double multiplyerIncreaseDomainTime)
{
	mMultiplyerIncreaseDomainTime = multiplyerIncreaseDomainTime;
}
double PeriodicStromalBoxBoundaryCondition3dIncreasingDomain::GetMultiplyerIncreaseDomainTime()
{
	return mMultiplyerIncreaseDomainTime;
}

void PeriodicStromalBoxBoundaryCondition3dIncreasingDomain::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<CellPopulationDepth>" << mCellPopulationDepth << "</CellPopulationDepth>\n";
    *rParamsFile << "\t\t<CellPopulationWidth>" << mCellPopulationWidth << "</CellPopulationWidth>\n";
    *rParamsFile << "\t\t<MaxHeightForPinnedCells>" << mMaxHeightForPinnedCells << "</MaxHeightForPinnedCells>\n";
	*rParamsFile << "\t\t<IncreaseDomainTime>" << mIncreaseDomainTime << "</IncreaseDomainTime>\n";
    *rParamsFile << "\t\t<EndIncreaseDomainTime>" << mEndIncreaseDomainTime << "</EndIncreaseDomainTime>\n";
    *rParamsFile << "\t\t<MultiplyerIncreaseDomainTime>" << mMultiplyerIncreaseDomainTime << "</MultiplyerIncreaseDomainTime>\n";

    
    // Call method on direct parent class
	AbstractCellPopulationBoundaryCondition<3>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PeriodicStromalBoxBoundaryCondition3dIncreasingDomain)
