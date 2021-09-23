#include "PeriodicBoxBoundaryCondition3dGhostsIncreasingDomain.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "StromalCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"

#include "Debug.hpp"

PeriodicBoxBoundaryCondition3dGhostsIncreasingDomain::PeriodicBoxBoundaryCondition3dGhostsIncreasingDomain(AbstractCellPopulation<3>* pCellPopulation)
    : AbstractCellPopulationBoundaryCondition<3>(pCellPopulation),
    mCellPopulationWidth(DOUBLE_UNSET),
    mCellPopulationDepth(DOUBLE_UNSET),
    mMaxHeightForPinnedCells(DOUBLE_UNSET),
	mIncreaseDomainTime(DBL_MAX),
	mEndIncreaseDomainTime(0.0),
	mMultiplyerIncreaseDomainTime(1.0)
{
}

PeriodicBoxBoundaryCondition3dGhostsIncreasingDomain::~PeriodicBoxBoundaryCondition3dGhostsIncreasingDomain()
{
} 


//void PeriodicBoxBoundaryCondition3d::ImposeBoundaryCondition(const std::vector< c_vector<double, 3> >& rOldLocations)
void PeriodicBoxBoundaryCondition3dGhostsIncreasingDomain::ImposeBoundaryCondition(const std::map<Node<3>*, c_vector<double, 3> >& rOldLocations)
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
	else if(sim_time >= mEndIncreaseDomainTime && mEndIncreaseDomainTime > 0.0)
	{
		gotCellPopulationWidth = GetCellPopulationWidth() + mMultiplyerIncreaseDomainTime*(mEndIncreaseDomainTime - mIncreaseDomainTime)/mEndIncreaseDomainTime;
		gotCellPopulationDepth = GetCellPopulationDepth() + sqrt(0.75)*mMultiplyerIncreaseDomainTime*(mEndIncreaseDomainTime - mIncreaseDomainTime)/mEndIncreaseDomainTime ;
	}


	// TRACE("Imposing BCs");
	// for (AbstractCellPopulation<3>::Iterator cell_iter = this->mpCellPopulation->Begin();
	// 	 cell_iter != this->mpCellPopulation->End();
	// 	 ++cell_iter)
	for(unsigned node_index=0; node_index<this->mpCellPopulation->GetNumNodes(); node_index++)
	{
		
		// assert(this->mpCellPopulation->IsCellAssociatedWithADeletedLocation(*cell_iter) == false);

		// Get index of node associated with cell
        // unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);

        // Get pointer to this node
        Node<3>* p_node = this->mpCellPopulation->GetNode(node_index);

		double ghost_length = 0.0;
		try
        {
			(this->mpCellPopulation->GetCellUsingLocationIndex(node_index));
			// TRACE("Got cell");
			ghost_length = 0.0;
        }
        catch(Exception&)
        {
            // TRACE("No cell");
			ghost_length = 0.5;
        }

		// if (this->mpCellPopulation->GetCellUsingLocationIndex(node_index)->GetMutationState()->IsType<WildTypeCellMutationState>() ||
		// 	this->mpCellPopulation->GetCellUsingLocationIndex(node_index)->GetMutationState()->IsType<StromalCellMutationState>() )
		// {
		// 	ghost_length = 0.0;
		// }
		// else
		// {
		// 	ghost_length = 0.1;
		// }

		// Pin the bottom cells (only stromal)
		// if ( (this->mpCellPopulation->GetCellUsingLocationIndex(node_index)->GetMutationState()->IsType<StromalCellMutationState>())
		// 		&& (p_node->rGetLocation()[2] < mMaxHeightForPinnedCells) )
		// {
		// 	// Get old node location
		// 	typename std::map<Node<3>*, c_vector<double, 3> >::const_iterator it = rOldLocations.find(p_node);
    	// 	c_vector<double, 3> old_node_location = it->second;
		// 	// c_vector<double, 3> old_node_location = rOldLocations[node_index];

		// 	// Return node to old location
		// 	p_node->rGetModifiableLocation()[0] = old_node_location[0];
		// 	p_node->rGetModifiableLocation()[1] = old_node_location[1];
		// 	// p_node->rGetModifiableLocation()[2] = 0.0;
		// 	p_node->rGetModifiableLocation()[2] = mMaxHeightForPinnedCells;
		// }

		// Want to enforce a periodic boundary box, such that cells which move outside the boundaries of the box get
		// moved to the opposite side (as in the cylindrical mesh)
		if (p_node->rGetLocation()[0] < 0.0 + ghost_length)
		{
			/* Move to the opposite edge */
			if(ghost_length > 0)
			{
				p_node->rGetModifiableLocation()[0] = ghost_length;
			}
			else
			{
				p_node->rGetModifiableLocation()[0] = p_node->rGetLocation()[0] + gotCellPopulationWidth;
			}
		}
		else if (p_node->rGetLocation()[0] > gotCellPopulationWidth - ghost_length)
		{
			/* Move to the opposite edge */
			if(ghost_length > 0)
			{
				p_node->rGetModifiableLocation()[0] = gotCellPopulationWidth - ghost_length;
			}
			else
			{
				p_node->rGetModifiableLocation()[0] = p_node->rGetLocation()[0] - gotCellPopulationWidth;
			}
		}
		if (p_node->rGetLocation()[1] < 0.0 + ghost_length)
		{
			/* Move to the opposite edge */
			if(ghost_length > 0)
			{
				p_node->rGetModifiableLocation()[1] = ghost_length;
			}
			else
			{
				p_node->rGetModifiableLocation()[1] = p_node->rGetLocation()[1] + gotCellPopulationDepth;
			}
		
		}
		else if (p_node->rGetLocation()[1] > gotCellPopulationDepth - ghost_length)
		{
			/* Move to the opposite edge */
			if(ghost_length > 0)
			{
				p_node->rGetModifiableLocation()[1] = gotCellPopulationDepth - ghost_length;
			}
			else
			{
				p_node->rGetModifiableLocation()[1] = p_node->rGetLocation()[1] - gotCellPopulationDepth;
			}
		}
	}
}

bool PeriodicBoxBoundaryCondition3dGhostsIncreasingDomain::VerifyBoundaryCondition()
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
    // for (AbstractCellPopulation<3>::Iterator cell_iter = this->mpCellPopulation->Begin();
    //      cell_iter != this->mpCellPopulation->End();
    //      ++cell_iter)
	for(unsigned node_index=0; node_index<this->mpCellPopulation->GetNumNodes(); node_index++)
    {
        // Get index of node associated with cell
        // unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);

        // Get pointer to this node
        Node<3>* p_node = this->mpCellPopulation->GetNode(node_index);

		double ghost_length = 0.0;
		try
        {
			(this->mpCellPopulation->GetCellUsingLocationIndex(node_index));
			// TRACE("Got cell");
			ghost_length = 0.0;
        }
        catch(Exception&)
        {
            // TRACE("No cell");
			ghost_length = 0.5;
        }
		// if (this->mpCellPopulation->GetCellUsingLocationIndex(node_index)->GetMutationState()->IsType<WildTypeCellMutationState>() ||
		// 	this->mpCellPopulation->GetCellUsingLocationIndex(node_index)->GetMutationState()->IsType<StromalCellMutationState>() )
		// {
		// 	ghost_length = 0.0;
		// }
		// else
		// {
		// 	ghost_length = 0.1;
		// }

        // If this node lies below the z=0 boundary (NO! - Dom), or outside of the periodic boundaries, break and return false
        // If this node lies below the z=0 boundary, or outside of the periodic boundaries, break and return false
        // if ( (p_node->rGetLocation()[2] < 0.0) || (p_node->rGetLocation()[0] < 0.0) || (p_node->rGetLocation()[0] > mCellPopulationWidth)
        // 		|| (p_node->rGetLocation()[1] < 0.0) || (p_node->rGetLocation()[1] > mCellPopulationDepth) )
        if ( (p_node->rGetLocation()[0] < 0.0 + ghost_length) || (p_node->rGetLocation()[0] > gotCellPopulationWidth - ghost_length)
          || (p_node->rGetLocation()[1] < 0.0 + ghost_length) || (p_node->rGetLocation()[1] > gotCellPopulationDepth - ghost_length) )
        {
        	PRINT_2_VARIABLES(node_index, SimulationTime::Instance()->GetTime());
        	PRINT_VECTOR(p_node->rGetLocation());
			TRACE("Boundaries are:");
			PRINT_3_VARIABLES(ghost_length, gotCellPopulationWidth-ghost_length, gotCellPopulationDepth-ghost_length);
        	boundary_condition_satisfied = false;
        	break;
        }
    }

    return boundary_condition_satisfied;
}

void PeriodicBoxBoundaryCondition3dGhostsIncreasingDomain::SetMaxHeightForPinnedCells(double maxHeightForPinnedCells)
{
	mMaxHeightForPinnedCells = maxHeightForPinnedCells;
}

double PeriodicBoxBoundaryCondition3dGhostsIncreasingDomain::GetMaxHeightForPinnedCells()
{
	return mMaxHeightForPinnedCells;
}

void PeriodicBoxBoundaryCondition3dGhostsIncreasingDomain::SetCellPopulationWidth(double cellPopulationWidth)
{
	mCellPopulationWidth = cellPopulationWidth;
}

double PeriodicBoxBoundaryCondition3dGhostsIncreasingDomain::GetCellPopulationWidth()
{
	return mCellPopulationWidth;
}

void PeriodicBoxBoundaryCondition3dGhostsIncreasingDomain::SetCellPopulationDepth(double cellPopulationDepth)
{
	mCellPopulationDepth = cellPopulationDepth;
}

double PeriodicBoxBoundaryCondition3dGhostsIncreasingDomain::GetCellPopulationDepth()
{
	return mCellPopulationDepth;
}

void PeriodicBoxBoundaryCondition3dGhostsIncreasingDomain::SetIncreaseDomainTime(double increaseDomainTime)
{
	mIncreaseDomainTime = increaseDomainTime;
}

double PeriodicBoxBoundaryCondition3dGhostsIncreasingDomain::GetIncreaseDomainTime()
{
	return mIncreaseDomainTime;
}

void PeriodicBoxBoundaryCondition3dGhostsIncreasingDomain::SetEndIncreaseDomainTime(double endIncreaseDomainTime)
{
	mEndIncreaseDomainTime = endIncreaseDomainTime;
}

double PeriodicBoxBoundaryCondition3dGhostsIncreasingDomain::GetEndIncreaseDomainTime()
{
	return mEndIncreaseDomainTime;
}

void PeriodicBoxBoundaryCondition3dGhostsIncreasingDomain::SetMultiplyerIncreaseDomainTime(double multiplyerIncreaseDomainTime)
{
	mMultiplyerIncreaseDomainTime = multiplyerIncreaseDomainTime;
}

double PeriodicBoxBoundaryCondition3dGhostsIncreasingDomain::GetMultiplyerIncreaseDomainTime()
{
	return mMultiplyerIncreaseDomainTime;
}

void PeriodicBoxBoundaryCondition3dGhostsIncreasingDomain::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
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
CHASTE_CLASS_EXPORT(PeriodicBoxBoundaryCondition3dGhostsIncreasingDomain)
