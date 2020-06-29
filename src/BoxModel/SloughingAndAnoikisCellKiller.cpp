#include "SloughingAndAnoikisCellKiller.hpp"
#include "StromalCellMutationState.hpp"
#include "Debug.hpp"

SloughingAndAnoikisCellKiller::SloughingAndAnoikisCellKiller(AbstractCellPopulation<2>* pCrypt, double cellPopulationWidth)
    : AbstractCellKiller<2>(pCrypt),
      mCellsRemovedByAnoikis(0),
      mCellsRemovedBySloughing(0),
      mCellPopulationWidth(cellPopulationWidth)
{
}
    
SloughingAndAnoikisCellKiller::~SloughingAndAnoikisCellKiller()
{
}

//void SloughingAndAnoikisCellKiller::SetCellPopulationWidth(double cellPopulationWidth) const
//{
//	mCellPopulationWidth = cellPopulationWidth;
//}

double SloughingAndAnoikisCellKiller::GetCellPopulationWidth() const
{
	return mCellPopulationWidth;
}


/* Method to get the neighbouring nodes (including ghost nodes) of a particular node
 * Can then be used to identify the type of cells that surround a particular cell */

std::set<unsigned> SloughingAndAnoikisCellKiller::GetNeighbouringNodeIndices(unsigned nodeIndex)
{
	// Create a set of neighbouring node indices
	std::set<unsigned> neighbouring_node_indices;
	 	   
	// Need access to the mesh but can't get to it because the cell killer only owns a 
	// pointer to an AbstractCellPopulation		
    MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*> (this->mpCellPopulation);  
	
	// Find the indices of the elements owned by this node
    std::set<unsigned> containing_elem_indices = p_tissue->rGetMesh().GetNode(nodeIndex)->rGetContainingElementIndices();
    	 	
    // Iterate over these elements
    for (std::set<unsigned>::iterator elem_iter=containing_elem_indices.begin();
         elem_iter != containing_elem_indices.end();
	     ++elem_iter)	
   	    {    	
		    // Get all the nodes contained in this element
		    unsigned neighbour_global_index;

		    for (unsigned local_index=0; local_index<3; local_index++)
		    {
		    	neighbour_global_index = p_tissue->rGetMesh().GetElement(*elem_iter)->GetNodeGlobalIndex(local_index);

		    	// Don't want to include the original node or ghost nodes
		    	if( (neighbour_global_index != nodeIndex) && (!p_tissue->IsGhostNode(neighbour_global_index)) )		
	            {
	            	neighbouring_node_indices.insert(neighbour_global_index);
	            }
		    }   
	    }

    return neighbouring_node_indices;
}

/** Method to determine if a transit cell has lost all contacts with the differentiated cells below
 * TRUE if cell has popped up
 * FALSE if cell remains in the monolayer
 */
bool SloughingAndAnoikisCellKiller::HasCellPoppedUp(unsigned nodeIndex)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*> (this->mpCellPopulation);

	bool has_cell_popped_up = false;	// Initialising

   	std::set<unsigned> neighbours = GetNeighbouringNodeIndices(nodeIndex);

   	unsigned num_differentiated_neighbours = 0;

   	// Iterate over the neighbouring cells to check the number of differentiated cell neighbours

   	for(std::set<unsigned>::iterator neighbour_iter=neighbours.begin();
   							neighbour_iter != neighbours.end();
   							++neighbour_iter)
	{
		if ( (!p_tissue->IsGhostNode(*neighbour_iter)) && (p_tissue->GetCellUsingLocationIndex(*neighbour_iter)->GetMutationState()->IsType<StromalCellMutationState>()==true) )
   		{
			num_differentiated_neighbours += 1;
		}
   	}

   	if(num_differentiated_neighbours < 1)
   	{
   		has_cell_popped_up = true;
   	}

	return has_cell_popped_up;
}

/** A method to return a vector that indicates which cells should be killed by anoikis
 *
 */
std::vector<c_vector<unsigned,3> > SloughingAndAnoikisCellKiller::RemoveByAnoikisOrSloughing()
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*> (this->mpCellPopulation);
    assert(p_tissue->GetVoronoiTessellation()!=NULL);	// This fails during archiving of a simulation as Voronoi stuff not archived yet

    std::vector<c_vector<unsigned,3> > cells_to_remove;
    c_vector<unsigned,3> individual_node_information;	// Will store the node index and whether to remove or not (1 or 0)

	for (AbstractCellPopulation<2>::Iterator cell_iter = p_tissue->Begin();
    	 cell_iter != p_tissue->End();
    	 ++cell_iter)
	{
		unsigned node_index = p_tissue->GetNodeCorrespondingToCell(*cell_iter)->GetIndex();
		assert((!p_tissue->IsGhostNode(node_index)));

		// Initialise
		individual_node_information[0] = node_index;
		individual_node_information[1] = 0;
		individual_node_information[2] = 0;

		// Examine each epithelial node to see if it should be removed by anoikis and then if it should be removed by sloughing

		if ( cell_iter->GetMutationState()->IsType<StromalCellMutationState>()==false )
		{
			// Determining whether to remove this cell by anoikis

			if(this->HasCellPoppedUp(node_index))
			{
				individual_node_information[1] = 1;
			}
		}

		double x = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter)[0];
		
		if ((cell_iter->GetMutationState()->IsType<StromalCellMutationState>()==false)
				&& ((x < 0.0) || (x>mCellPopulationWidth)) )
		{
			individual_node_information[2] = 1;
		}

		cells_to_remove.push_back(individual_node_information);
	}

	return cells_to_remove;
}

/* Cell Killer that kills transit cells that move beyond the walls of the box
 * and also any transit cells that pop upwards and become detached from the stromal
 * cells
*/
void SloughingAndAnoikisCellKiller::CheckAndLabelCellsForApoptosisOrDeath()
{	
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*> (this->mpCellPopulation);
    assert(p_tissue->GetVoronoiTessellation()!=NULL);	// This fails during archiving of a simulation as Voronoi stuff not archived yet
    
    // Get the information at this timestep for each node index that says whether to remove by anoikis or sloughing
    std::vector<c_vector<unsigned,3> > cells_to_remove = this->RemoveByAnoikisOrSloughing();

    // Keep a record of how many cells have been removed at this timestep
    this->SetNumberCellsRemoved(cells_to_remove);

    // Need to avoid trying to kill any cells twice (i.e. both by anoikis or sloughing)
    // Loop over these vectors individually and kill any cells that they tell you to

    for (unsigned i=0; i<cells_to_remove.size(); i++)
    {
    	if ( (cells_to_remove[i][1] == 1) || (cells_to_remove[i][2] == 1) )
    	{
    		// Get cell associated to this node
    		CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(cells_to_remove[i][0]);
    		p_cell->Kill();
    	}
    }
}

void SloughingAndAnoikisCellKiller::SetNumberCellsRemoved(std::vector<c_vector<unsigned,3> > cellsRemoved)
{
	unsigned num_removed_by_anoikis = 0;
	unsigned num_removed_by_sloughing = 0;

    for (unsigned i=0; i<cellsRemoved.size(); i++)
    {
    	if(cellsRemoved[i][1]==1)
    	{
    		num_removed_by_anoikis+=1;
    	}
    	if(cellsRemoved[i][2]==1)
    	{
    		num_removed_by_sloughing+=1;
    	}
    }

    mCellsRemovedByAnoikis += num_removed_by_anoikis;
    mCellsRemovedBySloughing += num_removed_by_sloughing;
}

c_vector<unsigned,2> SloughingAndAnoikisCellKiller::GetNumberCellsRemoved()
{
	c_vector<unsigned,2> number_cells_removed;
	number_cells_removed[0] = mCellsRemovedByAnoikis;
	number_cells_removed[1] = mCellsRemovedBySloughing;
	return number_cells_removed;
}


void SloughingAndAnoikisCellKiller::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<SloughSides>" << mSloughSides << "</SloughSides> \n";
    *rParamsFile << "\t\t\t<CellPopulationWidth>" << mCellPopulationWidth << "</CellPopulationWidth> \n";
    *rParamsFile << "\t\t\t<CellsRemovedByAnoikis>" << mCellsRemovedByAnoikis << "</CellsRemovedByAnoikis> \n";
    *rParamsFile << "\t\t\t<CellsRemovedBySloughing>" << mCellsRemovedBySloughing << "</CellsRemovedBySloughing> \n";

    // Call direct parent class
    AbstractCellKiller<2>::OutputCellKillerParameters(rParamsFile);
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(SloughingAndAnoikisCellKiller)
