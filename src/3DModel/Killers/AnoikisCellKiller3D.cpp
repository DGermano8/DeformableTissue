#include "AnoikisCellKiller3D.hpp"
#include "Debug.hpp"

/* Apoptosis for cells that are epithelial and lose contact with the basement membrane
 *
 */
AnoikisCellKiller3D::AnoikisCellKiller3D(AbstractCellPopulation<3>* pCrypt)
    : AbstractCellKiller<3>(pCrypt),
    mCellsRemovedByAnoikis(0)
{
    // Sets up output file
//	OutputFileHandler output_file_handler(mOutputDirectory + "AnoikisData/", false);
//	mAnoikisOutputFile = output_file_handler.OpenOutputFile("results.anoikis");
}

AnoikisCellKiller3D::~AnoikisCellKiller3D()
{
//    mAnoikisOutputFile->close();
}

void AnoikisCellKiller3D::SetOutputDirectory(std::string outputDirectory)
{
	mOutputDirectory = outputDirectory;
}

std::string AnoikisCellKiller3D::GetOutputDirectory()
{
	return mOutputDirectory;
}

/*
 * Method to get the neighbouring nodes (excluding ghost nodes) of a particular node
 * Can then be used to identify the type of cells that surround a particular cell.
 */
std::set<unsigned> AnoikisCellKiller3D::GetNeighbouringNodeIndices(unsigned nodeIndex)
{
	// Create a set of neighbouring node indices
	std::set<unsigned> neighbouring_node_indices;

	// Need access to the mesh but can't get to it because the cell killer only owns a
	// pointer to an AbstractCellPopulation
    MeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<MeshBasedCellPopulationWithGhostNodes<3>*> (this->mpCellPopulation);

	// Find the indices of the elements owned by this node
    std::set<unsigned> containing_elem_indices = p_tissue->rGetMesh().GetNode(nodeIndex)->rGetContainingElementIndices();

    // Iterate over these elements
    for (std::set<unsigned>::iterator elem_iter=containing_elem_indices.begin();
         elem_iter != containing_elem_indices.end();
	     ++elem_iter)
    {
	    // Get all the nodes contained in this element
	    unsigned neighbour_global_index;

	    for (unsigned local_index=0; local_index<4; local_index++)
	    {
	    	neighbour_global_index = p_tissue->rGetMesh().GetElement(*elem_iter)->GetNodeGlobalIndex(local_index);
	    	// Don't want to include the original node or ghost nodes
	    	if( (neighbour_global_index != nodeIndex) && (!p_tissue->IsGhostNode(neighbour_global_index)) )
            {
            	neighbouring_node_indices.insert(neighbour_global_index);
            }
	    }
    }
    return neighbouring_node_indices;		// This will contain repeats, but it doesn't matter
}

/** Method to determine if an epithelial cell has lost all contacts with the stromal cells below
 * TRUE if cell has popped up
 * FALSE if cell remains in the monolayer
 */
bool AnoikisCellKiller3D::HasCellPoppedUp(unsigned nodeIndex)
{
	MeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<MeshBasedCellPopulationWithGhostNodes<3>*> (this->mpCellPopulation);

	bool has_cell_popped_up = false;	// Initialising

   	std::set<unsigned> neighbours = GetNeighbouringNodeIndices(nodeIndex);

   	unsigned num_stromal_neighbours = 0;

   	// Iterate over the neighbouring cells to check the number of stromal cell neighbours

   	for(std::set<unsigned>::iterator neighbour_iter=neighbours.begin();
   							neighbour_iter != neighbours.end();
   							++neighbour_iter)
	{
		if ( (!p_tissue->IsGhostNode(*neighbour_iter))
				&& (p_tissue->GetCellUsingLocationIndex(*neighbour_iter)->GetMutationState()->IsType<StromalCellMutationState>()==true) )
   		{
			num_stromal_neighbours += 1;
		}
   	}

   	if(num_stromal_neighbours < 1)
   	{
//   		PRINT_4_VARIABLES(SimulationTime::Instance()->GetTime(),nodeIndex, neighbours.size(), num_stromal_neighbours);
   		has_cell_popped_up = true;
   	}
   
	return has_cell_popped_up;
}

/** A method to return a vector that indicates which cells should be killed by anoikis
 */
std::vector<c_vector<unsigned,2> > AnoikisCellKiller3D::RemoveByAnoikis()
{
	MeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<MeshBasedCellPopulationWithGhostNodes<3>*> (this->mpCellPopulation);

//This needs fixing
//   assert(p_tissue->GetVoronoiTessellation()!=NULL);	// This fails during archiving of a simulation as Voronoi stuff not archived yet

    std::vector<c_vector<unsigned,2> > cells_to_remove;
    c_vector<unsigned,2> individual_node_information;	// Will store the node index and whether to remove or not (1 or 0)

	for (AbstractCellPopulation<3>::Iterator cell_iter = p_tissue->Begin();
    	 cell_iter != p_tissue->End();
    	 ++cell_iter)
	{
		unsigned node_index = p_tissue->GetNodeCorrespondingToCell(*cell_iter)->GetIndex();
		assert((!p_tissue->IsGhostNode(node_index)));
		
		// Initialise
		individual_node_information[0] = node_index;
		individual_node_information[1] = 0;

		// Examine each epithelial node to see if it should be removed by anoikis
		if (cell_iter->GetMutationState()->IsType<StromalCellMutationState>()==false)
		{
			assert(cell_iter->GetMutationState()->IsType<WildTypeCellMutationState>()==true);
			
			// Determining whether to remove this cell by anoikis
			if((!cell_iter->IsDead()) && (this->HasCellPoppedUp(node_index)==true) )
			{
				individual_node_information[1] = 1;
			}
		}

		cells_to_remove.push_back(individual_node_information);
	}

	return cells_to_remove;
}

/* Cell Killer that kills transit cells that move beyond the walls of the box
 * and also any transit cells that pop upwards and become detached from the stromal
 * cells
*/
void AnoikisCellKiller3D::CheckAndLabelCellsForApoptosisOrDeath()
{
	MeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<MeshBasedCellPopulationWithGhostNodes<3>*> (this->mpCellPopulation);

//This needs fixing
//    assert(p_tissue->GetVoronoiTessellation()!=NULL);	// This fails during archiving of a simulation as Voronoi stuff not archived yet

    // Get the information at this timestep for each node index that says whether to remove by anoikis 
    std::vector<c_vector<unsigned,2> > cells_to_remove = this->RemoveByAnoikis();

    // Keep a record of how many cells have been removed at this timestep
    this->SetNumberCellsRemovedByAnoikis(cells_to_remove);
    this->SetLocationsOfCellsRemovedByAnoikis(cells_to_remove);

    // Need to avoid trying to kill any cells twice (i.e. both by anoikis or random apoptosis)
    // Loop over these vectors individually and kill any cells that they tell you to

    for (unsigned i=0; i<cells_to_remove.size(); i++)
    {
    	if (cells_to_remove[i][1] == 1)
    	{
    		// Get cell associated to this node
    		CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(cells_to_remove[i][0]);
    		p_cell->Kill();
    	}
    }
}


void AnoikisCellKiller3D::SetNumberCellsRemovedByAnoikis(std::vector<c_vector<unsigned,2> > cellsRemoved)
{
	unsigned num_removed_by_anoikis = 0;
	
    for (unsigned i=0; i<cellsRemoved.size(); i++)
    {
    	if(cellsRemoved[i][1]==1)
    	{
    		num_removed_by_anoikis+=1;
    	}
    }

    mCellsRemovedByAnoikis += num_removed_by_anoikis;
}

unsigned AnoikisCellKiller3D::GetNumberCellsRemovedByAnoikis()
{
	return mCellsRemovedByAnoikis;
}

/* Data stored: time - node_index - x - y - z */
void AnoikisCellKiller3D::SetLocationsOfCellsRemovedByAnoikis(std::vector<c_vector<unsigned,2> > cellsRemoved)
{
	MeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<MeshBasedCellPopulationWithGhostNodes<3>*> (this->mpCellPopulation);
	double x_location, y_location, z_location;
	c_vector<double, 5> node_time_and_location;

	// Need to use the node indices to store the locations of where cells are removed
    for (unsigned i=0; i<cellsRemoved.size(); i++)
    {
    	if (cellsRemoved[i][1] == 1)		// This cell has been removed by anoikis
    	{
    		node_time_and_location[0] = SimulationTime::Instance()->GetTime();
    		node_time_and_location[1] = cellsRemoved[i][0];
    		
			unsigned node_index = cellsRemoved[i][0];

			CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(node_index);
			x_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[0];
			y_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[1];
			z_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[2];

			node_time_and_location[2] = x_location;
			node_time_and_location[3] = y_location;
			node_time_and_location[4] = z_location;

			mLocationsOfAnoikisCells.push_back(node_time_and_location);
    	}
    }
}

std::vector<c_vector<double,5> > AnoikisCellKiller3D::GetLocationsOfCellsRemovedByAnoikis()
{
	return mLocationsOfAnoikisCells;
}

void AnoikisCellKiller3D::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CellsRemovedByAnoikis>" << mCellsRemovedByAnoikis << "</CellsRemovedByAnoikis> \n";
//    *rParamsFile << "\t\t\t<XLocationsOfAnoikisCells>" << mXLocationsOfAnoikisCells << "</XLocationsOfAnoikisCells> \n";

    // Call direct parent class
    AbstractCellKiller<3>::OutputCellKillerParameters(rParamsFile);
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(AnoikisCellKiller3D)
