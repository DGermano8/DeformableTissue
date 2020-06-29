#include "TissueSlabSimulation.hpp"
#include "Debug.hpp"

bool SortVectorAccordingToXCoordinate(const c_vector<double, 2> lhs, const c_vector<double, 2> rhs)
{
    return lhs[0] < rhs[0];
}

bool SortVectorAccordingToYCoordinate(const c_vector<double, 2> lhs, const c_vector<double, 2> rhs)
{
    return lhs[1] < rhs[1];
}

//bool TissueSlabSimulation::StoppingEventHasOccurred()
//{
//	// /todo: Static cast the mesh in the constructor
//	Cylindrical2dMesh* p_mesh = dynamic_cast<Cylindrical2dMesh*>(&(mpStaticCastCellPopulation->rGetMesh()));
//
//	bool mismatched_elements = p_mesh->GetInstanceOfMismatchedBoundaryNodes();
//
//	return mismatched_elements;
//}

TissueSlabSimulation::TissueSlabSimulation(AbstractCellPopulation<2>& rCellPopulation,
                  bool deleteCellPopulationAndForceCollection,
                  bool initialiseCells)
    : OffLatticeSimulation<2>(rCellPopulation,
                          deleteCellPopulationAndForceCollection,
                          initialiseCells)
{
    mpStaticCastCellPopulation = static_cast<MeshBasedCellPopulationWithGhostNodes<2>*>(&mrCellPopulation);
}

TissueSlabSimulation::~TissueSlabSimulation()
{
}

//// Method to find the width of the mesh of real cells (so that it doesn't take ghost nodes into account)
//
//double TissueSlabSimulation::GetWidth(const unsigned& rDimension) const
//{
//    // We are in two dimensions
//    assert(rDimension < 2);
//
//    // We must have at least one cell
//    assert(mpStaticCastCellPopulation->GetNumRealCells() > 0u);
//
//    double max = -1e200;
//    double min = 1e200;
//
//    // Iterate over cells, not nodes
//    for (AbstractCellPopulation<2>::Iterator cell_iter = mpStaticCastCellPopulation->Begin();
//         cell_iter != mpStaticCastCellPopulation->End();
//         ++cell_iter)
//    {
//    	double this_node_value = mpStaticCastCellPopulation->GetNodeCorrespondingToCell(&(*cell_iter))->rGetLocation()[rDimension];	
//            if (this_node_value>max)
//            {
//                max = this_node_value;
//                
//            }
//            if (this_node_value < min)
//            {
//                min = this_node_value;
//            }
//    }
//
//    return max - min;    
//}

// Defining the node locations when a cell divides - for the transit cells on the top layer of the box
// model, we define them to divide horizontally

c_vector<double, 2> TissueSlabSimulation::CalculateCellDivisionVector(CellPtr pParentCell)
{
    double separation = mpStaticCastCellPopulation->GetMeinekeDivisionSeparation();
    c_vector<double, 2> parent_coords = mpStaticCastCellPopulation->GetLocationOfCellCentre(pParentCell);
    c_vector<double, 2> daughter_coords;

    // Move the parent cell backwards by 0.5*sep in the x-direction and return the position of 
    // the daughter cell (0.5*sep forwards in x-direction)

    // Make a direction vector of the required length
    c_vector<double, 2> direction_vector;
           
    direction_vector(0) = 0.5*separation;
    direction_vector(1) = 0.0;		// Horizontal division only 
   
    daughter_coords = parent_coords + direction_vector;  
    parent_coords -= direction_vector;

    assert(daughter_coords(1)>=0.0); // to make sure dividing cells stay in the tissue
    assert(parent_coords(1)>=0.0);   // to make sure dividing cells stay in the tissue

    // Set the parent to use this location
    ChastePoint<2> parent_coords_point(parent_coords);

    unsigned node_index = mpStaticCastCellPopulation->GetLocationIndexUsingCell(pParentCell);
    mrCellPopulation.SetNode(node_index, parent_coords_point);

    return daughter_coords;
}

/* The following option is used to have random direction division - used to investigate the effect of altered
 * mitotic spindle orientation.
 */

//c_vector<double, 2> TissueSlabSimulation::CalculateCellDivisionVector(CellPtr pParentCell)
//{
//    // Location of parent and daughter cells
//    c_vector<double, 2> parent_coords = mpStaticCastCellPopulation->GetLocationOfCellCentre(pParentCell);
//    c_vector<double, 2> daughter_coords;
//
//    // Get separation parameter
//	double separation = mpStaticCastCellPopulation->GetMeinekeDivisionSeparation();
//
//    // Make a random direction vector of the required length
//    c_vector<double, 2> random_vector;
//
//    /*
//     * Pick a random direction and move the parent cell backwards by 0.5*separation
//     * in that direction and return the position of the daughter cell 0.5*separation
//     * forwards in that direction.
//     */
//
//    double random_angle = RandomNumberGenerator::Instance()->ranf();
//    random_angle *= 2.0*M_PI;
//
//    random_vector(0) = 0.5*separation*cos(random_angle);
//    random_vector(1) = 0.5*separation*sin(random_angle);
//
//    c_vector<double, 2> proposed_new_parent_coords = parent_coords - random_vector;
//    c_vector<double, 2> proposed_new_daughter_coords = parent_coords + random_vector;
//
//    if (   (proposed_new_parent_coords(1) >= 0.0)
//        && (proposed_new_daughter_coords(1) >= 0.0))
//    {
//        // We are not too close to the bottom of the tissue, so move parent
//        parent_coords = proposed_new_parent_coords;
//        daughter_coords = proposed_new_daughter_coords;
//    }
//    else
//    {
//        proposed_new_daughter_coords = parent_coords + 2.0*random_vector;
//        while (proposed_new_daughter_coords(1) < 0.0)
//        {
//            random_angle = RandomNumberGenerator::Instance()->ranf();
//            random_angle *= 2.0*M_PI;
//
//            random_vector(0) = separation*cos(random_angle);
//            random_vector(1) = separation*sin(random_angle);
//            proposed_new_daughter_coords = parent_coords + random_vector;
//        }
//        daughter_coords = proposed_new_daughter_coords;
//    }
//
//    assert(daughter_coords(1) >= 0.0); // to make sure dividing cells stay in the tissue
//    assert(parent_coords(1) >= 0.0);   // to make sure dividing cells stay in the tissue
//
//    // Set the parent to use this location
//    ChastePoint<2> parent_coords_point(parent_coords);
//
//    unsigned node_index = mpStaticCastCellPopulation->GetLocationIndexUsingCell(pParentCell);
//    mrCellPopulation.SetNode(node_index, parent_coords_point);
//
//    return daughter_coords;
//}

std::set<unsigned> TissueSlabSimulation::GetNeighbouringNodeIndices(unsigned nodeIndex)
{
    // Create a set of neighbouring node indices
	std::set<unsigned> neighbouring_node_indices; 
    
    // Find the indices of the elements (triangles!) owned by this node
	std::set<unsigned> containing_elem_indices = mpStaticCastCellPopulation->GetNode(nodeIndex)->rGetContainingElementIndices();

    // Iterate over these elements
    for (std::set<unsigned>::iterator elem_iter=containing_elem_indices.begin();
         elem_iter != containing_elem_indices.end();
         ++elem_iter)
    {
        // Get all the nodes contained in this element
        // (note that we've also included the original node, but this doesn't matter too much)
        unsigned neighbour_global_index;
        for (unsigned local_index=0; local_index<3; local_index++)
        {
            neighbour_global_index = mpStaticCastCellPopulation->rGetMesh().GetElement(*elem_iter)->GetNodeGlobalIndex(local_index);

            // Don't want to include the original node or ghost nodes
	    	if( (neighbour_global_index != nodeIndex) && (!mpStaticCastCellPopulation->IsGhostNode(neighbour_global_index)) )
            {
            	neighbouring_node_indices.insert(neighbour_global_index);
            }
        }
    }
    
    return neighbouring_node_indices;
}

/** Method to calculate the mean horizontal and vertical gap between pairs of transit cells in the 
 * monolayer. This ignores any transit cells that have popped up.
 */

c_vector<double,2> TissueSlabSimulation::CalculateMeanGapBetweenTransitCells()
{		
	double total_horizontal_gap = 0.0;	// Initialise
	double total_vertical_gap = 0.0;
	
	std::vector<c_vector<double, 2> > cell_coordinates;	
	
	// We use an iterator provided by the tissue to loop over cells 
	for (AbstractCellPopulation<2>::Iterator cell_iter = mrCellPopulation.Begin();
									 cell_iter != mrCellPopulation.End();
									 ++cell_iter)
	{    	  
		Node<2>* p_node = mpStaticCastCellPopulation->GetNodeCorrespondingToCell(*cell_iter);	// Pointer to node
		unsigned node_index = p_node->GetIndex();	
		
		// If the cell is a transit cell, and it has not popped up
		
	   	if((cell_iter->GetCellCycleModel()->GetCellProliferativeType() == TRANSIT) && !(HasCellPoppedUp(node_index)))
	   	{       		   		   		
	   		// Getting a vector of x coords of the transit cells which we will then sort	   		
	   		c_vector<double,2> coord_of_cell = mrCellPopulation.GetLocationOfCellCentre(*cell_iter);
	   		cell_coordinates.push_back(coord_of_cell);	   		    			
	   	}        	
	}
	
	// Now sort this vector of node coordinates in increasing order according to their x coordinate
	sort(cell_coordinates.begin(),cell_coordinates.end(),SortVectorAccordingToXCoordinate);
	
	for(unsigned i=0; i<cell_coordinates.size()-1; i++)
	{
		assert(cell_coordinates[i+1][0] >=  cell_coordinates[i][0]);		// Until you fix test for this
	}
	
	// Now we loop over the coordinates vector to find the distance between neighbouring cells
	// both in the horizontal and vertical direction
		
	for (unsigned i = 0; i<(cell_coordinates.size()-1); i++) 
	{			
		assert(cell_coordinates[i+1][0] >= cell_coordinates[i][0]);
			
		total_horizontal_gap +=  (cell_coordinates[i+1][0] - cell_coordinates[i][0]);
		total_vertical_gap +=  fabs(cell_coordinates[i+1][1] - cell_coordinates[i][1]);
	}
			
	// Now use this to find the average horizontal gap and vertical gap
	
	c_vector<double,2> average_gap;
	
	double average_horizontal_gap = total_horizontal_gap/(cell_coordinates.size()-1);
	average_gap(0) = average_horizontal_gap;
	double average_vertical_gap = total_vertical_gap/(cell_coordinates.size()-1);
	average_gap(1) = average_vertical_gap;
	
	return average_gap;
}

/** Method to calculate the maximum vertical gap between any pair of transit cells in the 
 * monolayer. This ignores any transit cells that have popped up.
 */

double TissueSlabSimulation::CalculateMaxVerticalGapBetweenTransitCells()
{			
	std::vector<c_vector<double, 2> > cell_coordinates;	
	
	// We use an iterator provided by the tissue to loop over cells 
	for (AbstractCellPopulation<2>::Iterator cell_iter = mrCellPopulation.Begin();
									 cell_iter != mrCellPopulation.End();
									 ++cell_iter)
	{    	  
		Node<2>* p_node = mpStaticCastCellPopulation->GetNodeCorrespondingToCell(*cell_iter);	// Pointer to node
		unsigned node_index = p_node->GetIndex();	
		
		// If the cell is a transit cell, and it has not popped up (don't think the HasCellPoppedUp method works..)
		
	   	if((cell_iter->GetCellCycleModel()->GetCellProliferativeType() == TRANSIT) && !(HasCellPoppedUp(node_index)))
	   	{       		   		   		
	   		// Getting a vector of x coords of the transit cells which we will then sort	   		
	   		c_vector<double,2> coord_of_cell = mrCellPopulation.GetLocationOfCellCentre(*cell_iter);
	   		cell_coordinates.push_back(coord_of_cell);	   		    			
	   	}        	
	}
	
	// Now sort this vector of node coordinates in increasing order according to their y coordinate
	sort(cell_coordinates.begin(),cell_coordinates.end(),SortVectorAccordingToYCoordinate);
	
	for(unsigned i=0; i<cell_coordinates.size()-1; i++)
	{
		assert(cell_coordinates[i+1][1] >=  cell_coordinates[i][1]);		// Until you fix test for this
	}
	
	// Now we subtract the smallest y coordinate from the largest to get the max vertical gap

	unsigned end_node_index = cell_coordinates.size()-1;
	
	double max_vertical_gap = cell_coordinates[end_node_index][1] - cell_coordinates[0][1];
	
	return max_vertical_gap;
}


/** Method to determine if a transit cell has lost all contacts with the differentiated cells below
 * TRUE if cell has popped up
 * FALSE if cell remains in the monolayer
 */
bool TissueSlabSimulation::HasCellPoppedUp(unsigned nodeIndex)
{
	bool has_cell_popped_up = false;	// Initialising

   	std::set<unsigned> neighbours = GetNeighbouringNodeIndices(nodeIndex);

   	unsigned num_differentiated_neighbours = 0;

   	// Iterate over the neighbouring cells to check the number of differentiated cell neighbours

   	for(std::set<unsigned>::iterator neighbour_iter=neighbours.begin();
   							neighbour_iter != neighbours.end();
   							++neighbour_iter)
	{
		if (!mpStaticCastCellPopulation->IsGhostNode(*neighbour_iter))
   		{
	   		if(mrCellPopulation.GetCellUsingLocationIndex(*neighbour_iter)->GetCellCycleModel()->GetCellProliferativeType()==DIFFERENTIATED)
	   		{
	   			num_differentiated_neighbours += 1;
	   		}
		}
   	}

   	if(num_differentiated_neighbours < 1)
   	{
   		has_cell_popped_up = true;
   	}

	return has_cell_popped_up;
}

//c_vector<unsigned,2> TissueSlabSimulation::GetNumberCellsRemoved()
//{
//	TRACE("CALCULATING NUMBER OF CELLS REMOVED");
//	c_vector<unsigned,2> num_cells_removed;
//
//	unsigned removed_by_anoikis = 0;
//	unsigned removed_by_sloughing = 0;
//
//	// Loop over all cells to see if they have been removed
//	for (AbstractCellPopulation<2>::Iterator cell_iter = mrCellPopulation.Begin();
//	                            cell_iter != mrCellPopulation.End();
//	                            ++cell_iter)
//	{
//		if (cell_iter->GetCellCycleModel()->GetCellProliferativeType()==TRANSIT)
//		{
//			unsigned node_index = mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);
//
//			bool has_cell_popped_up = this->HasCellPoppedUp(node_index);
//
//			if (this->HasCellPoppedUp(node_index))
//			{
//				TRACE("REMOVING A CELL BY ANOIKIS");
//				removed_by_anoikis += 1;
//			}
//
//            double x = mrCellPopulation.GetLocationOfCellCentre(*cell_iter)[0];
//
//            if( (mrCellPopulation.GetLocationOfCellCentre(*cell_iter)[0] < 0.0) || (mrCellPopulation.GetLocationOfCellCentre(*cell_iter)[0] > mCellPopulationWidth) )
//			{
//            	TRACE("REMOVING A CELL BY SLOUGHING");
//            	removed_by_sloughing += 1;
//			}
//		}
//	}
//
//	num_cells_removed[0] = removed_by_anoikis;
//	num_cells_removed[1] = removed_by_sloughing;
//	return num_cells_removed;
//}

void TissueSlabSimulation::WriteVisualizerSetupFile()
{
    *mpVizSetupFile << "MeshWidth\t" << mpStaticCastCellPopulation->rGetMesh().GetWidth(0u) << "\n";// get furthest distance between nodes in the x-direction
}

std::string TissueSlabSimulation::GetDataOutputFile()
{
	return mDataOutputFile;
}

void TissueSlabSimulation::AfterSolve()
{
    if (   ( mrCellPopulation.Begin() != mrCellPopulation.End() )  // there are any cells
        && ( dynamic_cast<AbstractVanLeeuwen2009WntSwatCellCycleModel*>(mrCellPopulation.Begin()->GetCellCycleModel())) ) // assume all the cells are the same
    {
        mBetaCatResultsFile->close();
    }

    OffLatticeSimulation<2>::AfterSolve();
}

void TissueSlabSimulation::WriteBasalLaminaResultsToDirectory(double basalLaminaParameter)
{
    // Open the output file
    OutputFileHandler output_file_handler(this->mSimulationOutputDirectory+"/", false);
    std::stringstream file_name;
    file_name << "basal_lamina_" << basalLaminaParameter << "_results.dat";
    out_stream output_file = output_file_handler.OpenOutputFile(file_name.str());

    // Outputting the total number and average height of the transit cells in the simulation

    double num_transit_cells = 0.0;
    double total_of_all_heights = 0.0;

    // We use an iterator provided by the tissue to loop over cells
    for (AbstractCellPopulation<2>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
        if (cell_iter->GetCellCycleModel()->GetCellProliferativeType() == TRANSIT)
        {
            num_transit_cells++;

            double cell_height = mrCellPopulation.GetLocationOfCellCentre(*cell_iter)[1];

            total_of_all_heights += cell_height;
        }
    }

    // The following variable is unused, triggered a warning/error on some compilers
//    double average_height = total_of_all_heights/num_transit_cells;
       
    // Now want to output the average horizontal distance between transit nodes in the
    // monolayer
    
    c_vector<double, 2> average_gap = CalculateMeanGapBetweenTransitCells();
    
    double average_horizontal_gap = average_gap[0];
    double average_vertical_gap = average_gap[1];
    double average_gap_ratio = average_vertical_gap/average_horizontal_gap;		
    
    // Data outputted in rows so that can be collated to plot more easily
        
    *output_file << basalLaminaParameter << "\t" << num_transit_cells << "\t" << average_horizontal_gap << "\t" <<
    		average_vertical_gap << "\t" << average_gap_ratio << "\n";
    //*output_file << "End time of simulation is " << SimulationTime::Instance()->GetTime() << "\n";
       
    output_file->close();
}

c_vector<double,6> TissueSlabSimulation::OutputVectorOfBasalLaminaResults(double basalLaminaParameter)
{
    // Outputting the average height of the transit cells in the simulation

    double num_transit_cells = 0.0;
    double total_of_all_heights = 0.0;

    // We use an iterator provided by the tissue to loop over cells
    for (AbstractCellPopulation<2>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
        if (cell_iter->GetCellCycleModel()->GetCellProliferativeType() == TRANSIT)
        {
            num_transit_cells++;

            double cell_height = mrCellPopulation.GetLocationOfCellCentre(*cell_iter)[1];

            total_of_all_heights += cell_height;
        }
    }
        
//    double average_height = total_of_all_heights/num_transit_cells;
    
    // Now want to find the average horizontal distance between transit nodes in the
    // monolayer
    
    c_vector<double, 2> average_gap = CalculateMeanGapBetweenTransitCells();
    
    double average_horizontal_gap = average_gap[0];
    double average_vertical_gap = average_gap[1];
    double average_gap_ratio = average_vertical_gap/average_horizontal_gap;		
    
    // Now, in the event that that major buckling of the layer occurs, it is useful to have the maximum 
    // vertical difference between any pair of transit cells - would indicate if the layer has moved upwards
    // at one end
    
    double max_height_difference = CalculateMaxVerticalGapBetweenTransitCells();     
    
    c_vector<double,6> results;
    
    results[0] = basalLaminaParameter;
    results[1] = num_transit_cells;
    results[2] = average_horizontal_gap;
    results[3] = average_vertical_gap;
    results[4] = average_gap_ratio;
    results[5] = max_height_difference;

    return results;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(TissueSlabSimulation)
