#include "CryptSimulation3dGhosts.hpp"
#include "Debug.hpp"

CryptSimulation3dGhosts::CryptSimulation3dGhosts(AbstractCellPopulation<3>& rCellPopulation,
                  bool deleteCellPopulationAndForceCollection,
                  bool initialiseCells)
    : OffLatticeSimulation<3>(rCellPopulation,
                          deleteCellPopulationAndForceCollection,
                          initialiseCells)
{
    mpStaticCastCellPopulation = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*>(&mrCellPopulation);
    mFoundNonExistentNode = false;
}


bool CryptSimulation3dGhosts::StoppingEventHasOccurred()
{
	bool incompatible_bc = false;//this->GetInstanceOfIncompatibleBoundaryConditions();
	bool non_existent_node = GetInstanceOfNonExistentNode();
	bool stopping_event = false;
	
	if ( (incompatible_bc) || (non_existent_node) )
	{
		stopping_event = true;
	}

	return stopping_event;
}

/*
 * Remove repeated entries in a 1D vector
 */
void CryptSimulation3dGhosts::RemoveDuplicates1D(std::vector<unsigned>& rVectorWithDuplicates)
{
    std::sort(rVectorWithDuplicates.begin(), rVectorWithDuplicates.end());
    rVectorWithDuplicates.erase(std::unique(rVectorWithDuplicates.begin(), rVectorWithDuplicates.end()), rVectorWithDuplicates.end());
}

/* Defining the node locations when a cell divides. For the 3D surface, we want cell division to occur
 * within the plane, and so parallel to the basement membrane. So we'll find the average of all of the 
 * stromal nodes that the parent epithelial cell is connected to, and instruct it to divide along the 
 * normal vector to that connecting this average to the parent node. So it's kind of like the node is 
 * dividing parallel to the basement membrane...
 */
c_vector<double,3> CryptSimulation3dGhosts::CalculateCellDivisionVector(CellPtr pParentCell)
{
	// Check that the parent cell is an epithelial cell
    assert(pParentCell->GetMutationState()->IsType<StromalCellMutationState>()== false);
    assert(pParentCell->GetMutationState()->IsType<WildTypeCellMutationState>()== true);

	double separation = mpStaticCastCellPopulation->GetMeinekeDivisionSeparation();

    unsigned parent_node_index = mpStaticCastCellPopulation->GetLocationIndexUsingCell(pParentCell);

    c_vector<double, 3> parent_coords = mpStaticCastCellPopulation->GetLocationOfCellCentre(pParentCell); //p_node->rGetLocation();
    c_vector<double, 3> daughter_coords, division_direction;

    // Now going to try cell division where it averages the positions of all the neighbouring epithelial nodes and divides along
    // the vector connecting the parent to this averaged position

    // Get all neighbouring cells
    std::set<unsigned> neighbours = GetNeighbouringNodeIndices(parent_node_index);
    
    // Want to isolate only the epithelial cells, and average their locations
    std::set<unsigned> epithelial_neighbour_indices;
    c_vector<double, 3> sum_of_locations; 	// Vector of coordinates of only the epithelial nodes it's connected to

    sum_of_locations[0] = 0.0;
	sum_of_locations[1] = 0.0;
	sum_of_locations[2] = 0.0;
	unsigned node_index;

    for(std::set<unsigned>::iterator neighbour_iter = neighbours.begin();
			         neighbour_iter != neighbours.end();
				     ++neighbour_iter)
    {
    	// Want to check if the node is deleted, as it should be if that cell has been killed, and not consider that node in this case
    	// This is the point at which you generate the error that index > num_nodes
    	if (!mpStaticCastCellPopulation->rGetMesh().GetNode(*neighbour_iter)->IsDeleted())
    	{
			if(!(mpStaticCastCellPopulation->IsGhostNode(*neighbour_iter)))
			{	
				CellPtr p_neighbour_cell = mpStaticCastCellPopulation->GetCellUsingLocationIndex(*neighbour_iter);

				// Only want neighbouring epithelial nodes that are not associated to dead cells
				if( (p_neighbour_cell->GetMutationState()->IsType<StromalCellMutationState>()== false)
						&& (!(*p_neighbour_cell).IsDead()) )
				{
					node_index = mpStaticCastCellPopulation->rGetMesh().GetNode(*neighbour_iter)->GetIndex();
					epithelial_neighbour_indices.insert(node_index);
					sum_of_locations[0] += mpStaticCastCellPopulation->rGetMesh().GetNode(*neighbour_iter)->rGetLocation()[0];
					sum_of_locations[1] += mpStaticCastCellPopulation->rGetMesh().GetNode(*neighbour_iter)->rGetLocation()[1];
					sum_of_locations[2] += mpStaticCastCellPopulation->rGetMesh().GetNode(*neighbour_iter)->rGetLocation()[2];
				}
			}
    	}
    }
    
	assert(!isnan(epithelial_neighbour_indices.size()));
	assert(epithelial_neighbour_indices.size() != 0);
	
	c_vector<double, 3> random_vector;
	random_vector[0] = 0.0;
	random_vector[1] = 0.0;
	random_vector[2] = 0.0;
	
    // If it only has one epithelial neighbour just divide randomly
	if (epithelial_neighbour_indices.size() == 1)
	{
        double random_zenith_angle = M_PI*RandomNumberGenerator::Instance()->ranf(); // phi
        double random_azimuth_angle = 2*M_PI*RandomNumberGenerator::Instance()->ranf(); // theta

        random_vector[0] = 0.5*separation*cos(random_azimuth_angle)*sin(random_zenith_angle);
        random_vector[1] = 0.5*separation*sin(random_azimuth_angle)*sin(random_zenith_angle);
        random_vector[2] = 0.5*separation*cos(random_zenith_angle);

        parent_coords = parent_coords - random_vector;
        daughter_coords = parent_coords + random_vector;

        // Set the parent to use this location
        ChastePoint<3> parent_coords_point(parent_coords);
        unsigned temp_node_index = mrCellPopulation.GetLocationIndexUsingCell(pParentCell);
        mrCellPopulation.SetNode(temp_node_index, parent_coords_point);		
	}
	else 
	{
		// Work out the average location of these stromal nodes
		c_vector<double, 3> average_location;
		average_location[0] = sum_of_locations[0]/(double)epithelial_neighbour_indices.size();
		average_location[1] = sum_of_locations[1]/(double)epithelial_neighbour_indices.size();
		average_location[2] = sum_of_locations[2]/(double)epithelial_neighbour_indices.size();

//		PRINT_VARIABLE(epithelial_neighbour_indices.size());
//		PRINT_VECTOR(parent_coords);
//		PRINT_VECTOR(average_location);
		// Now get the vector from the parent node to this average position
		c_vector<double, 3> vector_to_average_epithelial_node = mpStaticCastCellPopulation->rGetMesh().GetVectorFromAtoB(parent_coords, average_location);

		if( (vector_to_average_epithelial_node[0]==0) && (vector_to_average_epithelial_node[1]==0) && (vector_to_average_epithelial_node[2]==0) )
		{
			// Divide along vector connecting this to one of the epithelial nodes
			vector_to_average_epithelial_node = mpStaticCastCellPopulation->rGetMesh().GetVectorFromAtoB(parent_coords,mpStaticCastCellPopulation->rGetMesh().GetNode(node_index)->rGetLocation());
		}
		
//		PRINT_VECTOR(vector_to_average_epithelial_node);
		
		// Get this unit vector because this will be the division direction
		double length = norm_2(vector_to_average_epithelial_node);
		vector_to_average_epithelial_node /= length;
				
		/* Using this direction, move the parent cell backwards by 0.5*separation
		 * in that direction and return the position of the daughter cell 0.5*separation
		 * forwards in that direction.
		 */
		vector_to_average_epithelial_node *= 0.5*separation;

		parent_coords = parent_coords - vector_to_average_epithelial_node;
		daughter_coords = parent_coords + vector_to_average_epithelial_node;
		

		// Set the parent to use this location
		ChastePoint<3> parent_coords_point(parent_coords);
		mrCellPopulation.SetNode(parent_node_index, parent_coords_point);
		
	}
    
    
//    // This averages the stromal nodes instead
//    
//    // Get all neighbouring cells
//    std::set<unsigned> neighbours = GetNeighbouringNodeIndices(parent_node_index);
//        
//    // Want to isolate only the stromal cells, and average their locations
//    std::set<unsigned> stromal_neighbour_indices;
//    c_vector<double, 3> sum_of_locations; 	// Vector of coordinates of only the stromal nodes it's connected to
//
//    sum_of_locations[0] = 0.0;
//	sum_of_locations[1] = 0.0;
//	sum_of_locations[2] = 0.0;
//	unsigned node_index;
//
//    for(std::set<unsigned>::iterator neighbour_iter = neighbours.begin();
//			         neighbour_iter != neighbours.end();
//				     ++neighbour_iter)
//    {    	
//    	// Want to check if the node is deleted, as it should be if that cell has been killed, and not consider that node in this case
//    	// This is the point at which you generate the error that index > num_nodes
//    	if (!mpStaticCastCellPopulation->rGetMesh().GetNode(*neighbour_iter)->IsDeleted())
//    	{
//	    	CellPtr p_neighbour_cell = mpStaticCastCellPopulation->GetCellUsingLocationIndex(*neighbour_iter);
//
//	    	// Only want neighbouring stromal nodes that are not associated to dead cells
//	    	if( (!mpStaticCastCellPopulation->IsGhostNode(*neighbour_iter))
//	    			&& (p_neighbour_cell->GetMutationState()->IsType<StromalCellMutationState>()== true)
//	    			&& (!(*p_neighbour_cell).IsDead()) )
//	    	{
//		   		node_index = mpStaticCastCellPopulation->rGetMesh().GetNode(*neighbour_iter)->GetIndex();
//		   		stromal_neighbour_indices.insert(node_index);
//		   		sum_of_locations[0] += mpStaticCastCellPopulation->rGetMesh().GetNode(*neighbour_iter)->rGetLocation()[0];
//		   		sum_of_locations[1] += mpStaticCastCellPopulation->rGetMesh().GetNode(*neighbour_iter)->rGetLocation()[1];
//		   		sum_of_locations[2] += mpStaticCastCellPopulation->rGetMesh().GetNode(*neighbour_iter)->rGetLocation()[2];
//	    	}
//    	}
//    }
//
//	assert(!isnan(stromal_neighbour_indices.size()));
//
//	// Should include the case to cope with anoikis being switch off, and epithelial cells
//	// legitimately losing contact with the membrane. In this event, cell division happens in a random
//	// direction
//	c_vector<double, 3> random_vector;
//
//	if (stromal_neighbour_indices.size() == 0)
//	{
//        double random_zenith_angle = M_PI*RandomNumberGenerator::Instance()->ranf(); // phi
//        double random_azimuth_angle = 2*M_PI*RandomNumberGenerator::Instance()->ranf(); // theta
//
//        random_vector[0] = 0.5*separation*cos(random_azimuth_angle)*sin(random_zenith_angle);
//        random_vector[1] = 0.5*separation*sin(random_azimuth_angle)*sin(random_zenith_angle);
//        random_vector[2] = 0.5*separation*cos(random_zenith_angle);
//
//        parent_coords = parent_coords - random_vector;
//        daughter_coords = parent_coords + random_vector;
//
//        // Set the parent to use this location
//        ChastePoint<3> parent_coords_point(parent_coords);
//        unsigned temp_node_index = mrCellPopulation.GetLocationIndexUsingCell(pParentCell);
//        mrCellPopulation.SetNode(temp_node_index, parent_coords_point);
//	}
//
//	else
//	{
//		assert(stromal_neighbour_indices.size() != 0);
//
//		// Work out the average location of these stromal nodes
//		c_vector<double, 3> average_location;
//		average_location[0] = sum_of_locations[0]/(double)stromal_neighbour_indices.size();
//		average_location[1] = sum_of_locations[1]/(double)stromal_neighbour_indices.size();
//		average_location[2] = sum_of_locations[2]/(double)stromal_neighbour_indices.size();
//
//		// Now get the vector from the parent node to this average position
//		c_vector<double, 3> vector_to_average_stromal_node = mpStaticCastCellPopulation->rGetMesh().GetVectorFromAtoB(parent_coords, average_location);
//
//		// Location of one of the neighbouring stromal nodes
//		c_vector<double, 3> stromal_node_location = mpStaticCastCellPopulation->rGetMesh().GetNode(node_index)->rGetLocation();
//
//		// Now need the perpendicular vector to this and one other random vector
//		double random_zenith_angle = M_PI*RandomNumberGenerator::Instance()->ranf(); // phi
//		double random_azimuth_angle = 2*M_PI*RandomNumberGenerator::Instance()->ranf(); // theta
//
//		random_vector[0] = 0.5*separation*cos(random_azimuth_angle)*sin(random_zenith_angle);
//		random_vector[1] = 0.5*separation*sin(random_azimuth_angle)*sin(random_zenith_angle);
//		random_vector[2] = 0.5*separation*cos(random_zenith_angle);
//
//		// Now find the normal to these two vectors using the cross product
//		c_vector<double, 3> normal;
//		normal[0] = vector_to_average_stromal_node[1]*random_vector[2] - vector_to_average_stromal_node[2]*random_vector[1];
//		normal[1] = - vector_to_average_stromal_node[0]*random_vector[2] + vector_to_average_stromal_node[2]*random_vector[0];
//		normal[2] = vector_to_average_stromal_node[0]*random_vector[1] - vector_to_average_stromal_node[1]*random_vector[0];
//
//		// Get unit vector
//		double length = norm_2(normal);
//		normal /= length;
//
//		/* Using this direction, move the parent cell backwards by 0.5*separation
//		 * in that direction and return the position of the daughter cell 0.5*separation
//		 * forwards in that direction.
//		 */
//		normal *= 0.5*separation;
//
//		parent_coords = parent_coords - normal;
//		daughter_coords = parent_coords + normal;
//
//		// Set the parent to use this location
//		ChastePoint<3> parent_coords_point(parent_coords);
//		mrCellPopulation.SetNode(parent_node_index, parent_coords_point);
//	}


	/* Instead, can just define the cells to divide along the same z-coord as the parent
	 */

//	// Location of parent and daughter cells
//    c_vector<double, 3> parent_coords = mrCellPopulation.GetLocationOfCellCentre(pParentCell);
//    c_vector<double, 3> daughter_coords;
//
//    // Get separation parameter
//    double separation = static_cast<AbstractCentreBasedCellPopulation<3>*>(&mrCellPopulation)->GetMeinekeDivisionSeparation();
//
//    // Make a random direction vector of the required length
//    c_vector<double, 3> random_vector;
//
//    /*
//     * Pick a random direction and move the parent cell backwards by 0.5*separation
//     * in that direction and return the position of the daughter cell 0.5*separation
//     * forwards in that direction. Keep within the epithelial plane, so have z=0 in the division vector.
//     */
//
//	double random_angle = 2.0*M_PI*RandomNumberGenerator::Instance()->ranf();
//
//	random_vector[0] = 0.5*separation*cos(random_angle);
//	random_vector[1] = 0.5*separation*sin(random_angle);
//	random_vector[2] = 0.0;
//
//    parent_coords = parent_coords - random_vector;
//    daughter_coords = parent_coords + random_vector;
//
//    assert(daughter_coords[2] == parent_coords[2]);
//
//    // Set the parent to use this location
//    ChastePoint<3> parent_coords_point(parent_coords);
//    unsigned node_index = mrCellPopulation.GetLocationIndexUsingCell(pParentCell);
//    mrCellPopulation.SetNode(node_index, parent_coords_point);
//
	PRINT_VECTOR(daughter_coords);
	
    return daughter_coords;
}
 
/*
 * Method to return the nodes connected to a particular node via the Delaunay
 * triangulation, excluding ghost nodes.
 */
std::set<unsigned> CryptSimulation3dGhosts::GetNeighbouringNodeIndices(unsigned nodeIndex)
{
    // Get pointer to this node
    Node<3>* p_node = mpStaticCastCellPopulation->GetNode(nodeIndex);

    // Loop over containing elements
    std::set<unsigned> neighbouring_node_indices;
    
    for (Node<3>::ContainingElementIterator elem_iter = p_node->ContainingElementsBegin();
         elem_iter != p_node->ContainingElementsEnd();
         ++elem_iter)
    {
        // Get pointer to this containing element
        Element<3,3>* p_element = mpStaticCastCellPopulation->rGetMesh().GetElement(*elem_iter);
       
        // Loop over nodes contained in this element
        for (unsigned i=0; i<p_element->GetNumNodes(); i++)
        {
            // Get index of this node and add its index to the set if not the original node
            unsigned node_index = p_element->GetNodeGlobalIndex(i);
            
            if (node_index > 15000000)
            {
            	TRACE("It did seem to find some randomly massive node index as a global index");
            	PRINT_VARIABLE(node_index);
            	mFoundNonExistentNode = true;
            }
            
            if (node_index != nodeIndex)
            {
                neighbouring_node_indices.insert(node_index);
            }
        }
    }
    return neighbouring_node_indices;
	
		
//	// Create a set of neighbouring node indices
//	std::set<unsigned> neighbouring_node_indices;
//
//    // Find the indices of the elements owned by this node
//	std::set<unsigned> containing_elem_indices = mpStaticCastCellPopulation->GetNode(nodeIndex)->rGetContainingElementIndices();
//
//    // Iterate over these elements
//    for (std::set<unsigned>::iterator elem_iter = containing_elem_indices.begin();
//         elem_iter != containing_elem_indices.end();
//         ++elem_iter)
//    {
//        // Get all the nodes contained in this element
//        // Don't want to include the current node
//        unsigned neighbour_global_index;
//
//        for (unsigned local_index=0; local_index<4; local_index++)
//        {
//            neighbour_global_index = mpStaticCastCellPopulation->rGetMesh().GetElement(*elem_iter)->GetNodeGlobalIndex(local_index);
//
//            if( (neighbour_global_index != nodeIndex) && (!mpStaticCastCellPopulation->IsGhostNode(neighbour_global_index)) )
//            {
//            	neighbouring_node_indices.insert(neighbour_global_index);
//            }
//        }
//    }
//
//    return neighbouring_node_indices;
}

bool CryptSimulation3dGhosts::GetInstanceOfNonExistentNode()
{
	return mFoundNonExistentNode;
}

std::vector<c_vector<unsigned, 2> > CryptSimulation3dGhosts::GetNeighbouringEpithelialPairs(AbstractCellPopulation<3>& rCellPopulation)
{
    DomMeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*>(&rCellPopulation);

    // Create a vector to record the pairs of nodes corresponding to *joined* epithelial nodes
    std::vector<c_vector<unsigned, 2> > node_pairs;
    c_vector<double, 2> pair;

    // We iterate over all cells in the tissue, and deal only with those that are epithelial cells
    for (AbstractCellPopulation<3>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
    	boost::shared_ptr<AbstractCellMutationState> p_state = cell_iter->GetMutationState();

		// Need these to not be stromal cells (and not dead)
		if ( (p_state->IsType<StromalCellMutationState>()==false) && (!cell_iter->IsDead()) )	// an epithelial cell
		{
			Node<3>* p_node = p_tissue->GetNodeCorrespondingToCell(*cell_iter);	// Pointer to node
			unsigned node_index = p_node->GetIndex();

			assert(!(p_tissue->IsGhostNode(node_index)));  // bit unnecessary at this stage but paranoia demands it

			std::vector<unsigned> neighbouring_epithelial_nodes;

			for (Node<3>::ContainingElementIterator iter = p_node->ContainingElementsBegin();
				 iter != p_node->ContainingElementsEnd();
				 ++iter)
			{
				// Get a pointer to the element
				Element<3,3>* p_element = p_tissue->rGetMesh().GetElement(*iter);
                    
				// Iterate over nodes owned by this element
                for (unsigned local_index=0; local_index<4; local_index++)
                {
                    unsigned nodeBGlobalIndex = p_element->GetNodeGlobalIndex(local_index);

					if(!(p_tissue->IsGhostNode(nodeBGlobalIndex)))
					{
						CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

						if (p_cell->GetMutationState()->IsType<StromalCellMutationState>()==false)
						{
							// Store the index of each tissue node that is attached to the epithelial
							// node. There will be repetitions due to iterating over neighbouring elements
							if (nodeBGlobalIndex != node_index)
							{
								neighbouring_epithelial_nodes.push_back(nodeBGlobalIndex);
							}
						}
					}
                }
			}

			// Remove any nodes that have been found twice
			RemoveDuplicates1D(neighbouring_epithelial_nodes);
	
			// Now construct the vector of node pairs
			for (unsigned i=0; i<neighbouring_epithelial_nodes.size(); i++)
			{
				pair[0] = node_index;
				pair[1] = neighbouring_epithelial_nodes[i];
				node_pairs.push_back(pair);
			}
		}
    }
		
    return node_pairs;
}

void CryptSimulation3dGhosts::WriteVisualizerSetupFile()
{
    *mpVizSetupFile << "MeshWidth\t" << mpStaticCastCellPopulation->rGetMesh().GetWidth(0u) << "\n";// get farthest distance between nodes in the x-direction
}

void CryptSimulation3dGhosts::SetupWriteEpithelialCoordinateData()
{
    OutputFileHandler output_file_handler(this->mSimulationOutputDirectory + "/", false);
    mEpithelialCoordinateDataResultsFile = output_file_handler.OpenOutputFile("epithelial_coordinates.dat");
    *mpVizSetupFile << "EpithelialCoordinates\n";
}

/* Outputs: time - average_z_difference - x - y - z - x - y - z - .... */
void CryptSimulation3dGhosts::WriteEpithelialCoordinateData(double time)
{
    DomMeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*>(&mrCellPopulation);
	
	// Going to hijack this to output some data I want
//	
//	PRINT_VARIABLE(mrCellPopulation.GetNumNodes());
//	
//	c_vector<double, 3> cell_location_1, cell_location_2;
//	unsigned node_index_1, node_index_2;
//	double small_radius = 0.1;
//	
//	// Loop over all cells and find those that are really, really close....
//    for (AbstractCellPopulation<3>::Iterator cell_iter = mrCellPopulation.Begin();
//         cell_iter != mrCellPopulation.End();
//         ++cell_iter)
//    {   
//    	cell_location_1 = p_tissue->GetNodeCorrespondingToCell(*cell_iter)->rGetLocation();    	
//    	node_index_1 = p_tissue->GetNodeCorrespondingToCell(*cell_iter)->GetIndex();
//        
//    	for (AbstractCellPopulation<3>::Iterator cell_iter = mrCellPopulation.Begin();
//             cell_iter != mrCellPopulation.End();
//             ++cell_iter)
//        {
//    		cell_location_2 = p_tissue->GetNodeCorrespondingToCell(*cell_iter)->rGetLocation();    		
//    		node_index_2 = p_tissue->GetNodeCorrespondingToCell(*cell_iter)->GetIndex();
//    		
//    		if ( (node_index_1 != node_index_2) && 
//    				( (cell_location_2[0] - cell_location_1[0])*(cell_location_2[0] - cell_location_1[0]) 
//    						+  (cell_location_2[1] - cell_location_1[1])*(cell_location_2[1] - cell_location_1[1])
//    						+  (cell_location_2[2] - cell_location_1[2])*(cell_location_2[2] - cell_location_1[2]) < small_radius) )
//    		{
//    			TRACE("Two very close nodes:");
//    			PRINT_2_VARIABLES(node_index_1, node_index_2);
//    		}    			
//        }   
//    }
//	
	
	/////////// the actual method
    
    
    *mEpithelialCoordinateDataResultsFile <<  time << "\t";

    // Firstly we want to output alpha at each timestep
    
    // Get neighbouring epithelial pairs (keep a record of the number of pairs)
        
    std::vector<c_vector<unsigned, 2> > epithelial_neighbours = GetNeighbouringEpithelialPairs(mrCellPopulation);
    
    unsigned num_pairs = epithelial_neighbours.size();
    
    unsigned node_index_A = 0;
    unsigned node_index_B = 0;    
    double alpha = 0.0;    
    c_vector<double, 3> location_A, location_B;   
    double s = 0.0;
    
    for(unsigned i=0; i<epithelial_neighbours.size(); i++)
    {
    	node_index_A = epithelial_neighbours[i][0];
    	node_index_B = epithelial_neighbours[i][1];
    	
    	assert(node_index_A != node_index_B);

    	location_A = mrCellPopulation.GetNode(node_index_A)->rGetLocation();
    	location_B = mrCellPopulation.GetNode(node_index_B)->rGetLocation();
    	
    	// Calculate s_AB
    	s = pow( (location_A[0] - location_B[0])*(location_A[0] - location_B[0]) + (location_A[1] - location_B[1])*(location_A[1] - location_B[1]), 0.5);
    	
    	if (s==0.0)
    	{
    		alpha = 0.0;
    	}
    	else
    	{
    		alpha += fabs(location_A[2] - location_B[2])/s;
    	}    	    	
    }
    
    alpha /= num_pairs;
    
    *mEpithelialCoordinateDataResultsFile << alpha << "\t";

    // Now output the coordinates of the epithelial cells: x, y, z
    for (AbstractCellPopulation<3>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {    	        	
        if (cell_iter->GetMutationState()->IsType<StromalCellMutationState>()==false)
        {
        	*mEpithelialCoordinateDataResultsFile << p_tissue->GetNodeCorrespondingToCell(*cell_iter)->rGetLocation()[0] << "\t"
        	<< p_tissue->GetNodeCorrespondingToCell(*cell_iter)->rGetLocation()[1] << "\t" 
        	<< p_tissue->GetNodeCorrespondingToCell(*cell_iter)->rGetLocation()[2] << "\t";
        }
    }
    
    *mEpithelialCoordinateDataResultsFile << "\n";
}


std::string CryptSimulation3dGhosts::GetDataOutputFile()
{
	return mDataOutputFile;
}

void CryptSimulation3dGhosts::SetupSolve()
{
	OffLatticeSimulation<3>::SetupSolve();
	
    // Sets up output file - name differently
	OutputFileHandler output_file_handler(this->mSimulationOutputDirectory + "/", false);
	mNumCellsResultsFile = output_file_handler.OpenOutputFile("results.numcells");
    *mpVizSetupFile << "NumCells\n";

    SetupWriteEpithelialCoordinateData();
    double current_time = SimulationTime::Instance()->GetTime();
    WriteEpithelialCoordinateData(current_time);
}

void CryptSimulation3dGhosts::PostSolve()
{
    SimulationTime* p_time = SimulationTime::Instance();
    *mNumCellsResultsFile << p_time->GetTime() << " " << mrCellPopulation.GetNumRealCells() << "\n";

    if ((p_time->GetTimeStepsElapsed()+1)%mSamplingTimestepMultiple==0)
    {
        double time_next_step = p_time->GetTime() + p_time->GetTimeStep();
        WriteEpithelialCoordinateData(time_next_step);
    }
}

/*
void CryptSimulation3dGhosts::AfterSolve()
{	
	OffLatticeSimulation<3>::AfterSolve();

	if ( mrCellPopulation.Begin() != mrCellPopulation.End() )  // there are any cells
    {
    	mNumCellsResultsFile->close();
    }
    
    mEpithelialCoordinateDataResultsFile->close();
}
*/

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CryptSimulation3dGhosts)

