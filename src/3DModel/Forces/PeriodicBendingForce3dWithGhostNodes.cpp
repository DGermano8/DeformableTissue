#include "PeriodicBendingForce3dWithGhostNodes.hpp"
#include "SimulationTime.hpp"
#include "DomMeshBasedCellPopulationWithGhostNodes.hpp"

PeriodicBendingForce3dWithGhostNodes::PeriodicBendingForce3dWithGhostNodes()
    : AbstractForce<3>(),
      mBasementMembraneParameter(DOUBLE_UNSET),
	  mExponentParameter(DOUBLE_UNSET),
      mPeriodicDomainWidth(DOUBLE_UNSET),
      mPeriodicDomainDepth(DOUBLE_UNSET),
      mpExtendedMesh(NULL),
      mSetNonZeroTargetCurvatureRegion(false),
      mNonZeroTargetCurvature(DOUBLE_UNSET),
      mRadiusOfNonZeroTargetCurvatureRegion(DOUBLE_UNSET),
      mXCoordinateOfCentreOfNonZeroTargetCurvatureRegion(DOUBLE_UNSET),
      mYCoordinateOfCentreOfNonZeroTargetCurvatureRegion(DOUBLE_UNSET)
{
}

PeriodicBendingForce3dWithGhostNodes::~PeriodicBendingForce3dWithGhostNodes()
{
    // Avoid memory leaks
    if (mpExtendedMesh != NULL)
    {
        delete mpExtendedMesh;
    }
    
}

void PeriodicBendingForce3dWithGhostNodes::SetBasementMembraneParameter(double basementMembraneParameter)
{
	mBasementMembraneParameter = basementMembraneParameter;
}

double PeriodicBendingForce3dWithGhostNodes::GetBasementMembraneParameter()
{
	return mBasementMembraneParameter;
}

void PeriodicBendingForce3dWithGhostNodes::SetExponentParameter(double exponentParameter)
{
	mExponentParameter = exponentParameter;
}

double PeriodicBendingForce3dWithGhostNodes::GetExponentParameter()
{
	return mExponentParameter;
}


void PeriodicBendingForce3dWithGhostNodes::SetCircularNonZeroTargetCurvatureRegion(bool setNonZeroTargetCurvatureRegion, double nonZeroTargetCurvature,
		double radiusOfNonZeroTargetCurvatureRegion, double xCoordinateOfCentreOfNonZeroTargetCurvatureRegion,
		double yCoordinateOfCentreOfNonZeroTargetCurvatureRegion)
{
	mSetNonZeroTargetCurvatureRegion = setNonZeroTargetCurvatureRegion;
	mNonZeroTargetCurvature = nonZeroTargetCurvature;
	mRadiusOfNonZeroTargetCurvatureRegion = radiusOfNonZeroTargetCurvatureRegion;
	mXCoordinateOfCentreOfNonZeroTargetCurvatureRegion = xCoordinateOfCentreOfNonZeroTargetCurvatureRegion;
	mYCoordinateOfCentreOfNonZeroTargetCurvatureRegion = yCoordinateOfCentreOfNonZeroTargetCurvatureRegion;
}

bool PeriodicBendingForce3dWithGhostNodes::GetWhetherToSetNonZeroTargetCurvatureRegion()
{
	return mSetNonZeroTargetCurvatureRegion;
}

double PeriodicBendingForce3dWithGhostNodes::GetNonZeroTargetCurvature()
{
	return mNonZeroTargetCurvature;
}

double PeriodicBendingForce3dWithGhostNodes::GetRadiusOfNonZeroTargetCurvatureRegion()
{
	return mRadiusOfNonZeroTargetCurvatureRegion;
}

double PeriodicBendingForce3dWithGhostNodes::GetXCoordinateOfCentreOfNonZeroTargetCurvatureRegion()
{
	return mXCoordinateOfCentreOfNonZeroTargetCurvatureRegion;
}

double PeriodicBendingForce3dWithGhostNodes::GetYCoordinateOfCentreOfNonZeroTargetCurvatureRegion()
{
	return mYCoordinateOfCentreOfNonZeroTargetCurvatureRegion;
}

/*
 * Remove repeated entries in a 1D vector
 */
void PeriodicBendingForce3dWithGhostNodes::RemoveDuplicates1D(std::vector<unsigned>& rVectorWithDuplicates)
{
    std::sort(rVectorWithDuplicates.begin(), rVectorWithDuplicates.end());
    rVectorWithDuplicates.erase(std::unique(rVectorWithDuplicates.begin(), rVectorWithDuplicates.end()), rVectorWithDuplicates.end());
}

// Dom - This is used to find connecting pairs of cells - could modify it to find the neighbouring nodes. Otherwisem, won't need it at all
/*
 * A method to find all the pairs of connections between healthy epithelial cells and labelled tissue cells.
 * Returns a vector of node pairings, without repeats. The first of each pair is the epithelial node index,
 * and the second is the tissue/stromal node index.
 */
std::vector<c_vector<unsigned, 2> > PeriodicBendingForce3dWithGhostNodes::GetEpithelialTissuePairs(AbstractCellPopulation<3>& rCellPopulation)
{
    // Create a vector to record the pairs of nodes corresponding to *joined* epithelial and tissue nodes
    std::vector<c_vector<unsigned, 2> > node_pairs;
    c_vector<double, 2> pair;

    // Loop over nodes of mpExtendedMesh
    for (unsigned extended_node_index=0; extended_node_index<mpExtendedMesh->GetNumNodes(); extended_node_index++)
    {
        // Get a pointer to this node in mpExtendedMesh
        Node<3>* p_node = mpExtendedMesh->GetNode(extended_node_index);

        // Get the corresponding node index in rCellPopulation
        unsigned node_index = mExtendedMeshNodeIndexMap[extended_node_index];

        // Get the cell corresponding to this node
        CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(node_index);

        // Get mutation state of cell
    	boost::shared_ptr<AbstractCellMutationState> p_state = p_cell->GetMutationState();

    	// Get whether cell is dead
    	bool cell_is_dead = p_cell->IsDead();

    	// Get whether this cell is a live epithelial cell
//    	bool is_live_epithelial_cell = (p_state->IsType<StromalCellMutationState>()==false) && !cell_is_dead;
    	bool is_live_epithelial_cell = (p_state->IsType<WildTypeCellMutationState>()==true) && !cell_is_dead;

    	if (is_live_epithelial_cell)
    	{
    	    // Iterate over elements of mpExtendedMesh containing this node and no ghost nodes
    		std::vector<unsigned> tissue_nodes;

    		for (Node<3>::ContainingElementIterator elem_iter = p_node->ContainingElementsBegin();
		         elem_iter != p_node->ContainingElementsEnd();
		         ++elem_iter)
    		{
				// Get a pointer to the element
				Element<3,3>* p_element = mpExtendedMesh->GetElement(*elem_iter);

				// ITERATE OVER NODES owned by this element
				for (unsigned local_index=0; local_index<4; local_index++)
				{
					unsigned nodeBGlobalIndex = p_element->GetNodeGlobalIndex(local_index);

					// Get the corresponding node index in rCellPopulation
					unsigned neighbour_index = mExtendedMeshNodeIndexMap[nodeBGlobalIndex];
					CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(neighbour_index);

					if (p_cell->GetMutationState()->IsType<StromalCellMutationState>())
					{
						// Store the index of each tissue node that is attached to the epithelial
						// node. There will be repetitions due to iterating over neighbouring elements
						tissue_nodes.push_back(nodeBGlobalIndex);
					}
				}
			}

			// Remove any nodes that have been found twice
			RemoveDuplicates1D(tissue_nodes);

			// Now construct the vector of node pairs
			for (unsigned i=0; i<tissue_nodes.size(); i++)
			{
				pair[0] = extended_node_index;
				pair[1] = tissue_nodes[i];
				node_pairs.push_back(pair);

				///\todo Consider reimplementing this check
//				// Check that these node share a common element
//				bool has_common_element = false;
//
//				// The elements that contain this epithelial node:
//				std::set<unsigned> epithelial_elements = rCellPopulation.GetNode(node_index)->rGetContainingElementIndices();
//				assert(epithelial_elements.size() != 0);
//
//				// The elements that contain the tissue node:
//				std::set<unsigned> tissue_elements = rCellPopulation.GetNode(tissue_nodes[i])->rGetContainingElementIndices();
//				assert(tissue_elements.size() != 0);
//
//				// Loop over all elements that contain the tissue node
//				for (Node<3>::ContainingElementIterator elt_it = rCellPopulation.GetNode(tissue_nodes[i])->ContainingElementsBegin();
//					 elt_it != rCellPopulation.GetNode(tissue_nodes[i])->ContainingElementsEnd();
//					 ++elt_it)
//				{
//					unsigned elt_index = *elt_it;
//
//					bool elt_contains_ghost_nodes = DoesElementContainGhostNodes(rCellPopulation, elt_index);
//
//					// Keep only those elements that also contain the epithelial node, but do not have ghost nodes
//					if ( (elt_contains_ghost_nodes == false) && (epithelial_elements.find(elt_index) != epithelial_elements.end()) )
//					{
//						// Common element
//						has_common_element = true;
//						break;
//					}
//				}
//
//				if (!has_common_element)
//				{
//					PRINT_2_VARIABLES(node_index,tissue_nodes[i]);
//				}
//				assert(has_common_element);
			}
		}
    }

    return node_pairs;
	
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<c_vector<unsigned, 10> > PeriodicBendingForce3dWithGhostNodes::GetEpithelialNeighbours(std::vector<c_vector<unsigned, 3> > rEpithelialMeshVector, int number_of_cells)
{
	std::vector<c_vector<unsigned, 10> > node_neighbours;

	int temp_neighbours[4*number_of_cells][10];
	for(unsigned i=0; i<4*number_of_cells; i++)
	{
		for(unsigned j=0; j<10; j++)
		{
			temp_neighbours[i][j] = 0;
		}
	}

	for(unsigned i=0; i<rEpithelialMeshVector.size(); i++)
	{
		unsigned node_A = rEpithelialMeshVector[i][0];
		unsigned node_B = rEpithelialMeshVector[i][1];
		unsigned node_C = rEpithelialMeshVector[i][2];

		// Do node_A first
		int row_neighs_A[10];
		for(unsigned j=0; j<10; j++)
		{
			row_neighs_A[j] = temp_neighbours[node_A][j];
		}
		if(std::find(std::begin(row_neighs_A), std::end(row_neighs_A), node_B) == std::end(row_neighs_A))
		{
			int iter_A = 0;
			while(iter_A < 10)
			{
				if(temp_neighbours[node_A][iter_A] == 0)
				{
					temp_neighbours[node_A][iter_A] = node_B;
					iter_A = 10;
				}
				iter_A++;
			}
		}
		for(unsigned j=0; j<10; j++)
		{
			row_neighs_A[j] = temp_neighbours[node_A][j];
		}
		if(std::find(std::begin(row_neighs_A), std::end(row_neighs_A), node_C) == std::end(row_neighs_A))
		{
			int iter_A = 0;
			while(iter_A < 10)
			{
				if(temp_neighbours[node_A][iter_A] == 0)
				{
					temp_neighbours[node_A][iter_A] = node_C;
					iter_A = 10;
				}
				iter_A++;
			}
		}



		// Then do node_B
		int row_neighs_B[10];
		for(unsigned j=0; j<10; j++)
		{
			row_neighs_B[j] = temp_neighbours[node_B][j];
		}
		if(std::find(std::begin(row_neighs_B), std::end(row_neighs_B), node_A) == std::end(row_neighs_B))
		{
			int iter_B = 0;
			while(iter_B < 10)
			{
				if(temp_neighbours[node_B][iter_B] == 0)
				{
					temp_neighbours[node_B][iter_B] = node_A;
					iter_B = 10;
				}
				iter_B++;
			}
		}
		for(unsigned j=0; j<10; j++)
		{
			row_neighs_B[j] = temp_neighbours[node_B][j];
		}
		if(std::find(std::begin(row_neighs_B), std::end(row_neighs_B), node_C) == std::end(row_neighs_B))
		{
			int iter_B = 0;
			while(iter_B < 10)
			{
				if(temp_neighbours[node_B][iter_B] == 0)
				{
					temp_neighbours[node_B][iter_B] = node_C;
					iter_B = 10;
				}
				iter_B++;
			}
		}

		// Then do node_C
		int row_neighs_C[10];
		for(unsigned j=0; j<10; j++)
		{
			row_neighs_C[j] = temp_neighbours[node_C][j];
		}
		if(std::find(std::begin(row_neighs_C), std::end(row_neighs_C), node_A) == std::end(row_neighs_C))
		{
			int iter_B = 0;
			while(iter_B < 10)
			{
				if(temp_neighbours[node_C][iter_B] == 0)
				{
					temp_neighbours[node_C][iter_B] = node_A;
					iter_B = 10;
				}
				iter_B++;
			}
		}
		for(unsigned j=0; j<10; j++)
		{
			row_neighs_C[j] = temp_neighbours[node_C][j];
		}
		if(std::find(std::begin(row_neighs_C), std::end(row_neighs_C), node_B) == std::end(row_neighs_C))
		{
			int iter_B = 0;
			while(iter_B < 10)
			{
				if(temp_neighbours[node_C][iter_B] == 0)
				{
					temp_neighbours[node_C][iter_B] = node_B;
					iter_B = 10;
				}
				iter_B++;
			}
		}

	}
	
	for(unsigned i=0; i<4*number_of_cells; i++)
	{
		std::vector<unsigned> holder_neighbour;
		c_vector<unsigned, 10> another_holder_neighbour;
		//std::cout << i << "    ";
		for(unsigned j=0; j<10; j++)
		{
			holder_neighbour.push_back(temp_neighbours[i][j]);

		//	std::cout<< temp_neighbours[i][j] << " ";
		}
		//std::cout << "\n";
		
		another_holder_neighbour[0] = holder_neighbour[0];
		another_holder_neighbour[1] = holder_neighbour[1];
		another_holder_neighbour[2] = holder_neighbour[2];
		another_holder_neighbour[3] = holder_neighbour[3];
		another_holder_neighbour[4] = holder_neighbour[4];
		another_holder_neighbour[5] = holder_neighbour[5];
		another_holder_neighbour[6] = holder_neighbour[6];
		another_holder_neighbour[7] = holder_neighbour[7];
		another_holder_neighbour[8] = holder_neighbour[8];
		another_holder_neighbour[9] = holder_neighbour[9];

		node_neighbours.push_back(another_holder_neighbour);
	}
	return node_neighbours;

}



std::vector<c_vector<unsigned, 3> > PeriodicBendingForce3dWithGhostNodes::GetEpithelialMesh(AbstractCellPopulation<3>& rCellPopulation)
{
// Get a pointer to this node in mpExtendedMesh
std::vector<c_vector<unsigned, 3> > epithelial_triangulation;

for (unsigned extended_node_index=0; extended_node_index<mpExtendedMesh->GetNumNodes(); extended_node_index++)
    {
        Node<3>* p_node = mpExtendedMesh->GetNode(extended_node_index);

        // Get the corresponding node index in rCellPopulation
        unsigned node_index = mExtendedMeshNodeIndexMap[extended_node_index];

        // Get the cell corresponding to this node
        CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(node_index);

        // Get mutation state of cell
    	boost::shared_ptr<AbstractCellMutationState> p_state = p_cell->GetMutationState();

    	// Get whether cell is dead
    	bool cell_is_dead = p_cell->IsDead();

    	// Get whether this cell is a live epithelial cell
    	bool is_live_epithelial_cell = (p_state->IsType<WildTypeCellMutationState>()==true) && !cell_is_dead;
		
		DomMeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*>(&rCellPopulation);

    	if (is_live_epithelial_cell)
    	{
    	    // Iterate over elements of mpExtendedMesh containing this node and no ghost nodes
    		
			c_vector<unsigned, 3> tri_el;

			// This needs both stromal and epithelial cells... so cant do monolayer..
    		for (Node<3>::ContainingElementIterator elem_iter = p_node->ContainingElementsBegin();
		         elem_iter != p_node->ContainingElementsEnd();
		         ++elem_iter)
    		{
				// Get a pointer to the element
				Element<3,3>* p_element = mpExtendedMesh->GetElement(*elem_iter);
				
				// ITERATE OVER NODES owned by this element
				std::vector<unsigned> temp_triangular_element;
				bool is_element_connected_to_tissue = false;
				for (unsigned local_index=0; local_index<4; local_index++)
				{
					unsigned nodeBGlobalIndex = p_element->GetNodeGlobalIndex(local_index);

					// Get the corresponding node index in rCellPopulation
					unsigned neighbour_index = mExtendedMeshNodeIndexMap[nodeBGlobalIndex];
					CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(neighbour_index);

					//if(neighbour_index != nodeBGlobalIndex)
					//{
					//	std::cout<< "neighbour_index = " << neighbour_index << "  nodeBGlobalIndex = " << nodeBGlobalIndex << "\n";
					//}


					if (p_cell->GetMutationState()->IsType<WildTypeCellMutationState>())
					{
						//std::cout<< "nodeBGlobalIndex = " << nodeBGlobalIndex << " neighbour_index = " << neighbour_index << "\n";
						// nodeBGlobalIndex -> from extended mesh
						// neighbour_index -> mapped nodes from original mesh
						temp_triangular_element.push_back(nodeBGlobalIndex);
						//temp_triangular_element.push_back(neighbour_index);
					}
					else if(p_cell->GetMutationState()->IsType<StromalCellMutationState>() || p_tissue->IsGhostNode(nodeBGlobalIndex))
					{
						is_element_connected_to_tissue = true;
					}

				}

				// This means, we only consider adding if theres 3 epithelial cells AND 1 stromal cell
				if(temp_triangular_element.size() == 3 && is_element_connected_to_tissue)
				{
					std::sort(temp_triangular_element.begin(), temp_triangular_element.end());
					tri_el[0] = temp_triangular_element[0];
					tri_el[1] = temp_triangular_element[1];
					tri_el[2] = temp_triangular_element[2];

					// Check for duplicates
					bool put_element_in_tri = 1;
					for(unsigned i=0; i<epithelial_triangulation.size(); i++)
					{
						if(epithelial_triangulation[i][0] == tri_el[0])
						{
							if(epithelial_triangulation[i][1] == tri_el[1])
							{
								if(epithelial_triangulation[i][2] == tri_el[2])
								{
									put_element_in_tri = 0;
								}
							}
						}
					}

					if(put_element_in_tri == 1)
					{
						epithelial_triangulation.push_back(tri_el);
					}

				}
			}
			
		}
    }
	return epithelial_triangulation;
}

std::vector<c_vector<double, 3> > PeriodicBendingForce3dWithGhostNodes::FitPlaneAndFindImage(AbstractCellPopulation<3>& rCellPopulation, std::vector<unsigned> second_order_neighs,  unsigned cell_i)
{
	c_vector<double, 3> normal_vector;
	
	std::vector<c_vector<double, 3> > image_location_per_second_order_neighbours;

	//	a00, a10, a20, a01, a11, a21, a02, a12, a22
	c_vector<double, 9> ATA;
	c_vector<double, 3> ATz;
	c_vector<double, 9> iATA;

	for (unsigned i=0; i<second_order_neighs.size(); i++)
	{
	//	unsigned cell_i_ext = mExtendedMeshNodeIndexMap[second_order_neighs[i]];
		unsigned cell_i_ext =second_order_neighs[i];
		//int epithelialNodeIndex = second_order_neighs[i];
		c_vector<double, 3> epithelial_location = mpExtendedMesh->GetNode(cell_i_ext)->rGetLocation();

		// Find Transpose(A)*A
		ATA[0] += epithelial_location[0]*epithelial_location[0];
		ATA[1] += epithelial_location[0]*epithelial_location[1];
		ATA[2] += epithelial_location[0];

		ATA[3] += epithelial_location[0]*epithelial_location[1];
		ATA[4] += epithelial_location[1]*epithelial_location[1];
		ATA[5] += epithelial_location[1];

		ATA[6] += epithelial_location[0];
		ATA[7] += epithelial_location[1];
		ATA[8] += 1.0;

		// Calculate Transpose(A)*z
		ATz[0] += epithelial_location[0]*epithelial_location[2];
		ATz[1] += epithelial_location[1]*epithelial_location[2];
		ATz[2] += epithelial_location[2];
	}

	// determinant of Transpose(A)*A
	double ATA_det = ATA[0]*(ATA[4]*ATA[8]-ATA[5]*ATA[7]) - ATA[1]*(ATA[3]*ATA[8]-ATA[5]*ATA[6]) + ATA[2]*(ATA[3]*ATA[7]-ATA[6]*ATA[4]);

	if (ATA_det != 0)
	{
		// Calculate the inverse of Transpose(A)*A
		iATA[0] = (1.0/ATA_det)*(ATA[4]*ATA[8] - ATA[7]*ATA[5]);
		iATA[1] = (1.0/ATA_det)*(ATA[2]*ATA[7] - ATA[8]*ATA[1]);
		iATA[2] = (1.0/ATA_det)*(ATA[1]*ATA[5] - ATA[4]*ATA[2]);

		iATA[3] = (1.0/ATA_det)*(ATA[5]*ATA[6] - ATA[8]*ATA[3]);
		iATA[4] = (1.0/ATA_det)*(ATA[0]*ATA[8] - ATA[6]*ATA[2]);
		iATA[5] = (1.0/ATA_det)*(ATA[2]*ATA[3] - ATA[5]*ATA[0]);

		iATA[6] = (1.0/ATA_det)*(ATA[3]*ATA[7] - ATA[6]*ATA[4]);
		iATA[7] = (1.0/ATA_det)*(ATA[1]*ATA[6] - ATA[7]*ATA[0]);
		iATA[8] = (1.0/ATA_det)*(ATA[0]*ATA[4] - ATA[3]*ATA[1]);
		
		// Calculate  normal = inverse(Transpose(A)*A)*(Transpose(A)*z)
		normal_vector[0] = iATA[0]*ATz[0] + iATA[1]*ATz[1] + iATA[2]*ATz[2];
		normal_vector[1] = iATA[3]*ATz[0] + iATA[4]*ATz[1] + iATA[5]*ATz[2];
		normal_vector[2] = iATA[6]*ATz[0] + iATA[7]*ATz[1] + iATA[8]*ATz[2];
		//std::cout << "case 1" << "\n"; 

		normal_vector[0] = -1.0*normal_vector[0];
		normal_vector[1] = -1.0*normal_vector[1];
	}
	else if (ATA_det == 0)
	{
		normal_vector[0] = 0;
		normal_vector[1] = 0;
		normal_vector[2] = 1;
		//std::cout << "case 2" << "\n"; 
	}

	//PRINT_VECTOR(normal_vector);

	double sign_of_curvature = (mNonZeroTargetCurvature > 0) - (mNonZeroTargetCurvature < 0);

	c_vector<double, 3> unit_normal;
	unit_normal[0] = sign_of_curvature*normal_vector[0]/(sqrt(pow(normal_vector[0],2)+pow(normal_vector[1],2)+1));
	unit_normal[1] = sign_of_curvature*normal_vector[1]/(sqrt(pow(normal_vector[0],2)+pow(normal_vector[1],2)+1));
	unit_normal[2] = sign_of_curvature*1.0/(sqrt(pow(normal_vector[0],2)+pow(normal_vector[1],2)+1));

	double radius_from_curvature = 1/sqrt(mNonZeroTargetCurvature*sign_of_curvature);

	c_vector<double, 3> circle_centre;
	unsigned cell_ii_ext = mExtendedMeshNodeIndexMap[cell_i];
	c_vector<double, 3> epithelial_cell_i = mpExtendedMesh->GetNode(cell_ii_ext)->rGetLocation();
	circle_centre[0]  = radius_from_curvature*unit_normal[0] + epithelial_cell_i[0];
	circle_centre[1]  = radius_from_curvature*unit_normal[1] + epithelial_cell_i[1];
	circle_centre[2]  = radius_from_curvature*unit_normal[2] + epithelial_cell_i[2];

	//  Note: when we use unit_normal in this loop, we may need to be using normal_vector instead...?
	for (unsigned i=0; i<second_order_neighs.size(); i++)
	{
		//unsigned cell_i_ext = mExtendedMeshNodeIndexMap[second_order_neighs[i]];
		unsigned cell_i_ext = second_order_neighs[i];
		//int epithelialNodeIndex = second_order_neighs[i];
		c_vector<double, 3> epithelial_location = mpExtendedMesh->GetNode(cell_i_ext)->rGetLocation();
		
		//PRINT_VECTOR(unit_normal);
		//std::cout << pow(unit_normal[0],2) + pow(unit_normal[1],2) + pow(unit_normal[2],2) << "\n";
		double alpha = pow(normal_vector[0],2) + pow(normal_vector[1],2) + 1.0;
		double beta  = 2*normal_vector[0]*(epithelial_location[0]-circle_centre[0]) + 2*normal_vector[1]*(epithelial_location[1]-circle_centre[1]) + 2*1.0*(epithelial_location[2]-circle_centre[2]);
		double delta = pow(epithelial_location[0]-circle_centre[0],2) + pow(epithelial_location[1]-circle_centre[1],2) + pow(epithelial_location[2]-circle_centre[2],2) - pow(radius_from_curvature,2);

		double param_T = (-1.0*beta - sign_of_curvature*sqrt(pow(beta,2) - 4*alpha*delta))/(2.0*alpha);
		c_vector<double, 3> point_on_sphere;
		
		if (isnan(param_T) || pow(beta,2) - 4*alpha*delta < 0 )
		{
			//param_T = 0.0;
			// i.e. don't move the cell
			point_on_sphere[0] = epithelial_location[0];
			point_on_sphere[1] = epithelial_location[1];
			point_on_sphere[2] = epithelial_location[2];
		}
		else
		{
			point_on_sphere[0] = normal_vector[0]*param_T + epithelial_location[0];
			point_on_sphere[1] = normal_vector[1]*param_T + epithelial_location[1];
			point_on_sphere[2] = 1.0*param_T + epithelial_location[2];
		}
		
		
		

		double k_point = (normal_vector[0]*point_on_sphere[0] + normal_vector[1]*point_on_sphere[1] + 1.0*point_on_sphere[2] - normal_vector[2])/alpha;
		
		c_vector<double, 3> point_on_plane;
		point_on_plane[0] = point_on_sphere[0] - k_point*normal_vector[0];
		point_on_plane[1] = point_on_sphere[1] - k_point*normal_vector[1];
		point_on_plane[2] = point_on_sphere[2] - k_point*1.0;

		double point_on_plane_DOT_point_on_sphere = sqrt(pow(point_on_plane[0]-point_on_sphere[0],2) + pow(point_on_plane[1]-point_on_sphere[1],2) + pow(point_on_plane[2]-point_on_sphere[2],2));

		c_vector<double, 3> temp_cell_location;
		temp_cell_location[0] = epithelial_location[0] - unit_normal[0]*point_on_plane_DOT_point_on_sphere;
		temp_cell_location[1] = epithelial_location[1] - unit_normal[1]*point_on_plane_DOT_point_on_sphere;
		temp_cell_location[2] = epithelial_location[2] - unit_normal[2]*point_on_plane_DOT_point_on_sphere;		
		
		//c_vector<double, 3> temp_cell_location;
		//temp_cell_location[0] = epithelial_location[0] - unit_normal[0]*point_on_plane_DOT_point_on_sphere;
		//temp_cell_location[1] = epithelial_location[1] - unit_normal[1]*point_on_plane_DOT_point_on_sphere;
		//temp_cell_location[2] = epithelial_location[2] - unit_normal[2]*point_on_plane_DOT_point_on_sphere;
		
		

		//std::cout << "T = "<< param_T << " desc = " << pow(beta,2) - 4*alpha*delta << "\n";
		
		//PRINT_VECTOR(epithelial_location);
		//if (isnan(param_T))
		//{
		//	PRINT_2_VARIABLES(param_T, (pow(beta,2) - 4*alpha*delta));
		//	PRINT_VECTOR(normal_vector);
		//}
		
		image_location_per_second_order_neighbours.push_back(temp_cell_location);
	}

	return image_location_per_second_order_neighbours;
}



c_vector<double, 3> PeriodicBendingForce3dWithGhostNodes::GetForceDueToDiscreteCurvature(AbstractCellPopulation<3>& rCellPopulation, std::vector<c_vector<unsigned, 3> > rEpithelialMeshVector, std::vector<unsigned> first_order_neighs,  std::vector<unsigned> second_order_neighs, std::vector<c_vector<double, 3> > image_location, unsigned cell_i, unsigned num_cells)
{
	
	/* 
	Need:
	Traingulation		    ->  rEpithelialMeshVector
	First order neighbours  ->  first_order_neighs 
	image node positions    ->  image_location 
	  needs (for node IDs)  ->  second_order_neighs 
	Cell_i					->  cell_i
	*/
	//std::cout<<"\n\n";
	double exponent_parameter = GetExponentParameter();
	c_vector<double, 3> force_due_to_curvature;
	force_due_to_curvature[0] = 0.0;
	force_due_to_curvature[1] = 0.0;
	force_due_to_curvature[2] = 0.0;
	//PRINT_VECTOR(force_due_to_curvature);

	int inverse_sec_neighs[4*num_cells];
	for(unsigned i=0; i<4*num_cells; i++)
	{
		inverse_sec_neighs[i] = 0;

		//std::cout << "i " << i << " s " << second_order_neighs[i] << "\n";
	}

	for(unsigned i=0; i<second_order_neighs.size(); i++)
	{
		inverse_sec_neighs[second_order_neighs[i]] = i;

		//std::cout << "i " << i << " s " << second_order_neighs[i] << "\n";
	}


	// Calculate angles about the first order neighbours in the mesh
	std::vector<double> angle_sum;
	std::vector<double> angle_sum_indicats;
	// std::cout << "\n\n";
	for(unsigned j=0; j<first_order_neighs.size(); j++)
	{
		double ang_sum = 0.0;
		c_vector<double, 3> grad_i_phi_j_el;
		grad_i_phi_j_el[0] = 0.0;
		grad_i_phi_j_el[1] = 0.0;
		grad_i_phi_j_el[2] = 0.0;

		int cell_a = first_order_neighs[j];
		int inv_cell_a = inverse_sec_neighs[cell_a];
		c_vector<double, 3> cell_a_loc, cell_b_loc, cell_c_loc;
		cell_a_loc[0] = image_location[inv_cell_a][0];
		cell_a_loc[1] = image_location[inv_cell_a][1];
		cell_a_loc[2] = image_location[inv_cell_a][2];

		

		for(unsigned i=0; i<rEpithelialMeshVector.size(); i++)
		{
			int cell_b, cell_c, inv_cell_b, inv_cell_c;
			bool have_triangle = false;

			if(rEpithelialMeshVector[i][0] == cell_a)
			{
				cell_b = rEpithelialMeshVector[i][1];
				cell_c = rEpithelialMeshVector[i][2];
				have_triangle = true;
			}
			else if(rEpithelialMeshVector[i][1] == cell_a)
			{
				cell_b = rEpithelialMeshVector[i][0];
				cell_c = rEpithelialMeshVector[i][2];
				have_triangle = true;
			}
			else if(rEpithelialMeshVector[i][2] == cell_a)
			{
				cell_b = rEpithelialMeshVector[i][0];
				cell_c = rEpithelialMeshVector[i][1];
				have_triangle = true;
			}
			// have_triangle = bool ( (cell_a == cell_i) + (cell_b == cell_i) + (cell_c == cell_i) );

			if(have_triangle)
			{
				inv_cell_b = inverse_sec_neighs[cell_b];
				cell_b_loc[0] = image_location[inv_cell_b][0];
				cell_b_loc[1] = image_location[inv_cell_b][1];
				cell_b_loc[2] = image_location[inv_cell_b][2];

				inv_cell_c = inverse_sec_neighs[cell_c];
				cell_c_loc[0] = image_location[inv_cell_c][0];
				cell_c_loc[1] = image_location[inv_cell_c][1];
				cell_c_loc[2] = image_location[inv_cell_c][2];

				c_vector<double, 3> vect_ab = cell_b_loc - cell_a_loc;
				c_vector<double, 3> vect_ac = cell_c_loc - cell_a_loc;

				double mag_ab = sqrt(pow(vect_ab[0],2) + pow(vect_ab[1],2) + pow(vect_ab[2],2));
				double mag_ac = sqrt(pow(vect_ac[0],2) + pow(vect_ac[1],2) + pow(vect_ac[2],2));

				double cos_abc = (vect_ab[0]*vect_ac[0] + vect_ab[1]*vect_ac[1] + vect_ab[2]*vect_ac[2])/(mag_ab*mag_ac);

				// Think these edge cases, they shouldn't be happening, but they do...
				if (!isnan(ang_sum))
				{
					if (cos_abc >= 0.99999999 || cos_abc == 1)
					{
						//PRINT_VARIABLE("+1");
						ang_sum += 0.0;
						cos_abc = 0.9999;
					}
					else if (cos_abc <= -0.99999999 || cos_abc == -1)
					{
						//PRINT_VARIABLE("-1");
						ang_sum += 3.141592653589793;
						cos_abc = -0.9999;
					}
					else
					{
						ang_sum += acos(cos_abc);
					}
				}
				else if (isnan(ang_sum))
				{
					cos_abc = 0;
				}

				c_vector<double, 3> grad_hold;

				// Calculate grad_i phi_i
				if(cell_a == cell_i)
				{
					// TRACE("one");
					//PRINT_VECTOR(cell_a_loc);
					//PRINT_VECTOR(cell_b_loc);
					//PRINT_VECTOR(cell_c_loc);
					grad_hold = (1.0/(sqrt(1.0 - cos_abc*cos_abc)))*(1.0/(mag_ab*mag_ab*mag_ac*mag_ac))*(mag_ab*(mag_ac-mag_ab*cos_abc)*vect_ac + mag_ac*(mag_ab-mag_ac*cos_abc)*vect_ab);

					grad_i_phi_j_el += grad_hold;
					
					//PRINT_VECTOR(grad_hold);
					//std::cout << "case 1 " << cell_a << "\n";
				}
				
				else if(cell_b == cell_i)
				{
					// TRACE("two");
					//PRINT_VECTOR(cell_a_loc);
					//PRINT_VECTOR(cell_b_loc);
					//PRINT_VECTOR(cell_c_loc);

					grad_hold = (1.0/(sqrt(1.0 - cos_abc*cos_abc))) * (vect_ac/(mag_ab*mag_ac) - cos_abc*vect_ab/(mag_ab*mag_ab));
					//grad_hold = (1.0/(sqrt(1.0 - cos_abc*cos_abc)))* (1.0/(mag_ac)) * (vect_ab/mag_ab - cos_abc*vect_ac/mag_ac);
					grad_i_phi_j_el += grad_hold;
					
					//PRINT_VECTOR(grad_hold);
					//std::cout << "case 2 " << cell_a << "\n";
				}
				else if(cell_c == cell_i)
				{
					// TRACE("three");
					//PRINT_VECTOR(cell_a_loc);
					//PRINT_VECTOR(cell_b_loc);
					//PRINT_VECTOR(cell_c_loc);
					
					grad_hold = (1.0/(sqrt(1.0 - cos_abc*cos_abc))) * (vect_ab/(mag_ab*mag_ac) - cos_abc*vect_ac/(mag_ac*mag_ac));
					//grad_hold = (1.0/(sqrt(1.0 - cos_abc*cos_abc)))* (1.0/(mag_ab)) * (vect_ac/mag_ac - cos_abc*vect_ab/mag_ab);
					grad_i_phi_j_el += grad_hold;

					//PRINT_VECTOR(grad_hold);
					//std::cout << "case 3 " << cell_a << "\n";
				}
				else
				{
					grad_hold[0] = 0.0;
					grad_hold[1] = 0.0;
					grad_hold[2] = 0.0;
					grad_i_phi_j_el += grad_hold;
				}

				// if (isnan(grad_i_phi_j_el[0]) || isnan(grad_i_phi_j_el[1]) || isnan(grad_i_phi_j_el[2]))
				// {
				// 	PRINT_VECTOR(grad_i_phi_j_el);
				// }

				// if ( cell_i == 154 )
				// {
				// 	PRINT_3_VARIABLES(cell_a,cell_b,cell_c);
				// 	PRINT_VECTOR(grad_hold);
				// }
			}
			
		}
		
		//if(cell_a == cell_i)
		//{
		//	PRINT_VECTOR(grad_i_phi_j_el);
		//}
		//std::cout << "Angle Sum = " << ang_sum << "\n";
		//PRINT_VARIABLE(ang_sum);

		angle_sum.push_back(ang_sum);
		angle_sum_indicats.push_back(cell_a);

		//PRINT_VECTOR(force_due_to_curvature);
		// Using alpha = 1.5 and beta = 2.0 for now, will fix this later though
		double sign_of_anlge = (round((2*3.141592653589793 - ang_sum)*pow(10.0,14))/pow(10.0,14) >= 0) - (round((2*3.141592653589793 - ang_sum)*pow(10.0,14))/pow(10.0,14) < 0);
		// if ( cell_i == 154 )
		// {
		// 	double hold = (exponent_parameter)*(sign_of_anlge)*pow(sqrt((2*3.141592653589793 - ang_sum)*(2*3.141592653589793 - ang_sum)),(exponent_parameter-1.0));
		// 	PRINT_2_VARIABLES(sign_of_anlge, hold);
		// }
		c_vector<double, 3> force_due_to_curvature_i = (exponent_parameter)*(sign_of_anlge)*pow(sqrt((2*3.141592653589793 - ang_sum)*(2*3.141592653589793 - ang_sum)),(exponent_parameter-1.0))*grad_i_phi_j_el;

		force_due_to_curvature +=  force_due_to_curvature_i;

		//PRINT_VECTOR(force_due_to_curvature_i);
		
		

		//std::cout << round((2*3.141592653589793 - ang_sum)*pow(10.0,14))/pow(10.0,14) << "\n";
		//std::cout<< "sign_of_anlge = " << sign_of_anlge << "\n";
	}
	

	//PRINT_VECTOR(force_due_to_curvature);

	//std::cout << "cell_i = " << cell_i << "\n"; 
	//PRINT_VECTOR(angle_sum);
	//PRINT_VECTOR(angle_sum_indicats);
	//PRINT_VECTOR(first_order_neighs);
	//PRINT_VECTOR(force_due_to_curvature);
//	std::cout << "\n";
	
	
	
	return force_due_to_curvature;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Dom - looks like some legacy code...
/*
 * Method to determine whether an element contains ghost nodes
 */
bool PeriodicBendingForce3dWithGhostNodes::DoesElementContainGhostNodes(AbstractCellPopulation<3>& rCellPopulation, unsigned elementIndex)
{
	DomMeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*>(&rCellPopulation);

	bool element_contains_ghost_nodes = false;

	// Get a pointer to the element
	Element<3,3>* p_element = p_tissue->rGetMesh().GetElement(elementIndex);

	// ITERATE OVER NODES owned by this element
	for (unsigned local_index=0; local_index<4; local_index++)
	{
		if (p_tissue->IsGhostNode(p_element->GetNodeGlobalIndex(local_index)) == true)
		{
			element_contains_ghost_nodes = true;
		}
	}

	return element_contains_ghost_nodes;
}

// Dom - used within the curvature calculation - probably wont need it later on
/*
 * Method to determine whether an element contains a long edge
 */
bool PeriodicBendingForce3dWithGhostNodes::DoesElementContainLongEdge(AbstractCellPopulation<3>& rCellPopulation, unsigned elementIndex, double maxEdgeLength)
{
	bool element_contains_long_edge = false;

	// Get a pointer to the element
	Element<3,3>* p_element = mpExtendedMesh->GetElement(elementIndex);

	c_vector<double, 3> location0 = mpExtendedMesh->GetNode(p_element->GetNodeGlobalIndex(0))->rGetLocation();
	c_vector<double, 3> location1 = mpExtendedMesh->GetNode(p_element->GetNodeGlobalIndex(1))->rGetLocation();
	c_vector<double, 3> location2 = mpExtendedMesh->GetNode(p_element->GetNodeGlobalIndex(2))->rGetLocation();
	c_vector<double, 3> location3 = mpExtendedMesh->GetNode(p_element->GetNodeGlobalIndex(3))->rGetLocation();

    double edge_length_01 = norm_2(location0-location1);
    double edge_length_02 = norm_2(location0-location2);
    double edge_length_03 = norm_2(location0-location3);
    double edge_length_12 = norm_2(location1-location2);
    double edge_length_13 = norm_2(location1-location3);
    double edge_length_32 = norm_2(location3-location2);

	if ( (edge_length_01>maxEdgeLength) || (edge_length_02>maxEdgeLength) || (edge_length_03>maxEdgeLength)
			|| (edge_length_12>maxEdgeLength) || (edge_length_13>maxEdgeLength) || (edge_length_32>maxEdgeLength) )
	{
		element_contains_long_edge = true;
	}

	return element_contains_long_edge;
}

// Dom - adds the force to the cell - will use this
void PeriodicBendingForce3dWithGhostNodes::AddForceContribution(AbstractCellPopulation<3>& rCellPopulation)
//void PeriodicBendingForce3dWithGhostNodes::AddForceContribution(std::vector<c_vector<double, 3> >& rForces,
//                                                                   AbstractCellPopulation<3>& rCellPopulation)
{
	DomMeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*>(&rCellPopulation);

	unsigned num_cells = rCellPopulation.GetNumRealCells();
	//unsigned num_nodes = rCellPopulation.GetNumNodes();
	
	//PRINT_2_VARIABLES(num_cells,num_nodes);

	//std::vector<c_vector<double, 3>> rForces(3*num_cells); // Legacy code
    
	// If the width of the periodic domain has not been specified, use the initial width of the cell population
    if (mPeriodicDomainWidth == DOUBLE_UNSET)
    {
        mPeriodicDomainWidth = rCellPopulation.GetWidth(0);
    }

	// If the width of the periodic domain has not been specified, use the initial width of the cell population
    if (mPeriodicDomainDepth == DOUBLE_UNSET)
    {
    	mPeriodicDomainDepth = rCellPopulation.GetWidth(1);
    }

    mExtendedMeshNodeIndexMap.clear();

    // Create a vector of nodes for use in constructing mpExtendedMesh
    std::vector<Node<3>*> extended_nodes(4*num_cells);

    // We iterate over all cells in the population
    unsigned count = 0;

	// Dom - Create a copy of original mesh
	for (AbstractCellPopulation<3>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // First, create and store a copy of this real node and cell
        unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        c_vector<double, 3> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

        // Create a copy of the node corresponding to this cell and store it
        Node<3>* p_real_node = new Node<3>(real_node_index, real_node_location);
        extended_nodes[count] = p_real_node;

        // Populate mExtendedMeshNodeIndexMap
        mExtendedMeshNodeIndexMap[count] = real_node_index;


        count++;
    }


    // First, extend the mesh in the x-direction
    for (AbstractCellPopulation<3>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // First, create and store a copy of this real node and cell
        unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        c_vector<double, 3> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

        // Create a copy of the node corresponding to this cell and store it
        Node<3>* p_real_node = new Node<3>(real_node_index, real_node_location);
	 //        extended_nodes[count] = p_real_node;

        // Compute the location of the image node corresponding to this node
        c_vector<double,3> image_node_location = real_node_location;
        if (real_node_location[0] >= mPeriodicDomainWidth*0.5)
        {
            image_node_location[0] -= mPeriodicDomainWidth;
        }
        else if (real_node_location[0] <  mPeriodicDomainWidth*0.5)
        {
            image_node_location[0] += mPeriodicDomainWidth;
        }

        // Create a copy of the node corresponding to this cell, suitable translated, and store it
	 //        Node<3>* p_image_node = new Node<3>(num_cells+count, image_node_location);
        Node<3>* p_image_node = new Node<3>(count, image_node_location);

	 //        extended_nodes[num_cells+count] = p_image_node;
        extended_nodes[count] = p_image_node;

        // Populate mExtendedMeshNodeIndexMap
        mExtendedMeshNodeIndexMap[count] = real_node_index;

        count++;
    }

    // Second, extend this extended mesh in the y-direction
    // (We don't need to store the real nodes anymore)
    for (AbstractCellPopulation<3>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // First, create and store a copy of this real node and cell
        unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        c_vector<double, 3> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

        // Compute the location of the image node corresponding to this node
        c_vector<double,3> image_node_location = real_node_location;

        if (real_node_location[1] >= mPeriodicDomainDepth*0.5)
        {
            image_node_location[1] -= mPeriodicDomainDepth;
        }
        else if (real_node_location[1] <  mPeriodicDomainDepth*0.5)
        {
            image_node_location[1] += mPeriodicDomainDepth;
        }

        // Create a copy of the node corresponding to this cell, suitable translated, and store it
	 //        Node<3>* p_image_node = new Node<3>(num_cells+count, image_node_location);
        Node<3>* p_image_node = new Node<3>(count, image_node_location);
	 //        extended_nodes[num_cells+count] = p_image_node;
        extended_nodes[count] = p_image_node;

        // Populate mExtendedMeshNodeIndexMap
	 //        mExtendedMeshNodeIndexMap[num_cells+count] = real_node_index;
        mExtendedMeshNodeIndexMap[count] = real_node_index;

        count++;
    }

	// Thirdly, extend this extended mesh so that we cover the corners too
    // (We don't need to store the real nodes anymore)
    for (AbstractCellPopulation<3>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // First, create and store a copy of this real node and cell
        unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        c_vector<double, 3> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

        // Compute the location of the image node corresponding to this node
        c_vector<double,3> image_node_location = real_node_location;

        if (real_node_location[1] >= mPeriodicDomainDepth*0.5)
        {
			if (real_node_location[0] >= mPeriodicDomainWidth*0.5)
			{
				image_node_location[0] -= mPeriodicDomainWidth;
			}
			else if (real_node_location[0] <  mPeriodicDomainWidth*0.5)
			{
				image_node_location[0] += mPeriodicDomainWidth;
			}
            image_node_location[1] -= mPeriodicDomainDepth;
        }
        else if (real_node_location[1] <  mPeriodicDomainDepth*0.5)
        {
			if (real_node_location[0] >= mPeriodicDomainWidth*0.5)
			{
				image_node_location[0] -= mPeriodicDomainWidth;
			}
			else if (real_node_location[0] <  mPeriodicDomainWidth*0.5)
			{
				image_node_location[0] += mPeriodicDomainWidth;
			}
            image_node_location[1] += mPeriodicDomainDepth;
        }

        // Create a copy of the node corresponding to this cell, suitable translated, and store it
	 //        Node<3>* p_image_node = new Node<3>(num_cells+count, image_node_location);
        Node<3>* p_image_node = new Node<3>(count, image_node_location);
	 //        extended_nodes[num_cells+count] = p_image_node;
        extended_nodes[count] = p_image_node;

        // Populate mExtendedMeshNodeIndexMap
	 //        mExtendedMeshNodeIndexMap[num_cells+count] = real_node_index;
        mExtendedMeshNodeIndexMap[count] = real_node_index;

        count++;
    }

	

    // We now construct mpExtendedMesh using extended_nodes
    if (mpExtendedMesh != NULL)
    {
    	delete mpExtendedMesh;
    }
    mpExtendedMesh = new MutableMesh<3,3>(extended_nodes);
    
	///////////////////////////////////////////////////////////////////////////////////////
	// Dom - everything above this is just creating the extended mesh - leave it as it is
	///////////////////////////////////////////////////////////////////////////////////////
    
	///\todo We may need to modify rCellPopulation to ensure thatmMarkedSprings is correct at this point

	// Dom - dont need it
	//// Identify epithelial-stromal cell pairs
	//std::vector<c_vector<unsigned, 2> > node_pairs = GetEpithelialTissuePairs(rCellPopulation);
	//std::cout<<"\n";
	//PRINT_VARIABLE("Start");

	std::vector<c_vector<unsigned, 3> > epithelial_triangulation = GetEpithelialMesh(rCellPopulation);

	std::vector<c_vector<unsigned, 10> > epithelial_neighbours = GetEpithelialNeighbours(epithelial_triangulation, num_cells);
	//	int second_order_neighs[100];
	//std::cout << "num_cells = " << num_cells << "\n";


	double basement_membrane_parameter = GetBasementMembraneParameter();

	
	//std::cout<<"\n\n";
	for(unsigned i=0; i<num_cells; i++)
	{
		
		
		//std::cout<<"\n";
        //PRINT_VARIABLE(i);
        //std::cout<<"\n";


        unsigned cell_i_ext = mExtendedMeshNodeIndexMap[i];
		//PRINT_VARIABLE(cell_i_ext);
		//std::cout << "cell_i_ext = " << cell_i_ext << "  i = " << i << "\n";
		CellPtr p_cell_i_ext = rCellPopulation.GetCellUsingLocationIndex(cell_i_ext);
    	bool cell_is_dead = p_cell_i_ext->IsDead();
		bool is_ghost = p_tissue->IsGhostNode(cell_i_ext);
		
		//PRINT_VARIABLE(cell_i_ext);
		//bool cell_is_ghost = p_cell_i_ext->IsDead();

		//Element<3,3>* p_element = p_tissue->rGetMesh().GetElement(elementIndex);

		//		if (p_tissue->IsGhostNode(p_element->GetNodeGlobalIndex(local_index)) == true)


		//bool apply_force = (p_cell_i_ext->GetMutationState()->IsType<WildTypeCellMutationState>() == true && cell_is_dead == false );
		//PRINT_3_VARIABLES(cell_i_ext,apply_force,is_ghost);

        if(p_cell_i_ext->GetMutationState()->IsType<WildTypeCellMutationState>() == true && cell_is_dead == false && is_ghost == false)
		{
			std::vector<unsigned> second_order_neighs;

			int first_order_neighs[10];
			
			for(unsigned j=0; j<10; j++)
			{
				first_order_neighs[j] = epithelial_neighbours[cell_i_ext][j];
				
						
				for(unsigned k=0; k<10; k++)
				{

					second_order_neighs.push_back(epithelial_neighbours[first_order_neighs[j]][k]);

				}
				
			}

			std::sort(second_order_neighs.begin(), second_order_neighs.end());
			second_order_neighs.erase(std::unique(second_order_neighs.begin(), second_order_neighs.end()), second_order_neighs.end());

			// Remove the value zero
			second_order_neighs.erase(std::remove(second_order_neighs.begin(), second_order_neighs.end(), 0), second_order_neighs.end());
			// Remove the value i
			//second_order_neighs.erase(std::remove(second_order_neighs.begin(), second_order_neighs.end(), i), second_order_neighs.end());
			// May aswell include i since we need it for the linear regression part.

			//std::cout<< second_order_neighs.size() << " ";

			c_vector<double, 3> epithelial_location = mpExtendedMesh->GetNode(cell_i_ext)->rGetLocation();
			
			std::vector<c_vector<double, 3> > image_location_per_second_order_neighbours;


			bool epithelial_cell_in_non_zero_target_curvature_region = (pow((epithelial_location[0] - mXCoordinateOfCentreOfNonZeroTargetCurvatureRegion),2)
																	  + pow((epithelial_location[1] - mYCoordinateOfCentreOfNonZeroTargetCurvatureRegion),2) 
																	 <= pow(mRadiusOfNonZeroTargetCurvatureRegion,2) );

			// here we will find
			if(epithelial_cell_in_non_zero_target_curvature_region == true)
			{
				// PRINT_VECTOR(epithelial_location);
				// PRINT_VARIABLE(cell_i_ext);
				// Are i and cell_i_ext always the same?
				//std::cout  << "cell_i_ext = " << cell_i_ext << " i = " << i << "\n";
				image_location_per_second_order_neighbours = FitPlaneAndFindImage(rCellPopulation, second_order_neighs, cell_i_ext);
				
				 
				basement_membrane_parameter = mBasementMembraneParameter;

				// Troubleshooting
				/*
				for(unsigned j=0; j<image_location_per_second_order_neighbours.size(); j++)
				{
					c_vector<double, 3> image_location_j;
					image_location_j[0] = image_location_per_second_order_neighbours[j][0];
					image_location_j[1] = image_location_per_second_order_neighbours[j][1];
					image_location_j[2] = image_location_per_second_order_neighbours[j][2];
					PRINT_VECTOR(image_location_j);
				}
				std::cout << "\n \n";
				*/
			}
			else if(epithelial_cell_in_non_zero_target_curvature_region == false)
			{
				basement_membrane_parameter = +1.0*mBasementMembraneParameter;

				for (unsigned j=0; j<second_order_neighs.size(); j++)
				{
					c_vector<double, 3> location_of_neighbour_j = mpExtendedMesh->GetNode(second_order_neighs[j])->rGetLocation();
					image_location_per_second_order_neighbours.push_back(location_of_neighbour_j);
				}
			}

			std::vector<unsigned> first_order_neighs_vect;
			for(unsigned j=0; j<10; j++)
			{
				first_order_neighs_vect.push_back(epithelial_neighbours[cell_i_ext][j]);
			}
			first_order_neighs_vect.push_back(cell_i_ext);
			// order and remove doubles
			std::sort(first_order_neighs_vect.begin(), first_order_neighs_vect.end());
			first_order_neighs_vect.erase(std::unique(first_order_neighs_vect.begin(), first_order_neighs_vect.end()), first_order_neighs_vect.end());
			// Remove zero from vector
			first_order_neighs_vect.erase(std::remove(first_order_neighs_vect.begin(), first_order_neighs_vect.end(), 0), first_order_neighs_vect.end());

			/* 
			Need:
			First order neighbours  ->  first_order_neighs_vect
			Traingulation		    ->  epithelial_triangulation
			image node positions    ->  image_location_per_second_order_neighbours 
			  needs (for node IDs)  ->  second_order_neighs 
			Cell_i					->  cell_i_ext
			*/
			//std::cout << "Implementing cell ";
			//PRINT_VARIABLE(cell_i_ext);

			c_vector<double, 3> force_due_to_curvature;

			if(first_order_neighs_vect.size() >= 3)
			{
				//PRINT_VECTOR(first_order_neighs_vect);
				force_due_to_curvature = GetForceDueToDiscreteCurvature(rCellPopulation, epithelial_triangulation, first_order_neighs_vect, second_order_neighs, image_location_per_second_order_neighbours, cell_i_ext, num_cells);
			}
			else if(first_order_neighs_vect.size() < 3)
			{
				//PRINT_VECTOR(first_order_neighs_vect);
				force_due_to_curvature[0] = 0.0;
				force_due_to_curvature[1] = 0.0;
				force_due_to_curvature[2] = 0.0;
			}
			
			//std::cout << " fin ";
			//PRINT_VARIABLE(cell_i_ext);

			//std::cout << "cell i = " << cell_i_ext << "\n";
			//PRINT_VECTOR(force_due_to_curvature);
			if (isnan(force_due_to_curvature[0]) || isnan(force_due_to_curvature[1]) || isnan(force_due_to_curvature[2]))
			{
				PRINT_VECTOR(force_due_to_curvature);
				force_due_to_curvature[0] = 0.0;
				force_due_to_curvature[1] = 0.0;
				force_due_to_curvature[2] = 0.0;
			}
			
			rCellPopulation.GetNode(cell_i_ext)->AddAppliedForceContribution(basement_membrane_parameter*force_due_to_curvature);
			//HOW_MANY_TIMES_HERE("inside for loop");

		
			//c_vector<double, 3> location_of_neighbour_j = mpExtendedMesh->GetNode(cell_i_ext)->rGetLocation();
			//std::cout << "cell_i = " << cell_i_ext << " ";
			//PRINT_VECTOR(location_of_neighbour_j);

			
			//c_vector<double, 3> node_i_location = p_tissue->GetNode(cell_i_ext)->rGetLocation();
			//PRINT_VECTOR(node_i_location);
			//PRINT_VARIABLE(basement_membrane_parameter);
			//PRINT_VECTOR(force_due_to_curvature);
		}
		
	}
	//PRINT_VARIABLE("Finish");
	
}


double PeriodicBendingForce3dWithGhostNodes::GetPeriodicDomainWidth()
{
	return mPeriodicDomainWidth;
}

void PeriodicBendingForce3dWithGhostNodes::SetPeriodicDomainWidth(double periodicDomainWidth)
{
	mPeriodicDomainWidth = periodicDomainWidth;
}

double PeriodicBendingForce3dWithGhostNodes::GetPeriodicDomainDepth()
{
	return mPeriodicDomainDepth;
}

void PeriodicBendingForce3dWithGhostNodes::SetPeriodicDomainDepth(double periodicDomainDepth)
{
	mPeriodicDomainDepth = periodicDomainDepth;
}


void PeriodicBendingForce3dWithGhostNodes::OutputForceParameters(out_stream& rParamsFile)
{
//	*rParamsFile <<  "\t\t\t<UsePositionDependentMembraneForce>"<<  mUsePositionDependentMembraneForce << "</UsePositionDependentMembraneForce> \n" ;
//	*rParamsFile <<  "\t\t\t<MembraneForceMultiplier>"<<  mMembraneForceMultiplier << "</MembraneForceMultiplier> \n" ;
//	*rParamsFile <<  "\t\t\t<UseOneWaySprings>"<<  mUseOneWaySprings << "</mUseOneWaySprings> \n" ;
//	*rParamsFile <<  "\t\t\t<CryptBaseCurvature>"<<  mCryptBaseCurvature << "</CryptBaseCurvature> \n" ;
//	*rParamsFile <<  "\t\t\t<CryptBaseLevel>"<<  mCryptBaseLevel << "</CryptBaseLevel> \n" ;
//	*rParamsFile <<  "\t\t\t<LeftBoundary>"<<  mLeftBoundary << "</LeftBoundary> \n" ;
//	*rParamsFile <<  "\t\t\t<RightBoundary>"<<  mRightBoundary << "</RightBoundary> \n" ;
	*rParamsFile <<  "\t\t\t<BasementMembraneParameter>"<<  mBasementMembraneParameter << "</BasementMembraneParameter> \n" ;
	*rParamsFile <<  "\t\t\t<ExponentParameter>"<<  mExponentParameter << "</ExponentParameter> \n" ;
	*rParamsFile <<  "\t\t\t<PeriodicDomainWidth>"<<  mPeriodicDomainWidth << "</PeriodicDomainWidth> \n" ;
	*rParamsFile <<  "\t\t\t<PeriodicDomainDepth>"<<  mPeriodicDomainDepth << "</PeriodicDomainDepth> \n" ;
	*rParamsFile <<  "\t\t\t<NonZeroTargetCurvature>"<<  mNonZeroTargetCurvature << "<NonZeroTargetCurvature> \n" ;
	*rParamsFile <<  "\t\t\t<RadiusOfCurvature>"<<  mRadiusOfNonZeroTargetCurvatureRegion << "<RadiusOfCurvature> \n" ;
	*rParamsFile <<  "\t\t\t<XCoordinateOfCurvature>"<<  mXCoordinateOfCentreOfNonZeroTargetCurvatureRegion << "<XCoordinateOfCurvature> \n" ;
	*rParamsFile <<  "\t\t\t<YCoordinateOfCurvature>"<<  mYCoordinateOfCentreOfNonZeroTargetCurvatureRegion << "<YCoordinateOfCurvature> \n" ;

	// Call direct parent class
	AbstractForce<3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PeriodicBendingForce3dWithGhostNodes)


