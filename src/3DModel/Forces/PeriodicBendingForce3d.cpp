#include "PeriodicBendingForce3d.hpp"
#include "SimulationTime.hpp"

PeriodicBendingForce3d::PeriodicBendingForce3d()
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

PeriodicBendingForce3d::~PeriodicBendingForce3d()
{
    // Avoid memory leaks
    if (mpExtendedMesh != NULL)
    {
        delete mpExtendedMesh;
    }
    
}

void PeriodicBendingForce3d::SetBasementMembraneParameter(double basementMembraneParameter)
{
	mBasementMembraneParameter = basementMembraneParameter;
}

double PeriodicBendingForce3d::GetBasementMembraneParameter()
{
	return mBasementMembraneParameter;
}

void PeriodicBendingForce3d::SetExponentParameter(double exponentParameter)
{
	mExponentParameter = exponentParameter;
}

double PeriodicBendingForce3d::GetExponentParameter()
{
	return mExponentParameter;
}


void PeriodicBendingForce3d::SetCircularNonZeroTargetCurvatureRegion(bool setNonZeroTargetCurvatureRegion, double nonZeroTargetCurvature,
		double radiusOfNonZeroTargetCurvatureRegion, double xCoordinateOfCentreOfNonZeroTargetCurvatureRegion,
		double yCoordinateOfCentreOfNonZeroTargetCurvatureRegion)
{
	mSetNonZeroTargetCurvatureRegion = setNonZeroTargetCurvatureRegion;
	mNonZeroTargetCurvature = nonZeroTargetCurvature;
	mRadiusOfNonZeroTargetCurvatureRegion = radiusOfNonZeroTargetCurvatureRegion;
	mXCoordinateOfCentreOfNonZeroTargetCurvatureRegion = xCoordinateOfCentreOfNonZeroTargetCurvatureRegion;
	mYCoordinateOfCentreOfNonZeroTargetCurvatureRegion = yCoordinateOfCentreOfNonZeroTargetCurvatureRegion;
}

bool PeriodicBendingForce3d::GetWhetherToSetNonZeroTargetCurvatureRegion()
{
	return mSetNonZeroTargetCurvatureRegion;
}

double PeriodicBendingForce3d::GetNonZeroTargetCurvature()
{
	return mNonZeroTargetCurvature;
}

double PeriodicBendingForce3d::GetRadiusOfNonZeroTargetCurvatureRegion()
{
	return mRadiusOfNonZeroTargetCurvatureRegion;
}

double PeriodicBendingForce3d::GetXCoordinateOfCentreOfNonZeroTargetCurvatureRegion()
{
	return mXCoordinateOfCentreOfNonZeroTargetCurvatureRegion;
}

double PeriodicBendingForce3d::GetYCoordinateOfCentreOfNonZeroTargetCurvatureRegion()
{
	return mYCoordinateOfCentreOfNonZeroTargetCurvatureRegion;
}

/*
 * Remove repeated entries in a 1D vector
 */
void PeriodicBendingForce3d::RemoveDuplicates1D(std::vector<unsigned>& rVectorWithDuplicates)
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
std::vector<c_vector<unsigned, 2> > PeriodicBendingForce3d::GetEpithelialTissuePairs(AbstractCellPopulation<3>& rCellPopulation)
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

std::vector<c_vector<unsigned, 10> > PeriodicBendingForce3d::GetEpithelialNeighbours(std::vector<c_vector<unsigned, 3> > rEpithelialMeshVector, int number_of_cells)
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

		/*
		if(std::find(temp_neighbours[node_A].begin(), temp_neighbours[node_A].end(), node_B) != temp_neighbours[node_A].end())
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
		if(std::find(temp_neighbours[node_A].begin(), temp_neighbours[node_A].end(), node_C) != temp_neighbours[node_A].end())
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

		}
		
		if(std::find(temp_neighbours[node_B].begin(), temp_neighbours[node_B].end(), node_A) != temp_neighbours[node_B].end())
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
		if(std::find(temp_neighbours[node_B].begin(), temp_neighbours[node_B].end(), node_C) != temp_neighbours[node_B].end())
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
		
		if(std::find(temp_neighbours[node_C].begin(), temp_neighbours[node_C].end(), node_A) != temp_neighbours[node_C].end())
		{
			int iter_C = 0;
			while(iter_C < 10)
			{
				if(temp_neighbours[node_C][iter_C] == 0)
				{
					temp_neighbours[node_C][iter_C] = node_A;
					iter_C = 10;
				}
				iter_C++;
			}
		}
		if(std::find(temp_neighbours[node_C].begin(), temp_neighbours[node_C].end(), node_B) != temp_neighbours[node_C].end())
		{
			int iter_C = 0;
			while(iter_C < 10)
			{
				if(temp_neighbours[node_C][iter_C] == 0)
				{
					temp_neighbours[node_C][iter_C] = node_B;
					iter_C = 10;
				}
				iter_C++;
			}
		}*/

		/*int iter_A = 0;
		while(iter_A < 10)
		{
			if(temp_neighbours[node_A][iter_A] == 0)
			{
				temp_neighbours[node_A][iter_A] = node_B;
				temp_neighbours[node_A][iter_A+1] = node_C;
				iter_A = 10;
			}
			iter_A++;
		}
		int iter_B = 0;
		while(iter_B < 10)
		{
			if(temp_neighbours[node_B][iter_B] == 0)
			{
				temp_neighbours[node_B][iter_B] = node_A;
				temp_neighbours[node_B][iter_B+1] = node_C;
				iter_B = 10;
			}
			iter_B++;
		}
		int iter_C = 0;
		while(iter_C < 10)
		{
			if(temp_neighbours[node_C][iter_C] == 0)
			{
				temp_neighbours[node_C][iter_C] = node_A;
				temp_neighbours[node_C][iter_C+1] = node_B;
				iter_C = 10;
			}
			iter_C++;
		}*/

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

// Dom - attempt to write a function to find epithelial node neighbours

//std::vector<c_vector<unsigned, 10> > PeriodicBendingForce3d::GetEpithelialNeighbours(AbstractCellPopulation<3>& rCellPopulation);

std::vector<c_vector<unsigned, 3> > PeriodicBendingForce3d::GetEpithelialMesh(AbstractCellPopulation<3>& rCellPopulation)
{
	//int num_cells = rCellPopulation.GetNumRealCells();

	//std::array<unsigned, 10> node_neighbours = {};
	//This explicitly relys on the fact that cell 0 is not an epithelial cell - which is NOT true for a monolayer model
	// This uses a magic number of 10, meaning that we assume epithelial cells have at most 10 neighbouring epithelial cells
	// This may need to be reviewed later on.
	//unsigned node_neighbours[9*num_cells][10] = {};
	

	//std::vector<c_vector<unsigned, 2> > node_pairs; Delete this
	
	//for (typename MeshBasedCellPopulation<DIM>::SpringIterator spring_iterator=(static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))->SpringsBegin();
    //    spring_iterator!=(static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))->SpringsEnd();
    //    ++spring_iterator)
	/*
	for (typename MeshBasedCellPopulation<3>::SpringIterator spring_iterator=mpExtendedMesh->SpringsBegin();
        spring_iterator!=mpExtendedMesh->SpringsEnd();
        ++spring_iterator)	
    {
        unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
        unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();

		unsigned real_nodeA = mExtendedMeshNodeIndexMap[nodeA_global_index];
		unsigned real_nodeB = mExtendedMeshNodeIndexMap[nodeB_global_index];
        
		CellPtr p_real_nodeA = rCellPopulation.GetCellUsingLocationIndex(real_nodeA);
 		CellPtr p_real_nodeB = rCellPopulation.GetCellUsingLocationIndex(real_nodeB);

		if((p_real_nodeA->GetMutationState()->IsType<WildTypeCellMutationState>()) && (p_real_nodeB->GetMutationState()->IsType<WildTypeCellMutationState>()) )
		{
			 //Its going to looks something like:
			 //
			 //node_neighbours[real_nodeA][iterator_for_A_somehow] = real_nodeB;
			 //node_neighbours[real_nodeB][iterator_for_B_somehow] = real_nodeA;
			 //
			 //but probably using push_back, i.e.:
			 //
			 int iter_A = 0;
			 while(iter_A < 10)
			 {
				 if(node_neighbours[real_nodeA][iter_A] == 0)
				{
					node_neighbours[real_nodeA][iter_A] = real_nodeB;
					iter_A = 10;
				}
				 iter_A++;
			 }
			 int iter_B = 0;
			 while(iter_B < 10)
			 {
				 if(node_neighbours[real_nodeB][iter_B] == 0)
				{
					node_neighbours[real_nodeB][iter_B] = real_nodeA;
					iter_B = 10;
				}
				 iter_B++;
			 }
		}

    }*/


	/////////
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

    	if (is_live_epithelial_cell)
    	{
    	    // Iterate over elements of mpExtendedMesh containing this node and no ghost nodes
    		
			c_vector<unsigned, 3> tri_el;

    		for (Node<3>::ContainingElementIterator elem_iter = p_node->ContainingElementsBegin();
		         elem_iter != p_node->ContainingElementsEnd();
		         ++elem_iter)
    		{
				// Get a pointer to the element
				Element<3,3>* p_element = mpExtendedMesh->GetElement(*elem_iter);

				// ITERATE OVER NODES owned by this element
				std::vector<unsigned> temp_triangular_element;
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

				}

				if(temp_triangular_element.size() == 3)
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

	/*
	std::cout << epithelial_triangulation.size() << "\n";
	for(unsigned i=0; i<epithelial_triangulation.size(); i++)
	{
		if(1)//epithelial_triangulation[i][0] == 36)
		{
			std::cout << epithelial_triangulation[i][0] << " " << epithelial_triangulation[i][1] << " " << epithelial_triangulation[i][2] << "\n";
		}
	}
	*/


	return epithelial_triangulation;
}

std::vector<c_vector<double, 3> > PeriodicBendingForce3d::FitPlaneAndFindImage(AbstractCellPopulation<3>& rCellPopulation, std::vector<unsigned> second_order_neighs,  unsigned cell_i)
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

	double radius_from_curvature = 1/sqrt(mNonZeroTargetCurvature);

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
		point_on_sphere[0] = normal_vector[0]*param_T + epithelial_location[0];
		point_on_sphere[1] = normal_vector[1]*param_T + epithelial_location[1];
		point_on_sphere[2] = 1.0*param_T + epithelial_location[2];

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
		
		//std::cout << "T = "<< param_T << " desc = " << pow(beta,2) - 4*alpha*delta << "\n";
		//PRINT_VECTOR(epithelial_location);
		image_location_per_second_order_neighbours.push_back(temp_cell_location);
	}

	/*for (unsigned i=0; i<second_order_neighs.size(); i++)
	{
		//int epithelialNodeIndex = second_order_neighs[i];
		c_vector<double, 3> epithelial_location = mpExtendedMesh->GetNode(second_order_neighs[i])->rGetLocation();
		
		//PRINT_VECTOR(unit_normal);
		//std::cout << pow(unit_normal[0],2) + pow(unit_normal[1],2) + pow(unit_normal[2],2) << "\n";
		double alpha = pow(unit_normal[0],2) + pow(unit_normal[1],2) + pow(unit_normal[2],2);
		double beta  = 2*unit_normal[0]*(epithelial_location[0]-circle_centre[0]) + 2*unit_normal[1]*(epithelial_location[1]-circle_centre[1]) + 2*unit_normal[2]*(epithelial_location[2]-circle_centre[2]);
		double delta = pow(epithelial_location[0]-circle_centre[0],2) + pow(epithelial_location[1]-circle_centre[1],2) + pow(epithelial_location[2]-circle_centre[2],2) - pow(radius_from_curvature,2);

		double param_T = (-1.0*beta - sign_of_curvature*sqrt(pow(beta,2) - 4*alpha*delta))/(2.0*alpha);
		
		c_vector<double, 3> point_on_sphere;
		point_on_sphere[0] = unit_normal[0]*param_T + epithelial_location[0];
		point_on_sphere[1] = unit_normal[1]*param_T + epithelial_location[1];
		point_on_sphere[2] = unit_normal[2]*param_T + epithelial_location[2];

		double k_point = (unit_normal[0]*point_on_sphere[0] + unit_normal[1]*point_on_sphere[1] + unit_normal[2]*point_on_sphere[2] - normal_vector[2])/alpha;
		
		c_vector<double, 3> point_on_plane;
		point_on_plane[0] = point_on_sphere[0] - k_point*unit_normal[0];
		point_on_plane[1] = point_on_sphere[1] - k_point*unit_normal[1];
		point_on_plane[2] = point_on_sphere[2] - k_point*unit_normal[2];

		double point_on_plane_DOT_point_on_sphere = sqrt(pow(point_on_plane[0]-point_on_sphere[0],2) + pow(point_on_plane[1]-point_on_sphere[1],2) + pow(point_on_plane[2]-point_on_sphere[2],2));

		c_vector<double, 3> temp_cell_location;
		temp_cell_location[0] = epithelial_location[0] - unit_normal[0]*point_on_plane_DOT_point_on_sphere;
		temp_cell_location[1] = epithelial_location[1] - unit_normal[1]*point_on_plane_DOT_point_on_sphere;
		temp_cell_location[2] = epithelial_location[2] - unit_normal[2]*point_on_plane_DOT_point_on_sphere;
		
		std::cout << pow(beta,2) - 4*alpha*delta << "\n";
		PRINT_VECTOR(temp_cell_location);
		image_location_per_second_order_neighbours.push_back(temp_cell_location);
	}*/
	//std::cout<< "\n";

	// for(unsigned i=0; i<image_location_per_second_order_neighbours.size(); i++)
	//{
	//	c_vector<double, 3> temp_cell_loc;
	//	temp_cell_loc[0] = image_location_per_second_order_neighbours[i][0];
	//	temp_cell_loc[1] = image_location_per_second_order_neighbours[i][1];
	//	temp_cell_loc[2] = image_location_per_second_order_neighbours[i][2];

	//	PRINT_VECTOR(temp_cell_loc);
	//}
	//std::cout << "\n";
	return image_location_per_second_order_neighbours;
}



c_vector<double, 3> PeriodicBendingForce3d::GetForceDueToDiscreteCurvature(AbstractCellPopulation<3>& rCellPopulation, std::vector<c_vector<unsigned, 3> > rEpithelialMeshVector, std::vector<unsigned> first_order_neighs,  std::vector<unsigned> second_order_neighs, std::vector<c_vector<double, 3> > image_location, unsigned cell_i, unsigned num_cells)
{
	
	/* 
	Need:
	Traingulation		    ->  rEpithelialMeshVector
	First order neighbours  ->  first_order_neighs 
	image node positions    ->  image_location 
	  needs (for node IDs)  ->  second_order_neighs 
	Cell_i					->  cell_i
	*/
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

				ang_sum += acos(cos_abc);

				// Calculate grad_i phi_i
				if(cell_a == cell_i)
				{
					grad_i_phi_j_el += (1.0/(sqrt(1.0 - cos_abc*cos_abc)))*(1.0/(mag_ab*mag_ab*mag_ac*mag_ac))*(mag_ab*(mag_ac-mag_ab*cos_abc)*vect_ac + mag_ac*(mag_ab-mag_ac*cos_abc)*vect_ab);
					//std::cout << "case 1 " << cell_a << "\n";
				}
				
				else if(cell_b == cell_i)
				{
					grad_i_phi_j_el += (1.0/(sqrt(1.0 - cos_abc*cos_abc)))* (1.0/(mag_ab)) * (vect_ac/mag_ac - cos_abc*vect_ab/mag_ab);
					//std::cout << "case 2 " << cell_a << "\n";
				}

				else if(cell_c == cell_i)
				{
					grad_i_phi_j_el += (1.0/(sqrt(1.0 - cos_abc*cos_abc)))* (1.0/(mag_ac)) * (vect_ab/mag_ab - cos_abc*vect_ac/mag_ac);
					//std::cout << "case 3 " << cell_a << "\n";
				}
			}
			
		}
		//if(cell_a == cell_i)
		//{
		//	PRINT_VECTOR(grad_i_phi_j_el);
		//}
		//std::cout << "Angle Sum = " << ang_sum << "\n";

		angle_sum.push_back(ang_sum);
		angle_sum_indicats.push_back(cell_a);

		//PRINT_VECTOR(force_due_to_curvature);
		// Using alpha = 1.5 and beta = 2.0 for now, will fix this later though
		double sign_of_anlge = (round((2*3.141592653589793 - ang_sum)*pow(10.0,14))/pow(10.0,14) > 0) - (round((2*3.141592653589793 - ang_sum)*pow(10.0,14))/pow(10.0,14) < 0);
		force_due_to_curvature +=  (exponent_parameter)*sign_of_anlge*pow(sqrt((2*3.141592653589793 - ang_sum)*(2*3.141592653589793 - ang_sum)),exponent_parameter-1.0)*grad_i_phi_j_el;
		
		//std::cout << round((2*3.141592653589793 - ang_sum)*pow(10.0,14))/pow(10.0,14) << "\n";
		//std::cout<< "sign_of_anlge = " << sign_of_anlge << "\n";
	}

	

	//std::cout << "cell_i = " << cell_i << "\n"; 
	//PRINT_VECTOR(angle_sum);
	//PRINT_VECTOR(angle_sum_indicats);
	//PRINT_VECTOR(first_order_neighs);
	//PRINT_VECTOR(force_due_to_curvature);
//	std::cout << "\n";
	
	
	
	return force_due_to_curvature;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// Dom - here is where the curvature is calculated (externally), which seems like a nice idea.
/* A method to take in an Epithelial-Tissue node pairing and return the curvature value
 * This takes in indices of nodes in the extended mesh, so to get the corresponding cells in the cell population, we use the
 * extended node map
 */
double PeriodicBendingForce3d::GetCurvatureDueToSpringMidpoints(AbstractCellPopulation<3>& rCellPopulation, unsigned epithelialNodeIndex, unsigned tissueNodeIndex)
{

	// Get the node indices that correspond to the real cells in the cell population that have been imaged
	unsigned real_epithelial_node = mExtendedMeshNodeIndexMap[epithelialNodeIndex];
	unsigned real_tissue_node = mExtendedMeshNodeIndexMap[tissueNodeIndex];
	
	CellPtr p_epithelial_cell = rCellPopulation.GetCellUsingLocationIndex(real_epithelial_node);
 	CellPtr p_tissue_cell = rCellPopulation.GetCellUsingLocationIndex(real_tissue_node);

	assert(p_epithelial_cell->GetMutationState()->IsType<WildTypeCellMutationState>());
 	assert(p_tissue_cell->GetMutationState()->IsType<StromalCellMutationState>());

    // Seeking the common elements that contain both the epithelial node and the tissue node
	// Note: Don't want to count any elements that have ghost nodes (this shouldn't matter if you just are using cutoff lengths)

	std::vector<unsigned> common_elements;	// Initialising
    std::vector<double> local_dAds;	// The value of dA/ds for a particular common element (these will all be summed)
    c_vector<double, 3> cross_product, cross_product_derivative;

	// The elements that contain this node in the extended mesh
    std::set<unsigned> epithelial_elements = mpExtendedMesh->GetNode(epithelialNodeIndex)->rGetContainingElementIndices();
    assert(!epithelial_elements.empty());

    // Check there are elements that contain the tissue node:
    assert(mpExtendedMesh->GetNode(tissueNodeIndex)->GetNumContainingElements() != 0);

    // Loop over all elements that contain the tissue node
    for (Node<3>::ContainingElementIterator elt_it = mpExtendedMesh->GetNode(tissueNodeIndex)->ContainingElementsBegin();
         elt_it != mpExtendedMesh->GetNode(tissueNodeIndex)->ContainingElementsEnd();
         ++elt_it)
    {
    	unsigned elt_index = *elt_it;

    	// We don't want to consider any elements with really long edges
    	bool elt_contains_long_edge = DoesElementContainLongEdge(rCellPopulation, elt_index, 1.5);

    	// Keep only those elements that also contain the epithelial node (there shouldn't be any ghost nodes in the extended mesh)
        if ( (elt_contains_long_edge == false) && (epithelial_elements.find(elt_index) != epithelial_elements.end()) )
        {
        	// Common element
           	common_elements.push_back(elt_index);
        }
    }

	//	assert(!common_elements.empty());		// This is bad - but may happen if we eliminate all cases with long edge...

	// We iterate over these common elements to find the midpoints of the other springs that
	// connect epithelial and  tissue nodes

	c_vector<double, 3> epithelial_location = mpExtendedMesh->GetNode(epithelialNodeIndex)->rGetLocation();
	c_vector<double, 3> tissue_location = mpExtendedMesh->GetNode(tissueNodeIndex)->rGetLocation();

	// The original spring midpoint, between the original epithelial and tissue nodes, which corresponds to v1 in your workings
	c_vector<double, 3> central_spring_midpoint = epithelial_location + 0.5*(mpExtendedMesh->GetVectorFromAtoB(epithelial_location, tissue_location));

	// We now have to consider each common tetrahedral element to calculate the other midpoints (v2 and v3) and work out each
	// triangle individually
	for (std::vector<unsigned>::iterator elem_iter = common_elements.begin();
										 elem_iter != common_elements.end();
										 ++elem_iter)
	{
		double local_area = 0.0;
		double sum_terms = 0.0;

		// Need to determine the cell type of each local node in the element as we want to find
		// the midpoint between epithelial and tissue pairs

		// Looking at the four nodes which form the vertices of the element

		unsigned global_index_A = mpExtendedMesh->GetElement(*elem_iter)->GetNodeGlobalIndex(0);
		unsigned global_index_B = mpExtendedMesh->GetElement(*elem_iter)->GetNodeGlobalIndex(1);
		unsigned global_index_C = mpExtendedMesh->GetElement(*elem_iter)->GetNodeGlobalIndex(2);
		unsigned global_index_D = mpExtendedMesh->GetElement(*elem_iter)->GetNodeGlobalIndex(3);

		// We will assign labels to the local nodes, retaining the epithelial node and the tissue node of interest
		unsigned E=UINT_MAX;  // E - the epithelial node we're interested in,
		unsigned T=UINT_MAX;  // T - the tissue node it's connected to, and
		unsigned P=UINT_MAX;  // P - one of the other nodes in the element,
		unsigned Q=UINT_MAX;  // Q - one of the other nodes in the element,

		// Need to assign these nodes correctly, so that we know which are the epithelial node and tissue node of interest
		if (global_index_A == epithelialNodeIndex)
		{
			E = global_index_A;

			if (global_index_B == tissueNodeIndex)
			{
				T = global_index_B;
				P = global_index_C;
				Q = global_index_D;
			}
			else if (global_index_C == tissueNodeIndex)
			{
				T = global_index_C;
				P = global_index_B;
				Q = global_index_D;
			}
			else if (global_index_D == tissueNodeIndex)
			{
				T = global_index_D;
				P = global_index_B;
				Q = global_index_C;
			}
		}
		else if (global_index_B == epithelialNodeIndex)
		{
			E = global_index_B;

			if (global_index_A == tissueNodeIndex)
			{
				T = global_index_A;
				P = global_index_C;
				Q = global_index_D;
			}
			else if (global_index_C == tissueNodeIndex)
			{
				T = global_index_C;
				P = global_index_A;
				Q = global_index_D;
			}
			else if (global_index_D == tissueNodeIndex)
			{
				T = global_index_D;
				P = global_index_A;
				Q = global_index_C;
			}
		}
		else if (global_index_C == epithelialNodeIndex)
		{
			E = global_index_C;

			if (global_index_A == tissueNodeIndex)
			{
				T = global_index_A;
				P = global_index_B;
				Q = global_index_D;
			}
			else if (global_index_B == tissueNodeIndex)
			{
				T = global_index_B;
				P = global_index_A;
				Q = global_index_D;
			}
			else if (global_index_D == tissueNodeIndex)
			{
				T = global_index_D;
				P = global_index_A;
				Q = global_index_B;
			}
		}
		else if (global_index_D == epithelialNodeIndex)
		{
			E = global_index_D;

			if (global_index_A == tissueNodeIndex)
			{
				T = global_index_A;
				P = global_index_B;
				Q = global_index_C;
			}
			else if (global_index_B == tissueNodeIndex)
			{
				T = global_index_B;
				P = global_index_A;
				Q = global_index_C;
			}
			else if (global_index_C == tissueNodeIndex)
			{
				T = global_index_C;
				P = global_index_A;
				Q = global_index_B;
			}
		}

		assert(E<UINT_MAX);
		assert(T<UINT_MAX);
		assert(P<UINT_MAX);
		assert(Q<UINT_MAX);

		// Check that we've identified the original epithelial and tissue pair correctly
		assert(E == epithelialNodeIndex);
		assert(T == tissueNodeIndex);

		/* Now we work with E (the epithelial node), T (the tissue node it's connected to) and P and Q, the other nodes.
		 * We need to work out if P and Q are epithelial or tissue nodes to then find the midpoint we require in this element.
		 */
		c_vector<double, 3> vector_E_to_P = mpExtendedMesh->GetVectorFromAtoB(mpExtendedMesh->GetNode(E)->rGetLocation(),mpExtendedMesh->GetNode(P)->rGetLocation());
		c_vector<double, 3> vector_P_to_T = mpExtendedMesh->GetVectorFromAtoB(mpExtendedMesh->GetNode(P)->rGetLocation(),mpExtendedMesh->GetNode(T)->rGetLocation());
		c_vector<double, 3> vector_E_to_Q = mpExtendedMesh->GetVectorFromAtoB(mpExtendedMesh->GetNode(E)->rGetLocation(),mpExtendedMesh->GetNode(Q)->rGetLocation());
		c_vector<double, 3> vector_Q_to_T = mpExtendedMesh->GetVectorFromAtoB(mpExtendedMesh->GetNode(Q)->rGetLocation(),mpExtendedMesh->GetNode(T)->rGetLocation());

	 // Dom - Check E is Epithelial cell
	 //		unsigned real_other_node_E = mExtendedMeshNodeIndexMap[E];
	 //		boost::shared_ptr<AbstractCellMutationState> E_mutation_state = rCellPopulation.GetCellUsingLocationIndex(real_other_node_E)->GetMutationState();
	 //		std::cout << "Check if E is Epithelial" << (E_mutation_state->IsType<WildTypeCellMutationState>()==true) << '\n';

		// Now check the mutation state of P and Q to see if they are tissue cells or epithelial cells
		unsigned real_other_node_P = mExtendedMeshNodeIndexMap[P];
		boost::shared_ptr<AbstractCellMutationState> P_mutation_state = rCellPopulation.GetCellUsingLocationIndex(real_other_node_P)->GetMutationState();

		unsigned real_other_node_Q = mExtendedMeshNodeIndexMap[Q];
		boost::shared_ptr<AbstractCellMutationState> Q_mutation_state = rCellPopulation.GetCellUsingLocationIndex(real_other_node_Q)->GetMutationState();

		// Now choose the appropriate midpoints to use to form the triangle. Note that we don't want to choose the midpoint of PQ,
		// so we always choose the midpoints of EP, EQ, TP or TQ, depending on whether P and Q are epithelial or tissue nodes
		c_vector<double, 3> midpoint_2, midpoint_3;

		if (P_mutation_state->IsType<StromalCellMutationState>()==false)		// P = Epithelial, i.e not labelled
		{
			midpoint_3 = mpExtendedMesh->GetNode(P)->rGetLocation() + 0.5*vector_P_to_T;

			if (Q_mutation_state->IsType<StromalCellMutationState>()==false)	// Q = Epithelial
			{
				midpoint_2 = mpExtendedMesh->GetNode(Q)->rGetLocation() + 0.5*vector_Q_to_T;
			}
			else if (Q_mutation_state->IsType<StromalCellMutationState>()==true)	// Q = Tissue, so labelled
			{
				midpoint_2 = mpExtendedMesh->GetNode(E)->rGetLocation() + 0.5*vector_E_to_Q;
			}
		}
		else if (P_mutation_state->IsType<StromalCellMutationState>()==true)	// P = Tissue
		{
			midpoint_3 = mpExtendedMesh->GetNode(E)->rGetLocation() + 0.5*vector_E_to_P;

			if(Q_mutation_state->IsType<StromalCellMutationState>()==true)	// Q = Tissue
			{
				midpoint_2 = mpExtendedMesh->GetNode(E)->rGetLocation() + 0.5*vector_E_to_Q;
			}
			else if (Q_mutation_state->IsType<StromalCellMutationState>()==false)	// Q = Epithelial
			{
				midpoint_2 = mpExtendedMesh->GetNode(Q)->rGetLocation() + 0.5*vector_Q_to_T;
			}
		}

		// Now we want the area of the triangle connecting these two midpoints to the central midpoint

		c_vector<double, 3> vector_to_midpoint_2 = mpExtendedMesh->GetVectorFromAtoB(central_spring_midpoint, midpoint_2);
		c_vector<double, 3> vector_to_midpoint_3 = mpExtendedMesh->GetVectorFromAtoB(central_spring_midpoint, midpoint_3);

		// Consider the vertices of the triangle as labelled A (central midpoint), B (midpoint 1) and C (midpoint 2)
		local_area = GetAreaOfTriangle(vector_to_midpoint_2, vector_to_midpoint_3);		// Checked

		if (isnan(local_area))
		{
			PRINT_VECTOR(vector_to_midpoint_2);
			PRINT_VECTOR(vector_to_midpoint_3);
			// NB. These indices correspond to the cell population, not the extended mesh
			PRINT_4_VARIABLES(SimulationTime::Instance()->GetTime(), epithelialNodeIndex, tissueNodeIndex, local_area);
		}
		assert(!isnan(local_area));

		// New definition of the central midpoint dependent on s, for use in the curvature calculation
		// TODO: WILL THIS CHANGE WHEN YOU HAVE TO THINK ABOUT PERIODICITY?
        //central_spring_midpoint = tissue_location + 0.5*(epithelial_location - tissue_location);

		// In your notes: v_5 = epithelial node, v_4 = tissue node, v_3 and v_2 are the spring midpoints from the common element
		// (and v_1 is the midpoint between the epithelial and tissue node)

		cross_product[0] = midpoint_3[1]*midpoint_2[2] - 0.5*tissue_location[1]*midpoint_2[2] - 0.5*epithelial_location[1]*midpoint_2[2]
						   - midpoint_2[1]*midpoint_3[2] + 0.5*tissue_location[1]*midpoint_3[2] + 0.5*epithelial_location[1]*midpoint_3[2]
						   + 0.5*midpoint_2[1]*tissue_location[2] - 0.5*midpoint_3[1]*tissue_location[2]
						   + 0.5*midpoint_2[1]*epithelial_location[2] - 0.5*midpoint_3[1]*epithelial_location[2]; // Double checked

		cross_product[1] = -midpoint_3[0]*midpoint_2[2] + 0.5*tissue_location[0]*midpoint_2[2] + 0.5*epithelial_location[0]*midpoint_2[2]
						   + midpoint_2[0]*midpoint_3[2] - 0.5*tissue_location[0]*midpoint_3[2] - 0.5*epithelial_location[0]*midpoint_3[2]
						   - 0.5*midpoint_2[0]*tissue_location[2] + 0.5*midpoint_3[0]*tissue_location[2] - 0.5*midpoint_2[0]*epithelial_location[2]
						   + 0.5*midpoint_3[0]*epithelial_location[2];  // Double checked

		cross_product[2] = midpoint_3[0]*midpoint_2[1] - 0.5*tissue_location[0]*midpoint_2[1] - 0.5*epithelial_location[0]*midpoint_2[1]
						   - midpoint_2[0]*midpoint_3[1] + 0.5*tissue_location[0]*midpoint_3[1] + 0.5*epithelial_location[0]*midpoint_3[1]
						   + 0.5*midpoint_2[0]*tissue_location[1] - 0.5*midpoint_3[0]*tissue_location[1] + 0.5*midpoint_2[0]*epithelial_location[1]
						   - 0.5*midpoint_3[0]*epithelial_location[1];  // Double checked

		cross_product_derivative[0] = tissue_location[1]*midpoint_2[2] - epithelial_location[1]*midpoint_2[2] - tissue_location[1]*midpoint_3[2]
									  + epithelial_location[1]*midpoint_3[2] - midpoint_2[1]*tissue_location[2] + midpoint_3[1]*tissue_location[2]
									  + midpoint_2[1]*epithelial_location[2] - midpoint_3[1]*epithelial_location[2];  // Double checked

		cross_product_derivative[1] = -tissue_location[0]*midpoint_2[2] + epithelial_location[0]*midpoint_2[2] + tissue_location[0]*midpoint_3[2]
									  - epithelial_location[0]*midpoint_3[2] + midpoint_2[0]*tissue_location[2] - midpoint_3[0]*tissue_location[2]
									  - midpoint_2[0]*epithelial_location[2] + midpoint_3[0]*epithelial_location[2];  // Double checked

		cross_product_derivative[2] = tissue_location[0]*midpoint_2[1] - epithelial_location[0]*midpoint_2[1] - tissue_location[0]*midpoint_3[1]
									  + epithelial_location[0]*midpoint_3[1] - midpoint_2[0]*tissue_location[1] + midpoint_3[0]*tissue_location[1]
									  + midpoint_2[0]*epithelial_location[1] - midpoint_3[0]*epithelial_location[1];  // Double checked

		sum_terms = cross_product[0]*cross_product_derivative[0] + cross_product[1]*cross_product_derivative[1] + cross_product[2]*cross_product_derivative[2];

		assert(!isnan(sum_terms));

		double dAds;

		if (local_area == 0.0)
		{
			TRACE("triangle area zero");
			dAds = 0.0;
		}
		else
		{
			dAds = 0.5*(1.0/local_area)*sum_terms;		// for this triangle
		}

		assert(!isnan(dAds));

		local_dAds.push_back(dAds);
		// At this point you want to have dA/ds for this common element (push back into a vector which can then be summed
		// just after this loop ends)
	}

	assert(local_dAds.size() == common_elements.size());    // Need a triangle for every common tetrahedra

	double total_dAds = 0.0;		// This should be summed for each epithelial-tissue pair, and so should be the sum of each triangle found

	for (unsigned k=0; k<local_dAds.size(); k++)
	{
		total_dAds += local_dAds[k];
	}

	double curvature = total_dAds;

	if (mSetNonZeroTargetCurvatureRegion)	// The values will only be set if this is true
	{
		bool epithelial_cell_in_non_zero_target_curvature_region = (pow(epithelial_location[0] - mXCoordinateOfCentreOfNonZeroTargetCurvatureRegion,2)
				+ pow(epithelial_location[1] - mYCoordinateOfCentreOfNonZeroTargetCurvatureRegion,2) <= mRadiusOfNonZeroTargetCurvatureRegion*mRadiusOfNonZeroTargetCurvatureRegion);

		if(epithelial_cell_in_non_zero_target_curvature_region)
		{
			curvature += mNonZeroTargetCurvature;
		}
	}

	assert(!isnan(curvature));
	return curvature;
}

// Dom - Probalby wont need this once my force is implemented
/* Method to return the area of the triangle defined by it's three vertices: A, B and C
 *
 * area = 0.5 * ( |AB|^2*|AC|^2 - (AB.AC)^2 ) ^0.5
 *
 */
double PeriodicBendingForce3d::GetAreaOfTriangle(c_vector<double,3> vectorOneToTwo, c_vector<double,3> vectorOneToThree)
{
	// Need to find |AB|^2 and |AC|^2, i.e. the dot products of these vectors with themselves

	double vector_one_to_two_squared = vectorOneToTwo[0]*vectorOneToTwo[0] + vectorOneToTwo[1]*vectorOneToTwo[1] + vectorOneToTwo[2]*vectorOneToTwo[2];
	double vector_one_to_three_squared = vectorOneToThree[0]*vectorOneToThree[0] + vectorOneToThree[1]*vectorOneToThree[1] + vectorOneToThree[2]*vectorOneToThree[2];

	// Now need the dot product of the two vectors
	double dot_product = vectorOneToTwo[0]*vectorOneToThree[0] + vectorOneToTwo[1]*vectorOneToThree[1] + vectorOneToTwo[2]*vectorOneToThree[2];

	// Now the area can be calculated:
	double area = 0.5*sqrt( vector_one_to_two_squared*vector_one_to_three_squared - pow(dot_product,2) );

	return area;
}

// Dom - looks like some legacy code...
/*
 * Method to determine whether an element contains ghost nodes
 */
bool PeriodicBendingForce3d::DoesElementContainGhostNodes(AbstractCellPopulation<3>& rCellPopulation, unsigned elementIndex)
{
	MeshBasedCellPopulation<3>* p_tissue = static_cast<MeshBasedCellPopulation<3>*>(&rCellPopulation);

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
bool PeriodicBendingForce3d::DoesElementContainLongEdge(AbstractCellPopulation<3>& rCellPopulation, unsigned elementIndex, double maxEdgeLength)
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
void PeriodicBendingForce3d::AddForceContribution(AbstractCellPopulation<3>& rCellPopulation)
//void PeriodicBendingForce3d::AddForceContribution(std::vector<c_vector<double, 3> >& rForces,
//                                                                   AbstractCellPopulation<3>& rCellPopulation)
{
	MeshBasedCellPopulation<3>* p_tissue = static_cast<MeshBasedCellPopulation<3>*>(&rCellPopulation);

	unsigned num_cells = rCellPopulation.GetNumRealCells();
	
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
	std::vector<c_vector<unsigned, 3> > epithelial_triangulation = GetEpithelialMesh(rCellPopulation);
	std::vector<c_vector<unsigned, 10> > epithelial_neighbours = GetEpithelialNeighbours(epithelial_triangulation, num_cells);
//	int second_order_neighs[100];
	//std::cout << "num_cells = " << num_cells << "\n";
	
	double basement_membrane_parameter = GetBasementMembraneParameter();

	for(unsigned i=0; i<num_cells; i++)
	{
        unsigned cell_i_ext = mExtendedMeshNodeIndexMap[i];
		//PRINT_VARIABLE(cell_i_ext);
		//std::cout << "cell_i_ext = " << cell_i_ext << "  i = " << i << "\n";
		CellPtr p_cell_i_ext = rCellPopulation.GetCellUsingLocationIndex(cell_i_ext);
    	bool cell_is_dead = p_cell_i_ext->IsDead();

        if(p_cell_i_ext->GetMutationState()->IsType<WildTypeCellMutationState>() == true && cell_is_dead == false )
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


			bool epithelial_cell_in_non_zero_target_curvature_region = (pow(epithelial_location[0] - mXCoordinateOfCentreOfNonZeroTargetCurvatureRegion,2)
			+ pow(epithelial_location[1] - mYCoordinateOfCentreOfNonZeroTargetCurvatureRegion,2) <= mRadiusOfNonZeroTargetCurvatureRegion*mRadiusOfNonZeroTargetCurvatureRegion);

			// here we will find
			if(epithelial_cell_in_non_zero_target_curvature_region == true)
			{
				// Are i and cell_i_ext always the same?
				//std::cout  << "cell_i_ext = " << cell_i_ext << " i = " << i << "\n";
				image_location_per_second_order_neighbours = FitPlaneAndFindImage(rCellPopulation, second_order_neighs, cell_i_ext);
				
				/* Troubleshooting
				for(unsigned j=0; j<image_location_per_second_order_neighbours.size(); j++)
				{
					c_vector<double, 3> image_location_j;
					image_location_j[0] = image_location_per_second_order_neighbours[j][0];
					image_location_j[1] = image_location_per_second_order_neighbours[j][1];
					image_location_j[2] = image_location_per_second_order_neighbours[j][2];
					//PRINT_VECTOR(image_location_j);
				}
				std::cout << "\n \n";
				*/
			}
			else if(epithelial_cell_in_non_zero_target_curvature_region == false)
			{
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

			rCellPopulation.GetNode(cell_i_ext)->AddAppliedForceContribution(basement_membrane_parameter*force_due_to_curvature);
			
			//HOW_MANY_TIMES_HERE("inside for loop");

		
			//c_vector<double, 3> location_of_neighbour_j = mpExtendedMesh->GetNode(cell_i_ext)->rGetLocation();
			//std::cout << "cell_i = " << cell_i_ext << " ";
			//PRINT_VECTOR(location_of_neighbour_j);

			
			//c_vector<double, 3> node_i_location = p_tissue->GetNode(cell_i_ext)->rGetLocation();
			//PRINT_VECTOR(node_i_location);
		}
		
	}

	////
/*
	for (unsigned i=0; i<num_cells; i++)
	{

		unsigned extended_epithelial_node_index = node_pairs[i][0];
		unsigned extended_tissue_node_index = node_pairs[i][1];
		// Get the corresponding node index in rCellPopulation
        unsigned epithelial_node_index = mExtendedMeshNodeIndexMap[extended_epithelial_node_index];
        unsigned tissue_node_index = mExtendedMeshNodeIndexMap[extended_tissue_node_index];

        // Get the cells corresponding to these nodes to check the cell types
        CellPtr p_epithelial_cell = rCellPopulation.GetCellUsingLocationIndex(epithelial_node_index);
        assert(p_epithelial_cell->GetMutationState()->IsType<StromalCellMutationState>() == false);
        assert(p_epithelial_cell->GetMutationState()->IsType<WildTypeCellMutationState>() == true);

		CellPtr p_tissue_cell = rCellPopulation.GetCellUsingLocationIndex(tissue_node_index);
		assert(p_tissue_cell->GetMutationState()->IsType<StromalCellMutationState>() == true);		

		// Get the locations in the extended mesh
		c_vector<double, 3> epithelial_location = mpExtendedMesh->GetNode(extended_epithelial_node_index)->rGetLocation();
		c_vector<double, 3> tissue_location = mpExtendedMesh->GetNode(extended_tissue_node_index)->rGetLocation();

		// The force due to the basal lamina acts along the spring connecting the epithelial and tissue nodes, T->E direction
		c_vector<double, 3> curvature_force_direction = mpExtendedMesh->GetVectorFromAtoB(tissue_location, epithelial_location);

		double distance_between_nodes = norm_2(curvature_force_direction);
		assert(distance_between_nodes > 0);
		assert(!isnan(distance_between_nodes));

		curvature_force_direction /= distance_between_nodes;

		//////
		// Safety check
		unsigned real_epithelial_node = mExtendedMeshNodeIndexMap[extended_epithelial_node_index];
		unsigned real_tissue_node = mExtendedMeshNodeIndexMap[extended_tissue_node_index];
		
		CellPtr p_ext_epithelial_cell = rCellPopulation.GetCellUsingLocationIndex(real_epithelial_node);
		CellPtr p_ext_tissue_cell = rCellPopulation.GetCellUsingLocationIndex(real_tissue_node);
	
		assert(p_ext_epithelial_cell->GetMutationState()->IsType<WildTypeCellMutationState>());
		assert(p_ext_tissue_cell->GetMutationState()->IsType<StromalCellMutationState>());

		boost::shared_ptr<AbstractCellMutationState> p_extended_epithelial_node_mutation_state = p_ext_epithelial_cell->GetMutationState();
    	boost::shared_ptr<AbstractCellMutationState> p_extended_tissue_node_mutation_state = p_ext_tissue_cell->GetMutationState();
		////

		double curvature = GetCurvatureDueToSpringMidpoints(rCellPopulation, extended_epithelial_node_index, extended_tissue_node_index);

	 // Dom
	 //		if (extended_epithelial_node_index == 11)
	 //		{
	 //			std::cout<< "curvature at node 11 = " << curvature << "\n";
	 //		}
			

		double basement_membrane_parameter = GetBasementMembraneParameter();

		c_vector<double, 3> force_due_to_basement_membrane = -basement_membrane_parameter*curvature*curvature_force_direction;
		// Had this (v) original version with a minus sign...
		//c_vector<double, 3> force_due_to_basement_membrane = -basement_membrane_parameter*curvature*curvature_force_direction;

        // Now we make sure that we only apply the force to the real node and not the image node
        if (extended_epithelial_node_index < num_cells)
        {
			//rForces[epithelial_node_index] += force_due_to_basement_membrane;
			rCellPopulation.GetNode(epithelial_node_index)->AddAppliedForceContribution(-force_due_to_basement_membrane);
        }
		
	}
*/
}


double PeriodicBendingForce3d::GetPeriodicDomainWidth()
{
	return mPeriodicDomainWidth;
}

void PeriodicBendingForce3d::SetPeriodicDomainWidth(double periodicDomainWidth)
{
	mPeriodicDomainWidth = periodicDomainWidth;
}

double PeriodicBendingForce3d::GetPeriodicDomainDepth()
{
	return mPeriodicDomainDepth;
}

void PeriodicBendingForce3d::SetPeriodicDomainDepth(double periodicDomainDepth)
{
	mPeriodicDomainDepth = periodicDomainDepth;
}


void PeriodicBendingForce3d::OutputForceParameters(out_stream& rParamsFile)
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

	// Call direct parent class
	AbstractForce<3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PeriodicBendingForce3d)


