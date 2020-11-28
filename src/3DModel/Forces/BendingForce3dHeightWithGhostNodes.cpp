#include "BendingForce3dHeightWithGhostNodes.hpp"
#include "SimulationTime.hpp"
#include <cmath>
#include "DomMeshBasedCellPopulationWithGhostNodes.hpp"

BendingForce3dHeightWithGhostNodes::BendingForce3dHeightWithGhostNodes()
    : AbstractForce<3>(),
      mBasementMembraneParameter(DOUBLE_UNSET),
	  mExponentParameter(DOUBLE_UNSET),
	  mHeightDependantCurvatureParameter(1.0),
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

BendingForce3dHeightWithGhostNodes::~BendingForce3dHeightWithGhostNodes()
{
    // Avoid memory leaks
    if (mpExtendedMesh != NULL)
    {
        delete mpExtendedMesh;
    }
    
}

void BendingForce3dHeightWithGhostNodes::SetBasementMembraneParameter(double basementMembraneParameter)
{
	mBasementMembraneParameter = basementMembraneParameter;
}

double BendingForce3dHeightWithGhostNodes::GetBasementMembraneParameter()
{
	return mBasementMembraneParameter;
}

void BendingForce3dHeightWithGhostNodes::SetExponentParameter(double exponentParameter)
{
	mExponentParameter = exponentParameter;
}

double BendingForce3dHeightWithGhostNodes::GetExponentParameter()
{
	return mExponentParameter;
}

void BendingForce3dHeightWithGhostNodes::SetHeightDependantCurvatureParameter(double heightparameter)
{
	mHeightDependantCurvatureParameter = heightparameter;
}

double BendingForce3dHeightWithGhostNodes::GetHeightDependantCurvatureParameter()
{
	return mHeightDependantCurvatureParameter;
}

void BendingForce3dHeightWithGhostNodes::SetCircularNonZeroTargetCurvatureRegion(bool setNonZeroTargetCurvatureRegion, double nonZeroTargetCurvature,
		double radiusOfNonZeroTargetCurvatureRegion, double xCoordinateOfCentreOfNonZeroTargetCurvatureRegion,
		double yCoordinateOfCentreOfNonZeroTargetCurvatureRegion)
{
	mSetNonZeroTargetCurvatureRegion = setNonZeroTargetCurvatureRegion;
	mNonZeroTargetCurvature = nonZeroTargetCurvature;
	mRadiusOfNonZeroTargetCurvatureRegion = radiusOfNonZeroTargetCurvatureRegion;
	mXCoordinateOfCentreOfNonZeroTargetCurvatureRegion = xCoordinateOfCentreOfNonZeroTargetCurvatureRegion;
	mYCoordinateOfCentreOfNonZeroTargetCurvatureRegion = yCoordinateOfCentreOfNonZeroTargetCurvatureRegion;
}

bool BendingForce3dHeightWithGhostNodes::GetWhetherToSetNonZeroTargetCurvatureRegion()
{
	return mSetNonZeroTargetCurvatureRegion;
}

double BendingForce3dHeightWithGhostNodes::GetNonZeroTargetCurvature()
{
	return mNonZeroTargetCurvature;
}

double BendingForce3dHeightWithGhostNodes::GetRadiusOfNonZeroTargetCurvatureRegion()
{
	return mRadiusOfNonZeroTargetCurvatureRegion;
}

double BendingForce3dHeightWithGhostNodes::GetXCoordinateOfCentreOfNonZeroTargetCurvatureRegion()
{
	return mXCoordinateOfCentreOfNonZeroTargetCurvatureRegion;
}

double BendingForce3dHeightWithGhostNodes::GetYCoordinateOfCentreOfNonZeroTargetCurvatureRegion()
{
	return mYCoordinateOfCentreOfNonZeroTargetCurvatureRegion;
}

/*
 * Remove repeated entries in a 1D vector
 */
void BendingForce3dHeightWithGhostNodes::RemoveDuplicates1D(std::vector<unsigned>& rVectorWithDuplicates)
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
std::vector<c_vector<unsigned, 2> > BendingForce3dHeightWithGhostNodes::GetEpithelialTissuePairs(AbstractCellPopulation<3>& rCellPopulation)
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

std::vector<c_vector<unsigned, 10> > BendingForce3dHeightWithGhostNodes::GetEpithelialNeighbours(std::vector<c_vector<unsigned, 3> > rEpithelialMeshVector, int number_of_cells)
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
		c_vector<unsigned, 10> another_holder_neighbour = zero_vector<double>(10);
		for(unsigned j=0; j<10; j++)
		{
			holder_neighbour.push_back(temp_neighbours[i][j]);

		}
		
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



std::vector<c_vector<unsigned, 3> > BendingForce3dHeightWithGhostNodes::GetEpithelialMesh(AbstractCellPopulation<3>& rCellPopulation)
{
// Get a pointer to this node in mpExtendedMesh
std::vector<c_vector<unsigned, 3> > epithelial_triangulation;


DomMeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*>(&rCellPopulation);
for (unsigned elem_index=0; elem_index<mpExtendedMesh->GetNumElements(); elem_index++) 
{ 

    // Get a pointer to the element
    Element<3,3>* p_element = mpExtendedMesh->GetElement(elem_index);
        
    int Node0Index = p_element->GetNodeGlobalIndex(0);
    int Node1Index = p_element->GetNodeGlobalIndex(1);
    int Node2Index = p_element->GetNodeGlobalIndex(2);
    int Node3Index = p_element->GetNodeGlobalIndex(3);
    int node_index[4] = {Node0Index, Node1Index, Node2Index, Node3Index};

	int node0GlobalIndex = mExtendedMeshNodeIndexMap[Node0Index];
	int node1GlobalIndex = mExtendedMeshNodeIndexMap[Node1Index];
	int node2GlobalIndex = mExtendedMeshNodeIndexMap[Node2Index];
	int node3GlobalIndex = mExtendedMeshNodeIndexMap[Node3Index];
	int node_global_intex[4] = {node0GlobalIndex, node1GlobalIndex, node2GlobalIndex, node3GlobalIndex};

	for(int j=0; j<4; j++)
    {
		c_vector<unsigned, 3> tri_el = zero_vector<double>(3);

        CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(node_global_intex[(0+j)%4]);
        boost::shared_ptr<AbstractCellMutationState> p_state = p_cell->GetMutationState();

        bool is_0_stromal_cell = (p_state->IsType<StromalCellMutationState>()==true);

        if(is_0_stromal_cell) 
        {
            CellPtr p_cell_1 = rCellPopulation.GetCellUsingLocationIndex(node_global_intex[(1+j)%4]);
            boost::shared_ptr<AbstractCellMutationState> p_state_1 = p_cell_1->GetMutationState();

            CellPtr p_cell_2 = rCellPopulation.GetCellUsingLocationIndex(node_global_intex[(2+j)%4]);
            boost::shared_ptr<AbstractCellMutationState> p_state_2 = p_cell_2->GetMutationState();
                    
            CellPtr p_cell_3 = rCellPopulation.GetCellUsingLocationIndex(node_global_intex[(3+j)%4]);
            boost::shared_ptr<AbstractCellMutationState> p_state_3 = p_cell_3->GetMutationState();
                    
            int number_of_epithelial_cell = (p_state_1->IsType<WildTypeCellMutationState>()==true) + (p_state_2->IsType<WildTypeCellMutationState>()==true) + (p_state_3->IsType<WildTypeCellMutationState>()==true);

            if (number_of_epithelial_cell == 3)
            {
				tri_el[0] = node_index[(1+j)%4];
				tri_el[1] = node_index[(2+j)%4];
				tri_el[2] = node_index[(3+j)%4];

				epithelial_triangulation.push_back(tri_el);
            }
                    
        }
                
                
    }


}

	return epithelial_triangulation;
}

std::vector<c_vector<double, 3> > BendingForce3dHeightWithGhostNodes::FitPlaneAndFindImage(AbstractCellPopulation<3>& rCellPopulation, std::vector<unsigned> second_order_neighs,  unsigned cell_i)
{
	c_vector<double, 3> normal_vector;
	
	std::vector<c_vector<double, 3> > image_location_per_second_order_neighbours;

	c_vector<double, 9> ATA = zero_vector<double>(9);
	c_vector<double, 3> ATz = zero_vector<double>(3);
	c_vector<double, 9> iATA = zero_vector<double>(9);

	for (unsigned i=0; i<second_order_neighs.size(); i++)
	{
		unsigned cell_i_ext =second_order_neighs[i];
		c_vector<double, 3> epithelial_location = this->mpExtendedMesh->GetNode(cell_i_ext)->rGetLocation();


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

	if ( (ATA_det >= pow(10,-10)) || (ATA_det <= -1.0*pow(10,-10)) )
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
	else if ( (ATA_det <= pow(10,-10)) && (ATA_det >= -1.0*pow(10,-10)) )
	{
		normal_vector[0] = 0;
		normal_vector[1] = 1;
		normal_vector[2] = 0;
		//std::cout << "case 2" << "\n";


		// THIS NEEDS TO BE DONE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// c_vector<double, 9> ATA2;
		// c_vector<double, 3> ATy;
		// c_vector<double, 9> iATA2;

		// for (unsigned i=0; i<second_order_neighs.size(); i++)
		// {
		// //	unsigned cell_i_ext = mExtendedMeshNodeIndexMap[second_order_neighs[i]];
		// 	unsigned cell_i_ext =second_order_neighs[i];
		// 	//int epithelialNodeIndex = second_order_neighs[i];
		// 	// c_vector<double, 3> epithelial_location = mpExtendedMesh->GetNode(cell_i_ext)->rGetLocation();
		// 	c_vector<double, 3> epithelial_location = this->mpExtendedMesh->GetNode(cell_i_ext)->rGetLocation();


		// 	// Find Transpose(A)*A
		// 	ATA2[0] += epithelial_location[0]*epithelial_location[0];
		// 	ATA2[1] += epithelial_location[0]*epithelial_location[2];
		// 	ATA2[2] += epithelial_location[0];

		// 	ATA2[3] += epithelial_location[0]*epithelial_location[2];
		// 	ATA2[4] += epithelial_location[2]*epithelial_location[2];
		// 	ATA2[5] += epithelial_location[2];

		// 	ATA2[6] += epithelial_location[0];
		// 	ATA2[7] += epithelial_location[2];
		// 	ATA2[8] += 1.0;

		// 	// Calculate Transpose(A)*z
		// 	ATy[0] += epithelial_location[0]*epithelial_location[1];
		// 	ATy[1] += epithelial_location[2]*epithelial_location[1];
		// 	ATy[2] += epithelial_location[1];
		// }

		// // determinant of Transpose(A)*A
		// double ATA_det2 = ATA2[0]*(ATA2[4]*ATA2[8]-ATA2[5]*ATA2[7]) - ATA2[1]*(ATA2[3]*ATA2[8]-ATA2[5]*ATA2[6]) + ATA2[2]*(ATA2[3]*ATA2[7]-ATA2[6]*ATA2[4]);

		// if ( ATA_det2 >= pow(10,-10) || ATA_det2 <= -1.0*pow(10,-10) )
		// {
		// 	// Calculate the inverse of Transpose(A)*A
		// 	iATA2[0] = (1.0/ATA_det2)*(ATA2[4]*ATA2[8] - ATA2[7]*ATA2[5]);
		// 	iATA2[1] = (1.0/ATA_det2)*(ATA2[2]*ATA2[7] - ATA2[8]*ATA2[1]);
		// 	iATA2[2] = (1.0/ATA_det2)*(ATA2[1]*ATA2[5] - ATA2[4]*ATA2[2]);

		// 	iATA2[3] = (1.0/ATA_det2)*(ATA2[5]*ATA2[6] - ATA2[8]*ATA2[3]);
		// 	iATA2[4] = (1.0/ATA_det2)*(ATA2[0]*ATA2[8] - ATA2[6]*ATA2[2]);
		// 	iATA2[5] = (1.0/ATA_det2)*(ATA2[2]*ATA2[3] - ATA2[5]*ATA2[0]);

		// 	iATA2[6] = (1.0/ATA_det2)*(ATA2[3]*ATA2[7] - ATA2[6]*ATA2[4]);
		// 	iATA2[7] = (1.0/ATA_det2)*(ATA2[1]*ATA2[6] - ATA2[7]*ATA2[0]);
		// 	iATA2[8] = (1.0/ATA_det2)*(ATA2[0]*ATA2[4] - ATA2[3]*ATA2[1]);
			
		// 	// Calculate  normal = inverse(Transpose(A)*A)*(Transpose(A)*z)
		// 	normal_vector[0] = iATA2[0]*ATy[0] + iATA2[1]*ATy[1] + iATA2[2]*ATy[2];
		// 	normal_vector[1] = iATA2[3]*ATy[0] + iATA2[4]*ATy[1] + iATA2[5]*ATy[2];
		// 	normal_vector[2] = iATA2[6]*ATy[0] + iATA2[7]*ATy[1] + iATA2[8]*ATy[2];
		// 	//std::cout << "case 1" << "\n"; 

		// 	normal_vector[0] = -1.0*normal_vector[0];
		// 	normal_vector[1] = -1.0*normal_vector[1];
		// }

		// else if ( ATA_det2 <= pow(10,-10) & ATA_det2 >= -1.0*pow(10,-10) )
		// {
		// 	c_vector<double, 9> ATA3;
		// 	c_vector<double, 3> ATx;
		// 	c_vector<double, 9> iATA3;

		// 	for (unsigned i=0; i<second_order_neighs.size(); i++)
		// 	{
		// 	//	unsigned cell_i_ext = mExtendedMeshNodeIndexMap[second_order_neighs[i]];
		// 		unsigned cell_i_ext =second_order_neighs[i];
		// 		//int epithelialNodeIndex = second_order_neighs[i];
		// 		// c_vector<double, 3> epithelial_location = mpExtendedMesh->GetNode(cell_i_ext)->rGetLocation();
		// 		c_vector<double, 3> epithelial_location = this->mpExtendedMesh->GetNode(cell_i_ext)->rGetLocation();


		// 		// Find Transpose(A)*A
		// 		ATA3[0] += epithelial_location[1]*epithelial_location[1];
		// 		ATA3[1] += epithelial_location[1]*epithelial_location[2];
		// 		ATA3[2] += epithelial_location[1];

		// 		ATA3[3] += epithelial_location[1]*epithelial_location[2];
		// 		ATA3[4] += epithelial_location[2]*epithelial_location[2];
		// 		ATA3[5] += epithelial_location[2];

		// 		ATA3[6] += epithelial_location[1];
		// 		ATA3[7] += epithelial_location[2];
		// 		ATA3[8] += 1.0;

		// 		// Calculate Transpose(A)*z
		// 		ATx[0] += epithelial_location[1]*epithelial_location[0];
		// 		ATx[1] += epithelial_location[2]*epithelial_location[0];
		// 		ATx[2] += epithelial_location[0];
		// 	}

		// 	// determinant of Transpose(A)*A
		// 	double ATA_det3 = ATA3[0]*(ATA3[4]*ATA3[8]-ATA3[5]*ATA3[7]) - ATA3[1]*(ATA3[3]*ATA3[8]-ATA3[5]*ATA3[6]) + ATA3[2]*(ATA3[3]*ATA3[7]-ATA3[6]*ATA3[4]);

		// 	if ( ATA_det3 >= pow(10,-10) || ATA_det3 <= -1.0*pow(10,-10) )
		// 	{
		// 		// Calculate the inverse of Transpose(A)*A
		// 		iATA3[0] = (1.0/ATA_det3)*(ATA3[4]*ATA3[8] - ATA3[7]*ATA3[5]);
		// 		iATA3[1] = (1.0/ATA_det3)*(ATA3[2]*ATA3[7] - ATA3[8]*ATA3[1]);
		// 		iATA3[2] = (1.0/ATA_det3)*(ATA3[1]*ATA3[5] - ATA3[4]*ATA3[2]);

		// 		iATA3[3] = (1.0/ATA_det3)*(ATA3[5]*ATA3[6] - ATA3[8]*ATA3[3]);
		// 		iATA3[4] = (1.0/ATA_det3)*(ATA3[0]*ATA3[8] - ATA3[6]*ATA3[2]);
		// 		iATA3[5] = (1.0/ATA_det3)*(ATA3[2]*ATA3[3] - ATA3[5]*ATA3[0]);

		// 		iATA3[6] = (1.0/ATA_det3)*(ATA3[3]*ATA3[7] - ATA3[6]*ATA3[4]);
		// 		iATA3[7] = (1.0/ATA_det3)*(ATA3[1]*ATA3[6] - ATA3[7]*ATA3[0]);
		// 		iATA3[8] = (1.0/ATA_det3)*(ATA3[0]*ATA3[4] - ATA3[3]*ATA3[1]);
				
		// 		// Calculate  normal = inverse(Transpose(A)*A)*(Transpose(A)*z)
		// 		normal_vector[0] = iATA3[0]*ATx[0] + iATA3[1]*ATx[1] + iATA3[2]*ATx[2];
		// 		normal_vector[1] = iATA3[3]*ATx[0] + iATA3[4]*ATx[1] + iATA3[5]*ATx[2];
		// 		normal_vector[2] = iATA3[6]*ATx[0] + iATA3[7]*ATx[1] + iATA3[8]*ATx[2];
		// 		//std::cout << "case 1" << "\n"; 

		// 		normal_vector[0] = -1.0*normal_vector[0];
		// 		normal_vector[1] = -1.0*normal_vector[1];
		// 	}
		// }

	}

	// PRINT_VECTOR(normal_vector);

	double sign_of_curvature = (mNonZeroTargetCurvature > 0) - (mNonZeroTargetCurvature < 0);

	c_vector<double, 3> unit_normal = zero_vector<double>(3);
	unit_normal[0] = sign_of_curvature*normal_vector[0]/(sqrt(pow(normal_vector[0],2)+pow(normal_vector[1],2)+1));
	unit_normal[1] = sign_of_curvature*normal_vector[1]/(sqrt(pow(normal_vector[0],2)+pow(normal_vector[1],2)+1));
	unit_normal[2] = sign_of_curvature*1.0/(sqrt(pow(normal_vector[0],2)+pow(normal_vector[1],2)+1));
	// PRINT_VECTOR(unit_normal);

	double radius_from_curvature = 1/sqrt(mNonZeroTargetCurvature*sign_of_curvature);

	c_vector<double, 3> circle_centre = zero_vector<double>(3);
	unsigned cell_ii_ext = mExtendedMeshNodeIndexMap[cell_i];
	
	// c_vector<double, 3> epithelial_cell_i = mpExtendedMesh->GetNode(cell_ii_ext)->rGetLocation();
	c_vector<double, 3> epithelial_cell_i = this->mpExtendedMesh->GetNode(cell_ii_ext)->rGetLocation();


	circle_centre[0]  = radius_from_curvature*unit_normal[0] + epithelial_cell_i[0];
	circle_centre[1]  = radius_from_curvature*unit_normal[1] + epithelial_cell_i[1];
	circle_centre[2]  = radius_from_curvature*unit_normal[2] + epithelial_cell_i[2];
	// PRINT_VECTOR(circle_centre);
	// PRINT_VECTOR(epithelial_cell_i);
	
	double height_of_plane = normal_vector[0]*epithelial_cell_i[0] + normal_vector[1]*epithelial_cell_i[1] + epithelial_cell_i[2];

	//  Note: when we use unit_normal in this loop, we may need to be using normal_vector instead...?
	for (unsigned i=0; i<second_order_neighs.size(); i++)
	{
		//unsigned cell_i_ext = mExtendedMeshNodeIndexMap[second_order_neighs[i]];
		unsigned cell_i_ext = second_order_neighs[i];
		//int epithelialNodeIndex = second_order_neighs[i];
		// c_vector<double, 3> epithelial_location = mpExtendedMesh->GetNode(cell_i_ext)->rGetLocation();
		c_vector<double, 3> epithelial_location = this->mpExtendedMesh->GetNode(cell_i_ext)->rGetLocation();

		
		//PRINT_VECTOR(unit_normal);
		//std::cout << pow(unit_normal[0],2) + pow(unit_normal[1],2) + pow(unit_normal[2],2) << "\n";
		double alpha = pow(normal_vector[0],2) + pow(normal_vector[1],2) + 1.0;
		double beta  = 2*normal_vector[0]*(epithelial_location[0]-circle_centre[0]) + 2*normal_vector[1]*(epithelial_location[1]-circle_centre[1]) + 2*1.0*(epithelial_location[2]-circle_centre[2]);
		double delta = pow(epithelial_location[0]-circle_centre[0],2) + pow(epithelial_location[1]-circle_centre[1],2) + pow(epithelial_location[2]-circle_centre[2],2) - pow(radius_from_curvature,2);

		// Original!!!!!
		//
		// double param_T = (-1.0*beta - sign_of_curvature*sqrt(pow(beta,2) - 4*alpha*delta))/(2.0*alpha);
		//
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
		// PRINT_VECTOR(point_on_sphere);
		
		

		double k_point = (normal_vector[0]*point_on_sphere[0] + normal_vector[1]*point_on_sphere[1] + 1.0*point_on_sphere[2] - normal_vector[2])/alpha;
		
		c_vector<double, 3> point_on_plane;
		point_on_plane[0] = point_on_sphere[0] - k_point*normal_vector[0];
		point_on_plane[1] = point_on_sphere[1] - k_point*normal_vector[1];
		// point_on_plane[2] = point_on_sphere[2] - k_point*1.0;
		point_on_plane[2] = (-1.0*normal_vector[0])*point_on_plane[0] + (-1.0*normal_vector[1])*point_on_plane[1] + height_of_plane;

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
		// PRINT_VECTOR(temp_cell_location);
		
		image_location_per_second_order_neighbours.push_back(temp_cell_location);
	}
	// std::cout<<"_1_\n";
	return image_location_per_second_order_neighbours;
}


c_vector<double, 3> BendingForce3dHeightWithGhostNodes::GetForceDueToDiscreteCurvature(AbstractCellPopulation<3>& rCellPopulation, std::vector<c_vector<unsigned, 3> > rEpithelialMeshVector, std::vector<unsigned> first_order_neighs,  std::vector<unsigned> second_order_neighs, std::vector<c_vector<double, 3> > image_location, unsigned cell_i, unsigned num_cells)
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
	c_vector<double, 3> force_due_to_curvature = zero_vector<double>(3);
	// force_due_to_curvature[0] = 0.0;
	// force_due_to_curvature[1] = 0.0;
	// force_due_to_curvature[2] = 0.0;

	int inverse_sec_neighs[4*num_cells];
	for(unsigned i=0; i<4*num_cells; i++)
	{
		inverse_sec_neighs[i] = 0;

	}

	for(unsigned i=0; i<second_order_neighs.size(); i++)
	{
		inverse_sec_neighs[second_order_neighs[i]] = i;

	}


	// Calculate angles about the first order neighbours in the mesh
	// std::vector<double> angle_sum;
	// std::vector<double> angle_sum_indicats;

	for(unsigned j=0; j<first_order_neighs.size(); j++)
	{
		double ang_sum = 0.0;
		c_vector<double, 3> grad_i_phi_j_el = zero_vector<double>(3);;
		// grad_i_phi_j_el[0] = 0.0;
		// grad_i_phi_j_el[1] = 0.0;
		// grad_i_phi_j_el[2] = 0.0;

		int cell_a = first_order_neighs[j];
		int inv_cell_a = inverse_sec_neighs[cell_a];
		c_vector<double, 3> a_loc, b_loc, c_loc;
		a_loc[0] = image_location[inv_cell_a][0];
		a_loc[1] = image_location[inv_cell_a][1];
		a_loc[2] = image_location[inv_cell_a][2];

		// if(cell_i == 154)
		// {
		// 	PRINT_VARIABLE(cell_a);
		// }

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
				b_loc[0] = image_location[inv_cell_b][0];
				b_loc[1] = image_location[inv_cell_b][1];
				b_loc[2] = image_location[inv_cell_b][2];

				inv_cell_c = inverse_sec_neighs[cell_c];
				c_loc[0] = image_location[inv_cell_c][0];
				c_loc[1] = image_location[inv_cell_c][1];
				c_loc[2] = image_location[inv_cell_c][2];

				c_vector<double, 3> vect_ab = b_loc - a_loc;
				c_vector<double, 3> vect_ac = c_loc - a_loc;

				// double mag_ab = sqrt(pow(vect_ab[0],2) + pow(vect_ab[1],2) + pow(vect_ab[2],2));
				// double mag_ac = sqrt(pow(vect_ac[0],2) + pow(vect_ac[1],2) + pow(vect_ac[2],2));
				double mag_ab = sqrt(pow(b_loc[0] - a_loc[0],2) + pow(b_loc[1] - a_loc[1],2) + pow(b_loc[2] - a_loc[2],2));
				double mag_ac = sqrt(pow(c_loc[0] - a_loc[0],2) + pow(c_loc[1] - a_loc[1],2) + pow(c_loc[2] - a_loc[2],2));

				double abDotac = (b_loc[0] - a_loc[0])*(c_loc[0] - a_loc[0]) + (b_loc[1] - a_loc[1])*(c_loc[1] - a_loc[1]) + (b_loc[2] - a_loc[2])*(c_loc[2] - a_loc[2]);

				// double cos_abc = (vect_ab[0]*vect_ac[0] + vect_ab[1]*vect_ac[1] + vect_ab[2]*vect_ac[2])/(mag_ab*mag_ac);
				double cos_abc = (abDotac)/(mag_ab*mag_ac);

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

				c_vector<double, 3> grad_hold = zero_vector<double>(3);;

				// Calculate grad_i phi_i
				if(cell_a == mCellPopulationNodeIndexMap[cell_i])
				{
					// grad_hold = ( -1.0/sqrt(1.0 - cos_abc*cos_abc) )*(1.0/(mag_ab*mag_ac))*(2.0*a_loc - b_loc - c_loc + (vect_ab[0]*vect_ac[0] + vect_ab[1]*vect_ac[1] + vect_ab[2]*vect_ac[2])*( (b_loc-a_loc)/pow(mag_ab,2) +(c_loc-a_loc)/pow(mag_ac,2) ) );
					grad_hold = ( -1.0/sqrt(1.0 - pow(cos_abc,2)) )*(1.0/(mag_ab*mag_ac))*(2.0*a_loc - b_loc - c_loc + (abDotac)*( (b_loc-a_loc)/pow(mag_ab,2) +(c_loc-a_loc)/pow(mag_ac,2) ) );


					grad_i_phi_j_el += grad_hold;

					// if(cell_i == 154)
					// {
					// 	std::cout<< " a {" << grad_hold[0] << ", " << grad_hold[1] << ", " << grad_hold[2] << "} \n";
					// 	// PRINT_VECTOR(grad_hold);
					// }

				}
				
				else if(cell_b == mCellPopulationNodeIndexMap[cell_i])
				{

					// grad_hold = (-1.0/(sqrt(1.0 - cos_abc*cos_abc)))*(1.0/(mag_ab*mag_ac))*( c_loc-a_loc + (b_loc-a_loc)*((vect_ab[0]*vect_ac[0] + vect_ab[1]*vect_ac[1] + vect_ab[2]*vect_ac[2]))/pow(mag_ab,2) ) ;
					grad_hold = (-1.0/(sqrt(1.0 - pow(cos_abc,2))))*(1.0/(mag_ab*mag_ac))*( (c_loc-a_loc) - (b_loc-a_loc)*((abDotac))/pow(mag_ab,2) ) ;

					grad_i_phi_j_el += grad_hold;

					// if(cell_i == 154)
					// {
					// 	std::cout<< " b {" << grad_hold[0] << ", " << grad_hold[1] << ", " << grad_hold[2] << "} \n";
					// 	// PRINT_VECTOR(grad_hold);
					// }

				}
				else if(cell_c == mCellPopulationNodeIndexMap[cell_i])
				{
					
					// grad_hold = (-1.0/(sqrt(1.0 - cos_abc*cos_abc)))*(1.0/(mag_ab*mag_ac))*( b_loc-a_loc + (c_loc-a_loc)*((vect_ab[0]*vect_ac[0] + vect_ab[1]*vect_ac[1] + vect_ab[2]*vect_ac[2]))/pow(mag_ac,2) ) ;
					grad_hold = (-1.0/(sqrt(1.0 -  pow(cos_abc,2))))*(1.0/(mag_ab*mag_ac))*( (b_loc-a_loc) - (c_loc-a_loc)*((abDotac))/pow(mag_ac,2) ) ;

					grad_i_phi_j_el += grad_hold;

					// if(cell_i == 154)
					// {
					// 	std::cout<< " c {" << grad_hold[0] << ", " << grad_hold[1] << ", " << grad_hold[2] << "} \n";
					// 	// PRINT_VECTOR(grad_hold);
					// }

				}
				// else
				// {
				// 	grad_hold[0] = 0.0;
				// 	grad_hold[1] = 0.0;
				// 	grad_hold[2] = 0.0;
				// 	grad_i_phi_j_el += grad_hold;
				// }

			}
			
		}


		// angle_sum.push_back(ang_sum);
		// angle_sum_indicats.push_back(cell_a);
		// if(cell_i == 154)
		// {
		// 	PRINT_VECTOR(grad_i_phi_j_el);
		// 	std::cout<<"\n";
		// 	// PRINT_VECTOR(a_loc);
		// }

		double target_angle  = 0.0;
		if(first_order_neighs.size() > 5)
		{
			target_angle = 2*3.141592653589793;
		}
		else if (first_order_neighs.size() > 3 && first_order_neighs.size() < 6)
		{
			target_angle = 3.141592653589793;
		}
		else
		{
			target_angle = 0.5*3.141592653589793;
		}

		double sign_of_anlge = (round((target_angle - ang_sum)*pow(10.0,14))/pow(10.0,14) >= 0) - (round((target_angle - ang_sum)*pow(10.0,14))/pow(10.0,14) < 0);
		c_vector<double, 3> force_due_to_curvature_i = 1.0*(exponent_parameter)*sign_of_anlge*pow(sqrt(pow( (target_angle - ang_sum) ,2)),(exponent_parameter-1.0))*grad_i_phi_j_el;

		// c_vector<double, 3> force_due_to_curvature_i = 1.0*(exponent_parameter)*( pow((2*3.141592653589793 - ang_sum),(exponent_parameter-1.0)) )*grad_i_phi_j_el;

		// c_vector<double, 3> force_due_to_curvature_i = 	(1.0/(1.0+exp(-1.0*exponent_parameter*(ang_sum - 2*3.141592653589793)) ) - 0.5 ) *grad_i_phi_j_el;

		force_due_to_curvature +=  force_due_to_curvature_i;

		// if(cell_a == 154)
		// {
		// 	PRINT_VARIABLE(ang_sum);
		// }

	}

	// if(cell_i == 154)
	// 	{
	// 		// PRINT_VECTOR(force_due_to_curvature);
	// 		// PRINT_VECTOR(a_loc);
	// 	}
	// PRINT_VECTOR(force_due_to_curvature);
	// std::cout<<"_2_\n";

	
	return force_due_to_curvature;
}

// ORIGINAL
// c_vector<double, 3> BendingForce3dHeightWithGhostNodes::GetForceDueToDiscreteCurvature(AbstractCellPopulation<3>& rCellPopulation, std::vector<c_vector<unsigned, 3> > rEpithelialMeshVector, std::vector<unsigned> first_order_neighs,  std::vector<unsigned> second_order_neighs, std::vector<c_vector<double, 3> > image_location, unsigned cell_i, unsigned num_cells)
// {
	
// 	/* 
// 	Need:
// 	Traingulation		    ->  rEpithelialMeshVector
// 	First order neighbours  ->  first_order_neighs 
// 	image node positions    ->  image_location 
// 	  needs (for node IDs)  ->  second_order_neighs 
// 	Cell_i					->  cell_i
// 	*/
// 	//std::cout<<"\n\n";
// 	double exponent_parameter = GetExponentParameter();
// 	c_vector<double, 3> force_due_to_curvature;
// 	force_due_to_curvature[0] = 0.0;
// 	force_due_to_curvature[1] = 0.0;
// 	force_due_to_curvature[2] = 0.0;
// 	//PRINT_VECTOR(force_due_to_curvature);

// 	int inverse_sec_neighs[4*num_cells];
// 	for(unsigned i=0; i<4*num_cells; i++)
// 	{
// 		inverse_sec_neighs[i] = 0;

// 		//std::cout << "i " << i << " s " << second_order_neighs[i] << "\n";
// 	}

// 	for(unsigned i=0; i<second_order_neighs.size(); i++)
// 	{
// 		inverse_sec_neighs[second_order_neighs[i]] = i;

// 		//std::cout << "i " << i << " s " << second_order_neighs[i] << "\n";
// 	}


// 	// Calculate angles about the first order neighbours in the mesh
// 	std::vector<double> angle_sum;
// 	std::vector<double> angle_sum_indicats;
// 	// std::cout << "\n\n";
// 	for(unsigned j=0; j<first_order_neighs.size(); j++)
// 	{
// 		double ang_sum = 0.0;
// 		c_vector<double, 3> grad_i_phi_j_el;
// 		grad_i_phi_j_el[0] = 0.0;
// 		grad_i_phi_j_el[1] = 0.0;
// 		grad_i_phi_j_el[2] = 0.0;

// 		int cell_a = first_order_neighs[j];
// 		int inv_cell_a = inverse_sec_neighs[cell_a];
// 		c_vector<double, 3> cell_a_loc, cell_b_loc, cell_c_loc;
// 		cell_a_loc[0] = image_location[inv_cell_a][0];
// 		cell_a_loc[1] = image_location[inv_cell_a][1];
// 		cell_a_loc[2] = image_location[inv_cell_a][2];

		

// 		for(unsigned i=0; i<rEpithelialMeshVector.size(); i++)
// 		{
// 			int cell_b, cell_c, inv_cell_b, inv_cell_c;
// 			bool have_triangle = false;

// 			if(rEpithelialMeshVector[i][0] == cell_a)
// 			{
// 				cell_b = rEpithelialMeshVector[i][1];
// 				cell_c = rEpithelialMeshVector[i][2];
// 				have_triangle = true;
// 			}
// 			else if(rEpithelialMeshVector[i][1] == cell_a)
// 			{
// 				cell_b = rEpithelialMeshVector[i][0];
// 				cell_c = rEpithelialMeshVector[i][2];
// 				have_triangle = true;
// 			}
// 			else if(rEpithelialMeshVector[i][2] == cell_a)
// 			{
// 				cell_b = rEpithelialMeshVector[i][0];
// 				cell_c = rEpithelialMeshVector[i][1];
// 				have_triangle = true;
// 			}
// 			// have_triangle = bool ( (cell_a == cell_i) + (cell_b == cell_i) + (cell_c == cell_i) );

// 			if(have_triangle)
// 			{
// 				inv_cell_b = inverse_sec_neighs[cell_b];
// 				cell_b_loc[0] = image_location[inv_cell_b][0];
// 				cell_b_loc[1] = image_location[inv_cell_b][1];
// 				cell_b_loc[2] = image_location[inv_cell_b][2];

// 				inv_cell_c = inverse_sec_neighs[cell_c];
// 				cell_c_loc[0] = image_location[inv_cell_c][0];
// 				cell_c_loc[1] = image_location[inv_cell_c][1];
// 				cell_c_loc[2] = image_location[inv_cell_c][2];

// 				c_vector<double, 3> vect_ab = cell_b_loc - cell_a_loc;
// 				c_vector<double, 3> vect_ac = cell_c_loc - cell_a_loc;

// 				double mag_ab = sqrt(pow(vect_ab[0],2) + pow(vect_ab[1],2) + pow(vect_ab[2],2));
// 				double mag_ac = sqrt(pow(vect_ac[0],2) + pow(vect_ac[1],2) + pow(vect_ac[2],2));

// 				double cos_abc = (vect_ab[0]*vect_ac[0] + vect_ab[1]*vect_ac[1] + vect_ab[2]*vect_ac[2])/(mag_ab*mag_ac);

// 				// Think these edge cases, they shouldn't be happening, but they do...
// 				if (!isnan(ang_sum))
// 				{
// 					if (cos_abc >= 0.99999999 || cos_abc == 1)
// 					{
// 						//PRINT_VARIABLE("+1");
// 						ang_sum += 0.0;
// 						cos_abc = 0.9999;
// 					}
// 					else if (cos_abc <= -0.99999999 || cos_abc == -1)
// 					{
// 						//PRINT_VARIABLE("-1");
// 						ang_sum += 3.141592653589793;
// 						cos_abc = -0.9999;
// 					}
// 					else
// 					{
// 						ang_sum += acos(cos_abc);
// 					}
// 				}
// 				else if (isnan(ang_sum))
// 				{
// 					cos_abc = 0;
// 				}

// 				c_vector<double, 3> grad_hold;

// 				// Calculate grad_i phi_i
// 				if(cell_a == mCellPopulationNodeIndexMap[cell_i])
// 				{
// 					// TRACE("one");
// 					//PRINT_VECTOR(cell_a_loc);
// 					//PRINT_VECTOR(cell_b_loc);
// 					//PRINT_VECTOR(cell_c_loc);
// 					grad_hold = (1.0/(sqrt(1.0 - cos_abc*cos_abc)))*(1.0/(mag_ab*mag_ab*mag_ac*mag_ac))*(mag_ab*(mag_ac-mag_ab*cos_abc)*vect_ac + mag_ac*(mag_ab-mag_ac*cos_abc)*vect_ab);

// 					grad_i_phi_j_el += grad_hold;
					
// 					//PRINT_VECTOR(grad_hold);
// 					//std::cout << "case 1 " << cell_a << "\n";
// 				}
				
// 				else if(cell_b == mCellPopulationNodeIndexMap[cell_i])
// 				{
// 					// TRACE("two");
// 					//PRINT_VECTOR(cell_a_loc);
// 					//PRINT_VECTOR(cell_b_loc);
// 					//PRINT_VECTOR(cell_c_loc);

// 					grad_hold = (1.0/(sqrt(1.0 - cos_abc*cos_abc))) * (vect_ac/(mag_ab*mag_ac) - cos_abc*vect_ab/(mag_ab*mag_ab));
// 					//grad_hold = (1.0/(sqrt(1.0 - cos_abc*cos_abc)))* (1.0/(mag_ac)) * (vect_ab/mag_ab - cos_abc*vect_ac/mag_ac);
// 					grad_i_phi_j_el += grad_hold;
					
// 					//PRINT_VECTOR(grad_hold);
// 					//std::cout << "case 2 " << cell_a << "\n";
// 				}
// 				else if(cell_c == mCellPopulationNodeIndexMap[cell_i])
// 				{
// 					// TRACE("three");
// 					//PRINT_VECTOR(cell_a_loc);
// 					//PRINT_VECTOR(cell_b_loc);
// 					//PRINT_VECTOR(cell_c_loc);
					
// 					grad_hold = (1.0/(sqrt(1.0 - cos_abc*cos_abc))) * (vect_ab/(mag_ab*mag_ac) - cos_abc*vect_ac/(mag_ac*mag_ac));
// 					//grad_hold = (1.0/(sqrt(1.0 - cos_abc*cos_abc)))* (1.0/(mag_ab)) * (vect_ac/mag_ac - cos_abc*vect_ab/mag_ab);
// 					grad_i_phi_j_el += grad_hold;

// 					//PRINT_VECTOR(grad_hold);
// 					//std::cout << "case 3 " << cell_a << "\n";
// 				}
// 				else
// 				{
// 					grad_hold[0] = 0.0;
// 					grad_hold[1] = 0.0;
// 					grad_hold[2] = 0.0;
// 					grad_i_phi_j_el += grad_hold;
// 				}

// 				// if (isnan(grad_i_phi_j_el[0]) || isnan(grad_i_phi_j_el[1]) || isnan(grad_i_phi_j_el[2]))
// 				// {
// 				// 	PRINT_VECTOR(grad_i_phi_j_el);
// 				// }

// 				// if ( cell_i == 154 )
// 				// {
// 				// 	PRINT_3_VARIABLES(cell_a,cell_b,cell_c);
// 				// 	PRINT_VECTOR(grad_hold);
// 				// }
// 			}
			
// 		}
		
// 		//if(cell_a == cell_i)
// 		//{
// 		//	PRINT_VECTOR(grad_i_phi_j_el);
// 		//}
// 		//std::cout << "Angle Sum = " << ang_sum << "\n";
// 		//PRINT_VARIABLE(ang_sum);

// 		angle_sum.push_back(ang_sum);
// 		angle_sum_indicats.push_back(cell_a);

// 		//PRINT_VECTOR(force_due_to_curvature);
// 		// Using alpha = 1.5 and beta = 2.0 for now, will fix this later though
// 		double sign_of_anlge = (round((2*3.141592653589793 - ang_sum)*pow(10.0,14))/pow(10.0,14) >= 0) - (round((2*3.141592653589793 - ang_sum)*pow(10.0,14))/pow(10.0,14) < 0);
// 		// if ( cell_i == 154 )
// 		// {
// 		// 	double hold = (exponent_parameter)*(sign_of_anlge)*pow(sqrt((2*3.141592653589793 - ang_sum)*(2*3.141592653589793 - ang_sum)),(exponent_parameter-1.0));
// 		// 	PRINT_2_VARIABLES(sign_of_anlge, hold);
// 		// }
// 		// c_vector<double, 3> force_due_to_curvature_i = (exponent_parameter)*pow((2*3.141592653589793 - ang_sum),(exponent_parameter-1.0))*grad_i_phi_j_el;// +
// 		// c_vector<double, 3> force_due_to_curvature_i = (exponent_parameter)*(sign_of_anlge)*pow(sqrt((2*3.141592653589793 - ang_sum)*(2*3.141592653589793 - ang_sum)),(exponent_parameter-1.0))*grad_i_phi_j_el;// +
// 														//  (1.1)*(sign_of_anlge)*pow(sqrt((2*3.141592653589793 - ang_sum)*(2*3.141592653589793 - ang_sum)),(1.1-1.0))*grad_i_phi_j_el + 
// 														//  (2.0)*(sign_of_anlge)*pow(sqrt((2*3.141592653589793 - ang_sum)*(2*3.141592653589793 - ang_sum)),(2.0-1.0))*grad_i_phi_j_el;

// 		// c_vector<double, 3> force_due_to_curvature_i = 1.0*(exponent_parameter)*( pow((2*3.141592653589793 - ang_sum),(exponent_parameter-1.0))*(1.0+pow((2*3.141592653589793 - ang_sum),(exponent_parameter)))*exp(pow((2*3.141592653589793 - ang_sum),(exponent_parameter)))  )*grad_i_phi_j_el;
		
// 		// c_vector<double, 3> force_due_to_curvature_i = 1.0*(exponent_parameter)*( (1.0+exponent_parameter*pow((2*3.141592653589793 - ang_sum),(exponent_parameter)))*exp(pow((2*3.141592653589793 - ang_sum),(exponent_parameter)))  )*grad_i_phi_j_el;
		
// 		// c_vector<double, 3> force_due_to_curvature_i = 1.0*(exponent_parameter)*exp(pow((2*3.141592653589793 - ang_sum),(exponent_parameter))) * (1.0+pow((2*3.141592653589793 - ang_sum),(exponent_parameter))) * pow((2*3.141592653589793 - ang_sum),(exponent_parameter-1.0)) *grad_i_phi_j_el;
		
// 		c_vector<double, 3> force_due_to_curvature_i = 1.0*(exponent_parameter)*( pow((2*3.141592653589793 - ang_sum),(exponent_parameter-1.0)) )*grad_i_phi_j_el;
// 		// c_vector<double, 3> force_due_to_curvature_i = sign_of_anlge*1.0*(exponent_parameter)*( pow((2*3.141592653589793 - ang_sum),(exponent_parameter-1.0)) )*grad_i_phi_j_el;
		
// 		force_due_to_curvature +=  force_due_to_curvature_i;

// 		//PRINT_VECTOR(force_due_to_curvature_i);
		
		

// 		//std::cout << round((2*3.141592653589793 - ang_sum)*pow(10.0,14))/pow(10.0,14) << "\n";
// 		//std::cout<< "sign_of_anlge = " << sign_of_anlge << "\n";
// 	}
	

// 	//PRINT_VECTOR(force_due_to_curvature);

// 	//std::cout << "cell_i = " << cell_i << "\n"; 
// 	//PRINT_VECTOR(angle_sum);
// 	//PRINT_VECTOR(angle_sum_indicats);
// 	//PRINT_VECTOR(first_order_neighs);
// 	//PRINT_VECTOR(force_due_to_curvature);
// //	std::cout << "\n";
	
	
	
// 	return force_due_to_curvature;
// }

// c_vector<unsigned, 2> BendingForce3dHeightWithGhostNodes::GetMaxMinEpithelialCells(AbstractCellPopulation<3>& rCellPopulation)
// {
// 	unsigned num_cells = rCellPopulation.GetNumRealCells();
// 	DomMeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*>(&rCellPopulation);

// 	c_vector<double, 2> max_min_ep_cells;

// 	double max_height = -10.0;
// 	double min_height = 10.0;

// 	for(unsigned i=0; i<num_cells; i++)
// 	{

//         unsigned cell_i_ext = mExtendedMeshNodeIndexMap[i];

// 		CellPtr p_cell_i_ext = rCellPopulation.GetCellUsingLocationIndex(cell_i_ext);
//     	bool cell_is_dead = p_cell_i_ext->IsDead();
// 		bool is_ghost = p_tissue->IsGhostNode(cell_i_ext);

		

//         if((p_cell_i_ext->GetMutationState()->IsType<WildTypeCellMutationState>() == true) && 
// 		   (p_cell_i_ext->GetMutationState()->IsType<StromalCellMutationState>() == false) && 
// 		   cell_is_dead == false && is_ghost == false)
// 		{
// 			// PRINT_VARIABLE(cell_i_ext);
// 			// c_vector<double, 3> epithelial_location = mpExtendedMesh->GetNode(cell_i_ext)->rGetLocation();
// 			c_vector<double, 3> epithelial_location = this->mpExtendedMesh->GetNode(cell_i_ext)->rGetLocation();


// 			if (epithelial_location[2] > max_height)
// 			{
// 				max_height = epithelial_location[2];
// 			}

// 			if (epithelial_location[2] < min_height)
// 			{
// 				min_height = epithelial_location[2];
// 			}

// 		}
		
// 	}
// 	// std::cout<<"\n\n";
// 	max_min_ep_cells[0] = max_height;
// 	max_min_ep_cells[1] = min_height;

// 	return max_min_ep_cells;
// }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Dom - looks like some legacy code...
/*
 * Method to determine whether an element contains ghost nodes
 */
bool BendingForce3dHeightWithGhostNodes::DoesElementContainGhostNodes(AbstractCellPopulation<3>& rCellPopulation, unsigned elementIndex)
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
bool BendingForce3dHeightWithGhostNodes::DoesElementContainLongEdge(AbstractCellPopulation<3>& rCellPopulation, unsigned elementIndex, double maxEdgeLength)
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
void BendingForce3dHeightWithGhostNodes::AddForceContribution(AbstractCellPopulation<3>& rCellPopulation)
{
	DomMeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*>(&rCellPopulation);

	unsigned num_cells = rCellPopulation.GetNumRealCells();
	unsigned num_nodes = rCellPopulation.GetNumNodes();
	

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

	mCellPopulationNodeIndexMap.clear();

	for(int iter; iter<num_nodes; iter++)
	{
		mCellPopulationNodeIndexMap[iter] = -1;
	}

    // Create a vector of nodes for use in constructing mpExtendedMesh
    std::vector<Node<3>*> extended_nodes(num_cells);

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
		mCellPopulationNodeIndexMap[real_node_index] = count;
		

        count++;
    }


    // // First, extend the mesh in the x-direction
    // for (AbstractCellPopulation<3>::Iterator cell_iter = rCellPopulation.Begin();
    //      cell_iter != rCellPopulation.End();
    //      ++cell_iter)
    // {
    //     // First, create and store a copy of this real node and cell
    //     unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
    //     c_vector<double, 3> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

    //     // Create a copy of the node corresponding to this cell and store it
    //     Node<3>* p_real_node = new Node<3>(real_node_index, real_node_location);

    //     // Compute the location of the image node corresponding to this node
    //     c_vector<double,3> image_node_location = real_node_location;
    //     if (real_node_location[0] >= mPeriodicDomainWidth*0.5)
    //     {
    //         image_node_location[0] -= mPeriodicDomainWidth;
    //     }
    //     else if (real_node_location[0] <  mPeriodicDomainWidth*0.5)
    //     {
    //         image_node_location[0] += mPeriodicDomainWidth;
    //     }

    //     // Create a copy of the node corresponding to this cell, suitable translated, and store it
    //     Node<3>* p_image_node = new Node<3>(count, image_node_location);

    //     extended_nodes[count] = p_image_node;

    //     // Populate mExtendedMeshNodeIndexMap
    //     mExtendedMeshNodeIndexMap[count] = real_node_index;

    //     count++;
    // }

    // // Second, extend this extended mesh in the y-direction
    // // (We don't need to store the real nodes anymore)
    // for (AbstractCellPopulation<3>::Iterator cell_iter = rCellPopulation.Begin();
    //      cell_iter != rCellPopulation.End();
    //      ++cell_iter)
    // {
    //     // First, create and store a copy of this real node and cell
    //     unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
    //     c_vector<double, 3> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

    //     // Compute the location of the image node corresponding to this node
    //     c_vector<double,3> image_node_location = real_node_location;

    //     if (real_node_location[1] >= mPeriodicDomainDepth*0.5)
    //     {
    //         image_node_location[1] -= mPeriodicDomainDepth;
    //     }
    //     else if (real_node_location[1] <  mPeriodicDomainDepth*0.5)
    //     {
    //         image_node_location[1] += mPeriodicDomainDepth;
    //     }

    //     // Create a copy of the node corresponding to this cell, suitable translated, and store it
    //     Node<3>* p_image_node = new Node<3>(count, image_node_location);
    //     extended_nodes[count] = p_image_node;

    //     // Populate mExtendedMeshNodeIndexMap
    //     mExtendedMeshNodeIndexMap[count] = real_node_index;

    //     count++;
    // }

	// // Thirdly, extend this extended mesh so that we cover the corners too
    // // (We don't need to store the real nodes anymore)
    // for (AbstractCellPopulation<3>::Iterator cell_iter = rCellPopulation.Begin();
    //      cell_iter != rCellPopulation.End();
    //      ++cell_iter)
    // {
    //     // First, create and store a copy of this real node and cell
    //     unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
    //     c_vector<double, 3> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

    //     // Compute the location of the image node corresponding to this node
    //     c_vector<double,3> image_node_location = real_node_location;

    //     if (real_node_location[1] >= mPeriodicDomainDepth*0.5)
    //     {
	// 		if (real_node_location[0] >= mPeriodicDomainWidth*0.5)
	// 		{
	// 			image_node_location[0] -= mPeriodicDomainWidth;
	// 		}
	// 		else if (real_node_location[0] <  mPeriodicDomainWidth*0.5)
	// 		{
	// 			image_node_location[0] += mPeriodicDomainWidth;
	// 		}
    //         image_node_location[1] -= mPeriodicDomainDepth;
    //     }
    //     else if (real_node_location[1] <  mPeriodicDomainDepth*0.5)
    //     {
	// 		if (real_node_location[0] >= mPeriodicDomainWidth*0.5)
	// 		{
	// 			image_node_location[0] -= mPeriodicDomainWidth;
	// 		}
	// 		else if (real_node_location[0] <  mPeriodicDomainWidth*0.5)
	// 		{
	// 			image_node_location[0] += mPeriodicDomainWidth;
	// 		}
    //         image_node_location[1] += mPeriodicDomainDepth;
    //     }

    //     // Create a copy of the node corresponding to this cell, suitable translated, and store it
    //     Node<3>* p_image_node = new Node<3>(count, image_node_location);
    //     extended_nodes[count] = p_image_node;

    //     mExtendedMeshNodeIndexMap[count] = real_node_index;

    //     count++;
    // }

	

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

	std::vector<c_vector<unsigned, 3> > epithelial_triangulation = GetEpithelialMesh(rCellPopulation);

	std::vector<c_vector<unsigned, 10> > epithelial_neighbours = GetEpithelialNeighbours(epithelial_triangulation, num_cells);


	double get_basement_membrane_parameter = GetBasementMembraneParameter();

	// c_vector<unsigned, 2> CryptMaxMin = GetMaxMinEpithelialCells(rCellPopulation);

	////
	double max_height = -10.0;
	double min_height = 10.0;

	double tissue_max_height = -10.0;
	double tissue_min_height = 10.0;

	for (AbstractCellPopulation<3>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // First, create and store a copy of this real node and cell
        unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        c_vector<double, 3> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

		CellPtr p_cell_i_ext = rCellPopulation.GetCellUsingLocationIndex(real_node_index);

		bool is_wild = p_cell_i_ext->GetMutationState()->IsType<WildTypeCellMutationState>();
        bool is_stom = p_cell_i_ext->GetMutationState()->IsType<StromalCellMutationState>();

		if (is_wild && !is_stom)
		{
			if (real_node_location[2] > max_height)
			{
				max_height = real_node_location[2];
			}
			if (real_node_location[2] < min_height)
			{
				min_height = real_node_location[2];
			}
		}

	}
	
	double CryptHeight = max_height - min_height;
	double heightparameter_ratio = GetHeightDependantCurvatureParameter();
	double max_height_for_curvature = min_height + heightparameter_ratio*CryptHeight + 10.0;
	// std::cout<<"_____________________________________________________________\n\n";
	for (AbstractCellPopulation<3>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
		


        // First, create and store a copy of this real node and cell
        unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        c_vector<double, 3> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

		CellPtr p_cell_i_ext = rCellPopulation.GetCellUsingLocationIndex(real_node_index);



        unsigned cell_i_ext = real_node_index;

    	bool cell_is_dead = p_cell_i_ext->IsDead();
		bool is_ghost = p_tissue->IsGhostNode(cell_i_ext);
		
				
        if(p_cell_i_ext->GetMutationState()->IsType<WildTypeCellMutationState>() == true && cell_is_dead == false && is_ghost == false)
		{

			std::vector<unsigned> second_order_neighs;

			int first_order_neighs[10];
			
			for(unsigned j=0; j<10; j++)
			{
				unsigned extended_node_index = mCellPopulationNodeIndexMap[real_node_index];
				first_order_neighs[j] = epithelial_neighbours[extended_node_index][j];
						
				for(unsigned k=0; k<10; k++)
				{
					unsigned extended_node_index_j = first_order_neighs[j];
					second_order_neighs.push_back(epithelial_neighbours[extended_node_index_j][k]);

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

			c_vector<double, 3> epithelial_location = real_node_location;
			
			std::vector<c_vector<double, 3> > image_location_per_second_order_neighbours;

			bool epithelial_cell_in_disk_target_curvature_region = (pow((epithelial_location[0] - mXCoordinateOfCentreOfNonZeroTargetCurvatureRegion),2)
																  + pow((epithelial_location[1] - mYCoordinateOfCentreOfNonZeroTargetCurvatureRegion),2) 
																 <= pow(mRadiusOfNonZeroTargetCurvatureRegion,2) );

			bool epithelial_cell_in_height_target_curvature_region = (epithelial_location[2] < max_height_for_curvature);

			// here we will find
			double basement_membrane_parameter = 0.0;
			if(epithelial_cell_in_height_target_curvature_region == true && epithelial_cell_in_disk_target_curvature_region == true)
			{
				image_location_per_second_order_neighbours = FitPlaneAndFindImage(rCellPopulation, second_order_neighs, cell_i_ext);
				
				basement_membrane_parameter = 1.0*get_basement_membrane_parameter;

			}
			else 
			{
				// TRACE("no bending cell");
				basement_membrane_parameter = +1.0*get_basement_membrane_parameter;

				for (unsigned j=0; j<second_order_neighs.size(); j++)
				{
					c_vector<double, 3> location_of_neighbour_j = mpExtendedMesh->GetNode(second_order_neighs[j])->rGetLocation();
					image_location_per_second_order_neighbours.push_back(location_of_neighbour_j);
				}
			}

			std::vector<unsigned> first_order_neighs_vect;
			for(unsigned j=0; j<10; j++)
			{
				first_order_neighs_vect.push_back(epithelial_neighbours[mCellPopulationNodeIndexMap[cell_i_ext]][j]);
				// first_order_neighs_vect.push_back(epithelial_neighbours[cell_i_ext][j]);
			}
			first_order_neighs_vect.push_back(mCellPopulationNodeIndexMap[cell_i_ext]);
			// order and remove doubles
			std::sort(first_order_neighs_vect.begin(), first_order_neighs_vect.end());
			first_order_neighs_vect.erase(std::unique(first_order_neighs_vect.begin(), first_order_neighs_vect.end()), first_order_neighs_vect.end());
			// Remove zero from vector
			first_order_neighs_vect.erase(std::remove(first_order_neighs_vect.begin(), first_order_neighs_vect.end(), 0), first_order_neighs_vect.end());
			// PRINT_VECTOR(first_order_neighs_vect);
			/* 
			* Need:
			* First order neighbours  ->  first_order_neighs_vect
			* Traingulation		    ->  epithelial_triangulation
			* image node positions    ->  image_location_per_second_order_neighbours 
			*   needs (for node IDs)  ->  second_order_neighs 
			* Cell_i					->  cell_i_ext
			*/


			c_vector<double, 3> force_due_to_curvature = zero_vector<double>(3);;

			if(first_order_neighs_vect.size() >= 3)
			{
				force_due_to_curvature = GetForceDueToDiscreteCurvature(rCellPopulation, epithelial_triangulation, first_order_neighs_vect, second_order_neighs, image_location_per_second_order_neighbours, cell_i_ext, num_cells);

			}
			else if(first_order_neighs_vect.size() < 3)
			{
				// if(cell_i_ext > 199)
				// {
				// 	PRINT_VARIABLE(cell_i_ext);
					TRACE("Has no neighbours")
				// }
				force_due_to_curvature[0] = 0.0;
				force_due_to_curvature[1] = 0.0;
				force_due_to_curvature[2] = 0.0;
			}
			
			if (isnan(force_due_to_curvature[0]) || isnan(force_due_to_curvature[1]) || isnan(force_due_to_curvature[2]))
			{
				PRINT_VECTOR(force_due_to_curvature);
				force_due_to_curvature[0] = 0.0;
				force_due_to_curvature[1] = 0.0;
				force_due_to_curvature[2] = 0.0;
			}

			rCellPopulation.GetNode(cell_i_ext)->AddAppliedForceContribution(basement_membrane_parameter*force_due_to_curvature);
			int sim_time = (SimulationTime::Instance()->GetTime())/(SimulationTime::Instance()->GetTimeStep());
			// if(cell_i_ext == 154 && (sim_time%10 == 0) )
			// {
			// 	PRINT_VARIABLE(SimulationTime::Instance()->GetTime());
			// 	PRINT_VECTOR(force_due_to_curvature);
			// }
		}		
		
	}
	
}


double BendingForce3dHeightWithGhostNodes::GetPeriodicDomainWidth()
{
	return mPeriodicDomainWidth;
}

void BendingForce3dHeightWithGhostNodes::SetPeriodicDomainWidth(double periodicDomainWidth)
{
	mPeriodicDomainWidth = periodicDomainWidth;
}

double BendingForce3dHeightWithGhostNodes::GetPeriodicDomainDepth()
{
	return mPeriodicDomainDepth;
}

void BendingForce3dHeightWithGhostNodes::SetPeriodicDomainDepth(double periodicDomainDepth)
{
	mPeriodicDomainDepth = periodicDomainDepth;
}


void BendingForce3dHeightWithGhostNodes::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile <<  "\t\t\t<BasementMembraneParameter>"<<  mBasementMembraneParameter << "</BasementMembraneParameter> \n" ;
	*rParamsFile <<  "\t\t\t<ExponentParameter>"<<  mExponentParameter << "</ExponentParameter> \n" ;
	*rParamsFile <<  "\t\t\t<HeightDependantCurvatureParameter>"<<  mHeightDependantCurvatureParameter << "</HeightDependantCurvatureParameter> \n" ;
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
CHASTE_CLASS_EXPORT(BendingForce3dHeightWithGhostNodes)


