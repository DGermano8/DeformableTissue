#include "PeriodicBasementMembraneForce.hpp"
#include "Debug.hpp"

/**
 * To avoid warnings on some compilers, C++ style initialization of member
 * variables should be done in the order they are defined in the header file.
 */
PeriodicBasementMembraneForce::PeriodicBasementMembraneForce()
    : AbstractForce<2>(),
      mBasementMembraneParameter(DOUBLE_UNSET),
      mCryptBaseCurvature(DOUBLE_UNSET),
      mLeftBoundary(DOUBLE_UNSET),
      mRightBoundary(DOUBLE_UNSET),
      mUsePositionDependentMembraneForce(false),
      mMembraneForceMultiplier(DOUBLE_UNSET),
      mPeriodicDomainWidth(DOUBLE_UNSET),
      mpExtendedMesh(NULL)
{
    // Sets up output file
//	OutputFileHandler output_file_handler("CurvatureData/", false);
//	mMeinekeOutputFile = output_file_handler.OpenOutputFile("results.curvature");
}

PeriodicBasementMembraneForce::~PeriodicBasementMembraneForce()
{
    // Avoid memory leaks
    if (mpExtendedMesh != NULL)
    {
        delete mpExtendedMesh;
    }
}

void PeriodicBasementMembraneForce::SetBasementMembraneParameter(double basementMembraneParameter)
{
	mBasementMembraneParameter = basementMembraneParameter;
}

double PeriodicBasementMembraneForce::GetBasementMembraneParameter()
{
	return mBasementMembraneParameter;
}

void PeriodicBasementMembraneForce::SetCryptBaseCurvature(double cryptBaseCurvature, double leftBoundary, double rightBoundary)
{
	mCryptBaseCurvature = cryptBaseCurvature;
	mLeftBoundary = leftBoundary;
	mRightBoundary = rightBoundary;
}

double PeriodicBasementMembraneForce::GetCryptBaseCurvature()
{
	return mCryptBaseCurvature;
}

void PeriodicBasementMembraneForce::SetPositionDependentMultiplier(bool usePositionDependentMembraneForce, double membraneForceMultiplier)
{
	mUsePositionDependentMembraneForce = usePositionDependentMembraneForce;
	mMembraneForceMultiplier = membraneForceMultiplier;
}

double PeriodicBasementMembraneForce::GetPositionDependentMultiplier()
{
	return mMembraneForceMultiplier;
}

void PeriodicBasementMembraneForce::RemoveDuplicates1D(std::vector<unsigned>& rVectorWithDuplicates)
{
    std::sort(rVectorWithDuplicates.begin(), rVectorWithDuplicates.end());
    rVectorWithDuplicates.erase(std::unique(rVectorWithDuplicates.begin(), rVectorWithDuplicates.end()), rVectorWithDuplicates.end());
}

/*
 * A method to find all the pairs of connections between healthy epithelial cells and labelled tissue cells.
 * Returns a vector of node pairings, without repeats. The first of each pair is the epithelial node index,
 * and the second is the tissue node index. Updating so that it also returns mutant-labelled cell pairs.
 */
std::vector<c_vector<unsigned, 2> > PeriodicBasementMembraneForce::GetEpithelialTissuePairs(AbstractCellPopulation<2>& rCellPopulation)
{
    MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

    // Create a vector to record the pairs of nodes corresponding to *joined* epithelial and tissue nodes
    std::vector<c_vector<unsigned, 2> > node_pairs;
    c_vector<double, 2> pair;

    // Loop over nodes of mpExtendedMesh
    for (unsigned extended_node_index=0; extended_node_index<mpExtendedMesh->GetNumNodes(); extended_node_index++)
    {
        // Get a pointer to this node in mpExtendedMesh
        Node<2>* p_node = mpExtendedMesh->GetNode(extended_node_index);

        // Get the corresponding node index in rCellPopulation
        unsigned node_index = mExtendedMeshNodeIndexMap[extended_node_index];

        // Get the cell corresponding to this node
        CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(node_index);

        // Get mutation state of cell
    	boost::shared_ptr<AbstractCellMutationState> p_state = p_cell->GetMutationState();

    	// Get whether cell is dead
    	bool cell_is_dead = p_cell->IsDead();

    	// Get whether this cell is a live epithelial cell
    	bool is_live_epithelial_cell = (p_state->IsType<StromalCellMutationState>()==false) && !cell_is_dead;

    	if (is_live_epithelial_cell)
    	{
    	    // Iterate over elements of mpExtendedMesh containing this node and no ghost nodes
    		std::vector<unsigned> tissue_nodes;

    		for (Node<2>::ContainingElementIterator elem_iter = p_node->ContainingElementsBegin();
		         elem_iter != p_node->ContainingElementsEnd();
		         ++elem_iter)
    		{
    			bool element_contains_ghost_nodes = false;

    			// Get a pointer to the element
    			Element<2,2>* p_element = mpExtendedMesh->GetElement(*elem_iter);

    			// ITERATE OVER NODES owned by this element
    			for (unsigned local_index=0; local_index<3; local_index++)
    			{
    				unsigned nodeBGlobalIndex = p_element->GetNodeGlobalIndex(local_index);

    				// Get the corresponding node index in rCellPopulation
                    unsigned neighbour_index = mExtendedMeshNodeIndexMap[nodeBGlobalIndex];
    				
    				if (p_tissue->IsGhostNode(neighbour_index))
    				{
    					element_contains_ghost_nodes = true;
    					break; 				// This should break out of the inner for loop
    				}
    			}

				if (!element_contains_ghost_nodes)
				{
                    // ITERATE OVER NODES owned by this element
                    for (unsigned local_index=0; local_index<3; local_index++)
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
			}

    		// Remove any nodes that have been found twice
    		RemoveDuplicates1D(tissue_nodes);

    		// Now construct the vector of node pairs
    		for (unsigned i=0; i<tissue_nodes.size(); i++)
    		{
    			pair[0] = extended_node_index;
    			pair[1] = tissue_nodes[i];
    			node_pairs.push_back(pair);
    			    			
///\todo consider re-implementing this check
//    			// Check that these node share a common element
//				bool has_common_element = false;
//
//				// The elements that contain this epithelial node:
//				std::set<unsigned> epithelial_elements = rCellPopulation.GetNode(node_index)->rGetContainingElementIndices();
//				assert(!epithelial_elements.empty());
//
//				// The elements that contain the tissue node:
//				std::set<unsigned> tissue_elements = rCellPopulation.GetNode(tissue_nodes[i])->rGetContainingElementIndices();
//				assert(!tissue_elements.empty());
//
//				// Loop over all elements that contain the tissue node
//				for (  Node<2>::ContainingElementIterator elt_it = rCellPopulation.GetNode(tissue_nodes[i])->ContainingElementsBegin();
//				         elt_it != rCellPopulation.GetNode(tissue_nodes[i])->ContainingElementsEnd();
//				         ++elt_it)
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
//					TRACE("No common element between:");
//					PRINT_2_VARIABLES(node_index,tissue_nodes[i]);
//				}
//				assert(has_common_element);
    		}    		
    	}
    }

    return node_pairs;
}

/*
 * Method to determine whether an element contains ghost nodes
 */
bool PeriodicBasementMembraneForce::DoesElementContainGhostNodes(AbstractCellPopulation<2>& rCellPopulation, unsigned elementIndex)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

	bool element_contains_ghost_nodes = false;

	// Get a pointer to the element
	Element<2,2>* p_element = p_tissue->rGetMesh().GetElement(elementIndex);

	// ITERATE OVER NODES owned by this element
	for (unsigned local_index=0; local_index<3; local_index++)
	{
		if (p_tissue->IsGhostNode(p_element->GetNodeGlobalIndex(local_index)) == true)
		{
			element_contains_ghost_nodes = true;

		}
	}

	return element_contains_ghost_nodes;
}

/* Method to return the current coordinates of the crypt orifice and
 * crypt base - these can be used to accurately define the region of the
 * crypt base. (This will be the y-coordinate in 2D, or the z coordinate in 3D)
 * [0] - y or z-coordinate of orifice
 * [1] - y or z-coordinate of base
 */
c_vector<double,2> PeriodicBasementMembraneForce::GetCryptHeightExtremes(AbstractCellPopulation<2>& rCellPopulation)
{
    MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

    // Create a vector to store the y-coordinates of the lowest point of the crypt base and the highest point of the
    // crypt orifice
    c_vector<double,2> height_extremes;

    double max_height = 0.0;
    double min_height = DBL_MAX;

    double current_height_coordinate;

    // We iterate over all cells in the tissue, and deal only with those that are epithelial cells
    for (  AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
    	boost::shared_ptr<AbstractCellMutationState> p_state = cell_iter->GetMutationState();

	   	// Need these to not be labelled cells
	   	if ( (p_state->IsType<StromalCellMutationState>()==false) )
	   	{
	   		Node<2>* p_node = p_tissue->GetNodeCorrespondingToCell(*cell_iter);

   			current_height_coordinate = p_node->rGetLocation()[1];

	    	if (current_height_coordinate > max_height)
	    	{
	    		max_height = current_height_coordinate;
	    	}
	    	else if (current_height_coordinate < min_height)
	    	{
	    		min_height = current_height_coordinate;
	    	}
	    }
    }

    height_extremes[0] = max_height;
    height_extremes[1] = min_height;

    return height_extremes;
}

/*
 * Using the vector of node pairs found in GetEpithelialTissuePairs to determine the curvature
 * of the curve that passes through the midpoints of the neighbouring springs, given a epithelial-tissue
 * node pairing. (Note - it ignores the end pairs because one of the common elements will contain ghost nodes, but this
 * should only crop up if you don't have periodic boundary conditions)
 * Updating this so that it will still find the curvature if one of the epithelial cells is a mutant cell, eg. apc2 hit
 */
double PeriodicBasementMembraneForce::GetCurvatureFromNodePair(AbstractCellPopulation<2>& rCellPopulation,
                                                               unsigned epithelialNodeIndex,
															   unsigned tissueNodeIndex)
{    
    // Seeking the common elements that contain both the epithelial node and the tissue node
	// Note: Don't want to count any elements that have ghost nodes (but there shouldn't be any ghost nodes
	// in the extended mesh.

	std::vector<unsigned> common_elements;	// Initialising
	
	// The elements that contain this epithelial node:
	
    std::set<unsigned> epithelial_elements = mpExtendedMesh->GetNode(epithelialNodeIndex)->rGetContainingElementIndices();

    // Safety check
    assert(!epithelial_elements.empty());
    assert(mpExtendedMesh->GetNode(tissueNodeIndex)->GetNumContainingElements() != 0);
        
    // Loop over all elements that contain the tissue node
    for (Node<2>::ContainingElementIterator elt_it = mpExtendedMesh->GetNode(tissueNodeIndex)->ContainingElementsBegin();
         elt_it != mpExtendedMesh->GetNode(tissueNodeIndex)->ContainingElementsEnd();
         ++elt_it)
    {
    	unsigned elt_index = *elt_it;

///\todo there should not be any ghost nodes in mpExtendedMesh
//    	bool elt_contains_ghost_nodes = DoesElementContainGhostNodes(rCellPopulation, elt_index);

    	// Keep only those elements that also contain the epithelial node, but do not have ghost nodes

        if (epithelial_elements.find(elt_index) != epithelial_elements.end())
        {
        	// Common element
           	common_elements.push_back(elt_index);
        }
    }

	assert(!common_elements.empty()); // This is bad - the nodes should be connected in the first place...
	
	// We iterate over these common elements to find the midpoints of the springs that
	// connect epithelial and  tissue nodes

	c_vector<double, 2> spring_midpoint_a, spring_midpoint_b, spring_midpoint_c;

	spring_midpoint_b = mpExtendedMesh->GetNode(epithelialNodeIndex)->rGetLocation() + 0.5*(mpExtendedMesh->GetVectorFromAtoB(mpExtendedMesh->GetNode(epithelialNodeIndex)->rGetLocation(), mpExtendedMesh->GetNode(tissueNodeIndex)->rGetLocation()));

	// If there is only one such common element, then this epithelial node will be at either end and so we don't
    // consider the force along the very first / very last spring (only happens if you don't use a cylindrical
	// mesh!)
    if (common_elements.size() == 1)
    {
    	double curvature = 0.0;
    	return curvature;
    }
    else
    {
    	assert(common_elements.size() == 2);		// Should only be two common elements

    	for (std::vector<unsigned>::iterator elem_iter = common_elements.begin();
		   	 elem_iter != common_elements.end();
		   	 ++elem_iter)
    	{
    		// Need to determine the cell type of each local node in the element
    		// Want to find the midpoint between epithelial and tissue pairs

    		// Looking at the three nodes which form the vertices of the element
    	    unsigned global_index_A = mpExtendedMesh->GetElement(*elem_iter)->GetNodeGlobalIndex(0);
            unsigned global_index_B = mpExtendedMesh->GetElement(*elem_iter)->GetNodeGlobalIndex(1);
            unsigned global_index_C = mpExtendedMesh->GetElement(*elem_iter)->GetNodeGlobalIndex(2);
    	    	   		
	   		unsigned E=UINT_MAX;  // E - the epithelial node we're interested in,
	   		unsigned T=UINT_MAX;  // T - the tissue node it's connected to, and
	   		unsigned P=UINT_MAX;  // P - the other node in the element.

	   		if (global_index_A == epithelialNodeIndex)
	   		{
	   			E = global_index_A;

	   			if (global_index_B == tissueNodeIndex)
	   			{
	   				T = global_index_B;
	   				P = global_index_C;
	   			}
	   			else
	   			{
	   				P = global_index_B;
	   				T = global_index_C;
	   			}
	   		}
	   		else if (global_index_B == epithelialNodeIndex)
	   		{
	   			E = global_index_B;

	   			if (global_index_A == tissueNodeIndex)
	   			{
	   				T = global_index_A;
	   				P = global_index_C;
	   			}
	   			else
	   			{
	   				P = global_index_A;
	   				T = global_index_C;
	   			}
	   		}
	   		else if (global_index_C == epithelialNodeIndex)
	   		{
	   			E = global_index_C;

	   			if (global_index_A == tissueNodeIndex)
	   			{
	   				T = global_index_A;
	   				P = global_index_B;
	   			}
	   			else
	   			{
	   				P = global_index_A;
	   				T = global_index_B;
	   			}
	   		}

	   		assert(E<UINT_MAX);
            assert(E == epithelialNodeIndex);
	   		
            assert(T<UINT_MAX);
	   		assert(T == tissueNodeIndex);
	   		
            assert(P<UINT_MAX);
   		
	   		// Get the real node indices to check the associated cell types
	   		unsigned real_epithelial_node = mExtendedMeshNodeIndexMap[E];
	   		unsigned real_tissue_node = mExtendedMeshNodeIndexMap[T];
            unsigned real_other_node = mExtendedMeshNodeIndexMap[P];
	   		
	   		// Check that we have assigned the epithelial node correctly
	  
	   		CellPtr p_E_cell = rCellPopulation.GetCellUsingLocationIndex(real_epithelial_node);
	   		boost::shared_ptr<AbstractCellMutationState> p_E_state = p_E_cell->GetMutationState();
	   		assert( (p_E_state->IsType<WildTypeCellMutationState>()) || (p_E_state->IsType<ApcTwoHitCellMutationState>()) );
	   		assert(!p_E_state->IsType<StromalCellMutationState>());

	   		// Check that we have assigned the tissue node correctly
	   		CellPtr p_T_cell = rCellPopulation.GetCellUsingLocationIndex(real_tissue_node);
	   		boost::shared_ptr<AbstractCellMutationState> p_T_state = p_T_cell->GetMutationState();
	   		assert(p_T_state->IsType<StromalCellMutationState>());
	   		
	   		/*
	   		 * Now we work with E (the epithelial node), T (the tissue node it's connected to) and P, the other node,
	   		 * which we will now assign as either P1 (det < 0) or P2 (det > 0).
	   		 *
	   		 * We also need to determine whether P1 and P2 are epithelial or tissue, as this will affect which vector we
	   		 * choose to take the spring midpoint from.
	   		 */
	   		c_vector<double, 2> vector_T_to_E = mpExtendedMesh->GetVectorFromAtoB(mpExtendedMesh->GetNode(T)->rGetLocation(),mpExtendedMesh->GetNode(E)->rGetLocation());
	   		c_vector<double, 2> vector_E_to_T = mpExtendedMesh->GetVectorFromAtoB(mpExtendedMesh->GetNode(E)->rGetLocation(),mpExtendedMesh->GetNode(T)->rGetLocation());
	   		c_vector<double, 2> vector_E_to_P = mpExtendedMesh->GetVectorFromAtoB(mpExtendedMesh->GetNode(E)->rGetLocation(),mpExtendedMesh->GetNode(P)->rGetLocation());
	   		c_vector<double, 2> vector_T_to_P = vector_T_to_E + vector_E_to_P;
	   		c_vector<double, 2> vector_P_to_T = mpExtendedMesh->GetVectorFromAtoB(mpExtendedMesh->GetNode(P)->rGetLocation(),mpExtendedMesh->GetNode(T)->rGetLocation());

	   		// Now calculate the determinant (TE x EP)

	   		double det = vector_T_to_E[0]*vector_E_to_P[1] - vector_T_to_E[1]*vector_E_to_P[0];
	  		assert(!isnan(det));
	   		assert(det != 0.0);

	   		/*
	   		 * If det < 0 then P = P1 and we can assign spring_midpoint_c
	   		 * If det > 0 then P = P2 and we can assign spring_midpoint_a
	   		 * Also need to take into account whether P is a epithelial or tissue node
	   		 * to choose the right spring
	   		 */
	   		CellPtr p_P_cell = rCellPopulation.GetCellUsingLocationIndex(real_other_node);
	   		boost::shared_ptr<AbstractCellMutationState> p_state = p_P_cell->GetMutationState();
	   		
	   		if ( (det < 0) && (p_state->IsType<StromalCellMutationState>()==false) )	// P = Epithelial, not labelled
	   		{
	   			spring_midpoint_c = mpExtendedMesh->GetNode(E)->rGetLocation() + vector_E_to_P + 0.5*vector_P_to_T;
	   		}
	   		else if ((det < 0) && (p_state->IsType<StromalCellMutationState>()==true))	// P = Tissue
	   		{
	   			spring_midpoint_c = mpExtendedMesh->GetNode(E)->rGetLocation() + 0.5*vector_E_to_P;
	   		}

	   		else if ( (det > 0) && (p_state->IsType<StromalCellMutationState>()==false) )	// P = Epithelial, not labelled
	   		{
	   			spring_midpoint_a = mpExtendedMesh->GetNode(E)->rGetLocation() + vector_E_to_T + 0.5*vector_T_to_P;
	   		}

	   		else if ((det > 0) && (p_state->IsType<StromalCellMutationState>()==true))	// P = Tissue
	   		{
	   			spring_midpoint_a = mpExtendedMesh->GetNode(E)->rGetLocation() + 0.5*vector_E_to_P;
	   		}
    	}

    	double curvature = FindParametricCurvature(spring_midpoint_a, spring_midpoint_b, spring_midpoint_c);

    	// Need to subtract the target curvature if the epithelial node lies in the crypt base
    	c_vector<double, 2> epithelial_location = mpExtendedMesh->GetNode(epithelialNodeIndex)->rGetLocation();

		// We need to use the current crypt height and base levels to determine where to apply the non-zero curvature
		// [0] - y-coordinate of orifice
		// [1] - y-coordinate of base
//		c_vector<double, 2> height_extremes = GetCryptHeightExtremes(rCellPopulation);
//
//		double top_of_crypt_base = height_extremes(1) + (height_extremes(0) - height_extremes(1))*0.2;
//
////		// For cross section crypt model:
//		if ( (epithelial_location[1] <= top_of_crypt_base))
//		{
//			curvature -= mCryptBaseCurvature;
//		}

		// For sandwich box model:
//		if ( (epithelial_location[1] <= 20.0) && (epithelial_location[0] > mLeftBoundary) && (epithelial_location[0] < mRightBoundary) )
//		{
//			curvature -= mCryptBaseCurvature;
//		}

    	assert(!isnan(curvature));
    	return curvature;
    }
}

/*
 * Function to return the curvature between three midpoints parametrically - in this case, we find the normal
 * to the vector joining the left and right midpoints, and then find the perpendicular distance of the centre midpoint
 * from the left->right vector
 */

double PeriodicBasementMembraneForce::GetCurvatureFromMidpoints(AbstractCellPopulation<2>& rCellPopulation,
																c_vector<double, 2> leftMidpoint,
																c_vector<double, 2> centreMidpoint,
																c_vector<double, 2> rightMidpoint)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

	// Firstly find the normal to the vector joining the left and right midpoints
	c_vector<double, 2>	vector_left_to_right = p_tissue->rGetMesh().GetVectorFromAtoB(leftMidpoint, rightMidpoint);
	c_vector<double, 2> normal_vector;
	normal_vector[0] = vector_left_to_right[1];
	normal_vector[1] = -vector_left_to_right[0];

	c_vector<double, 2> vector_left_to_centre = p_tissue->rGetMesh().GetVectorFromAtoB(leftMidpoint, centreMidpoint);

	double curvature = normal_vector[0]*vector_left_to_centre[0] + normal_vector[1]*vector_left_to_centre[1];

	return curvature;
}

/*
859	 * Function to return the curvature between three points parametrically - the midpoints of the springs connecting the
860	 * transit cells to the differentiated cells. NB. The input arguments need to be in order from either left to right
861	 * or right to left. If they are wrongly arranged (eg. middle, left, right) then you get a different curvature,
862	 * but left->right = -(right-> left).
863	 */

double PeriodicBasementMembraneForce::FindParametricCurvature(c_vector<double, 2> leftMidpoint,
															c_vector<double, 2> centreMidpoint,
															c_vector<double, 2> rightMidpoint)
{
	// Firstly find the parametric intervals
	double left_s = sqrt(pow(centreMidpoint[0] - leftMidpoint[0],2) + pow(centreMidpoint[1] - leftMidpoint[1],2));
	double right_s = sqrt(pow(rightMidpoint[0] - centreMidpoint[0],2) + pow(rightMidpoint[1] - centreMidpoint[1],2));

	assert(left_s >= 0);
	assert(right_s >= 0);

	double sum_intervals = left_s + right_s;

	double x_prime = (rightMidpoint[0] - leftMidpoint[0])/sum_intervals;
	double y_prime = (rightMidpoint[1] - leftMidpoint[1])/sum_intervals;

	double x_double_prime = 2*(left_s*rightMidpoint[0] - sum_intervals*centreMidpoint[0] + right_s*leftMidpoint[0])/(left_s*right_s*sum_intervals);
	double y_double_prime = 2*(left_s*rightMidpoint[1] - sum_intervals*centreMidpoint[1] + right_s*leftMidpoint[1])/(left_s*right_s*sum_intervals);

	double curvature = (x_prime*y_double_prime - y_prime*x_double_prime)/pow((pow(x_prime,2) + pow(y_prime,2)),3/2);

	return curvature;
}


/*
 * A method to return the number of elements that contain a particular node,
 * excluding those elements that have ghost nodes
 */

unsigned PeriodicBasementMembraneForce::GetNumContainingElementsWithoutGhostNodes(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

    // Get pointer to the node
    Node<2>* p_node = p_tissue->GetNode(nodeIndex);
    assert(!(p_tissue->IsGhostNode(nodeIndex)));

    unsigned num_elements_with_no_ghost_nodes = 0;		// Initialise

    // Iterate over containing elements
    for (  Node<2>::ContainingElementIterator iter = p_node->ContainingElementsBegin();
         iter != p_node->ContainingElementsEnd(); ++iter)
    {
        bool element_contains_ghost_nodes = false;
        Element<2,2>* p_element = p_tissue->rGetMesh().GetElement(*iter);

        // Iterate over nodes owned by this element
        for (unsigned local_index=0; local_index<3; local_index++)
        {
            if (p_tissue->IsGhostNode(p_element->GetNodeGlobalIndex(local_index)) == true)
            {
                element_contains_ghost_nodes = true;
                break; // I think this should break out of the inner for loop
            }
        }

        if (element_contains_ghost_nodes==false)
        {
            // This element contains no ghost nodes
            num_elements_with_no_ghost_nodes++;
        }
    }

    return num_elements_with_no_ghost_nodes;
}

/*
 * Method to return the nodes connected to a particular node via the Delaunay
 * triangulation, excluding ghost nodes.
 */

std::set<unsigned> PeriodicBasementMembraneForce::GetNeighbouringNodeIndices(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

	assert(!(p_tissue->IsGhostNode(nodeIndex)));

	// Create a set of neighbouring node indices
	std::set<unsigned> neighbouring_node_indices;

    // Find the indices of the elements owned by this node
	std::set<unsigned> containing_elem_indices = p_tissue->GetNode(nodeIndex)->rGetContainingElementIndices();

    // Iterate over these elements
    for (std::set<unsigned>::iterator elem_iter = containing_elem_indices.begin();
         elem_iter != containing_elem_indices.end();
         ++elem_iter)
    {
        // Get all the nodes contained in this element
        // Don't want to include the current node
        unsigned neighbour_global_index;

        for (unsigned local_index=0; local_index<3; local_index++)
        {
            neighbour_global_index = p_tissue->rGetMesh().GetElement(*elem_iter)->GetNodeGlobalIndex(local_index);

            if( (neighbour_global_index != nodeIndex) && (!p_tissue->IsGhostNode(neighbour_global_index)) )
            {
            	neighbouring_node_indices.insert(neighbour_global_index);
            }
        }
    }

    return neighbouring_node_indices;
}

/** Method to determine if an epithelial cell has lost all contacts with the stromal cells below
 * TRUE if cell has popped up
 * FALSE if cell remains in the monolayer
 */

bool PeriodicBasementMembraneForce::HasEpithelialCellDetachedFromBasementMembrane(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex)
{
	MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);

	bool has_cell_detached = false;	// Initialising

   	std::set<unsigned> neighbours = GetNeighbouringNodeIndices(rCellPopulation, nodeIndex);

   	unsigned num_stromal_neighbours = 0;

   	// Iterate over the neighbouring cells to check the number of tissue cell neighbours

   	for(std::set<unsigned>::iterator neighbour_iter=neighbours.begin();
   							neighbour_iter != neighbours.end();
   							++neighbour_iter)
	{
   		boost::shared_ptr<AbstractCellMutationState> p_state = p_tissue->GetCellUsingLocationIndex(*neighbour_iter)->GetMutationState();
		if ( (!p_tissue->IsGhostNode(*neighbour_iter)) && (p_state->IsType<StromalCellMutationState>()==true))
   		{
			num_stromal_neighbours += 1;
		}
   	}

   	if(num_stromal_neighbours < 1)
   	{
   		has_cell_detached = true;
   	}

	return has_cell_detached;
}

void PeriodicBasementMembraneForce::AddForceContribution(std::vector<c_vector<double, 2> >& rForces,
                                                         AbstractCellPopulation<2>& rCellPopulation)
{
    mExtendedMeshNodeIndexMap.clear();

    // Create a vector of nodes for use in constructing mpExtendedMesh
    unsigned num_cells = rCellPopulation.GetNumRealCells();
    std::vector<Node<2>*> extended_nodes(2*num_cells);

	// The width of the extended mesh
	double extended_mesh_width = mPeriodicDomainWidth;

    // We iterate over all cells in the population
	unsigned count = 0;
    for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // First, create and store a copy of this real node and cell
        unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        c_vector<double, 2> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

        // Create a copy of the node corresponding to this cell and store it
        Node<2>* p_real_node = new Node<2>(real_node_index, real_node_location);
        extended_nodes[count] = p_real_node;
        
        /**
         * \todo The code block below would need to be amended for 3d to cope with
         * the z direction too and to make sure we copy from left to right and from
         * front to back.
         */

    	// Compute the location of the image node corresponding to this node
        c_vector<double,2> image_node_location = real_node_location;
		if (real_node_location[0] >= mPeriodicDomainWidth*0.5) //centroid(0)) // Right-hand boundary node
		{
			image_node_location[0] -= extended_mesh_width;
		}
		else if (real_node_location[0] <  mPeriodicDomainWidth*0.5) //centroid(0))
		{
			image_node_location[0] += extended_mesh_width;
		}

        // Create a copy of the node corresponding to this cell, suitable translated, and store it
        Node<2>* p_image_node = new Node<2>(num_cells+count, image_node_location);
        extended_nodes[num_cells+count] = p_image_node;

        // Populate mExtendedMeshNodeIndexMap
        mExtendedMeshNodeIndexMap[count] = real_node_index;
        mExtendedMeshNodeIndexMap[num_cells+count] = real_node_index;

		count++;
    }

    // We now construct mpExtendedMesh using extended_nodes
    if (mpExtendedMesh != NULL)
    {
    	delete mpExtendedMesh;
    }
    mpExtendedMesh = new MutableMesh<2,2>(extended_nodes);

    ///\todo We may need to modify rCellPopulation to ensure thatmMarkedSprings is correct at this point

	// Identify epithelial-stromal cell pairs
	std::vector<c_vector<unsigned, 2> > node_pairs = GetEpithelialTissuePairs(rCellPopulation);
	
	// We loop over the epithelial-tissue node pairs to find the force acting on that
	// epithelial node, and the direction in which it acts
	for (unsigned i=0; i<node_pairs.size(); i++)
	{
		unsigned extended_epithelial_node_index = node_pairs[i][0];
		unsigned extended_tissue_node_index = node_pairs[i][1];
				
		// Get the corresponding node index in rCellPopulation
        unsigned epithelial_node_index = mExtendedMeshNodeIndexMap[extended_epithelial_node_index];
        unsigned tissue_node_index = mExtendedMeshNodeIndexMap[extended_tissue_node_index];
                
        // Get the cells corresponding to these nodes to check the cell types
        CellPtr p_epithelial_cell = rCellPopulation.GetCellUsingLocationIndex(epithelial_node_index);
        assert(p_epithelial_cell->GetMutationState()->IsType<StromalCellMutationState>() == false);

		CellPtr p_tissue_cell = rCellPopulation.GetCellUsingLocationIndex(tissue_node_index);		
		assert(p_tissue_cell->GetMutationState()->IsType<StromalCellMutationState>() == true);

		// Get the locations in the extended mesh
		c_vector<double, 2> epithelial_location = mpExtendedMesh->GetNode(extended_epithelial_node_index)->rGetLocation();		   
		c_vector<double, 2> tissue_location = mpExtendedMesh->GetNode(extended_tissue_node_index)->rGetLocation();

		// The force due to the basal lamina acts along the spring connecting the epithelial and tissue nodes, T->E direction
		c_vector<double, 2> curvature_force_direction = mpExtendedMesh->GetVectorFromAtoB(tissue_location, epithelial_location);

		double distance_between_nodes = norm_2(curvature_force_direction);
		assert(distance_between_nodes > 0);
		assert(!isnan(distance_between_nodes));

		curvature_force_direction /= distance_between_nodes;

		double curvature = GetCurvatureFromNodePair(rCellPopulation, extended_epithelial_node_index, extended_tissue_node_index);

		double basement_membrane_parameter = GetBasementMembraneParameter();

		c_vector<double, 2> force_due_to_basement_membrane = basement_membrane_parameter*curvature*curvature_force_direction;
	
		// Only want to add the force to real cells

        // Now we make sure that we only apply the force to the real node and not the image node	
        if (extended_epithelial_node_index < num_cells)
        {
			rForces[epithelial_node_index] += force_due_to_basement_membrane;
        }
	}
}

double PeriodicBasementMembraneForce::GetPeriodicDomainWidth()
{
	return mPeriodicDomainWidth;
}

void PeriodicBasementMembraneForce::SetPeriodicDomainWidth(double periodicDomainWidth)
{
	mPeriodicDomainWidth = periodicDomainWidth;
}

void PeriodicBasementMembraneForce::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile <<  "\t\t\t<BasementMembraneParameter>"<<  mBasementMembraneParameter << "</BasementMembraneParameter> \n" ;
	*rParamsFile <<  "\t\t\t<CryptBaseCurvature>"<<  mCryptBaseCurvature << "</CryptBaseCurvature> \n" ;
	*rParamsFile <<  "\t\t\t<LeftBoundary>"<<  mLeftBoundary << "</LeftBoundary> \n" ;
	*rParamsFile <<  "\t\t\t<RightBoundary>"<<  mRightBoundary << "</RightBoundary> \n" ;
	*rParamsFile <<  "\t\t\t<UsePositionDependentMembraneForce>"<<  mUsePositionDependentMembraneForce << "</UsePositionDependentMembraneForce> \n" ;
	*rParamsFile <<  "\t\t\t<MembraneForceMultiplier>"<<  mMembraneForceMultiplier << "</MembraneForceMultiplier> \n" ;
	*rParamsFile <<  "\t\t\t<PeriodicDomainWidth>"<<  mPeriodicDomainWidth << "</PeriodicDomainWidth> \n" ;

	// Call direct parent class
	AbstractForce<2>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PeriodicBasementMembraneForce)
