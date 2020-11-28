#include "PeriodicCryptModelInteractionForceWithGhostNodes.hpp"
#include "Debug.hpp"
#include <cmath>
#include <list>
#include <fstream>
#include "AbstractSimplePhaseBasedCellCycleModel.hpp"
#include "AbstractCellProperty.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"

/**
 * To avoid warnings on some compilers, C++ style initialization of member
 * variables should be done in the order they are defined in the header file.
 */
template<unsigned DIM>
PeriodicCryptModelInteractionForceWithGhostNodes<DIM>::PeriodicCryptModelInteractionForceWithGhostNodes()
   : LinearSpringWithVariableSpringConstantsForce<DIM>(),
   mPeriodicDomainWidth(DOUBLE_UNSET),
   mPeriodicDomainDepth(DOUBLE_UNSET),
   mpExtendedMesh(NULL),
   mUseCellTypeDependentSprings(false),
   mTransitTransitMultiplier(DOUBLE_UNSET),
   mDifferentiatedDifferentiatedMultiplier(DOUBLE_UNSET),
   mTransitDifferentiatedMultiplier(DOUBLE_UNSET),
   mUseEpithelialStromalCellDependentSprings(false),
   mEpithelialEpithelialMultiplier(DOUBLE_UNSET),
   mStromalStromalMultiplier(DOUBLE_UNSET),
   mEpithelialStromalMultiplier(DOUBLE_UNSET),
   mApcTwoHitStromalMultiplier(DOUBLE_UNSET),
   mUseEdgeBasedSpringConstant(false),
   mUseOneWaySprings(false),
   mUsePositionDependentSpringConstants(false),
   mSpringConstantsMultiplier(DOUBLE_UNSET)
{
    // Sets up output file
//	OutputFileHandler output_file_handler("CurvatureData/", false);
//	mMeinekeOutputFile = output_file_handler.OpenOutputFile("results.curvature");
}

template<unsigned DIM>
PeriodicCryptModelInteractionForceWithGhostNodes<DIM>::~PeriodicCryptModelInteractionForceWithGhostNodes()
{
    // Avoid memory leaks
    if (mpExtendedMesh != NULL)
    {
        delete mpExtendedMesh;
    }
//    mMeinekeOutputFile->close();
}

template<unsigned DIM>
void PeriodicCryptModelInteractionForceWithGhostNodes<DIM>::SetCellTypeDependentSprings(bool useCellTypeDependentSprings,
		double transitTransitMultiplier,
		double differentiatedDifferentiatedMultiplier,
		double transitDifferentiatedMultiplier)
{
    mUseCellTypeDependentSprings = useCellTypeDependentSprings;
    mTransitTransitMultiplier = transitTransitMultiplier;
    mDifferentiatedDifferentiatedMultiplier = differentiatedDifferentiatedMultiplier;
    mTransitDifferentiatedMultiplier = transitDifferentiatedMultiplier;
}

template<unsigned DIM>
void PeriodicCryptModelInteractionForceWithGhostNodes<DIM>::SetEpithelialStromalCellDependentSprings(bool useEpithelialStromalCellDependentSprings,
		double epithelialEpithelialMultiplier,
		double stromalStromalMultiplier,
		double epithelialStromalMultiplier,
		double apcTwoHitStromalMultiplier)
{
	mUseEpithelialStromalCellDependentSprings = useEpithelialStromalCellDependentSprings;
	mEpithelialEpithelialMultiplier = epithelialEpithelialMultiplier;
    mStromalStromalMultiplier = stromalStromalMultiplier;
    mEpithelialStromalMultiplier = epithelialStromalMultiplier;
    mApcTwoHitStromalMultiplier = apcTwoHitStromalMultiplier;
}

template<unsigned DIM>
void PeriodicCryptModelInteractionForceWithGhostNodes<DIM>::SetEdgeBasedSpringConstant(bool useEdgeBasedSpringConstant)
{
    assert(DIM == 2);
    mUseEdgeBasedSpringConstant = useEdgeBasedSpringConstant;
}

template<unsigned DIM>
void PeriodicCryptModelInteractionForceWithGhostNodes<DIM>::SetUseOneWaySprings(bool useOneWaySprings)
{
	mUseOneWaySprings = useOneWaySprings;
}

template<unsigned DIM>
void PeriodicCryptModelInteractionForceWithGhostNodes<DIM>::SetPositionDependentSpringConstants(bool usePositionDependentSpringConstants, double springConstantsMultiplier)
{
	mUsePositionDependentSpringConstants = usePositionDependentSpringConstants;
	mSpringConstantsMultiplier = springConstantsMultiplier;
}

template<unsigned DIM>
double PeriodicCryptModelInteractionForceWithGhostNodes<DIM>::GetPositionDependentSpringConstants()
{
	return mSpringConstantsMultiplier;
}

template<unsigned DIM>
double PeriodicCryptModelInteractionForceWithGhostNodes<DIM>::VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex,
																							unsigned nodeBGlobalIndex,
																							AbstractCellPopulation<DIM>& rCellPopulation,
																							bool isCloserThanRestLength)
{
    double multiplication_factor = LinearSpringWithVariableSpringConstantsForce<DIM>::VariableSpringConstantMultiplicationFactor(nodeAGlobalIndex,
																																nodeBGlobalIndex,
																																rCellPopulation,
																																isCloserThanRestLength);

    DomMeshBasedCellPopulationWithGhostNodes<DIM>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);
    assert(!(p_tissue->IsGhostNode(nodeAGlobalIndex)));
    assert(!(p_tissue->IsGhostNode(nodeBGlobalIndex)));

    CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
    CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (mUseCellTypeDependentSprings)
    {
        boost::shared_ptr<AbstractCellProperty> cell_A_type = p_cell_A->GetCellProliferativeType();
        boost::shared_ptr<AbstractCellProperty> cell_B_type = p_cell_B->GetCellProliferativeType();

        if (cell_A_type == cell_B_type)
        {
            if (cell_A_type->IsType<TransitCellProliferativeType>())
            {
            	multiplication_factor *= mTransitTransitMultiplier;
            }

            if (cell_A_type->IsType<DifferentiatedCellProliferativeType>())
            {
            	multiplication_factor *= mDifferentiatedDifferentiatedMultiplier;
            }
        }
        else
        {
        	multiplication_factor *= mTransitDifferentiatedMultiplier;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (mUseEpithelialStromalCellDependentSprings)
    {
        if ( (p_cell_A->GetMutationState()->IsType<StromalCellMutationState>()==false)			// If both not stromal => epithelial
			&& (p_cell_B->GetMutationState()->IsType<StromalCellMutationState>()==false) )
		{
			multiplication_factor *= mEpithelialEpithelialMultiplier;
		}
		else if ( (p_cell_A->GetMutationState()->IsType<StromalCellMutationState>()==true)		// If both stromal
				&& (p_cell_B->GetMutationState()->IsType<StromalCellMutationState>()==true) )
		{
			multiplication_factor *= mStromalStromalMultiplier;
		}
        else if ( ( (p_cell_A->GetMutationState()->IsType<StromalCellMutationState>()==false) && (p_cell_B->GetMutationState()->IsType<StromalCellMutationState>()==true) )
        	|| ( (p_cell_A->GetMutationState()->IsType<StromalCellMutationState>()==true) && (p_cell_B->GetMutationState()->IsType<StromalCellMutationState>()==false) ) )
        {
        	multiplication_factor *= mEpithelialStromalMultiplier;
        }
        else if ( ( (p_cell_A->GetMutationState()->IsType<ApcTwoHitCellMutationState>()) && (p_cell_A->GetMutationState()->IsType<StromalCellMutationState>()==false) && (p_cell_B->GetMutationState()->IsType<StromalCellMutationState>()==true) )
        		||  ( (p_cell_A->GetMutationState()->IsType<StromalCellMutationState>()==true) && (p_cell_B->GetMutationState()->IsType<ApcTwoHitCellMutationState>()) && (p_cell_B->GetMutationState()->IsType<StromalCellMutationState>()==false) ) )
        {
        	multiplication_factor *= mApcTwoHitStromalMultiplier;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (mUsePositionDependentSpringConstants)
    {
    	// Only need to worry about the connections between epithelial cells, as there is zero
    	// attractive force between epithelial and stromal cells

    	// Get the y-coordinate for the top of the crypt base
    	c_vector<double, 2> height_extremes = GetCryptHeightExtremes(rCellPopulation);

		double top_of_crypt_base = height_extremes(1) + (height_extremes(0) - height_extremes(1))*0.2;

		c_vector<double, 2> node_A_location = p_tissue->GetNode(nodeAGlobalIndex)->rGetLocation();
		c_vector<double, 2> node_B_location = p_tissue->GetNode(nodeBGlobalIndex)->rGetLocation();

        if ( (p_cell_A->GetMutationState()->IsType<StromalCellMutationState>()==false)			// If both not labelled => healthy
			&& (p_cell_B->GetMutationState()->IsType<StromalCellMutationState>()==false)
			&& ( (node_A_location[1] < top_of_crypt_base) || (node_B_location[1] < top_of_crypt_base) ))
        {
        	multiplication_factor *= mSpringConstantsMultiplier;
        }
    }

    return multiplication_factor;
}

/* Method to return the current coordinates of the crypt orifice and
 * crypt base - these can be used to accurately define the region of the
 * crypt base. (This will be the y-coordinate in 2D, or the z coordinate in 3D)
 * [0] - y or z-coordinate of orifice
 * [1] - y or z-coordinate of base
 */
template<unsigned DIM>
c_vector<double,2> PeriodicCryptModelInteractionForceWithGhostNodes<DIM>::GetCryptHeightExtremes(AbstractCellPopulation<DIM>& rCellPopulation)
{
    DomMeshBasedCellPopulationWithGhostNodes<DIM>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

    // Create a vector to store the y-coordinates of the lowest point of the crypt base and the highest point of the
    // crypt orifice
    c_vector<double,2> height_extremes;

    double max_height = 0.0;
    double min_height = DBL_MAX;

    double current_height_coordinate;

    // We iterate over all cells in the tissue, and deal only with those that are epithelial cells
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
    	boost::shared_ptr<AbstractCellMutationState> p_state = cell_iter->GetMutationState();

	   	// Need these to not be labelled cells
	   	if ( (p_state->IsType<StromalCellMutationState>()==false) )
	   	{
	   		Node<DIM>* p_node = p_tissue->GetNodeCorrespondingToCell(*cell_iter);

   			current_height_coordinate = p_node->rGetLocation()[DIM-1];

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

template<unsigned DIM>
void PeriodicCryptModelInteractionForceWithGhostNodes<DIM>::RemoveDuplicates1D(std::vector<unsigned>& rVectorWithDuplicates)
{
    std::sort(rVectorWithDuplicates.begin(), rVectorWithDuplicates.end());
    rVectorWithDuplicates.erase(std::unique(rVectorWithDuplicates.begin(), rVectorWithDuplicates.end()), rVectorWithDuplicates.end());
}

/*
 * Method to determine whether an element contains ghost nodes
 */
template<unsigned DIM>
bool PeriodicCryptModelInteractionForceWithGhostNodes<DIM>::DoesElementContainGhostNodes(AbstractCellPopulation<DIM>& rCellPopulation, unsigned elementIndex)
{
	DomMeshBasedCellPopulationWithGhostNodes<DIM>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

	bool element_contains_ghost_nodes = false;

	// Get a pointer to the element
	Element<DIM,DIM>* p_element = p_tissue->rGetMesh().GetElement(elementIndex);

	// ITERATE OVER NODES owned by this element
	for (unsigned local_index=0; local_index<DIM+1; local_index++)
	{
		if (p_tissue->IsGhostNode(p_element->GetNodeGlobalIndex(local_index)) == true)
		{
			element_contains_ghost_nodes = true;

		}
	}

	return element_contains_ghost_nodes;
}

/*
 * A method to return the number of elements that contain a particular node,
 * excluding those elements that have ghost nodes
 */
template<unsigned DIM>
unsigned PeriodicCryptModelInteractionForceWithGhostNodes<DIM>::GetNumContainingElementsWithoutGhostNodes(AbstractCellPopulation<DIM>& rCellPopulation, unsigned nodeIndex)
{
	DomMeshBasedCellPopulationWithGhostNodes<DIM>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

    // Get pointer to the node
    Node<DIM>* p_node = p_tissue->GetNode(nodeIndex);
    assert(!(p_tissue->IsGhostNode(nodeIndex)));

    unsigned num_elements_with_no_ghost_nodes = 0;		// Initialise

    // Iterate over containing elements
    for (typename Node<DIM>::ContainingElementIterator iter = p_node->ContainingElementsBegin();
         iter != p_node->ContainingElementsEnd(); ++iter)
    {
        bool element_contains_ghost_nodes = false;
        Element<DIM,DIM>* p_element = p_tissue->rGetMesh().GetElement(*iter);

        // Iterate over nodes owned by this element
        for (unsigned local_index=0; local_index<DIM+1; local_index++)
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
template<unsigned DIM>
std::set<unsigned> PeriodicCryptModelInteractionForceWithGhostNodes<DIM>::GetNeighbouringNodeIndices(AbstractCellPopulation<DIM>& rCellPopulation, unsigned nodeIndex)
{
	DomMeshBasedCellPopulationWithGhostNodes<DIM>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

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
template<unsigned DIM>
bool PeriodicCryptModelInteractionForceWithGhostNodes<DIM>::HasEpithelialCellDetachedFromBasementMembrane(AbstractCellPopulation<DIM>& rCellPopulation, unsigned nodeIndex)
{
	DomMeshBasedCellPopulationWithGhostNodes<DIM>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

	bool has_cell_detached = false;	// Initialising

   	std::set<unsigned> neighbours = GetNeighbouringNodeIndices(rCellPopulation, nodeIndex);

   	unsigned num_stromal_neighbours = 0;

   	// Iterate over the neighbouring cells to check the number of differentiated cell neighbours

   	for(std::set<unsigned>::iterator neighbour_iter=neighbours.begin();
   							neighbour_iter != neighbours.end();
   							++neighbour_iter)
	{
   		boost::shared_ptr<AbstractCellMutationState> p_state = p_tissue->GetCellUsingLocationIndex(*neighbour_iter)->GetMutationState();
		if ( (!p_tissue->IsGhostNode(*neighbour_iter)) && (p_state->IsType<StromalCellMutationState>()==true) )
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

template<unsigned DIM> 
//void PeriodicCryptModelInteractionForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
//                                                                   AbstractCellPopulation<DIM>& rCellPopulation)
void PeriodicCryptModelInteractionForceWithGhostNodes<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // If the width of the periodic domain has not been specified, use the initial width of the cell population
    if (mPeriodicDomainWidth == DOUBLE_UNSET)
    {
        mPeriodicDomainWidth = rCellPopulation.GetWidth(0);
    }

    mExtendedMeshNodeIndexMap.clear();

    // Create a vector of nodes for use in constructing mpExtendedMesh
    unsigned num_cells = rCellPopulation.GetNumRealCells();
    std::vector<Node<DIM>*> extended_nodes(DIM*num_cells);

    // We iterate over all cells in the population
    unsigned count = 0;

    switch (DIM)
    {
        case 1:
        {
            break;
        }
        case 2:
        {
            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                 cell_iter != rCellPopulation.End();
                 ++cell_iter)
            {
                // First, create and store a copy of this real node and cell
                unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                c_vector<double, DIM> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

                // Create a copy of the node corresponding to this cell and store it
                Node<DIM>* p_real_node = new Node<DIM>(real_node_index, real_node_location);
                extended_nodes[count] = p_real_node;

                // Compute the location of the image node corresponding to this node
                c_vector<double,DIM> image_node_location = real_node_location;
                if (real_node_location[0] >= mPeriodicDomainWidth*0.5)
                {
                    image_node_location[0] -= mPeriodicDomainWidth;
                }
                else if (real_node_location[0] <  mPeriodicDomainWidth*0.5)
                {
                    image_node_location[0] += mPeriodicDomainWidth;
                }

                // Create a copy of the node corresponding to this cell, suitable translated, and store it
                Node<DIM>* p_image_node = new Node<DIM>(num_cells+count, image_node_location);
                extended_nodes[num_cells+count] = p_image_node;

                // Populate mExtendedMeshNodeIndexMap
                mExtendedMeshNodeIndexMap[count] = real_node_index;
                mExtendedMeshNodeIndexMap[num_cells+count] = real_node_index;

                count++;
            }

            // We now construct mpExtendedMesh using extended_nodes (and avoid memory leaks)
            if (mpExtendedMesh != NULL)
            {
            	delete mpExtendedMesh;
            }
            if (mpExtendedMesh != NULL)
            {
            	delete mpExtendedMesh;
            }
            mpExtendedMesh = new MutableMesh<DIM,DIM>(extended_nodes);

        	// Now loop over the extended mesh and calculate the force acting on real nodes
        	// (using the edge iterator ensures that each edge is visited exactly once)
            for (typename MutableMesh<DIM,DIM>::EdgeIterator edge_iterator = mpExtendedMesh->EdgesBegin();
                 edge_iterator != mpExtendedMesh->EdgesEnd();
                 ++edge_iterator)
            {
                unsigned nodeA_global_index = edge_iterator.GetNodeA()->GetIndex();
                unsigned nodeB_global_index = edge_iterator.GetNodeB()->GetIndex();

                c_vector<double, DIM> force = CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index, rCellPopulation);

                // Apply this force to any real nodes (i.e. nodes whose indices are less than num_real_nodes)
                if (nodeA_global_index < num_cells)
                {
                    unsigned real_node_index_A = mExtendedMeshNodeIndexMap[nodeA_global_index];
                    //rForces[real_node_index_A] += force;
                   // rCellPopulation.GetNode(real_node_index_A)->AddAppliedForceContribution(force);
                    rCellPopulation.GetNode(real_node_index_A)->AddAppliedForceContribution(force);
                }
                if (nodeB_global_index < num_cells)
                {
                    unsigned real_node_index_B = mExtendedMeshNodeIndexMap[nodeB_global_index];
                    //rForces[real_node_index_B] -= force;
                   // rCellPopulation.GetNode(real_node_index_B)->AddAppliedForceContribution(-force);
                    rCellPopulation.GetNode(real_node_index_B)->AddAppliedForceContribution(-force);
                }
            }

			break;
        }
        case 3:
        {
            // If the width of the periodic domain has not been specified, use the initial width of the cell population
            if (mPeriodicDomainDepth == DOUBLE_UNSET)
            {
            	mPeriodicDomainDepth = rCellPopulation.GetWidth(1);
            }

            // First, extend the mesh in the x-direction

            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                 cell_iter != rCellPopulation.End();
                 ++cell_iter)
            {
                //p_tissue->IsGhostNode(p_element->GetNodeGlobalIndex(local_index))
                // First, create and store a copy of this real node and cell
                unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                c_vector<double, DIM> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

                // Create a copy of the node corresponding to this cell and store it
                Node<DIM>* p_real_node = new Node<DIM>(real_node_index, real_node_location);
                extended_nodes[count] = p_real_node;

                // Compute the location of the image node corresponding to this node
                c_vector<double,DIM> image_node_location = real_node_location;
                if (real_node_location[0] >= mPeriodicDomainWidth*0.5)
                {
                    image_node_location[0] -= mPeriodicDomainWidth;
                }
                else if (real_node_location[0] <  mPeriodicDomainWidth*0.5)
                {
                    image_node_location[0] += mPeriodicDomainWidth;
                }

                // Create a copy of the node corresponding to this cell, suitable translated, and store it
                Node<DIM>* p_image_node = new Node<DIM>(num_cells+count, image_node_location);
                extended_nodes[num_cells+count] = p_image_node;

                // Populate mExtendedMeshNodeIndexMap
                mExtendedMeshNodeIndexMap[count] = real_node_index;
                mExtendedMeshNodeIndexMap[num_cells+count] = real_node_index;

                count++;
            }

            // Second, extend this extended mesh in the y-direction, so that we cover the corners too
            // (We don't need to store the real nodes anymore

            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                 cell_iter != rCellPopulation.End();
                 ++cell_iter)
            {
                // First, create and store a copy of this real node and cell
                unsigned real_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                c_vector<double, DIM> real_node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

                // Compute the location of the image node corresponding to this node
                c_vector<double,DIM> image_node_location = real_node_location;

                if (real_node_location[1] >= mPeriodicDomainDepth*0.5)
                {
                    image_node_location[1] -= mPeriodicDomainDepth;
                }
                else if (real_node_location[1] <  mPeriodicDomainDepth*0.5)
                {
                    image_node_location[1] += mPeriodicDomainDepth;
                }

                // Create a copy of the node corresponding to this cell, suitable translated, and store it
                Node<DIM>* p_image_node = new Node<DIM>(num_cells+count, image_node_location);
                extended_nodes[num_cells+count] = p_image_node;

                // Populate mExtendedMeshNodeIndexMap
                mExtendedMeshNodeIndexMap[num_cells+count] = real_node_index;

                count++;
            }

            // We now construct mpExtendedMesh using extended_nodes
            if (mpExtendedMesh != NULL)
            {
            	delete mpExtendedMesh;
            }
            mpExtendedMesh = new MutableMesh<DIM,DIM>(extended_nodes);

        	// Now loop over the extended mesh and calculate the force acting on real nodes
        	// (using the edge iterator ensures that each edge is visited exactly once)
            for (typename MutableMesh<DIM,DIM>::EdgeIterator edge_iterator = mpExtendedMesh->EdgesBegin();
                 edge_iterator != mpExtendedMesh->EdgesEnd();
                 ++edge_iterator)
            {
                unsigned nodeA_global_index = edge_iterator.GetNodeA()->GetIndex();
                unsigned nodeB_global_index = edge_iterator.GetNodeB()->GetIndex();

                c_vector<double, DIM> force = CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index, rCellPopulation);
                c_vector<double, DIM> neg_force = -1.0*force;

                // Apply this force to any real nodes (i.e. nodes whose indices are less than num_real_nodes)
                if (nodeA_global_index < num_cells)
                {
                    // easy check to print neighbours
                    /*
                    if (nodeA_global_index==19)
                    {
                        std::cout<< nodeB_global_index << " \n ";
                        PRINT_VECTOR(force);
                    }
                    */

                    unsigned real_node_index_A = mExtendedMeshNodeIndexMap[nodeA_global_index];
                    //rForces[real_node_index_A] += force;
                   // rCellPopulation.GetNode(real_node_index_A)->AddAppliedForceContribution(force);
                    rCellPopulation.GetNode(real_node_index_A)->AddAppliedForceContribution(-neg_force);
                }
                if (nodeB_global_index < num_cells)
                {
                    // easy check to print neighbours
                    /*
                    if (nodeB_global_index==19)
                    {
                        std::cout<< nodeA_global_index << " \n ";
                        PRINT_VECTOR(neg_force);
                    }
                    */
                    unsigned real_node_index_B = mExtendedMeshNodeIndexMap[nodeB_global_index];
                    //rForces[real_node_index_B] -= force;
                   // rCellPopulation.GetNode(real_node_index_B)->AddAppliedForceContribution(neg_force);
                    rCellPopulation.GetNode(real_node_index_B)->AddAppliedForceContribution(-force);
                    
                }
            }
           
			break;
        }

        default:
            // This can't happen
            NEVER_REACHED;
    }
    

}


template<unsigned DIM>
c_vector<double, DIM> PeriodicCryptModelInteractionForceWithGhostNodes<DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
																						 unsigned nodeBGlobalIndex,
                                                                                         AbstractCellPopulation<DIM>& rCellPopulation)
{
    assert(dynamic_cast<DomMeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation) != nullptr);
 //   assert(rCellPopulation.IsMeshBasedCellPopulation());
 //   assert(bool(dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation)))


    // We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex != nodeBGlobalIndex);

    // Get the node locations
    c_vector<double, DIM> node_a_location = this->mpExtendedMesh->GetNode(nodeAGlobalIndex)->rGetLocation();
    c_vector<double, DIM> node_b_location = this->mpExtendedMesh->GetNode(nodeBGlobalIndex)->rGetLocation();

    /*
     * Get the unit vector parallel to the line joining the two nodes.
     *
     * We use the mesh method GetVectorFromAtoB() to compute the direction of the
     * unit vector along the line joining the two nodes, rather than simply subtract
     * their positions, because this method can be overloaded (e.g. to enforce a
     * periodic boundary in Cylindrical2dMesh).
     */
    c_vector<double, DIM> unit_difference = this->mpExtendedMesh->GetVectorFromAtoB(node_a_location, node_b_location);

    // Calculate the distance between the two nodes
    double distance_between_nodes = norm_2(unit_difference);
    assert(distance_between_nodes > 0);
    assert(!std::isnan(distance_between_nodes));

    unit_difference /= distance_between_nodes;

    /*
     * If mUseCutoffLength has been set, then there is zero force between
     * two nodes located a distance apart greater than mUseCutoffPoint.
     */
    if (this->mUseCutOffLength)
    {
        if (distance_between_nodes >= this->GetCutOffLength())
        {
            return zero_vector<double>(DIM); // c_vector<double,DIM>() is not guaranteed to be fresh memory
        }
    }

    // Calculate the rest length of the spring connecting the two nodes

    double rest_length = 1.0;

    ///\todo Extend force class to cope with newly divided cells (#1856)

    double a_rest_length = rest_length*0.5;
    double b_rest_length = a_rest_length;

    ///\todo Extend force class to cope with apoptotic cells (#1856)

    rest_length = a_rest_length + b_rest_length;
    assert(rest_length <= 1.0+1e-12);

	// Get the corresponding node index in rCellPopulation
    unsigned real_A_node_index = mExtendedMeshNodeIndexMap[nodeAGlobalIndex];
    unsigned real_B_node_index = mExtendedMeshNodeIndexMap[nodeBGlobalIndex];

    // Get the corresponding cells
    CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(real_A_node_index);
    CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(real_B_node_index);
//
//    double ageA = p_cell_A->GetAge();
//    double ageB = p_cell_B->GetAge();
//
//    assert(!std::isnan(ageA));
//    assert(!std::isnan(ageB));
//
//    /*
//     * If the cells are both newly divided, then the rest length of the spring
//     * connecting them grows linearly with time, until 1 hour after division.
//     */
//    if (ageA < this->mMeinekeSpringGrowthDuration && ageB < this->mMeinekeSpringGrowthDuration)
//    {
//        if (rCellPopulation.IsMeshBasedCellPopulation())
//        {
//            MeshBasedCellPopulation<DIM>* p_static_cast_cell_population = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);
//
//            std::pair<CellPtr,CellPtr> cell_pair = p_static_cast_cell_population->CreateCellPair(p_cell_A, p_cell_B);
//
//            if (p_static_cast_cell_population->IsMarkedSpring(cell_pair))
//            {
//                // Spring rest length increases from a small value to the normal rest length over 1 hour
//                double lambda = this->GetMeinekeDivisionRestingSpringLength();
//                rest_length = lambda + (1.0 - lambda) * ageA/this->mMeinekeSpringGrowthDuration;
//            }
//            if (ageA + SimulationTime::Instance()->GetTimeStep() >= this->mMeinekeSpringGrowthDuration)
//            {
//                // This spring is about to go out of scope
//                p_static_cast_cell_population->UnmarkSpring(cell_pair);
//            }
//        }
//        else
//        {
//            // Spring rest length increases from mDivisionRestingSpringLength to normal rest length, 1.0, over 1 hour
//            double lambda = this->GetMeinekeDivisionRestingSpringLength();
//            rest_length = lambda + (1.0 - lambda) * ageA/this->mMeinekeSpringGrowthDuration;
//        }
//    }
//
//    double a_rest_length = rest_length*0.5;
//    double b_rest_length = a_rest_length;
//
//    /*
//     * If either of the cells has begun apoptosis, then the length of the spring
//     * connecting them decreases linearly with time.
//     */
//    if (p_cell_A->HasApoptosisBegun())
//    {
//        double time_until_death_a = p_cell_A->GetTimeUntilDeath();
//        a_rest_length = a_rest_length * time_until_death_a / p_cell_A->GetApoptosisTime();
//    }
//    if (p_cell_B->HasApoptosisBegun())
//    {
//        double time_until_death_b = p_cell_B->GetTimeUntilDeath();
//        b_rest_length = b_rest_length * time_until_death_b / p_cell_B->GetApoptosisTime();
//    }
//

    bool is_closer_than_rest_length = (distance_between_nodes - rest_length <= 0);

    // Although in this class the 'spring constant' is a constant parameter, in
    // subclasses it can depend on properties of each of the cells
    double multiplication_factor = VariableSpringConstantMultiplicationFactor(this->mExtendedMeshNodeIndexMap[nodeAGlobalIndex], this->mExtendedMeshNodeIndexMap[nodeBGlobalIndex], rCellPopulation, is_closer_than_rest_length);
    double spring_stiffness = this->GetMeinekeSpringStiffness();
    double overlap = distance_between_nodes - rest_length;


    /* Want to have one-way springs between epithelial and stromal nodes, so that there is only repulsion due to compression
     * of the spring, but no attraction due to extension
     */
    if ( (mUseOneWaySprings) && ( ( (p_cell_A->GetMutationState()->IsType<StromalCellMutationState>() == false) && (p_cell_B->GetMutationState()->IsType<StromalCellMutationState>() == true) )
    	    || ( (p_cell_A->GetMutationState()->IsType<StromalCellMutationState>() == true) && (p_cell_B->GetMutationState()->IsType<StromalCellMutationState>() == false) ) ) )
    {
        if (distance_between_nodes > rest_length)
        {
        	return zero_vector<double>(DIM); // c_vector<double,DIM>() is not guaranteed to be fresh memory
        }
    }

 //   if (rCellPopulation.IsMeshBasedCellPopulation())
 // if(bool(dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation) != nullptr))
 if (bool(dynamic_cast<DomMeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation)))
    {
        return multiplication_factor * spring_stiffness * unit_difference * overlap;
    }
    else
    {
        // A reasonably stable simple force law
        if (distance_between_nodes > rest_length)
        {
            double alpha = 5;
            c_vector<double, DIM> temp = spring_stiffness * unit_difference * overlap * exp(-alpha * overlap);
            for (unsigned i=0; i<DIM; i++)
            {
                assert(!std::isnan(temp[i]));
            }
            return temp;
        }
        else
        {
            c_vector<double, DIM> temp = spring_stiffness * unit_difference * log(1 + overlap);
            for (unsigned i=0; i<DIM; i++)
            {
                assert(!std::isnan(temp[i]));
            }
            return temp;
        }
    }
}


template<unsigned DIM>
double PeriodicCryptModelInteractionForceWithGhostNodes<DIM>::GetPeriodicDomainWidth()
{
	return mPeriodicDomainWidth;
}

template<unsigned DIM>
void PeriodicCryptModelInteractionForceWithGhostNodes<DIM>::SetPeriodicDomainWidth(double periodicDomainWidth)
{
	mPeriodicDomainWidth = periodicDomainWidth;
}

template<unsigned DIM>
double PeriodicCryptModelInteractionForceWithGhostNodes<DIM>::GetPeriodicDomainDepth()
{
	return mPeriodicDomainDepth;
}

template<unsigned DIM>
void PeriodicCryptModelInteractionForceWithGhostNodes<DIM>::SetPeriodicDomainDepth(double periodicDomainDepth)
{
	mPeriodicDomainDepth = periodicDomainDepth;
}


template<unsigned DIM>
void PeriodicCryptModelInteractionForceWithGhostNodes<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile <<  "\t\t\t<UseCellTypeDependentSprings>"<< mUseCellTypeDependentSprings << "</UseCellTypeDependentSprings> \n" ;
	*rParamsFile <<  "\t\t\t<TransitTransitMultiplier>"<< mTransitTransitMultiplier << "</TransitTransitMultiplier> \n" ;
	*rParamsFile <<  "\t\t\t<DifferentiatedDifferentiatedMultiplier>"<< mDifferentiatedDifferentiatedMultiplier << "</DifferentiatedDifferentiatedMultiplier> \n" ;
	*rParamsFile <<  "\t\t\t<TransitDifferentiatedMultiplier>"<< mTransitDifferentiatedMultiplier << "</TransitDifferentiatedMultiplier> \n" ;
	*rParamsFile <<  "\t\t\t<UseEpithelialStromalCellDependentSprings>"<< mUseEpithelialStromalCellDependentSprings << "</UseEpithelialStromalCellDependentSprings> \n" ;
	*rParamsFile <<  "\t\t\t<EpithelialEpithelialMultiplier>"<< mEpithelialEpithelialMultiplier << "</EpithelialEpithelialMultiplier> \n" ;
	*rParamsFile <<  "\t\t\t<StromalStromalMultiplier>"<< mStromalStromalMultiplier << "</StromalStromalMultiplier> \n" ;
	*rParamsFile <<  "\t\t\t<EpithelialStromalMultiplier>"<< mEpithelialStromalMultiplier << "</EpithelialStromalMultiplier> \n" ;
	*rParamsFile <<  "\t\t\t<ApcTwoHitStromalMultiplier>"<< mApcTwoHitStromalMultiplier << "</ApcTwoHitStromalMultiplier> \n" ;
	*rParamsFile <<  "\t\t\t<UseOneWaySprings>"<<  mUseOneWaySprings << "</mUseOneWaySprings> \n" ;
	*rParamsFile <<  "\t\t\t<UseEdgeBasedSpringConstant>"<<  mUseEdgeBasedSpringConstant << "</UseEdgeBasedSpringConstant> \n" ;
	*rParamsFile <<  "\t\t\t<PeriodicDomainWidth>"<<  mPeriodicDomainWidth << "</PeriodicDomainWidth> \n" ;
	*rParamsFile <<  "\t\t\t<PeriodicDomainDepth>"<<  mPeriodicDomainDepth << "</PeriodicDomainDepth> \n" ;

	// Call direct parent class
	LinearSpringWithVariableSpringConstantsForce<DIM>::OutputForceParameters(rParamsFile);
}



/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class PeriodicCryptModelInteractionForceWithGhostNodes<1>;
template class PeriodicCryptModelInteractionForceWithGhostNodes<2>;
template class PeriodicCryptModelInteractionForceWithGhostNodes<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PeriodicCryptModelInteractionForceWithGhostNodes)

