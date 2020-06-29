#include "CellCentreWithVariableCellInteraction.hpp"
#include "Debug.hpp"
#include <iostream>
#include <fstream>

/* This global function is to allow the vectors of coordinates to be compared in
 *  terms of their x-value
 */
bool SortAccordingToXCoordinate(const c_vector<double, 2> lhs, const c_vector<double, 2> rhs)
{
    return lhs[0] < rhs[0];
}


template<unsigned DIM>
CellCentreWithVariableCellInteraction<DIM>::CellCentreWithVariableCellInteraction()
   : GeneralisedLinearSpringForce<DIM>(),  
	mUseCellTypeDependentSprings(false),	// Will need to override this??
	mTransitTransitMultiplier(DOUBLE_UNSET),
	mDifferentiatedDifferentiatedMultiplier(DOUBLE_UNSET),
	mTransitDifferentiatedMultiplier(DOUBLE_UNSET),
    mUseEdgeBasedSpringConstant(false),
    mUseMutantSprings(false),
    mMutantMutantMultiplier(DOUBLE_UNSET),
    mNormalMutantMultiplier(DOUBLE_UNSET),
    mUseOneWaySprings(false),
    mUseBCatSprings(false),
    mBetaCatSpringScaler(18.14/6.0), // scale spring constant with beta-catenin level (divided by 6 for heaxagonal cells)
    mUseApoptoticSprings(false),
    mApoptoticSpringTensionStiffness(15.0*0.25),
    mApoptoticSpringCompressionStiffness(15.0*0.75)
    {
    }

template<unsigned DIM>
CellCentreWithVariableCellInteraction<DIM>::~CellCentreWithVariableCellInteraction()
{
}

template<unsigned DIM>
void CellCentreWithVariableCellInteraction<DIM>::SetCellTypeDependentSprings(bool useCellTypeDependentSprings,
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
void CellCentreWithVariableCellInteraction<DIM>::SetEdgeBasedSpringConstant(bool useEdgeBasedSpringConstant)
{
    assert(DIM == 2);
    mUseEdgeBasedSpringConstant = useEdgeBasedSpringConstant;
}

template<unsigned DIM>
void CellCentreWithVariableCellInteraction<DIM>::SetMutantSprings(bool useMutantSprings, double mutantMutantMultiplier, double normalMutantMultiplier)
{
    mUseMutantSprings = useMutantSprings;
    mMutantMutantMultiplier = mutantMutantMultiplier;
    mNormalMutantMultiplier = normalMutantMultiplier;
}

template<unsigned DIM>
void CellCentreWithVariableCellInteraction<DIM>::SetBetaCateninSprings(bool useBCatSprings)
{
    mUseBCatSprings = useBCatSprings;
}

template<unsigned DIM>
void CellCentreWithVariableCellInteraction<DIM>::SetApoptoticSprings(bool useApoptoticSprings)
{
    mUseApoptoticSprings = useApoptoticSprings;
}

template<unsigned DIM>
void CellCentreWithVariableCellInteraction<DIM>::SetBasalLaminaParameter(double basalLaminaParameter)
{
	mBasalLaminaParameter = basalLaminaParameter;
}

template<unsigned DIM>
double CellCentreWithVariableCellInteraction<DIM>::GetBasalLaminaParameter()
{
	return mBasalLaminaParameter;
}

template<unsigned DIM>
void CellCentreWithVariableCellInteraction<DIM>::SetUseOneWaySprings(bool useOneWaySprings)
{
	mUseOneWaySprings = useOneWaySprings;
}

template<unsigned DIM>
void CellCentreWithVariableCellInteraction<DIM>::SetTransitTransitMultiplier(double transitTransitMultiplier)
{
	mTransitTransitMultiplier = transitTransitMultiplier;
}

template<unsigned DIM>
double CellCentreWithVariableCellInteraction<DIM>::GetTransitTransitMultiplier()
{
	return mTransitTransitMultiplier;
}

template<unsigned DIM>
double CellCentreWithVariableCellInteraction<DIM>::VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex,
                                                           AbstractCellPopulation<DIM>& rCellPopulation, bool isCloserThanRestLength)
{	
    double multiplication_factor = GeneralisedLinearSpringForce<DIM>::VariableSpringConstantMultiplicationFactor(nodeAGlobalIndex,
                                                                                                            nodeBGlobalIndex,
                                                                                                            rCellPopulation,
                                                                                                            isCloserThanRestLength);
    
    CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
    CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);
    
    if (mUseEdgeBasedSpringConstant)
    {
        assert(rCellPopulation.IsMeshBasedCellPopulation());
        assert(!mUseBCatSprings);   // don't want to do both (both account for edge length)

        MeshBasedCellPopulation<DIM>* p_tissue = (static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation));

        multiplication_factor = p_tissue->GetVoronoiEdgeLength(nodeAGlobalIndex, nodeBGlobalIndex)*sqrt(3);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (mUseCellTypeDependentSprings)
    {
        CellProliferativeType cell_A_type = p_cell_A->GetCellCycleModel()->GetCellProliferativeType();
        CellProliferativeType cell_B_type = p_cell_B->GetCellCycleModel()->GetCellProliferativeType();

        if (cell_A_type == cell_B_type)
        {
            if (cell_A_type == TRANSIT)
            {
            	multiplication_factor *= mTransitTransitMultiplier;
            }

            if (cell_A_type == DIFFERENTIATED)
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
    
    if (mUseMutantSprings)
    {
        unsigned number_of_mutants = 0;

        if (p_cell_A->GetMutationState()->IsType<ApcTwoHitCellMutationState>()
         || p_cell_A->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())
        {
            // If cell A is mutant
            number_of_mutants++;
        }

        if (p_cell_B->GetMutationState()->IsType<ApcTwoHitCellMutationState>()
         || p_cell_B->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())
        {
            // If cell B is mutant
            number_of_mutants++;
        }

        switch (number_of_mutants)
        {
            case 1u:
            {
                multiplication_factor *= mNormalMutantMultiplier;
                break;
            }
            case 2u:
            {
                multiplication_factor *= mMutantMutantMultiplier;
                break;
            }
        }
    }

    if (mUseBCatSprings)
    {
        assert(rCellPopulation.IsMeshBasedCellPopulation());
        // If using beta-cat dependent springs, both cell-cycle models had better be VanLeeuwen2009WntSwatCellCycleModel
        AbstractVanLeeuwen2009WntSwatCellCycleModel* p_model_A = dynamic_cast<AbstractVanLeeuwen2009WntSwatCellCycleModel*>(p_cell_A->GetCellCycleModel());
        AbstractVanLeeuwen2009WntSwatCellCycleModel* p_model_B = dynamic_cast<AbstractVanLeeuwen2009WntSwatCellCycleModel*>(p_cell_B->GetCellCycleModel());

        assert(!mUseEdgeBasedSpringConstant);   // This already adapts for edge lengths - don't want to do it twice.
        double beta_cat_cell_1 = p_model_A->GetMembraneBoundBetaCateninLevel();
        double beta_cat_cell_2 = p_model_B->GetMembraneBoundBetaCateninLevel();

        MeshBasedCellPopulation<DIM>* p_static_cast_tissue = (static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation));

        double perim_cell_1 = p_static_cast_tissue->GetSurfaceAreaOfVoronoiElement(nodeAGlobalIndex);
        double perim_cell_2 = p_static_cast_tissue->GetSurfaceAreaOfVoronoiElement(nodeBGlobalIndex);
        double edge_length_between_1_and_2 = p_static_cast_tissue->GetVoronoiEdgeLength(nodeAGlobalIndex, nodeBGlobalIndex);

        double beta_cat_on_cell_1_edge = beta_cat_cell_1 *  edge_length_between_1_and_2 / perim_cell_1;
        double beta_cat_on_cell_2_edge = beta_cat_cell_2 *  edge_length_between_1_and_2 / perim_cell_2;

        double min_beta_Cat_of_two_cells = std::min(beta_cat_on_cell_1_edge, beta_cat_on_cell_2_edge);

        multiplication_factor *= min_beta_Cat_of_two_cells / mBetaCatSpringScaler;
    }

    if (mUseApoptoticSprings)
    {
        bool cell_A_is_apoptotic = p_cell_A->HasCellProperty<ApoptoticCellProperty>();
        bool cell_B_is_apoptotic = p_cell_B->HasCellProperty<ApoptoticCellProperty>();

        if (cell_A_is_apoptotic || cell_B_is_apoptotic)
        {
            double spring_a_stiffness = 2.0 * this->GetMeinekeSpringStiffness();
            double spring_b_stiffness = 2.0 * this->GetMeinekeSpringStiffness();

            if (cell_A_is_apoptotic)
            {
                if (!isCloserThanRestLength) // if under tension
                {
                    spring_a_stiffness = mApoptoticSpringTensionStiffness;
                }
                else // if under compression
                {
                    spring_a_stiffness = mApoptoticSpringCompressionStiffness;
                }
            }
            if (cell_B_is_apoptotic)
            {
                if (!isCloserThanRestLength) // if under tension
                {
                    spring_b_stiffness = mApoptoticSpringTensionStiffness;
                }
                else // if under compression
                {
                    spring_b_stiffness = mApoptoticSpringCompressionStiffness;
                }
            }

            multiplication_factor /= (1.0/spring_a_stiffness + 1.0/spring_b_stiffness)*this->GetMeinekeSpringStiffness();
        }
    }
    
    return multiplication_factor;
}

template<unsigned DIM>
void CellCentreWithVariableCellInteraction<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile <<  "\t\t<mUseEdgeBasedSpringConstant> " <<  mUseEdgeBasedSpringConstant << " </mUseEdgeBasedSpringConstant> \n" ;
	*rParamsFile <<  "\t\t<mUseMutantSprings> " <<  mUseMutantSprings << " </mUseMutantSprings> \n" ;
	*rParamsFile <<  "\t\t<mMutantMutantMultiplier> " <<  mMutantMutantMultiplier << " </mMutantMutantMultiplier> \n" ;
	*rParamsFile <<  "\t\t<mNormalMutantMultiplier> " <<  mNormalMutantMultiplier << " </mNormalMutantMultiplier> \n" ;
	*rParamsFile <<  "\t\t<mUseBCatSprings> " <<  mUseBCatSprings << " </mUseBCatSprings> \n" ;
	*rParamsFile <<  "\t\t<mUseApoptoticSprings> " <<  mUseApoptoticSprings << " </mUseApoptoticSprings> \n" ;
	*rParamsFile <<  "\t\t<mBetaCatSpringScaler> " <<  mBetaCatSpringScaler << " </mBetaCatSpringScaler> \n" ;
	*rParamsFile <<  "\t\t<mApoptoticSpringTensionStiffness> " <<  mApoptoticSpringTensionStiffness << " </mApoptoticSpringTensionStiffness> \n" ;
	*rParamsFile <<  "\t\t<mApoptoticSpringCompressionStiffness> " <<  mApoptoticSpringCompressionStiffness << " </mApoptoticSpringCompressionStiffness> \n" ;

	// Call direct parent class
	GeneralisedLinearSpringForce<DIM>::OutputForceParameters(rParamsFile);
}


template<unsigned DIM>
void CellCentreWithVariableCellInteraction<DIM>::RemoveDuplicates1D(std::vector<unsigned>& vectorWithDuplicates)
{
    std::sort(vectorWithDuplicates.begin(), vectorWithDuplicates.end());
    vectorWithDuplicates.erase(std::unique(vectorWithDuplicates.begin(), vectorWithDuplicates.end()), vectorWithDuplicates.end());
}

/* A method to find all the pairs of connections between transit cells and differentiated cells.
 * Returns a vector of node pairings, without repeats
 */

template<unsigned DIM>
std::vector<c_vector<unsigned, DIM> > CellCentreWithVariableCellInteraction<DIM>::GetTransitDifferentiatedPairs(AbstractCellPopulation<DIM>& rCellPopulation)
{
    MeshBasedCellPopulation<DIM>* p_tissue = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

    // Create a vector to record the pairs of nodes corresponding to joined transit and differentiated nodes
    std::vector<c_vector<unsigned, DIM> > node_pairs;
    c_vector<double, DIM> pair;

    // We iterate over all cells in the tissue, and deal only with those that are transit cells

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
    	CellProliferativeType cell_type = cell_iter->GetCellCycleModel()->GetCellProliferativeType();
    	  	    	
    	if(cell_type == TRANSIT)
    	{
    		Node<DIM>* p_node = p_tissue->GetNodeCorrespondingToCell(*cell_iter);	
    		unsigned node_index = p_node->GetIndex();	
    		    	    	
    		// ITERATE OVER CONTAINING ELEMENTS and only work with those that DO NOT contain ghost nodes
	
    		std::vector<unsigned> diff_nodes;
    		
    		for (typename Node<DIM>::ContainingElementIterator iter = p_node->ContainingElementsBegin();
		         						iter != p_node->ContainingElementsEnd();
		         						++iter)
    		{		
    			bool element_contains_ghost_nodes = false;
					             
    			// Get a pointer to the element
    			Element<DIM,DIM>* p_element = p_tissue->rGetMesh().GetElement(*iter);                
					            			   			
    			// ITERATE OVER NODES owned by this element
    			for (unsigned local_index=0; local_index<DIM+1; local_index++)
    			{        	           	        		
    				if (p_tissue->IsGhostNode(p_element->GetNodeGlobalIndex(local_index)) == true)		
    				{							
    					element_contains_ghost_nodes = true;			    
    					break; 				// This should break out of the inner for loop						   
    				}			        			
				
    				if (element_contains_ghost_nodes==false)
    				{
    					unsigned nodeBGlobalIndex = p_element->GetNodeGlobalIndex(local_index);
    					CellProliferativeType cellB_type = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex)->GetCellCycleModel()->GetCellProliferativeType();  								
    					    					
    					if(cellB_type == DIFFERENTIATED)
    					{
    						// Store the index of each differentiated node that is attached to the transit
    						// node. There will be repetitions due to iterating over neighbouring elements
    						diff_nodes.push_back(nodeBGlobalIndex);	
    					}															    				    					
    				}    				    				
    			}
    		}    	
    		
    		// Remove any nodes that have been found twice
    		RemoveDuplicates1D(diff_nodes);    		    	
    		
    		// Now construct the vector of node pairs
    		for (unsigned i=0; i<diff_nodes.size(); i++)
    		{
    			pair[0] = node_index;
    			pair[1] = diff_nodes[i];
    			node_pairs.push_back(pair);
    		}    		    	
    	}
    }	    
	return node_pairs;	
}


/* Method to determine whether an element contains ghost nodes
 */

template<unsigned DIM>
bool CellCentreWithVariableCellInteraction<DIM>::DoesElementContainGhostNodes(AbstractCellPopulation<DIM>& rCellPopulation, unsigned elementIndex)
{
	MeshBasedCellPopulation<DIM>* p_tissue = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

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


/* Using the vector of node pairs found in GetTransitDifferentiatedPairs to determine the curvature
 * of the curve that passes through the midpoints of the neighbouring springs, given a transit-diff
 * node pairing. Note - it ignores the end pairs because one of the common elements will contain ghost nodes.
 */

template<unsigned DIM>
double CellCentreWithVariableCellInteraction<DIM>::GetCurvatureFromNodePair(AbstractCellPopulation<DIM>& rCellPopulation, unsigned transitNodeIndex, unsigned diffNodeIndex)
{
	MeshBasedCellPopulation<DIM>* p_tissue = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);
		
    // Seeking the common elements that contain both the transit node and the differentiated node
	// Note: Don't want to count any elements that have ghost nodes 
	
	std::vector<unsigned> common_elements;	// Initialising
    
    std::set<unsigned> transit_elements = rCellPopulation.GetNode(transitNodeIndex)->rGetContainingElementIndices();
            
    // Loop over all elements that contain the differentiated node
    
    for (typename Node<DIM>::ContainingElementIterator elt_it = rCellPopulation.GetNode(diffNodeIndex)->ContainingElementsBegin();
                 elt_it != rCellPopulation.GetNode(diffNodeIndex)->ContainingElementsEnd();
                 ++elt_it)
    {    
    	unsigned elt_index = *elt_it;
    	
    	bool elt_contains_ghost_nodes = DoesElementContainGhostNodes(rCellPopulation, elt_index);
                
    	// Keep only those elements that also contain the transit node, but do not have ghost nodes
    	
        if ((transit_elements.find(elt_index) != transit_elements.end()) && (elt_contains_ghost_nodes == false))
        {     
        	// Common element
           	common_elements.push_back(elt_index);               
        }
    }
    
	if(common_elements.size() == 0)
	{
		TRACE(SimulationTime::Instance()->GetTime());
	}

    // If there is only one such common element, then this transit node will be at either end and so we don't 
    // consider the force along the very first / very last spring
    
    if(common_elements.size() == 1)
    {
    	double curvature = 0.0;
    	return curvature;
    }
    
    else
    {           
    	assert(common_elements.size()==2);		// Should only be two common elements
    
    	// We iterate over these common elements to find the midpoints of the springs that 
    	// connect transit and  differentiated nodes
	
    	std::vector<c_vector<double, DIM> > spring_midpoints;		// Vector to contain the spring midpoints
    
    	for (std::vector<unsigned>::iterator elem_iter = common_elements.begin();
    									   	 elem_iter != common_elements.end();
    									   	 ++elem_iter)
    	{	    		
    		// Need to determine the cell type of each local node in the element
    		// Want to find the midpoint between transit and differentiated pairs
    	   	
    		// Looking at the three nodes which form the vertices of the element
    		
    		unsigned global_index_A = p_tissue->rGetMesh().GetElement(*elem_iter)->GetNodeGlobalIndex(0);
	   		unsigned global_index_B = p_tissue->rGetMesh().GetElement(*elem_iter)->GetNodeGlobalIndex(1);
	   		unsigned global_index_C = p_tissue->rGetMesh().GetElement(*elem_iter)->GetNodeGlobalIndex(2);
	   		
	   		CellProliferativeType cell_type_A = p_tissue->GetCellUsingLocationIndex(global_index_A)->GetCellCycleModel()->GetCellProliferativeType();
	   		CellProliferativeType cell_type_B = p_tissue->GetCellUsingLocationIndex(global_index_B)->GetCellCycleModel()->GetCellProliferativeType();
	   		CellProliferativeType cell_type_C = p_tissue->GetCellUsingLocationIndex(global_index_C)->GetCellCycleModel()->GetCellProliferativeType();
	   		
	   		if((cell_type_A==TRANSIT && cell_type_B==DIFFERENTIATED) || (cell_type_B==TRANSIT && cell_type_A==DIFFERENTIATED))
	   		{
	   			c_vector<double, DIM> midpoint = 0.5*(p_tissue->GetNode(global_index_A)->rGetLocation() + p_tissue->GetNode(global_index_B)->rGetLocation());
	   			spring_midpoints.push_back(midpoint);
	   		}
	   		
	   		if((cell_type_A==TRANSIT && cell_type_C==DIFFERENTIATED) || (cell_type_C==TRANSIT && cell_type_A==DIFFERENTIATED))
	   		{
	   			c_vector<double, DIM> midpoint = 0.5*(p_tissue->GetNode(global_index_A)->rGetLocation() + p_tissue->GetNode(global_index_C)->rGetLocation());
	   			spring_midpoints.push_back(midpoint);
	   		}
	   		
	   		if((cell_type_B==TRANSIT && cell_type_C==DIFFERENTIATED) || (cell_type_C==TRANSIT && cell_type_B==DIFFERENTIATED))
	   		{
	   			c_vector<double, DIM> midpoint = 0.5*(p_tissue->GetNode(global_index_B)->rGetLocation() + p_tissue->GetNode(global_index_C)->rGetLocation());
	   			spring_midpoints.push_back(midpoint);
	   		}
	   			   	    	
    	}
    	
    	// Now sort the elements of this vector to be in order of increasing x coordinate
    	sort(spring_midpoints.begin(), spring_midpoints.end(), SortAccordingToXCoordinate); 	
    	
    	// After sorting, the middle two entries of the spring_midpoints vector will be the same (it finds
    	// the middle one twice) and so we remove this entry
    	  
    	spring_midpoints.erase(spring_midpoints.begin()+2);
    	
    	assert(spring_midpoints.size()==3); // Should only have 3 midpoints to consider  
        
    	// Now loop over the spring midpoints to find the curvature for each group of three
    		   		   		   			
    	double curvature = FindCurvature(spring_midpoints[0], spring_midpoints[1], spring_midpoints[2]);			
    	   	
    	// Now we're going to add a little to make the curvature non-zero to see what happens    	
//    	curvature = curvature + 0.1;
    	
    	return curvature;
  	
    }
  	
}



/* Gets an ordered vector of curvatures by finding the spring midpoints, sorting in order of 
 * x-coordinate and then looping over in groups of three
 */ 

template<unsigned DIM>
void CellCentreWithVariableCellInteraction<DIM>::GetOrderedCurvatureVector(AbstractCellPopulation<DIM>& rCellPopulation)
{
	MeshBasedCellPopulation<DIM>* p_tissue = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

	std::vector<c_vector<double, DIM> > spring_midpoints;		// Vector to contain the spring midpoints
	std::vector<unsigned> transit_nodes;
	
	unsigned num_transit_nodes = 0; // !!!DOUBLE CHECK THAT THIS IS THE BEST PLACE TO REINITIALISE THIS
									// Needs to be done everytime you work out the forces
	
	// Iterate over the springs
	
	for (typename MeshBasedCellPopulation<DIM>::SpringIterator spring_iterator=p_tissue->SpringsBegin();
	     spring_iterator!=p_tissue->SpringsEnd();
	     ++spring_iterator)
	{			
	    // Get the cell types
	    CellProliferativeType cellA_type = spring_iterator.GetCellA()->GetCellCycleModel()->GetCellProliferativeType();
	    CellProliferativeType cellB_type = spring_iterator.GetCellB()->GetCellCycleModel()->GetCellProliferativeType();

	    if (   (cellA_type==TRANSIT && cellB_type==DIFFERENTIATED) 
	        || (cellA_type==DIFFERENTIATED && cellB_type==TRANSIT) )
	    {	    	    	
	    	num_transit_nodes++;	// Keep track of the number of transit nodes	    		        		    
	    	    		    		    		    	    
	        // This spring connects a transit cell with a differentiated cell
	    	// so work out its midpoint
	    	c_vector<double, DIM> midpoint = 0.5*(spring_iterator.GetNodeA()->rGetLocation() + spring_iterator.GetNodeB()->rGetLocation());
	     	    		    
	    	spring_midpoints.push_back(midpoint);	    		    		    		    	    	
	    }	    	    
	    
	}
	
	// Sort the spring midpoints in order of increasing x coordinate (so that you move from left to right)
	
	sort(spring_midpoints.begin(), spring_midpoints.end(),SortAccordingToXCoordinate);	
	
	// Now loop over the spring midpoints to find the curvature for each group of three
	
	std::vector<double> curvature(spring_midpoints.size()-2);
	
	for (unsigned i = 1; i<(spring_midpoints.size()-1); i++)
	{
		c_vector<double, DIM> left_midpoint = spring_midpoints[i-1];
		c_vector<double, DIM> centre_midpoint = spring_midpoints[i];
		c_vector<double, DIM> right_midpoint = spring_midpoints[i+1];
		
	    curvature[i-1] = FindCurvature(left_midpoint, centre_midpoint, right_midpoint);			
	}

	// This isn't very helpful as the indexing of transit nodes changes when they divide - so
	// it isn't possible to loop over the transit nodes from left to right and apply what was found here

}


/* Function to return the curvature between three points  - the midpoints of the springs connecting the 
 * transit cells to the differentiated cells. NB. The input arguments need to be in order of increasing 
 * x-coordinate so that they correspond literally to 'left', 'middle' and 'right'.
 */

template<unsigned DIM>
double CellCentreWithVariableCellInteraction<DIM>::FindCurvature(c_vector<double, DIM> leftMidpoint,
															  c_vector<double, DIM> centreMidpoint, 
															  c_vector<double, DIM> rightMidpoint)
{
	// Firstly find the size of the intervals by subtracting the x coordinates
	double left_x_interval = centreMidpoint[0] - leftMidpoint[0];
	double left_y_interval = centreMidpoint[1] - leftMidpoint[1];
	double right_x_interval = rightMidpoint[0] - centreMidpoint[0];
	double right_y_interval = rightMidpoint[1] - centreMidpoint[1];
	
	double left_param = pow(pow(left_x_interval,2) + pow(left_y_interval,2),1/2);
	double right_param = pow(pow(right_x_interval,2) + pow(right_y_interval,2),1/2);	
	
	double sum_param = left_param + right_param;
	
	//TS_ASSERT_LESS_THAN(0,left_interval);		// Ensure that these are both positive
	//TS_ASSERT_LESS_THAN(0,right_interval);
	
	double x_prime = (rightMidpoint[0] - leftMidpoint[0])/sum_param;
	double y_prime = (rightMidpoint[1] - leftMidpoint[1])/sum_param;
	
	double x_double_prime = 2*(left_param*rightMidpoint[0] - sum_param*centreMidpoint[0] + right_param*leftMidpoint[0])/(left_param*right_param*sum_param);
	double y_double_prime = 2*(left_param*rightMidpoint[1] - sum_param*centreMidpoint[1] + right_param*leftMidpoint[1])/(left_param*right_param*sum_param);
		
	double curvature = (x_prime*y_double_prime - y_prime*x_double_prime)/pow((pow(x_prime,2) + pow(y_prime,2)),3/2);
	
	return curvature;
}

/* A method to return the number of elements that contain a particular node, 
 * excluding those elements that have ghost nodes 
 */

template<unsigned DIM>
unsigned CellCentreWithVariableCellInteraction<DIM>::GetNumContainingElementsWithoutGhostNodes(AbstractCellPopulation<DIM>& rCellPopulation, unsigned nodeIndex)
{
	MeshBasedCellPopulation<DIM>* p_tissue = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);
	
    // Get pointer to the node
    Node<DIM>* p_node = p_tissue->GetNode(nodeIndex);

    unsigned num_elements_with_no_ghost_nodes = 0;		// Initialise

    // Iterate over containing elements
    for (typename Node<DIM>::ContainingElementIterator iter = p_node->ContainingElementsBegin();
         iter != p_node->ContainingElementsEnd();
         ++iter)
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



template<unsigned DIM>
void CellCentreWithVariableCellInteraction<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                                           AbstractCellPopulation<DIM>& rCellPopulation)
{		
	// Firstly determine the force acting on each transit cell due to the basal lamina  
	// Start by identifying the transit-differentiated pairs
	
	std::vector<c_vector<unsigned, DIM> > node_pairs = GetTransitDifferentiatedPairs(rCellPopulation);
	
	// We loop over the transit-differentiated node pairs to find the force acting on that 
	// transit node, and the direction in which it acts
	
	for (unsigned i=0; i<node_pairs.size(); i++)
	{		
		unsigned transit_node_index = node_pairs[i][0];
		unsigned diff_node_index = node_pairs[i][1];
				
		c_vector<double, 2> transit_location = rCellPopulation.GetNode(transit_node_index)->rGetLocation();
		c_vector<double, 2> diff_location = rCellPopulation.GetNode(diff_node_index)->rGetLocation();
		
		// The force due to the basal lamina acts along the spring connecting the transit and diff nodes, in the 
		// direction D->T
		c_vector<double, 2> curvature_force_direction = transit_location - diff_location;
				
		double distance_between_nodes = norm_2(curvature_force_direction);
		assert(distance_between_nodes > 0);
		assert(!isnan(distance_between_nodes));
		
		curvature_force_direction /= distance_between_nodes;
		
		double curvature = GetCurvatureFromNodePair(rCellPopulation, transit_node_index, diff_node_index);	
		
		double basal_lamina_parameter = GetBasalLaminaParameter();

		c_vector<double, 2> force_due_to_basal_lamina = basal_lamina_parameter*curvature*curvature_force_direction;
				
		// Add the force due to the basal lamina to the forces acting on that transit node

//		std::cout << basal_lamina_parameter*curvature << "\t" << force_due_to_basal_lamina[0] << "\t" << force_due_to_basal_lamina[1] << "\n";

		rForces[transit_node_index] += force_due_to_basal_lamina;
	}
		
	for (typename MeshBasedCellPopulation<DIM>::SpringIterator spring_iterator=(static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))->SpringsBegin();
        spring_iterator!=(static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))->SpringsEnd();
        ++spring_iterator)
    {

        unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
        unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
              
        c_vector<double, DIM> force = CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index, rCellPopulation); 
        
//        std::cout << force[0] << "\t" << force[1] << "\n";

        rForces[nodeB_global_index] -= force;
        rForces[nodeA_global_index] += force;        
    }
			   		
}

template<unsigned DIM>
c_vector<double, DIM> CellCentreWithVariableCellInteraction<DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
																						 unsigned nodeBGlobalIndex,
                                                                                         AbstractCellPopulation<DIM>& rCellPopulation)
{
    // We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex != nodeBGlobalIndex);

    // Get the node locations
    c_vector<double, DIM> node_a_location = rCellPopulation.GetNode(nodeAGlobalIndex)->rGetLocation();
    c_vector<double, DIM> node_b_location = rCellPopulation.GetNode(nodeBGlobalIndex)->rGetLocation();

    // Get the unit vector parallel to the line joining the two nodes
    c_vector<double, DIM> unit_difference;
    
    if (rCellPopulation.IsMeshBasedCellPopulation())
    {
        /*
         * We use the mesh method GetVectorFromAtoB() to compute the direction of the
         * unit vector along the line joining the two nodes, rather than simply subtract
         * their positions, because this method can be overloaded (e.g. to enforce a
         * periodic boundary in Cylindrical2dMesh).
         */
        unit_difference = (static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))->rGetMesh().GetVectorFromAtoB(node_a_location, node_b_location);
    }
    else
    {
        unit_difference = node_b_location - node_a_location;
    }
   
    // Calculate the distance between the two nodes
    double distance_between_nodes = norm_2(unit_difference);
    assert(distance_between_nodes > 0);
    assert(!std::isnan(distance_between_nodes));

    unit_difference /= distance_between_nodes;

    /*
     * If mUseCutoffPoint has been set, then there is zero force between
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

    CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
    CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

    double ageA = p_cell_A->GetAge();
    double ageB = p_cell_B->GetAge();

    assert(!std::isnan(ageA));
    assert(!std::isnan(ageB));

    /*
     * If the cells are both newly divided, then the rest length of the spring
     * connecting them grows linearly with time, until 1 hour after division.
     */
    if (ageA < this->mMeinekeSpringGrowthDuration && ageB < this->mMeinekeSpringGrowthDuration)
    {
        if (rCellPopulation.IsMeshBasedCellPopulation())
        {
            MeshBasedCellPopulation<DIM>* p_static_cast_cell_population = static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation);

            std::pair<CellPtr,CellPtr> cell_pair = p_static_cast_cell_population->CreateCellPair(p_cell_A, p_cell_B);

            if (p_static_cast_cell_population->IsMarkedSpring(cell_pair))
            {
                // Spring rest length increases from a small value to the normal rest length over 1 hour
                double lambda = this->GetMeinekeDivisionRestingSpringLength();
                rest_length = lambda + (1.0 - lambda) * ageA/this->mMeinekeSpringGrowthDuration;
            }
            if (ageA + SimulationTime::Instance()->GetTimeStep() >= this->mMeinekeSpringGrowthDuration)
            {
                // This spring is about to go out of scope
                p_static_cast_cell_population->UnmarkSpring(cell_pair);
            }
        }
        else
        {
            // Spring rest length increases from mDivisionRestingSpringLength to normal rest length, 1.0, over 1 hour
            double lambda = this->GetMeinekeDivisionRestingSpringLength();
            rest_length = lambda + (1.0 - lambda) * ageA/this->mMeinekeSpringGrowthDuration;
        }
    }

    double a_rest_length = rest_length*0.5;
    double b_rest_length = a_rest_length;

    /*
     * If either of the cells has begun apoptosis, then the length of the spring
     * connecting them decreases linearly with time.
     */
    if (p_cell_A->HasApoptosisBegun())
    {
        double time_until_death_a = p_cell_A->GetTimeUntilDeath();
        a_rest_length = a_rest_length * time_until_death_a / p_cell_A->GetApoptosisTime();
    }
    if (p_cell_B->HasApoptosisBegun())
    {
        double time_until_death_b = p_cell_B->GetTimeUntilDeath();
        b_rest_length = b_rest_length * time_until_death_b / p_cell_B->GetApoptosisTime();
    }

    rest_length = a_rest_length + b_rest_length;
    assert(rest_length <= 1.0+1e-12);

    bool is_closer_than_rest_length = (distance_between_nodes - rest_length <= 0);

    // Although in this class the 'spring constant' is a constant parameter, in
    // subclasses it can depend on properties of each of the cells
    double multiplication_factor = VariableSpringConstantMultiplicationFactor(nodeAGlobalIndex, nodeBGlobalIndex, rCellPopulation, is_closer_than_rest_length);
    double spring_stiffness = this->GetMeinekeSpringStiffness();
    double overlap = distance_between_nodes - rest_length;
   	
    // Get the cell types
    CellProliferativeType cellA_type = p_cell_A->GetCellCycleModel()->GetCellProliferativeType();
    CellProliferativeType cellB_type = p_cell_B->GetCellCycleModel()->GetCellProliferativeType();
    
    /* Want to have one-way springs between epithelial and tissue nodes, so that there is only repulsion due to compression
     * of the spring, but no attraction due to extension
     */    
    if ( (cellA_type==TRANSIT && cellB_type==DIFFERENTIATED) || (cellA_type==DIFFERENTIATED && cellB_type==TRANSIT) )
    {
        if (distance_between_nodes > rest_length)
        {
            return zero_vector<double>(DIM); // c_vector<double,DIM>() is not guaranteed to be fresh memory
        }
    }    
    
    if (rCellPopulation.IsMeshBasedCellPopulation())
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

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class CellCentreWithVariableCellInteraction<1>;
template class CellCentreWithVariableCellInteraction<2>;
template class CellCentreWithVariableCellInteraction<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellCentreWithVariableCellInteraction)
