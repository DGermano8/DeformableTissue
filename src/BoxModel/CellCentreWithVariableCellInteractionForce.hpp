#ifndef CELLCENTREWITHVARIABLECELLINTERACTIONFORCE_HPP_
#define CELLCENTREWITHVARIABLECELLINTERACTIONFORCE_HPP_

#include "GeneralisedLinearSpringForce.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisOne.hpp"
#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo.hpp"
#include "Debug.hpp"
#include <cmath>
#include <list>
#include <fstream>

#include "ChasteSerialization.hpp"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>

/**
 * A subclass of MeinekeInteractionForce with variable spring constants.
 */
template<unsigned DIM>
class CellCentreWithVariableCellInteractionForce : public GeneralisedLinearSpringForce<DIM>
{
    friend class TestBoxModelSimulation;
    
private :

    /** Needed for serialization. */
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<GeneralisedLinearSpringForce<DIM> >(*this);
        archive & mUseCellTypeDependentSprings;
    	archive & mTransitTransitMultiplier;
    	archive & mDifferentiatedDifferentiatedMultiplier;
    	archive & mTransitDifferentiatedMultiplier;
        archive & mUseEdgeBasedSpringConstant;
        archive & mUseMutantSprings;
        archive & mMutantMutantMultiplier;
        archive & mNormalMutantMultiplier;
        archive & mUseBCatSprings;
        archive & mBetaCatSpringScaler;
        archive & mUseApoptoticSprings;
        archive & mApoptoticSpringTensionStiffness;
        archive & mApoptoticSpringCompressionStiffness;
        archive & mBasalLaminaParameter;
        archive & mUseOneWaySprings;
    }

protected :
	
	// Whether to use different spring stiffness dependent on the types of neighbouring cells
	bool mUseCellTypeDependentSprings;
    
	// Multiplier for spring stiffness if cells are both transit type
	double mTransitTransitMultiplier;

	// Multiplier for spring stiffness if cells are both differentiated type
	double mDifferentiatedDifferentiatedMultiplier;
	
	// Multiplier for spring stiffness if cells are of different types
	double mTransitDifferentiatedMultiplier;
		
    /** Whether to use spring constant proportional to cell-cell contact length/area (defaults to false) */
    bool mUseEdgeBasedSpringConstant;

    /** Whether to use different stiffnesses depending on whether either cell is a mutant */
    bool mUseMutantSprings;

    /** Multiplier for spring stiffness if mutant */
    double mMutantMutantMultiplier;

    /** Multiplier for spring stiffness if mutant */
    double mNormalMutantMultiplier;

    /** Use springs which are dependent on beta-catenin levels */
    bool mUseBCatSprings;
    
    /** Scaling factor for beta catenin to spring strength. */
    double mBetaCatSpringScaler;

    /** Use springs which are dependent on whether cells are apoptotic */
    bool mUseApoptoticSprings;
    
    /** Non-dimensionalized 'stiffness' of a apoptotic cell under tension. */
    double mApoptoticSpringTensionStiffness;

    /** Non-dimensionalized 'stiffness' of a apoptotic cell under compression. */
    double mApoptoticSpringCompressionStiffness;
    
    /** Parameter that multiplies the curvature to give the basal lamina force */
    double mBasalLaminaParameter;
    
    /** Use one way springs between epithelial and tissue nodes, that only act under repulsion */
    bool mUseOneWaySprings;

public :

    /**
     * Constructor.
     */
    CellCentreWithVariableCellInteractionForce();
    
    /**
     * Destructor.
     */
    ~CellCentreWithVariableCellInteractionForce();
    
    // Set whether to use spring constants dependent on the types of neighbouring cells:
    // Transit-transit, transit-differentiated, differentiated-differentiated
    
    // @param useCellTypeDependentSprings  whether to use mutant springs
    // @param transitTransitMultiplier  the multiplier for springs connecting two transit cells
    // @param transitDifferentiatedMultiplier  the multiplier for springs connecting a transit cell with a differentiated cell
    // @param differentiatedDifferentiatedMultiplier the multiplier for springs connecting two differentiated cells
    
    void SetCellTypeDependentSprings(bool useCellTypeDependentSprings, double transitTransitMultiplier=1, 
    		double differentiatedDifferentiatedMultiplier=1,
    		double transitDifferentiatedMultiplier=1);
    		
    /**
     * Set whether to use an edge-based spring constant.
     * 
     * @param useEdgeBasedSpringConstant
     */
    void SetEdgeBasedSpringConstant(bool useEdgeBasedSpringConstant);

    /**
     * Use different spring strengths depending on two cells:
     * Normal-normal, Normal-mutant, mutant-mutant
     * 
     * @param useMutantSprings  whether to use mutant springs
     * @param mutantMutantMultiplier  the multiplier for springs connecting two mutant cells
     * @param normalMutantMultiplier  the multiplier for springs connecting a mutant cell with a normal cell
     */
    void SetMutantSprings(bool useMutantSprings, double mutantMutantMultiplier=2, double normalMutantMultiplier=1.5);

    /**
     * Use the amount of beta-catenin on an edge to find spring constant.
     * 
     * @param useBCatSprings whether to use beta-catenin-dependent spring stiffness
     */
    void SetBetaCateninSprings(bool useBCatSprings);

    /**
     * Set spring stiffness to be dependent on whether cells are apoptotic
     * 
     * @param useApoptoticSprings whether to have apoptosis-dependent spring stiffness
     */
    void SetApoptoticSprings(bool useApoptoticSprings);  
        
    /* Set method for Basal Lamina Parameter
     */
    void SetBasalLaminaParameter(double basalLaminaParameter);
    
    /* Get method for Basal Lamina Parameter
     */
    double GetBasalLaminaParameter();
    
    /* Set method for using one way springs between epithelial and tissue nodes
     */
    void SetUseOneWaySprings(bool useOneWaySprings = false);

    /* Set method for Transit-Transit spring strength multiplier
     */
    void SetTransitTransitMultiplier(double transitTransitMultiplier);
    
    /* Get method for Transit-Transit spring strength multiplier
     */
    double GetTransitTransitMultiplier();  
    
    /**
     * @return mBetaCatSpringScaler
     */
    double GetBetaCatSpringScaler();

    /**
     * Set mBetaCatSpringScaler.
     * 
     * @param betaCatSpringScaler the new value of mBetaCatSpringScaler
     */
    void SetBetaCatSpringScaler(double betaCatSpringScaler);

    /**
     * @return mApoptoticSpringTensionStiffness
     */
    double GetApoptoticSpringTensionStiffness();

    /**
     * Set mApoptoticSpringTensionStiffness.
     * 
     * @param apoptoticSpringTensionStiffness the new value of mApoptoticSpringTensionStiffness
     */
    void SetApoptoticSpringTensionStiffness(double apoptoticSpringTensionStiffness);

    /**
     * @return mApoptoticSpringCompressionStiffness
     */
    double GetApoptoticSpringCompressionStiffness();

    /**
     * Set mApoptoticSpringCompressionStiffness.
     * 
     * @param apoptoticSpringCompressionStiffness the new value of mApoptoticSpringCompressionStiffness
     */
    void SetApoptoticSpringCompressionStiffness(double apoptoticSpringCompressionStiffness);
    
    /**
     * Outputs force Parameters to file
	 *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
    
    /* Removing duplicated entries of a vector 
     */    
    void RemoveDuplicates1D(std::vector<unsigned>& vec);
        
    /* Returns a boolean for whether the element contains ghost nodes
     */
    
    bool DoesElementContainGhostNodes(AbstractCellPopulation<DIM>& rCellPopulation, unsigned elementIndex);

    /* Finding the connected pairs of transit-differentiated nodes
     */
    std::vector<c_vector<unsigned, DIM> > GetTransitDifferentiatedPairs(AbstractCellPopulation<DIM>& rCellPopulation);
        
    /* Takes a transit node index and a differentiated node index and returns the curvature of 
     * the curve passing through the midpoints of the transit-differentiated springs of the 
     * common elements
     */
    
    double GetCurvatureFromNodePair(AbstractCellPopulation<DIM>& rCellPopulation, unsigned transitNode, unsigned diffNode);
    
    /* Finding an ordered vector of curvature, where the curvature is found from the midpoints
     * of the springs connecting transit and differentiated cells. Each curvature value is determined
     * from a set of three midpoints, moving in order from left to right.
     */
   
    void GetOrderedCurvatureVector(AbstractCellPopulation<DIM>& rCellPopulation);
    
    /* Finding the curvature between three points - taken to be the midpoints of three springs
     * connecting transit cells to differentiated cells
     */
    
    double FindCurvature(c_vector<double, DIM> leftMidpoint, 
			  			 c_vector<double, DIM> centreMidpoint, 
			  			 c_vector<double, DIM> rightMidpoint);
    
    /* Finding the number of elements that a node belongs to, which contain only real nodes
     * and not ghost nodes */
    
    unsigned GetNumContainingElementsWithoutGhostNodes(AbstractCellPopulation<DIM>& rCellPopulation, unsigned nodeIndex); 
    
    
    /**
     * Return a multiplication factor for the spring constant, which 
     * may depend on whether the given pair of neighbouring cells are 
     * e.g. undergoing apoptosis, have mutations, or experience variable 
     * levels of beta catenin.
     * 
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeAGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the tissue
     * @param isCloserThanRestLength whether the neighbouring nodes lie closer than the rest length of their connecting spring
     * 
     * @return the multiplication factor.
     */     
    double VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex, 
                                                      unsigned nodeBGlobalIndex, 
                                                      AbstractCellPopulation<DIM>& rCellPopulation, 
                                                      bool isCloserThanRestLength);
    
    /**
     * Overridden AddForceContribution method.
     * 
     * @param rForces reference to vector of forces on nodes
     * @param rCellPopulation reference to the tissue
     */
    void AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                              AbstractCellPopulation<DIM>& rCellPopulation);
    
    c_vector<double, DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractCellPopulation<DIM>& rCellPopulation);
 
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellCentreWithVariableCellInteractionForce)

#endif /*CELLCENTREWITHVARIABLECELLINTERACTIONFORCE_HPP_*/
