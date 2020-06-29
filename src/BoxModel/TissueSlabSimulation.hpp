#ifndef TISSUESLABSIMULATION_HPP_
#define TISSUESLABSIMULATION_HPP_

#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisOne.hpp"
#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo.hpp"
#include "AbstractMesh.hpp"
#include "OffLatticeSimulation.hpp"
#include "SimpleDataWriter.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "Cylindrical2dMesh.hpp"
#include "TrianglesMeshWriter.hpp"

// Needs to be included last
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>


// A 2D box model simulation object 

class TissueSlabSimulation : public OffLatticeSimulation<2>
{
    // Allow tests to access private members, in order to test computation of
    // private functions eg. DoCellBirth
    friend class TestTissueSlabSimulation;

private:

    /** Define a stopping event which says stop if you have mismatched boundary elements (mMismatchedBoundaryElements == true)*/
//    bool StoppingEventHasOccurred();

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & boost::serialization::base_object<OffLatticeSimulation<2> >(*this);
        archive & mBasalLaminaParameter;
    }

    // The file that the values of beta catenin is written out to
    out_stream mBetaCatResultsFile;

    // Helper member that is a static cast of the tissue
    MeshBasedCellPopulationWithGhostNodes<2>* mpStaticCastCellPopulation;

    // The parameter for the basal lamina force - multiplies the local curvature
    double mBasalLaminaParameter;

    // The output file directory for the simulation data - number of popped up cells, 
    // average distance between cells etc.
    std::string mDataOutputFile;

    /**
     * Calculates the new locations of a dividing cell's cell centres.
     * Moves the dividing node a bit and returns co-ordinates for the new node.
     * It does this by picking a random direction (0->2PI) and placing the parent
     * and daughter in opposing directions on this axis.
     *
     * @param pParentCell pointer to the parent cell
     *
     * @return daughter_coords the coordinates for the daughter cell.
     */
    c_vector<double, 2> CalculateCellDivisionVector(CellPtr pParentCell);

    // Sloughing off the top row of transit cells
    void SloughTransitCells();

    /**
     * Overridden WriteVisualizerSetupFile() method.
     * 
     * Writes out special information about the mesh to the visualizer.
     */
    void WriteVisualizerSetupFile();

    /**
     * Get the data output directory of the simulation.
     */
     std::string GetDataOutputFile();

    /**
     * Overridden AfterSolve() method.
     * 
     * Closes beta catenin results file if required, then calls 
     * the base class method.
     */
    void AfterSolve();

public :

    /**
     *  Constructor.
     *
     *  @param rCellPopulation A tissue facade class (contains a mesh and cells)
     *  @param forceCollection The mechanics to use in the simulation
     *  @param deleteCellPopulationAndForceCollection Whether to delete the tissue and force collection on destruction to free up memory
     *  @param initialiseCells whether to initialise cells (set to false when loading from an archive)
     *  @param tissueWidth Width of mesh of real nodes  
     *  @param sloughTransitCells whether to slough off those transit cells that move beyond the boundaries
     */
    TissueSlabSimulation(AbstractCellPopulation<2>& rCellPopulation,
                      bool deleteCellPopulationAndForceCollection=false,
                      bool initialiseCells=true);


    ~TissueSlabSimulation();

    /**
     * Method to calculate the average horizontal gap between transit cells - 
     * only those transit cells in the monolayer however, not considering transit cells
     * that have popped up
     */   
    c_vector<double,2> CalculateMeanGapBetweenTransitCells();    
    
    /** Method to calculate the maximum vertical gap between any pair of transit cells in the 
     * monolayer. This ignores any transit cells that have popped up.
     */    
    double CalculateMaxVerticalGapBetweenTransitCells();

    /**
     * Method to return true/false if a particular transit cell has lost connections
     * to the differentiated cells below, i.e. that it has popped up.
     */
    bool HasCellPoppedUp(unsigned nodeIndex);
//
//    /**
//     * Method to return the number of cells removed by anoikis ([0]) and sloughing ([1])
//     */
//    c_vector<unsigned,2> GetNumberCellsRemoved();

    /**
     * Given a node, find a set containing the indices of its neighbouring nodes.
     * 
     * @param nodeIndex global index of the node
     * @return its neighbouring node indices
     */
    std::set<unsigned> GetNeighbouringNodeIndices(unsigned nodeIndex);  

    /**
     * Output the results of the simulation (the number of transit cells 
     * which popped out during the simulation, and the mean horizontal and 
     * vertical distance between transit cells at the end of the simulation).
     * 
     * To be called after Solve() in tests.
     * 
     * @param basalLaminaParameter the basal lamina parameter 
     *                             (needs to be passed in as it is a member of the force law, not the tissue simulation)
     */
    void WriteBasalLaminaResultsToDirectory(double basalLaminaParameter);	// When directory is specified	
    
    c_vector<double,6> OutputVectorOfBasalLaminaResults(double basalLaminaParameter);
};


// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(TissueSlabSimulation)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a TissueSlabSimulation
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const TissueSlabSimulation * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<2> * p_tissue = &(t->rGetCellPopulation());
    ar & p_tissue;
}

/**
 * De-serialize constructor parameters and initialise Box Model.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, TissueSlabSimulation * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<2>* p_tissue;
    ar >> p_tissue;
    std::vector<AbstractForce<2>*> force_collection;
    ar >> force_collection;

    // Invoke inplace constructor to initialise instance
    ::new(t)TissueSlabSimulation(*p_tissue, true, false);
}
}
} // namespace



#endif /*BOXMODELSIMULATION_HPP_*/
