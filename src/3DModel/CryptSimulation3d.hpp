#ifndef CRYPTSIMULATION3D_HPP_
#define CRYPTSIMULATION3D_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>

#include "AbstractMesh.hpp"
#include "OffLatticeSimulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "AbstractVanLeeuwen2009WntSwatCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StromalCellMutationState.hpp"
#include "Cylindrical2dMesh.hpp"
#include "CryptModelInteractionForce.hpp"

class CryptSimulation3d : public OffLatticeSimulation<3>
{
    // Allow tests to access private members, in order to test computation of
    // private functions
    friend class Test3dBoxModel;

protected:
	
    /** Define a stopping event which says stop if you have mismatched boundary elements (mMismatchedBoundaryElements == true)*/
    bool StoppingEventHasOccurred();

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the simulation and member variable.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & boost::serialization::base_object<OffLatticeSimulation<3> >(*this);
    }

    // The file that the number of cells at each timestep are written out to
    out_stream mNumCellsResultsFile;

    // Helper member that is a static cast of the tissue
    MeshBasedCellPopulationWithGhostNodes<3>* mpStaticCastCellPopulation;

    // The output file directory for the simulation data - number of popped up cells,
    // average distance between cells etc.
    std::string mDataOutputFile;
    
    /* The output file for epithelial cell coordinate data
     */
    out_stream mEpithelialCoordinateDataResultsFile;
    
    /** Stopping condition */
    bool mFoundNonExistentNode;

    /**
     * Calculates the new locations of a dividing cell's cell centres.
     * There are two choices here: one returns a division vector that is parallel to
     * the vector connecting the neighbouring epithelial cell centres, and the second
     * chooses a random direction, as below.
     *
     * Moves the dividing node a bit and returns co-ordinates for the new node.
     * It does this by picking a random direction (0->2PI) and placing the parent
     * and daughter in opposing directions on this axis.
     *
     * @param pParentCell pointer to the parent cell
     *
     * @return daughter_coords the coordinates for the daughter cell.
     */
    c_vector<double, 3> CalculateCellDivisionVector(CellPtr pParentCell);
    
    /* A method to return a vector of epithelial node pairs, which is used to calulate the wiggliness 
     * of the epithelial layer
     */
    std::vector<c_vector<unsigned, 2> > GetNeighbouringEpithelialPairs(AbstractCellPopulation<3>& rCellPopulation);
    
    /**
     * Overridden WriteVisualizerSetupFile() method.
     *
     * Writes out special information about the mesh to the visualizer.
     */
    void WriteVisualizerSetupFile();

    /** 
     * Use an output file handler to create an epithelial cell coordinate results file
     */
    void SetupWriteEpithelialCoordinateData();
    
    /**
     * Write epithelial cell data to file
     * 
     * @param time the current time
     */
    void WriteEpithelialCoordinateData(double time);
    
    /**
     * Get the data output directory of the simulation.
     */
     std::string GetDataOutputFile();

     /**
      * Overridden SetupSolve() method.
      *
      * Write results of interest to file if required.
      */
     void SetupSolve();

     /**
      * Overridden PostSolve() method.
      *
      * Write results of interest to file if required.
      */
     void PostSolve();


    /**
     * Overridden AfterSolve() method.
     * 
     */
    void AfterSolve();

public :

    /**
     *  Constructor.
     *
     *  @param rCellPopulation A tissue facade class (contains a mesh and cells)
     *  @param deleteCellPopulationAndForceCollection Whether to delete the tissue and force collection on destruction to free up memory
     *  @param initialiseCells whether to initialise cells (set to false when loading from an archive)
     *  @param tissueWidth Width of mesh of real nodes
     *  @param sloughTransitCells whether to slough off those transit cells that move beyond the boundaries
     */
	CryptSimulation3d(AbstractCellPopulation<3>& rCellPopulation,
                      bool deleteCellPopulationAndForceCollection=false,
                      bool initialiseCells=true);
	
    /* Removing duplicated entries of a vector
     */
    void RemoveDuplicates1D(std::vector<unsigned>& rVectorWithDuplicates);

    /**
     * Given a node, find a set containing the indices of its neighbouring nodes.
     *
     * @param nodeIndex global index of the node
     * @return its neighbouring node indices
     */
    std::set<unsigned> GetNeighbouringNodeIndices(unsigned nodeIndex);
    
    /** SJ: For stopping condition */
    bool GetInstanceOfNonExistentNode();


};


// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CryptSimulation3d)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a BoxModelSimulation
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const CryptSimulation3d * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<3> * p_tissue = &(t->rGetCellPopulation());
    ar & p_tissue;
}

/**
 * De-serialize constructor parameters and initialise a CryptSimulation3d
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, CryptSimulation3d * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<3>* p_tissue;
    ar >> p_tissue;

    // Invoke inplace constructor to initialise instance
    ::new(t)CryptSimulation3d(*p_tissue, true, false);
}
}
} // namespace

#endif /* CRYPTSIMULATION3D_HPP_ */
