#ifndef SLOUGHINGANDANOIKISCELLKILLER_HPP_
#define SLOUGHINGANDANOIKISCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"
#include "MeshBasedCellPopulation.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 *  A cell killer that kills transit cells if they are outside the boundaries.
 *  
 *  The crypt is assumed to start at x=0 and y=0. By default only cells
 *  are sloughed if x < 0 or x > width. To slough the sides call the constructor
 *  with the appropriate parameter.
 */
class SloughingAndAnoikisCellKiller : public AbstractCellKiller<2>
{
private:

    /** Whether cells should be sloughed from the sides of the crypt. */
    bool mSloughSides;
    
    // Width of tissue box
    double mCellPopulationWidth;
    
    unsigned mCellsRemovedByAnoikis;

    unsigned mCellsRemovedBySloughing;

    /** Needed for serialization. */   
    friend class boost::serialization::access;
    
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<2> >(*this);
        
        archive & mSloughSides;
        archive & mCellPopulationWidth;
        archive & mCellsRemovedByAnoikis;
        archive & mCellsRemovedBySloughing;
    }

public:

    /**
     * Default constructor.
     * 
     * @param pCrypt pointer to a tissue
     * @param sloughSides whether to slough cells at the side of the crypt
     */
    SloughingAndAnoikisCellKiller(AbstractCellPopulation<2>* pCrypt, double cellPopulationWidth);

    // Destructor
    ~SloughingAndAnoikisCellKiller();
    
//    void SetCellPopulationWidth(double cellPopulationWidth) const;

    double GetCellPopulationWidth() const;
    
    std::set<unsigned> GetNeighbouringNodeIndices(unsigned nodeIndex);
    
    bool HasCellPoppedUp(unsigned nodeIndex);

    /** Method that returns a matrix of rows (nodeIndex / KillByAnoikis / KillBySloughing),
     * where the last two arguments are 1 if true, 0 if false. This is called by CheckAndLabelCellsForApoptosisOrDeath()
     * and used to determine which cells to kill
     */
    std::vector<c_vector<unsigned,3> > RemoveByAnoikisOrSloughing();

    /**
     *  Loops over cells and kills cells by anoikis or by sloughing
     */
    void CheckAndLabelCellsForApoptosisOrDeath();
    
    /* After each event of cell killing in CheckAndLabelCellsForApoptosisOrDeath(), the information of whether to kill each cell
     * or not is passed to this method which then increments the member variables corresponding to the total number of cells
     * killed by anoikis or sloughing
     */
    void SetNumberCellsRemoved(std::vector<c_vector<unsigned,3> > cellsRemoved);

    /* Returns the total number of cells removed by anoikis ([0]) and by sloughing ([1])
     *
     */
    c_vector<unsigned,2> GetNumberCellsRemoved();

    /**
     * Outputs cell killer parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellKillerParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(SloughingAndAnoikisCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a SloughingAndAnoikisCellKiller.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const SloughingAndAnoikisCellKiller * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<2>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;
    double width = t->GetCellPopulationWidth();
    ar << width;
}

/**
 * De-serialize constructor parameters and initialise Crypt.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, SloughingAndAnoikisCellKiller * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<2>* p_cell_population;
    ar >> p_cell_population;
    double width;
    ar >> width;
      
    // Invoke inplace constructor to initialise instance
    ::new(t)SloughingAndAnoikisCellKiller(p_cell_population, width);
}
}
} // namespace ...

#endif /*SLOUGHINGANDANOIKISCELLKILLER_HPP_*/
