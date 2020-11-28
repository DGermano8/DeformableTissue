#ifndef RANDOMRADIALCELLKILLER_HPP_
#define RANDOMRADIALCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "OutputFileHandler.hpp"
#include "StromalCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApoptoticCellProperty.hpp"
#include "SimpleWntCellCycleModel.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A cell killer that randomly kills cells in the outer rim of
 * the epithelial monolayer based on the user set probability.
 *
 * The probability passed into the constructor will be the probability
 * of any cell dying whenever CheckAndLabelCellsForApoptosis() is called.
 *
 * Note this does take into account timesteps - the input probability is the
 * probability that in an hour's worth of trying, the cell killer will have
 * successfully killed a given cell. In the method CheckAndLabelSingleCellForApoptosis()
 * this probability is used to calculate the probability that the cell is killed
 * at a given time step.
 */
class RandomRadialCellKiller3d : public AbstractCellKiller<3>
{
private:

    /**
     * Probability that in an hour's worth of trying, this cell killer
     * will have successfully killed a given cell.
     */
    double mProbabilityOfDeathInAnHour;

    /* The radius of the circle OUTSIDE of which cells are randomly killed */
    double mDeathRadius;

    /* x-coordinate of the centre of the tissue slab */
    double mCentroidXCoordinate;

    /* y-coordinate of the centre of the tissue slab */
    double mCentroidYCoordinate;

    // The output file directory for the simulation data that corresponds to the number of cells
    // killed by anoikis
    out_stream mRandomDeathOutputFile;

    std::string mOutputDirectory;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * Serialization of singleton objects must be done with care.
     * Before the object is serialized via a pointer, it *MUST* be
     * serialized directly, or an assertion will trip when a second
     * instance of the class is created on de-serialization.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<3> >(*this);

        archive & mProbabilityOfDeathInAnHour;
        archive & mDeathRadius;
        archive & mOutputDirectory;
    }

public:

    /**
     * Default constructor.
     *
     * @param pCellPopulation pointer to a tissue
     * @param sloughOrifice whether to slough compressed cells at crypt orifice
     */
    RandomRadialCellKiller3d(AbstractCellPopulation<3>* pCellPopulation, double probabilityOfDeathInAnHour = 0.5,
    		double deathRadius = 0.0, double centroidXCoordinate = 0.0, double centroidYCoordinate = 0.0);

	// Destructor
	~RandomRadialCellKiller3d();


    /*@ return mProbabilityOfDeathInAnHour */
    double GetDeathProbabilityInAnHour() const;

    /* Get the radius of the circle outside of which random cell death occurs */
    double GetDeathRadius() const;

    void SetOutputDirectory(std::string outputDirectory);

    std::string GetOutputDirectory();

    std::vector<c_vector<unsigned,2> > RemoveByRandomSelection();

    /**
     * Overridden method to test a given cell for apoptosis.
     *
     * @param pCell the cell to test for apoptosis
     */
    void CheckAndLabelSingleCellForApoptosis(CellPtr pCell);

    /**
     *  Loops over and kills cells by anoikis
     */
    void CheckAndLabelCellsForApoptosisOrDeath();

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
CHASTE_CLASS_EXPORT(RandomRadialCellKiller3d)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a RandomRadialCellKiller3d.
 */
template<class Archive >
inline void save_construct_data(
    Archive & ar, const RandomRadialCellKiller3d * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<3>* const p_tissue = t->GetCellPopulation();
    ar << p_tissue;
}

/**
 * De-serialize constructor parameters and initialise Crypt.
 */
template<class Archive >
inline void load_construct_data(
    Archive & ar, RandomRadialCellKiller3d * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<3>* p_tissue;
    ar >> p_tissue;

    // Invoke inplace constructor to initialise instance
    ::new(t)RandomRadialCellKiller3d(p_tissue);
}
}
}

#endif /* RANDOMRADIALCELLKILLER_HPP_ */
