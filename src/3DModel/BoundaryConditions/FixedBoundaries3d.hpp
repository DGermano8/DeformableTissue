#ifndef FIXEDBOUNDARIES3D_HPP_
#define FIXEDBOUNDARIES3D_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>


class FixedBoundaries3d : public AbstractCellPopulationBoundaryCondition<3>
{
private:

    // Parameter used to define the height below which cells should be pinned
    double mMaxHeightForPinnedCells;

    // The width of the original mesh of nodes (i.e. the maximum x value)
    double mCellPopulationWidth;

    // The depth of the original mesh of nodes (i.e. the maximum y value)
    double mCellPopulationDepth;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<3> >(*this);
    }

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param point a point on the boundary plane
     * @param normal the outward-facing unit normal vector to the boundary plane
     */
    FixedBoundaries3d(AbstractCellPopulation<3>* pCellPopulation, double maxHeightForPinnedCells,
    		double cellPopulationWidth, double cellPopulationDepth);


    double GetMaxHeightForPinnedCells();

    double GetCellPopulationWidth();

    double GetCellPopulationDepth();

    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     */
    void ImposeBoundaryCondition();

    /**
     * Overridden VerifyBoundaryCondition() method.
     * Verify the boundary conditions have been applied.
     * This is called after ImposeBoundaryCondition() to ensure the condition is still satisfied.
     *
     * @return whether the boundary conditions are satisfied.
     */
    bool VerifyBoundaryCondition();

    /**
     * Overridden OutputCellPopulationBoundaryConditionParameters() method.
     * Output cell population boundary condition parameters to file.
	 *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(FixedBoundaries3d)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a FixedBoundaries3d.
 */
template<class Archive >
inline void save_construct_data(
    Archive & ar, const FixedBoundaries3d * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<3>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;
    double mMaxHeightForPinnedCells;
    ar << mMaxHeightForPinnedCells;
    double mCellPopulationWidth;
    ar << mCellPopulationWidth;
    double mCellPopulationDepth;
    ar << mCellPopulationDepth;
}

/**
 * De-serialize constructor parameters and initialize a CryptBoundaryCondition3d.
 */
template<class Archive >
inline void load_construct_data(
    Archive & ar, FixedBoundaries3d * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<3>* p_cell_population;
    ar >> p_cell_population;
    double maxHeightForPinnedCells;
    ar >> maxHeightForPinnedCells;
    double cellPopulationWidth;
    ar >> cellPopulationWidth;
    double cellPopulationDepth;
    ar >> cellPopulationWidth;

    // Invoke inplace constructor to initialise instance
    ::new(t)FixedBoundaries3d(p_cell_population, maxHeightForPinnedCells, cellPopulationWidth, cellPopulationDepth);
}
}
} // namespace ...


#endif /* FIXEDBOUNDARIES3D_HPP_ */
