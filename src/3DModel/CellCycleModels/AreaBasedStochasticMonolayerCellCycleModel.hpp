#ifndef AREABASEDSTOCHASTICMONOLAYERCELLCYCLEMODEL_HPP_
#define AREABASEDSTOCHASTICMONOLAYERCELLCYCLEMODEL_HPP_

#include "StochasticMonolayerCellCycleModel.hpp"
#include "CellwiseData.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellLabel.hpp"
#include "SimulationTime.hpp"

/**
 * Simple stress-based cell cycle model.
 *
 * A simple stress-dependent cell cycle model that inherits from
 * StochasticMonolayerCellCycleModel. The duration of G1 phase depends
 * on the local stress. This model allows for quiescence imposed
 * by transient periods of high stress, followed by relaxation.
 */
class AreaBasedStochasticMonolayerCellCycleModel : public StochasticMonolayerCellCycleModel
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<StochasticMonolayerCellCycleModel>(*this);
        archive & mQuiescentAreaFraction;
    }

    /**
     * \todo Doxygen
     */
    double mQuiescentAreaFraction;

public:

    /**
     * Constructor.
     */
    AreaBasedStochasticMonolayerCellCycleModel();

    /**
     * Overridden UpdateCellCyclePhase() method.
     */
    void UpdateCellCyclePhase();

    /**
     * Overridden builder method to create new instances of
     * the cell cycle model.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Get the spatial dimension.
     *
     * @return mDimension
     */
    unsigned GetDimension();

    /**
     * Set #mQuiescentAreaFraction.
     */
    void SetQuiescentAreaFraction(double quiescentAreaFraction);

//    virtual bool ReadyToDivide();
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(AreaBasedStochasticMonolayerCellCycleModel)

namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a AreaBasedStochasticMonolayerCellCycleModel instance.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const AreaBasedStochasticMonolayerCellCycleModel * t, const unsigned int file_version)
{
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a AreaBasedStochasticMonolayerCellCycleModel instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, AreaBasedStochasticMonolayerCellCycleModel * t, const unsigned int file_version)
{
    /**
     * Invoke inplace constructor to initialise an instance of AreaBasedStochasticMonolayerCellCycleModel.
     * It doesn't actually matter what values we pass to our standard constructor,
     * provided they are valid parameter values, since the state loaded later
     * from the archive will overwrite their effect in this case.
     */

    ::new(t)AreaBasedStochasticMonolayerCellCycleModel;
}
}
} // namespace ...

#endif /*AREABASEDSTOCHASTICMONOLAYERCELLCYCLEMODEL_HPP_*/
