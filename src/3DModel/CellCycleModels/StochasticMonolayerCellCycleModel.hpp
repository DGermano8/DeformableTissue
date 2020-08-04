
#ifndef STOCHASTICMONOLAYERCELLCYCLEMODEL_HPP_
#define STOCHASTICMONOLAYERCELLCYCLEMODEL_HPP_

#include <cassert>
#include <iostream>

#include "AbstractSimpleGenerationBasedCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "Exception.hpp"

class StochasticMonolayerCellCycleModel : public AbstractSimpleGenerationBasedCellCycleModel
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell cycle model and random number generator, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleGenerationBasedCellCycleModel>(*this);

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        archive & *p_gen;
        archive & p_gen;
    }

    /**
     * Stochastically set the G1 duration.  Called on cell creation at
     * the start of a simulation, and for both parent and daughter
     * cells at cell division.
     */
    void SetG1Duration();

public:

    /**
     * Constructor - just a default, mBirthTime is now set in the AbstractCellCycleModel class.
     * mG1Duration is set very high, it is set for the individual cells when InitialiseDaughterCell is called
     */
    StochasticMonolayerCellCycleModel()
    {}

    /**
     * Overridden builder method to create new copies of
     * this cell cycle model.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(StochasticMonolayerCellCycleModel)

#endif /*STOCHASTICMONOLAYERCELLCYCLEMODEL_HPP_*/
