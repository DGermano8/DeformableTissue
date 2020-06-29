#ifndef STROMALCELLPROPERTY_HPP_
#define STROMALCELLPROPERTY_HPP_

#include "AbstractCellProperty.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/shared_ptr.hpp>

/**
 * Stromal cell property class.
 *
 * Each Cell owns a CellPropertyCollection, which may include a shared pointer
 * to an object of this type. When a Cell that is stromal divides, the daughter
 * cells are both stromal. But it shouldn't divide under normal circumstances...
 *
 * The StromalCellProperty object keeps track of the number of cells that have this property, as well
 * as what colour should be used by the visualizer to display cells with the property.
 */
class StromalCellProperty : public AbstractCellProperty
{
private:

    /**
     * Colour for use by visualizer.
     */
    unsigned mColour;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<AbstractCellProperty>(*this);
        archive & mColour;
    }

public:

    /**
     * Constructor.
     *
     * @param colour  what colour cells with this label should be in the visualizer (defaults to 5)
     */
    StromalCellProperty(unsigned colour=5);

    /**
     * Destructor.
     */
    virtual ~StromalCellProperty();

    /**
     * Get #mColour.
     */
    unsigned GetColour() const;
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(StromalCellProperty)

#endif /* STROMALCELLPROPERTY_HPP_ */
