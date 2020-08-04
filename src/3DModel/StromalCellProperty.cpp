#include "StromalCellProperty.hpp"

StromalCellProperty::StromalCellProperty(unsigned colour)
    : AbstractCellProperty(),
      mColour(2)
{
}

StromalCellProperty::~StromalCellProperty()
{
}

unsigned StromalCellProperty::GetColour() const
{
    return mColour;
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(StromalCellProperty)


