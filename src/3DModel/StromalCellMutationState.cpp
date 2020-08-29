#include "StromalCellMutationState.hpp"

StromalCellMutationState::StromalCellMutationState()
    : AbstractCellMutationState(1)
{}


#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(StromalCellMutationState)
