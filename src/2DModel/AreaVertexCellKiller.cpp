/*

Copyright (c) 2005-2018, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "AreaVertexCellKiller.hpp"
#include "Debug.hpp"
template<unsigned DIM>
AreaVertexCellKiller<DIM>::AreaVertexCellKiller(AbstractCellPopulation<DIM>* pCellPopulation, double probabilityOfDeathInAnHour)
        : AbstractCellKiller<DIM>(pCellPopulation),
          mProbabilityOfDeathInAnHour(probabilityOfDeathInAnHour)
{
    if ((mProbabilityOfDeathInAnHour<0) || (mProbabilityOfDeathInAnHour>1))
    {
        EXCEPTION("Probability of death must be between zero and one");
    }
}

template<unsigned DIM>
double AreaVertexCellKiller<DIM>::GetDeathProbabilityInAnHour() const
{
    return mProbabilityOfDeathInAnHour;
}

template<unsigned DIM>
void AreaVertexCellKiller<DIM>::CheckAndLabelSingleCellForApoptosis(CellPtr pCell)
{
    // std::cout<< "yo1\n";
    double cell_volume = this->mpCellPopulation->GetVolumeOfCell(pCell);
    
    // std::cout<< cell_volume << "\n";
    double volume_threshold = mProbabilityOfDeathInAnHour;

    if ( cell_volume < volume_threshold)
    {
        pCell->StartApoptosis();
        // std::cout<< "yo\n";
    }
}

template<unsigned DIM>
void AreaVertexCellKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{

    double domain_tollerance = 1.5;

    double mCellPopulationWidth = 20.0;
	double mCellPopulationDepth = 24.0*sqrt(0.75);
    // TRACE("IN");
    // double cell_volume_tot = 0.0;
    // double counter = 0.0;
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        // double cell_volume = this->mpCellPopulation->GetVolumeOfCell(*cell_iter);
        // cell_volume_tot = cell_volume_tot + cell_volume;
        // counter = counter + 1;
        if( !((*cell_iter)->IsDead()) &&  !((*cell_iter)->HasApoptosisBegun()) )
        {
            c_vector<double, DIM> rLocation = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
            double x_location = rLocation[0];
            double y_location = rLocation[1];
            
            // Make sure cell is well inside tissue and not just an edge
            if ( x_location <= domain_tollerance || x_location >= mCellPopulationWidth - domain_tollerance ||
            y_location <= domain_tollerance || y_location >= mCellPopulationDepth - domain_tollerance   )
            {			
                CheckAndLabelSingleCellForApoptosis(*cell_iter);
            }
        }
    }

    // cell_volume_tot = cell_volume_tot/counter;
    // std::cout<< " Average Cell Volume = " << cell_volume_tot << "\n";
	// for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
    //      cell_iter != this->mpCellPopulation->End();
    //      ++cell_iter)
    // {

    //     CellPtr p_cell = this->mpCellPopulation->GetCellUsingLocationIndex(cell_iter);
    //     if(!(p_cell->IsDead()))
    //     {
    //         if(p_cell->HasApoptosisBegun())
    //         {
    //             if (p_cell->GetTimeUntilDeath() <= 2*SimulationTime::Instance()->GetTimeStep())
    //             {
    //                 p_cell->Kill();
    //             }
    //         }
    //     }
    // }
    // TRACE("Out");
}

template<unsigned DIM>
void AreaVertexCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ProbabilityOfDeathInAnHour>" << mProbabilityOfDeathInAnHour << "</ProbabilityOfDeathInAnHour>\n";

    // Call method on direct parent class
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

// Explicit instantiation
template class AreaVertexCellKiller<1>;
template class AreaVertexCellKiller<2>;
template class AreaVertexCellKiller<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AreaVertexCellKiller)
