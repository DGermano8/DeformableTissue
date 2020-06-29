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

#include "CellPopulationNeighbourWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "ApoptoticCellProperty.hpp"
//#include "MeshBasedCellPopulationWithGhostNodes.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellPopulationNeighbourWriter<ELEMENT_DIM, SPACE_DIM>::CellPopulationNeighbourWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("cellneighbours.dat")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationNeighbourWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    assert(SPACE_DIM==2 || SPACE_DIM==3);

    //int cell_number = pCellPopulation->GetNumRealCells();
    //*this->mpOutStream << cell_number;

    
    std::set<unsigned> neighbouring_node_indices; 
    int numnNodeNeighbours = 0;
    
    for(int i=0; i<pCellPopulation->GetNumNodes(); i++)
    {
        neighbouring_node_indices = pCellPopulation->GetNeighbouringNodeIndices(i);

        *this->mpOutStream << neighbouring_node_indices.size()  << " , "; //Use this to see just how many neighbours there are
        
        /*
        if(!pCellPopulation->IsGhostNode(i))
        {
            
            for (auto elem : neighbouring_node_indices)
            {
                if(!pCellPopulation->IsGhostNode(elem))
                {
                    *this->mpOutStream << elem+1  << " "; //Plus 1 since arrays start at 0 here and using nnz in matlab for post
                    
                }
            }

            *this->mpOutStream  << "  ,  ";

        } 
        */
        
    }
    /*
    for(int i=0; i<pCellPopulation->GetNumNodes(); i++)
    {
        neighbouring_node_indices = pCellPopulation->GetNeighbouringNodeIndices(i);

        for (auto elem : neighbouring_node_indices)
        {
            if(!pCellPopulation->IsGhostNode(elem))
            {
                numnNodeNeighbours++;
                
            }
        }
        *this->mpOutStream << numnNodeNeighbours  << " ";

        numnNodeNeighbours = 0;
    }
    */
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationNeighbourWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("CellPopulationAreaWriter cannot be used with a CaBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationNeighbourWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("CellPopulationAreaWriter cannot be used with a NodeBasedCellPopulation");
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationNeighbourWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    assert(SPACE_DIM==2 || SPACE_DIM==3);

    int cell_number = pCellPopulation->GetNumElements();
    *this->mpOutStream << cell_number;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPopulationNeighbourWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("CellPopulationAreaWriter cannot be used with a VertexBasedCellPopulation");
}

// Explicit instantiation
template class CellPopulationNeighbourWriter<1,1>;
template class CellPopulationNeighbourWriter<1,2>;
template class CellPopulationNeighbourWriter<2,2>;
template class CellPopulationNeighbourWriter<1,3>;
template class CellPopulationNeighbourWriter<2,3>;
template class CellPopulationNeighbourWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellPopulationNeighbourWriter)
