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

#include "PeriodicNeighbourModifier.hpp"
#include "AbstractCellPopulation.hpp"

#include "DomMeshBasedCellPopulationWithGhostNodes.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PeriodicNeighbourModifier<ELEMENT_DIM, SPACE_DIM>::Visit(DomMeshBasedCellPopulationWithGhostNodes<SPACE_DIM>* pCellPopulation)
{
    assert(SPACE_DIM==3);

    DomMeshBasedCellPopulationWithGhostNodes<3>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<3>*> (pCellPopulation);

    
    unsigned num_nodes = pCellPopulation->GetNumNodes();
    std::vector<Node<3>*> extended_nodes(4*num_nodes);
	std::map<unsigned, unsigned> mExtendedMeshNodeIndexMap;
	MutableMesh<3,3>* mpExtendedMesh = nullptr;

	unsigned count = 0;
	// Dom - Create a copy of original mesh
    for (unsigned i=0; i<num_nodes; i++)
    {
        // First, create and store a copy of this real node and cell
        unsigned real_node_index = i;
        c_vector<double, 3> real_node_location = pCellPopulation->GetNode(real_node_index)->rGetLocation();

        // Create a copy of the node corresponding to this cell and store it
        Node<3>* p_real_node = new Node<3>(real_node_index, real_node_location);
        extended_nodes[count] = p_real_node;

        // Populate mExtendedMeshNodeIndexMap
        mExtendedMeshNodeIndexMap[count] = real_node_index;

        count++;
    }

    for (unsigned i=0; i<num_nodes; i++)
    {
        // First, create and store a copy of this real node and cell
        unsigned real_node_index = i;
        c_vector<double, 3> real_node_location = pCellPopulation->GetNode(real_node_index)->rGetLocation();

        // Compute the location of the image node corresponding to this node
        c_vector<double,3> image_node_location = real_node_location;
        if (real_node_location[0] >= mCellPopulationWidth*0.5)
        {
            image_node_location[0] -= mCellPopulationWidth;
        }
        else if (real_node_location[0] <  mCellPopulationWidth*0.5)
        {
            image_node_location[0] += mCellPopulationWidth;
        }

        // Create a copy of the node corresponding to this cell, suitable translated, and store it
        Node<3>* p_image_node = new Node<3>(count, image_node_location);
        extended_nodes[count] = p_image_node;

        // Populate mExtendedMeshNodeIndexMap
        mExtendedMeshNodeIndexMap[count] = real_node_index;

        count++;
    }

    for (unsigned i=0; i<num_nodes; i++)
    {
        // First, create and store a copy of this real node and cell
        unsigned real_node_index = i;
        c_vector<double, 3> real_node_location = pCellPopulation->GetNode(real_node_index)->rGetLocation();

        // Compute the location of the image node corresponding to this node
        c_vector<double,3> image_node_location = real_node_location;

        if (real_node_location[1] >= mCellPopulationDepth*0.5)
        {
            image_node_location[1] -= mCellPopulationDepth;
        }
        else if (real_node_location[1] <  mCellPopulationDepth*0.5)
        {
            image_node_location[1] += mCellPopulationDepth;
        }

        // Create a copy of the node corresponding to this cell, suitable translated, and store it
        Node<3>* p_image_node = new Node<3>(count, image_node_location);
        extended_nodes[count] = p_image_node;

        // Populate mExtendedMeshNodeIndexMap
        mExtendedMeshNodeIndexMap[count] = real_node_index;

        count++;
    }

    for (unsigned i=0; i<num_nodes; i++)
    {
        // First, create and store a copy of this real node and cell
        unsigned real_node_index = i;
        c_vector<double, 3> real_node_location = pCellPopulation->GetNode(real_node_index)->rGetLocation();

        // Compute the location of the image node corresponding to this node
        c_vector<double,3> image_node_location = real_node_location;

        if (real_node_location[1] >= mCellPopulationDepth*0.5)
        {
            image_node_location[1] -= mCellPopulationDepth;
        }
        else if (real_node_location[1] <  mCellPopulationDepth*0.5)
        {
            image_node_location[1] += mCellPopulationDepth;
        }
		if (real_node_location[0] >= mCellPopulationWidth*0.5)
        {
            image_node_location[0] -= mCellPopulationWidth;
        }
        else if (real_node_location[0] <  mCellPopulationWidth*0.5)
        {
            image_node_location[0] += mCellPopulationWidth;
        }

        // Create a copy of the node corresponding to this cell, suitable translated, and store it
        Node<3>* p_image_node = new Node<3>(count, image_node_location);
        extended_nodes[count] = p_image_node;

        // Populate mExtendedMeshNodeIndexMap
        mExtendedMeshNodeIndexMap[count] = real_node_index;

        count++;
    }

    if (mpExtendedMesh != NULL)
    {
    	delete mpExtendedMesh;
    }
    mpExtendedMesh = new MutableMesh<3,3>(extended_nodes);

    for (typename AbstractCellPopulation<3>::Iterator cell_iter = mpCellPopulation->Begin();
        cell_iter != mpCellPopulation->End();
        ++cell_iter)
	{
        double cell_size = 0.0;

		unsigned index = mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
        CellPtr p_cell = pCellPopulation->GetCellUsingLocationIndex(index);

        double x_location = pCellPopulation->GetLocationOfCellCentre(p_cell)[0];
		double y_location = pCellPopulation->GetLocationOfCellCentre(p_cell)[1];
		double z_location = pCellPopulation->GetLocationOfCellCentre(p_cell)[2];


        // Get pointer to this node
        Node<3>* p_node = mpExtendedMesh.GetNode(index);

        // Loop over containing elements
        std::set<unsigned> neighbouring_node_indices;

        for (typename Node<3>::ContainingElementIterator elem_iter = p_node->ContainingElementsBegin();
         elem_iter != p_node->ContainingElementsEnd();
         ++elem_iter)
        {
            // Get pointer to this containing element
            Element<3,3>* p_element = static_cast<MutableMesh<3,3>&>((mpExtendedMesh)).GetElement(*elem_iter);

            // Loop over nodes contained in this element
            for (unsigned i=0; i<p_element->GetNumNodes(); i++)
            {
                // Get index of this node and add its index to the set if not the original node
                unsigned node_index = p_element->GetNodeGlobalIndex(i);
                if (node_index != index && (!p_tissue->IsGhostNode(node_index)) )
                {
                    CellPtr p_cell_neigh = pCellPopulation->GetCellUsingLocationIndex(node_index);

                    double x_location_neigh = pCellPopulation->GetLocationOfCellCentre(p_cell_neigh)[0];
                    double y_location_neigh = pCellPopulation->GetLocationOfCellCentre(p_cell_neigh)[1];
                    double z_location_neigh = pCellPopulation->GetLocationOfCellCentre(p_cell_neigh)[2];

                    double dist = sqrt( pow(x_location - x_location_neigh,2) + pow(y_location - y_location_neigh,2) + pow(z_location - z_location_neigh,2));
                    cell_size = cell_size + dist;
                    neighbouring_node_indices.insert(node_index);
                }
            }
        }
        cell_size = cell_size/neighbouring_node_indices.size();
        *this->mpOutStream << index << ' ' << neighbouring_node_indices.size() << ' ' << cell_size;

    }


}


// Explicit instantiation
// template class PeriodicNeighbourModifier<1,1>;
// template class PeriodicNeighbourModifier<1,2>;
// template class PeriodicNeighbourModifier<2,2>;
// template class PeriodicNeighbourModifier<1,3>;
// template class PeriodicNeighbourModifier<2,3>;
template class PeriodicNeighbourModifier<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(PeriodicNeighbourModifier)
