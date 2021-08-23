/*

Copyright (c) 2005-2019, University of Oxford.
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
#include "DomMeshBasedCellPopulationWithGhostNodes.hpp"
#include "VtkMeshWriter.hpp"

// #include "AbstractCellMutationState.hpp"
// #include "WildTypeCellMutationState.hpp"
// #include "StromalCellMutationState.hpp"
#include "Debug.hpp"



template<unsigned DIM>
PeriodicNeighbourModifier<DIM>::PeriodicNeighbourModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
    mOutputDirectory(""),
    mCellPopulationWidth(0.0),
    mCellPopulationDepth(0.0),
    mCutoff(5.0)
{
}

template<unsigned DIM>
void PeriodicNeighbourModifier<DIM>::SetOutputDirectory(std::string outputDirectory)
{
	mOutputDirectory = outputDirectory;

    OutputFileHandler output_file_handler(mOutputDirectory+"/", false);
    out_stream locationFile = output_file_handler.OpenOutputFile("periodicNeighbours.dat");
    *locationFile << "time \t";
    *locationFile << "Cell ID" << "\t";
    *locationFile << "Neighbour Number " << "\t";
    *locationFile << "Cell Density " << "\n";
    locationFile->close();
}

template<unsigned DIM>
void PeriodicNeighbourModifier<DIM>::SetWidth(double width)
{
	mCellPopulationWidth = width;
}

template<unsigned DIM>
void PeriodicNeighbourModifier<DIM>::SetDepth(double depth)
{
	mCellPopulationDepth = depth;
}

template<unsigned DIM>
void PeriodicNeighbourModifier<DIM>::SetCutoff(double cutoff)
{
	mCutoff = cutoff;
}

// std::string MeshModifier::GetOutputDirectory()
// {
// 	return mOutputDirectory;
// }

template<unsigned DIM>
PeriodicNeighbourModifier<DIM>::~PeriodicNeighbourModifier()
{
}

template<unsigned DIM>
void PeriodicNeighbourModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM>& rCellPopulation)
{
}

template<unsigned DIM>
void PeriodicNeighbourModifier<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM>& rCellPopulation)
{
    
    DomMeshBasedCellPopulationWithGhostNodes<DIM>* p_tissue = static_cast<DomMeshBasedCellPopulationWithGhostNodes<DIM>*> (&rCellPopulation);

    
    unsigned num_nodes = rCellPopulation.GetNumNodes();
    std::vector<Node<DIM>*> extended_nodes(4*num_nodes);
	std::map<unsigned, unsigned> mExtendedMeshNodeIndexMap;
	MutableMesh<DIM,DIM>* mpExtendedMesh = nullptr;

	unsigned count = 0;
	// Dom - Create a copy of original mesh
    for (unsigned i=0; i<num_nodes; i++)
    {
        // First, create and store a copy of this real node and cell
        unsigned real_node_index = i;
        c_vector<double, DIM> real_node_location = rCellPopulation.GetNode(real_node_index)->rGetLocation();

        // Create a copy of the node corresponding to this cell and store it
        Node<DIM>* p_real_node = new Node<DIM>(real_node_index, real_node_location);
        extended_nodes[count] = p_real_node;

        // Populate mExtendedMeshNodeIndexMap
        mExtendedMeshNodeIndexMap[count] = real_node_index;

        count++;
    }

    for (unsigned i=0; i<num_nodes; i++)
    {
        // First, create and store a copy of this real node and cell
        unsigned real_node_index = i;
        c_vector<double, DIM> real_node_location = rCellPopulation.GetNode(real_node_index)->rGetLocation();

        // Compute the location of the image node corresponding to this node
        c_vector<double,DIM> image_node_location = real_node_location;
        if (real_node_location[0] >= mCellPopulationWidth*0.5)
        {
            image_node_location[0] -= mCellPopulationWidth;
        }
        else if (real_node_location[0] <  mCellPopulationWidth*0.5)
        {
            image_node_location[0] += mCellPopulationWidth;
        }

        // Create a copy of the node corresponding to this cell, suitable translated, and store it
        Node<DIM>* p_image_node = new Node<DIM>(count, image_node_location);
        extended_nodes[count] = p_image_node;

        // Populate mExtendedMeshNodeIndexMap
        mExtendedMeshNodeIndexMap[count] = real_node_index;

        count++;
    }

    for (unsigned i=0; i<num_nodes; i++)
    {
        // First, create and store a copy of this real node and cell
        unsigned real_node_index = i;
        c_vector<double, DIM> real_node_location = rCellPopulation.GetNode(real_node_index)->rGetLocation();

        // Compute the location of the image node corresponding to this node
        c_vector<double,DIM> image_node_location = real_node_location;

        if (real_node_location[1] >= mCellPopulationDepth*0.5)
        {
            image_node_location[1] -= mCellPopulationDepth;
        }
        else if (real_node_location[1] <  mCellPopulationDepth*0.5)
        {
            image_node_location[1] += mCellPopulationDepth;
        }

        // Create a copy of the node corresponding to this cell, suitable translated, and store it
        Node<DIM>* p_image_node = new Node<DIM>(count, image_node_location);
        extended_nodes[count] = p_image_node;

        // Populate mExtendedMeshNodeIndexMap
        mExtendedMeshNodeIndexMap[count] = real_node_index;

        count++;
    }

    for (unsigned i=0; i<num_nodes; i++)
    {
        // First, create and store a copy of this real node and cell
        unsigned real_node_index = i;
        c_vector<double, DIM> real_node_location = rCellPopulation.GetNode(real_node_index)->rGetLocation();

        // Compute the location of the image node corresponding to this node
        c_vector<double, DIM> image_node_location = real_node_location;

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
        Node<DIM>* p_image_node = new Node<DIM>(count, image_node_location);
        extended_nodes[count] = p_image_node;

        // Populate mExtendedMeshNodeIndexMap
        mExtendedMeshNodeIndexMap[count] = real_node_index;

        count++;
    }

    if (mpExtendedMesh != NULL)
    {
    	delete mpExtendedMesh;
    }
    mpExtendedMesh = new MutableMesh<DIM,DIM>(extended_nodes);

    OutputFileHandler output_file_handler(mOutputDirectory+"/", false);
	out_stream locationFile = output_file_handler.OpenOutputFile("periodicNeighbours.dat", std::ios::app);
    
    SimulationTime* p_time = SimulationTime::Instance();

	*locationFile << p_time->GetTime() << "\t";	

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
        cell_iter != rCellPopulation.End();
        ++cell_iter)
	{
        
        double cell_size = 0.0;

		unsigned index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(index);

        double x_location = rCellPopulation.GetLocationOfCellCentre(p_cell)[0];
		double y_location = rCellPopulation.GetLocationOfCellCentre(p_cell)[1];
		double z_location = rCellPopulation.GetLocationOfCellCentre(p_cell)[2];


        // Get pointer to this node
        Node<DIM>* p_node = mpExtendedMesh->GetNode(index);

        // Loop over containing elements
        std::vector<unsigned> neighbouring_node_indices;

        for (typename Node<DIM>::ContainingElementIterator elem_iter = p_node->ContainingElementsBegin();
         elem_iter != p_node->ContainingElementsEnd();
         ++elem_iter)
        {
            // Get pointer to this containing element
            // Element<3,3>* p_element = static_cast<MutableMesh<3,3>&>(mpExtendedMesh).GetElement(*elem_iter);
            Element<DIM,DIM>* p_element = mpExtendedMesh->GetElement(*elem_iter);

            // Loop over nodes contained in this element
            for (unsigned i=0; i<p_element->GetNumNodes(); i++)
            {
                // Get index of this node and add its index to the set if not the original node
                unsigned node_index = p_element->GetNodeGlobalIndex(i);

			    unsigned global_index = mExtendedMeshNodeIndexMap[node_index];

                if (node_index != index && !(p_tissue->IsGhostNode(global_index)) )
                {

                    if(std::find(neighbouring_node_indices.begin(), neighbouring_node_indices.end(), node_index) != neighbouring_node_indices.end()) 
                    {
                        /* neighbouring_node_indices contains node_index */
                        // Do nothing
                    } 
                    else 
                    {
                        /* neighbouring_node_indices does not contain node_index */
                        c_vector<double, DIM> cell_location = mpExtendedMesh->GetNode(node_index)->rGetLocation();
                        double dist = sqrt( pow(x_location - cell_location[0],2) + pow(y_location - cell_location[1],2) + pow(z_location - cell_location[2],2) );
                        
                        // PRINT_2_VARIABLES(node_index,dist);
                        
                        cell_size = cell_size + dist;

                        neighbouring_node_indices.push_back(node_index);
                    }

                    
                }
            }
        }
        
        int neigh_number = neighbouring_node_indices.size();
        // PRINT_VARIABLE(neigh_number);

        double cell_density = cell_size/(neigh_number);
        *locationFile << index << "\t";
        *locationFile << neigh_number << "\t";
        *locationFile << cell_density << "\t";

        

    }

    *locationFile << "\n";
    locationFile->close();

    if (mpExtendedMesh != NULL)
    {
    	delete mpExtendedMesh;
    }
}



template<unsigned DIM>
void PeriodicNeighbourModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateAtEndOfOutputTimeStep(rCellPopulation);


}


template<unsigned DIM>
void PeriodicNeighbourModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
// template class MeshModifier<1>;
// template class MeshModifier<2>;
template class PeriodicNeighbourModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PeriodicNeighbourModifier)

