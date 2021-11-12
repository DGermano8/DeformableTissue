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

#include "PointWriterModifier.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "VtkMeshWriter.hpp"
#include "NodesOnlyMesh.hpp"

#include "AbstractCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StromalCellMutationState.hpp"
#include "Debug.hpp"



template<unsigned DIM>
PointWriterModifier<DIM>::PointWriterModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
    mOutputDirectory(""),
    mWidth(0.0),
    mDepth(0.0),
    mCutoff(5.0)
{
}

template<unsigned DIM>
void PointWriterModifier<DIM>::SetOutputDirectory(std::string outputDirectory)
{
	mOutputDirectory = outputDirectory;
}

template<unsigned DIM>
void PointWriterModifier<DIM>::SetWidth(double width)
{
	mWidth = width;
}

template<unsigned DIM>
void PointWriterModifier<DIM>::SetDepth(double depth)
{
	mDepth = depth;
}

template<unsigned DIM>
void PointWriterModifier<DIM>::SetCutoff(double cutoff)
{
	mCutoff = cutoff;
}

// std::string MeshModifier::GetOutputDirectory()
// {
// 	return mOutputDirectory;
// }

template<unsigned DIM>
PointWriterModifier<DIM>::~PointWriterModifier()
{
}

template<unsigned DIM>
void PointWriterModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM>& rCellPopulation)
{
}

template<unsigned DIM>
void PointWriterModifier<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM>& rCellPopulation)
{
    MeshBasedCellPopulationWithGhostNodes<DIM>* p_tissue = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

    // if (true)
    {
        //std::cout<< "\n";
        //PRINT_VARIABLE("Start ");
        
        // Create mesh writer for VTK output
        std::ostringstream time_string;
        time_string << SimulationTime::Instance()->GetTimeStepsElapsed();
        VtkMeshWriter<DIM, DIM> cells_writer(mOutputDirectory, "results_"+time_string.str(), false);

        
        // Iterate over any cell writers that are present
        unsigned num_cells = rCellPopulation.GetNumAllCells();
        unsigned total_num_nodes = rCellPopulation.rGetMesh().GetNumNodes();
        unsigned num_cell_data_items = rCellPopulation.Begin()->GetCellData()->GetNumItems();
        std::vector<std::string> cell_data_names = rCellPopulation.Begin()->GetCellData()->GetKeys();

        //PRINT_2_VARIABLES(total_num_nodes, num_cell_data_items);
        
        unsigned cell_data[num_cell_data_items][total_num_nodes];

        std::vector<double> ghosts(total_num_nodes);

        // for (typename std::vector<boost::shared_ptr<AbstractCellWriter<DIM, DIM> > >::iterator cell_writer_iter = rCellPopulation.mCellWriters.begin();
        //      cell_writer_iter != rCellPopulation.mCellWriters.end();
        //      ++cell_writer_iter)
        // {
        //     // Create vector to store VTK cell data
        //     std::vector<double> vtk_cell_data(total_num_nodes);

        //     // Loop over cells
        //     for (unsigned node_index=0; node_index<total_num_nodes; node_index++)
        //     {
        //         // If this node corresponds to a ghost node, set any "cell" data to be -1.0
        //         if (p_tissue->IsGhostNode(node_index))
        //         {
        //             // Populate the vector of VTK cell data
        //             vtk_cell_data[node_index] = -1.0;
        //         }
        //         else
        //         {
        //             // Get the cell corresponding to this node
        //             CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(node_index);

        //             // Populate the vector of VTK cell data
        //             vtk_cell_data[node_index] = (*cell_writer_iter)->GetCellDataForVtkOutput(p_cell, rCellPopulation);
        //             //PRINT_VARIABLE((*cell_writer_iter)->GetCellDataForVtkOutput(p_cell, this));
        //         }
                
        //     }

        //     cells_writer.AddPointData((*cell_writer_iter)->GetVtkCellDataName(), vtk_cell_data);
        // }
        // Loop over cells
        for (unsigned node_index=0; node_index<total_num_nodes; node_index++)
        {

            // If this node corresponds to a ghost node, set any "cell" data to be -1.0
            if (p_tissue->IsGhostNode(node_index))
            {
                for (unsigned var=0; var<num_cell_data_items; var++)
                {
                    cell_data[var][node_index] = -1;
                    
                }
            }
            else
            {
                // Get the cell corresponding to this node
                CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(node_index);

                for (unsigned var=0; var<num_cell_data_items; var++)
                {
                    // TRACE("adding data");
                    cell_data[var][node_index] = p_cell->GetCellData()->GetItem(cell_data_names[var]);
                    // TRACE("got data");
                }
            }
            
        }
        
        // for (unsigned var=0; var<num_cell_data_items; var++)
        // {
        //     cells_writer.AddPointData(cell_data_names[var], cell_data[var]);
        // }
       
        // Make a copy of the nodes in a disposable mesh for writing
        {
            std::vector<Node<DIM>* > nodes;
            for (unsigned index=0; index<rCellPopulation.rGetMesh().GetNumNodes(); index++)
            {
            //
            //PRINT_VARIABLE((this->IsGhostNode(index)));
            //
            Node<DIM>* p_node = rCellPopulation.rGetMesh().GetNode(index);
            nodes.push_back(p_node);
            
            /*    if(!(this->IsGhostNode(index)))
                {
                    //TRACE(this->mrMesh.GetNode(index));
                    Node<DIM>* p_node = this->mrMesh.GetNode(index);
                    nodes.push_back(p_node);
                }
                else if (this->IsGhostNode(index))
                {
                    //Do something...
                }
                //PRINT_2_VARIABLES(index,(this->IsGhostNode(index)));
            */
            }
            NodesOnlyMesh<3> mesh;
            mesh.ConstructNodesWithoutMesh(nodes, 1.5); // Arbitrary cut off as connectivity not used.
            cells_writer.WriteFilesUsingMesh(mesh);
        }
        //PRINT_2_VARIABLES("Middle",num_timesteps);
        
        unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();

        
        // *(rCellPopulation.mpVtkMetaFile) << "        <DataSet timestep=\"";
        // *(rCellPopulation.mpVtkMetaFile) << num_timesteps;
        // *(rCellPopulation.mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
        // *(rCellPopulation.mpVtkMetaFile) << num_timesteps;
        // *(rCellPopulation.mpVtkMetaFile) << ".vtu\"/>\n";
        
        //PRINT_VARIABLE("Finish");
    }
}



template<unsigned DIM>
void PointWriterModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateAtEndOfOutputTimeStep(rCellPopulation);


}


template<unsigned DIM>
void PointWriterModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
// template class MeshModifier<1>;
// template class MeshModifier<2>;
template class PointWriterModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PointWriterModifier)

