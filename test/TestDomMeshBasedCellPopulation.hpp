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

#ifndef TESTMESHBASEDCELLPOPULATION_HPP_
#define TESTMESHBASEDCELLPOPULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "DomMeshBasedCellPopulationWithGhostNodes.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ArchiveOpener.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellLabel.hpp"
#include "CellAncestor.hpp"
#include "CellId.hpp"
#include "CellPropertyRegistry.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "ApoptoticCellProperty.hpp"
#include "FixedCentreBasedDivisionRule.hpp"

// Cell writers
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellLabelWriter.hpp"
#include "CellLocationIndexWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellVolumesWriter.hpp"

// Cell population writers
#include "CellPopulationAreaWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "VoronoiDataWriter.hpp"

#include "PetscSetupAndFinalize.hpp"
#include "Debug.hpp"

class TestMeshBasedCellPopulation : public AbstractCellBasedTestSuite
{
private:

    template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
    void TestSmallMeshBasedCellPopulation(std::string meshFilename)
    {
        // Create a simple mesh
        TrianglesMeshReader<ELEMENT_DIM,SPACE_DIM> mesh_reader(meshFilename);
        MutableMesh<ELEMENT_DIM,SPACE_DIM> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, SPACE_DIM> generator;
        generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create the cell population
        unsigned num_cells = cells.size();
        MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM> cell_population(mesh, cells);

        TS_ASSERT_EQUALS(cell_population.rGetMesh().GetNumNodes(), mesh.GetNumNodes());
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), num_cells);

        // Test set/get method of member variables
        TS_ASSERT_DELTA(cell_population.GetMeinekeDivisionSeparation(), 0.3, 1e-6);
        cell_population.SetMeinekeDivisionSeparation(0.5);
        TS_ASSERT_DELTA(cell_population.GetMeinekeDivisionSeparation(), 0.5, 1e-6);
        cell_population.SetMeinekeDivisionSeparation(0.3);

        TS_ASSERT_DELTA(cell_population.GetAbsoluteMovementThreshold(), 2.0, 1e-6);
        cell_population.SetAbsoluteMovementThreshold(1.5);
        TS_ASSERT_DELTA(cell_population.GetAbsoluteMovementThreshold(), 1.5, 1e-6);
        cell_population.SetMeinekeDivisionSeparation(0.5);

        unsigned counter = 0;
        for (typename AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Test operator* and that cells are in sync
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), counter);

            // Test operator-> and that cells are in sync
            TS_ASSERT_DELTA(cell_iter->GetAge(), (double)counter, 1e-12);

            counter++;
        }

        TS_ASSERT_EQUALS(counter, cell_population.GetNumRealCells());
        //Since no cells died, this should be all cells in the mesh.
        TS_ASSERT_EQUALS(counter, cell_population.GetNumAllCells());
    }

public:
    // This test checks that the cells and nodes are correctly archived.
    void TestArchivingMeshBasedCellPopulation()
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "dom_mesh_based_cell_population.arch";
        ArchiveLocationInfo::SetMeshFilename("dom_mesh_based_cell_population_mesh");

        std::vector<c_vector<double,2> > cell_locations;

        // Archive a cell population
        {
            // Need to set up time
            unsigned num_steps = 10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);

            // Create a simple mesh
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
            MutableMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            // Set up cells, one for each node. Give each a birth time of -node_index,
            // so the age = node_index
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

            std::vector<unsigned> cell_location_indices;
            for (unsigned i=0; i<cells.size(); i++)
            {
                cell_location_indices.push_back(i);
            }

            // Create the cell population
            DomMeshBasedCellPopulationWithGhostNodes<2>* const p_cell_population = new DomMeshBasedCellPopulationWithGhostNodes<2>(mesh, cells, cell_location_indices);

            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // loop over them to run to time 0.0;
            unsigned index_for_data = 0;
            for (AbstractCellPopulation<2>::Iterator cell_iter=p_cell_population->Begin();
                 cell_iter!=p_cell_population->End();
                 ++cell_iter)
            {
                cell_iter->ReadyToDivide();
                cell_locations.push_back(p_cell_population->GetLocationOfCellCentre(*cell_iter));
                // Add cell data
                cell_iter->GetCellData()->SetItem("data", (double) index_for_data);
                index_for_data++;
            }

            std::pair<CellPtr,CellPtr> cell_pair_0_1 = p_cell_population->CreateCellPair(p_cell_population->GetCellUsingLocationIndex(0), p_cell_population->GetCellUsingLocationIndex(1));
            p_cell_population->MarkSpring(cell_pair_0_1);

            // Set area-based viscosity
            p_cell_population->SetAreaBasedDampingConstant(true);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Write the cell population to the archive
            (*p_arch) << static_cast<const SimulationTime&> (*p_simulation_time);

            TRACE("Before Save");

            (*p_arch) << p_cell_population;

            TRACE("After Save");

            SimulationTime::Destroy();

            TRACE("After Destroy");

            delete p_cell_population;
        }

        // Restore cell population
        {
            TRACE("Before Loop");

            // Need to set up time
            unsigned num_steps=10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);
            p_simulation_time->IncrementTimeOneStep();

            DomMeshBasedCellPopulationWithGhostNodes<2>* p_cell_population;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();
            (*p_arch) >> *p_simulation_time;

            TRACE("Before Load");

            (*p_arch) >> p_cell_population;

            TRACE("After Load");

            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // this checks that individual cells and their models are archived.
            unsigned counter = 0;
            for (AbstractCellPopulation<2>::Iterator cell_iter=p_cell_population->Begin();
                 cell_iter!=p_cell_population->End();
                 ++cell_iter)
            {
                TS_ASSERT_DELTA(cell_iter->GetAge(),(double)(counter),1e-7);
                TS_ASSERT_DELTA(p_cell_population->GetLocationOfCellCentre(*cell_iter)[0], cell_locations[counter][0], 1e-9);
                TS_ASSERT_DELTA(p_cell_population->GetLocationOfCellCentre(*cell_iter)[1], cell_locations[counter][1], 1e-9);
                TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("data"), (double) p_cell_population->GetLocationIndexUsingCell(*cell_iter), 1e-12);
                counter++;
            }

            TRACE("After Loop");

            // TS_ASSERT_EQUALS(p_cell_population->rGetMesh().GetNumNodes(), 5u);

            // TRACE("After GetMesh().GetNumNodes");


            TS_ASSERT_EQUALS(p_cell_population->GetNode(0)->IsBoundaryNode(), true);
            TS_ASSERT_EQUALS(p_cell_population->GetNode(1)->IsBoundaryNode(), true);
            TS_ASSERT_EQUALS(p_cell_population->GetNode(2)->IsBoundaryNode(), true);
            TS_ASSERT_EQUALS(p_cell_population->GetNode(3)->IsBoundaryNode(), true);
            TS_ASSERT_EQUALS(p_cell_population->GetNode(4)->IsBoundaryNode(), false);

            // // Check the marked spring
            std::pair<CellPtr,CellPtr> cell_pair_0_1 = p_cell_population->CreateCellPair(p_cell_population->GetCellUsingLocationIndex(0), p_cell_population->GetCellUsingLocationIndex(1));
            TS_ASSERT_EQUALS(p_cell_population->IsMarkedSpring(cell_pair_0_1), true);

            TRACE("A");
            // Check the simulation time has been restored (through the cell)
            TS_ASSERT_EQUALS(p_simulation_time->GetTime(), 0.0);

            TRACE("B");

            // Check the cell population has been restored
            TS_ASSERT_EQUALS(p_cell_population->rGetCells().size(), 5u);

            TRACE("C");

            // Check area-based viscosity is still true
            TS_ASSERT_EQUALS(p_cell_population->UseAreaBasedDampingConstant(), true);

            TRACE("D");

            // TS_ASSERT_EQUALS((p_cell_population->GetGhostNodeIndicies()).size(), 0u);


            TS_ASSERT_EQUALS(p_cell_population->rGetMesh().GetNumNodes(), 5u);

            TRACE("Before Delete");

            delete p_cell_population;

            TRACE("After delete");
        }
    }

};

#endif /*TESTMESHBASEDCELLPOPULATION_HPP_*/
