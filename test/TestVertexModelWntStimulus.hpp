#include <cxxtest/TestSuite.h>
#include "CellBasedSimulationArchiver.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "CellsGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
// #include "ToroidalHoneycombMeshGenerator.hpp"
// #include "Toroidal2dMesh.hpp"
#include "ToroidalHoneycombVertexMeshGenerator.hpp"
#include "Toroidal2dVertexMesh.hpp"

#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "OnLatticeSimulation.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "RandomCellKiller.hpp"
#include "RepulsionForce.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "VoronoiDataWriter.hpp"
#include "DeltaNotchSrnModel.hpp"

#include "FakePetscSetup.hpp"

#include "TransitCellProliferativeType.hpp"

#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "DeltaNotchTrackingModifier.hpp"

#include "DomSimpleWntCellCycleModel.hpp"
#include "DomWntConcentration.hpp"

#include "AreaVertexCellKiller.hpp"
#include "T2SwapCellKiller.hpp"
#include "NodeVelocityWriter.hpp"
#include "CellIdWriter.hpp"

#include "Timer.hpp"
#include "RandomMotionForce.hpp"
#include "DomTargetAreaModifier.hpp"

class SimulationVT_1 : public AbstractCellBasedTestSuite
{
public:

    /* We next show how to modify the previous test to include 'ghost nodes', which do not
     * correspond to cells but are sometimes needed when using a Voronoi tessellation. We
     * will discuss ghost nodes in more detail in subsequent cell-based tutorials.
     */
    // void TestMeshBasedMonolayerWithGhostNodes()
    // {
    //     EXIT_IF_PARALLEL;

    //     // CylindricalHoneycombMeshGenerator generator(6, 8, 2);
    //     // Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

    //     ToroidalHoneycombVertexMeshGenerator generator(4, 4);
    //     Toroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();

    //     std::vector<CellPtr> cells;
    //     MAKE_PTR(TransitCellProliferativeType, p_transit_type);
    //     CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
    //     cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

    //     VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
    //     OffLatticeSimulation<2> simulator(cell_population);
        
    //     MAKE_PTR(NagaiHondaForce<2>, p_force);
    //     simulator.AddForce(p_force);

    //     simulator.SetOutputDirectory("CellBased_VT"); //**Changed**//
    //     simulator.SetSamplingTimestepMultiple(1);
    //     simulator.SetEndTime(1.0); //**Changed**//

    //     simulator.Solve();

    // }

    void TestVertexBasedMonolayer()
    {
        /* We include the next line because Vertex simulations cannot be run in parallel */
        EXIT_IF_PARALLEL;

        RandomNumberGenerator::Instance()->Reseed(2);

        Timer::Reset();

        /* First we create a regular vertex mesh. */
        unsigned width = 14;	   // x
        unsigned height = 18;      // y
        double width_space = 1.0;
        double height_space = 1.0*sqrt(0.75);
        double periodic_width = (double) (width+0.0)*width_space;
        double periodic_height = (double) (height+0.0)*height_space;

        ToroidalHoneycombVertexMeshGenerator generator(width, height);
        Toroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            // UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();

            DomSimpleWntCellCycleModel* p_cc_model = new DomSimpleWntCellCycleModel();

            p_cc_model->SetDimension(2);

            p_cc_model->SetTransitCellG1Duration(4);
            p_cc_model->SetSDuration(2);
            p_cc_model->SetG2Duration(1);
            p_cc_model->SetMDuration(1);

            /* We choose to initialise the concentrations to random levels in each cell. */
            // std::vector<double> initial_conditions;
            // initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            // DeltaNotchSrnModel* p_srn_model = new DeltaNotchSrnModel();
            // p_srn_model->SetInitialConditions(initial_conditions);

            // CellPtr p_cell(new Cell(p_state, p_cc_model, p_srn_model));
            CellPtr p_cell(new Cell(p_state, p_cc_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            
            double birth_time = -8-RandomNumberGenerator::Instance()->ranf()*(p_cc_model->GetTransitCellG1Duration() + p_cc_model->GetSG2MDuration());
            
            p_cell->SetApoptosisTime(1.5);

            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        DomWntConcentration<2>::Instance()->SetType(DomRADIAL);
        DomWntConcentration<2>::Instance()->SetCellPopulation(cell_population);
        DomWntConcentration<2>::Instance()->SetCryptLength(20.0);
        DomWntConcentration<2>::Instance()->SetCryptCentreX(0.5*periodic_width);
        DomWntConcentration<2>::Instance()->SetCryptCentreY(0.5*periodic_height);
        DomWntConcentration<2>::Instance()->SetCryptRadius(2);
        DomWntConcentration<2>::Instance()->SetWntConcentrationParameter(2.0);


        /* Using the vertex mesh and cells, we create a cell-based population object, and specify which results to
         * output to file. */
        // cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        // cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        // cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
        // cell_population.AddCellWriter<CellProliferativePhasesWriter>();
        // cell_population.AddCellWriter<CellAgesWriter>();
        // cell_population.AddCellWriter<CellVolumesWriter>();
        // cell_population.AddPopulationWriter<NodeVelocityWriter>();
        cell_population.AddCellWriter<CellIdWriter>();


        /* We are now in a position to create and configure the cell-based simulation object, pass a force law to it,
         * and run the simulation. We can make the simulation run for longer to see more patterning by increasing the end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVertexBasedMonolayer");
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetDt(0.01);
        simulator.SetEndTime(144.0);

        
        MAKE_PTR_ARGS(AreaVertexCellKiller<2>, p_killer, (&cell_population, 0.825, periodic_width, periodic_height));
        simulator.AddCellKiller(p_killer);

        // MAKE_PTR_ARGS(T2SwapCellKiller<2>, p_t2_killer, (&cell_population));
        // simulator.AddCellKiller(p_t2_killer);


        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);
        /* This modifier assigns target areas to each cell, which are required by the {{{NagaiHondaForce}}}.
         */
        MAKE_PTR(DomTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(0.05); //0.1 causes dissasociation, 0.001 is not enough
        simulator.AddForce(p_random_force);
        
        
        simulator.Solve();

        Timer::Print("Time Ellapsed");
    }

    /* To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dCentreCells /tmp/$USER/testoutput/CellBased_VT/results_from_time_0}}}.
     * We may have to do: {{{javac Visualize2dCentreCells.java}}} beforehand to create the
     * java executable.
     */

};
