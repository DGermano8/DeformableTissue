#ifndef TEST3DBOXMODEL_HPP_
#define TEST3DBOXMODEL_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "OffLatticeSimulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StromalCellMutationState.hpp"
#include "RandomMotionForce.hpp"
#include "PeriodicCryptModelInteractionForceWithGhostNodes.hpp"
#include "PeriodicBendingForce3dHeightWithGhostNodes.hpp"
#include "SloughingCellKiller3DWithGhostNodes.hpp"
#include "AnoikisCellKiller3DWithGhostNodes.hpp"
#include "UniformCellKiller3dWithGhostNodes.hpp"
#include "Debug.hpp"

#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "CellIdWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellPopulationEpithelialWriter.hpp"
#include "VoronoiDataWriter.hpp"

#include "PetscSetupAndFinalize.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "PlanarDivisionRule.hpp"
#include "DriftPreventForce.hpp"


#include "MeshRemeshModifier.hpp"


#include "CylindricalHoneycombMeshGenerator.hpp"

#include "BendingForce2DGhosts.hpp"



// Tests have been copied directly from TestOffLatticeSimulation3d.hpp - most likely that bits will get changed
// and broken by me along the way

class Test3dBoxModel : public AbstractCellBasedTestSuite
{
private:

public:
    
    /* A cube that consists of a block of stromal cells, and a single layer of epithelial cells.
     * Target curvature for the basement membrane is zero, and multiple division events occur.
     */
    void TestPeriodicCubeWithGhosts() throw (Exception)
    {
    	RandomNumberGenerator::Instance()->Reseed(2);


        unsigned width = 30;	   // x
        unsigned height = 10;      // y
        unsigned ghosts_top = 3;

        // Initialise the tissue in an equilibrum state
        double width_space = 1.0;
        double height_space = 1.0*sqrt(0.75);
        double ghost_sep = 1.0;
        unsigned cells_per_layer = width*height;
        unsigned cell_iter = 0;

        double periodic_width = (double) (width+0.0)*width_space;
        double periodic_height = (double) (height+0.0)*height_space;

        double tissue_base = 4.0; //Hieght of tissue to prevent drift
        double tissue_middle = 0.0;

        bool isGhost = false;
        double x_coordinate, y_coordinate;

        // Get the dimensions for the non-zero target curvature region
        double centre_x = periodic_width*0.5;
        double centre_y = periodic_height*0.5;
        // double centre_x = 20;
        // double centre_y = 20;

        double spring_strength = 25.0;
        bool add_springs = true;

        double radius =  4.5;//periodic_width+1.0;q
        double target_curvature = 0.45; // higher curvature means smaller sphere
        double beta_parameter = 0.05*spring_strength;
        double alpha_parameter = 1.5;

        double time_step = 0.005;
        double end_time = 12;

        bool include_springs = true;
        bool include_bending = true;

        int is_transit[height*width];
        int num_real_nodes = 0;

        
        //mesh.Translate(0.5, 0.5);

        CylindricalHoneycombMeshGenerator generator(width, height, ghosts_top);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        /* We create the cells, using the same method as before. Here, though, we use a {{{SimpleWntCellCycleModel}}}.*/
        // CryptCellsGenerator<SimpleWntCellCycleModel> cells_generator;
        // cells_generator.Generate(cells, p_mesh, location_indices, true);


        
        // Set up cells by iterating through the mesh nodes
        unsigned num_epithelial_cells = (width)*(height);
        unsigned num_tissue_cells = (width)*(height);

        boost::shared_ptr<AbstractCellMutationState> p_epithelial_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_stromal_state(new StromalCellMutationState);

        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        std::vector<CellPtr> cells;
        // Initialise Tissue cells (Stromal)
		for (unsigned i=0; i<(width)*(height-1); i++)
		{
			//StochasticDurationGenerationBasedCellCycleModel* p_model = new StochasticDurationGenerationBasedCellCycleModel();
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();

            CellPtr p_differentiated_cell(new Cell(p_stromal_state, p_model));
            p_differentiated_cell->SetCellProliferativeType(p_differentiated_type);
            p_model->SetDimension(2);

            double birth_time = - 10.0;
            p_differentiated_cell->SetBirthTime(birth_time);

			cells.push_back(p_differentiated_cell);
        }


        // Initialise Epithelial cells
        for (unsigned i=(width)*(height-1); i<(width)*(height); i++)
        {

        	// FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            // p_model->SetMaxTransitGenerations(1);
            // p_model->SetSDuration(2); 
            // p_model->SetG2Duration(2); 
            // p_model->SetMDuration(2); 
            // p_model->SetDimension(3);
            
            UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetMaxTransitGenerations(2);

            // p_model->SetTransitCellG1Duration(4);
            // p_model->SetSDuration(2);
            // p_model->SetG2Duration(2);
            // p_model->SetMDuration(1);

            // p_model->SetTransitCellG1Duration(12);
            // p_model->SetSDuration(6);
            // p_model->SetG2Duration(4);
            // p_model->SetMDuration(2);

            p_model->SetTransitCellG1Duration(2);
            p_model->SetSDuration(1);
            p_model->SetG2Duration(1);
            p_model->SetMDuration(1);
            
            CellPtr p_epithelial_cell(new Cell(p_epithelial_state, p_model));
            double birth_time = -8;

            // if(i == 154)
            // {
            //     p_epithelial_cell->SetCellProliferativeType(p_transit_type);
            //     birth_time = -(12+6+4+2) - 12;
            // }
            // else
            // {
            //     p_epithelial_cell->SetCellProliferativeType(p_differentiated_type);
            //     birth_time = -10;
            // }
            
            
            // if( i == 98 )
            // {
            //     p_epithelial_cell->SetCellProliferativeType(p_transit_type);
            //     birth_time = -RandomNumberGenerator::Instance()->ranf()*(p_model->GetTransitCellG1Duration() + p_model->GetSG2MDuration() );

            // }
            // else
            // {
            //     p_epithelial_cell->SetCellProliferativeType(p_differentiated_type);
            //     birth_time = -10.0;
            // }

            // if(i == 154)
            // {
            //     birth_time = -10;
            // }

            // p_epithelial_cell->SetCellProliferativeType(p_transit_type);
            p_epithelial_cell->SetCellProliferativeType(p_differentiated_type);

            // double birth_time = -RandomNumberGenerator::Instance()->ranf()*(p_model->GetTransitCellG1Duration() + p_model->GetSG2MDuration());
            
			
            p_epithelial_cell->SetBirthTime(birth_time);
            p_epithelial_cell->SetApoptosisTime(1);
            
            p_epithelial_cell->InitialiseCellCycleModel();
            
			cells.push_back(p_epithelial_cell);
        }


        

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        //MeshBasedCellPopulationWithGhostNodes<3> cell_population(mesh, cells, real_node_indices);
//        MeshBasedCellPopulation<3> cell_population(mesh, cells);
        assert(cell_population.GetNumRealCells() != 0);


        // Set the division rule for our population to be the random direction division rule
        // boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_division_rule_to_set(new PlanarDivisionRule());
        // cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

        OffLatticeSimulation<2> simulator(cell_population);

        // Make sure we have a Voronoi tessellation to begin with
        cell_population.CreateVoronoiTessellation();
        
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();
        cell_population.AddPopulationWriter<CellPopulationEpithelialWriter>();
        cell_population.AddPopulationWriter<VoronoiDataWriter>(); // paraview is pretty pointless at viewing this, worth looking into
        
        //cell_population.WriteVtkResultsToFile(output_directory);
        //To fix paraview
        cell_population.SetWriteVtkAsPoints(true);
        //std::cout<<cell_population.GetWriteVtkAsPoints() << "\n";
        //PRINT_VARIABLE(cell_population.GetWriteVtkAsPoints());
        cell_population.SetOutputMeshInVtk(true);


        std::map<Node<2>*, c_vector<double,2> > node_locations_before;
        for (std::list<CellPtr>::iterator cell_iter = cell_population.rGetCells().begin();
             cell_iter != cell_population.rGetCells().end();
             ++cell_iter)
        {
            Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
            node_locations_before[p_node] = p_node->rGetLocation();
            //PRINT_VECTOR(node_locations_before[p_node]);
        }

		// Create periodic spring force law
        // MAKE_PTR(PeriodicCryptModelInteractionForceWithGhostNodes<3>, periodic_spring_force);
        // periodic_spring_force->SetUseOneWaySprings(false); //turning this on makes the stromal cells act as ghosts..
        // periodic_spring_force->SetCutOffLength(1.2);
        // //                     SetEpithelialStromalCellDependentSprings(ind , Ep-Ep, Str-Str, Ep-Str, apcTwoHitStromalMultiplier);
        // periodic_spring_force->SetEpithelialStromalCellDependentSprings(true, 1.0,     0.5,     0.5,    1.0);
        // periodic_spring_force->SetPeriodicDomainWidth(periodic_width);
        // periodic_spring_force->SetPeriodicDomainDepth(periodic_height);
        // periodic_spring_force->SetMeinekeSpringStiffness(spring_strength);
        // if(include_springs)
        // {
        //     simulator.AddForce(periodic_spring_force);
        // }

	    // Create force law
	    MAKE_PTR(GeneralisedLinearSpringForce<2>, linear_force);
	    linear_force->SetCutOffLength(1.5);
        linear_force->SetMeinekeSpringStiffness(spring_strength);
        if(add_springs)
        {
            simulator.AddForce(linear_force);
        }

        MAKE_PTR(BendingForce2DGhosts<2>, bending_force);
        bending_force->SetBendingCoefficient(beta_parameter);
        bending_force->SetExponentParameter(alpha_parameter);
        bending_force->SetDomainWidth(periodic_width);
        bending_force->SetCircularNonZeroTargetCurvatureRegion(true, target_curvature, radius, centre_x);
        simulator.AddForce(bending_force);

        

		// Create periodic basement membrane force law
        // MAKE_PTR(PeriodicBendingForce3dHeightWithGhostNodes, periodic_bending_force);
        // periodic_bending_force->SetHeightDependantCurvatureParameter(1.2);
        // periodic_bending_force->SetBasementMembraneParameter(beta_parameter);
        // periodic_bending_force->SetExponentParameter(alpha_parameter);
        // // periodic_bending_force->SetCircularNonZeroTargetCurvatureRegion(true, target_curvature, radius, 20, 20);
        // periodic_bending_force->SetCircularNonZeroTargetCurvatureRegion(true, target_curvature, radius, centre_x, centre_y);
        // periodic_bending_force->SetPeriodicDomainWidth(periodic_width);
        // periodic_bending_force->SetPeriodicDomainDepth(periodic_height);
        // if(include_bending)
        // {
        //     simulator.AddForce(periodic_bending_force);
        // }


        // MAKE_PTR(DriftPreventForce<3>, p_drift_force);
        // p_drift_force->SetTissueMiddle(tissue_middle);
        // simulator.AddForce(p_drift_force);

        // Prevents getting stuck in a local minimums -> used to help break symmetry in cell anoikus
        // MAKE_PTR(RandomMotionForce<3>, p_random_force);
        // p_random_force->SetMovementParameter(0.001); //0.1 causes dissasociation, 0.001 is not enough
        // simulator.AddForce(p_random_force);

        


 
        // Add cell sloughing
        // MAKE_PTR_ARGS(SloughingCellKiller3DWithGhostNodes, sloughing, (&cell_population, periodic_width, periodic_height));
        // simulator.AddCellKiller(sloughing);

        double cut_off = 2.0;
        // Add anoikis cell killer
        // MAKE_PTR_ARGS(AnoikisCellKiller3DWithGhostNodes, anoikis, (&cell_population, cut_off, periodic_width, periodic_height));
        // simulator.AddCellKiller(anoikis);

        // std::string output_directory = "Test_WithSpringsBending_randz_alpha_2";
        std::string output_directory = "Test_WithSpringsBending_slave";
        simulator.SetOutputDirectory(output_directory);	 

        // Add random cell killer for death at the edges
        //                                                              ProbabilityOfDeathInAnHour, MinXBoundary, MaxXBoundary, MinYBoundary, MaxYBoundary
        // MAKE_PTR_ARGS(UniformCellKiller3dWithGhostNodes, random_cell_death, (&cell_population, 0.01, 1.5*width_space,  periodic_width-1.5*width_space, 1.5*height_space,  periodic_height-1.5*height_space, num_epithelial_cells+num_tissue_cells));
		// simulator.AddCellKiller(random_cell_death);

    
    //    for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
    //     {
    //         std::set<unsigned> neighbouring_node_indices = cell_population.GetNeighbouringNodeIndices(i);

    //         for (std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
    //          iter != neighbouring_node_indices.end();
    //          ++iter)
    //         {
    //             unsigned neighbour_node_global_index = *iter;
    //             std::cout<< neighbour_node_global_index  << " ";
    //         }
    //         std::cout<< "\n";

    //     }

        // //Test Forces
        // for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        // {
        //     cell_population.GetNode(i)->ClearAppliedForce();
        // }

        // SimulationTime::Destroy();
        // SimulationTime::Instance()->SetStartTime(0.0);
        // SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(0.001, 1);

        // periodic_bending_force->AddForceContribution(cell_population);
        // //periodic_force->AddForceContribution(cell_population);

        // SimulationTime::Destroy();

        // for (unsigned i=0; i<num_epithelial_cells+num_tissue_cells; i++)
        // {
        //     std::cout<< "Cell i = " << i << "\n";
        //     PRINT_VECTOR(cell_population.GetNode(i)->rGetAppliedForce());
        // }

        // ///////////////////
        // for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        // {
        //     cell_population.GetNode(i)->ClearAppliedForce();
        // }

        // SimulationTime::Destroy();
        // SimulationTime::Instance()->SetStartTime(0.0);
        // SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(0.001, 1);

        // periodic_spring_force->AddForceContribution(cell_population);
        // //periodic_force->AddForceContribution(cell_population);

        // SimulationTime::Destroy();

        // for (unsigned i=0; i<num_epithelial_cells+num_tissue_cells; i++)
        // {
        //     std::cout<< "Cell i = " << i << "\n";
        //     PRINT_VECTOR(cell_population.GetNode(i)->rGetAppliedForce());
        // }
        
               
    	
        simulator.SetSamplingTimestepMultiple(1);			// Every hour
		simulator.SetEndTime(end_time);
        simulator.SetDt(time_step);

        auto t1 = std::chrono::high_resolution_clock::now();
        simulator.Solve();
        auto t2 = std::chrono::high_resolution_clock::now();
        auto duration = pow(10.0,-6)*(std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count());

        std::cout << "\nTime taken = " << duration << " seconds\n\n";
        
        
        // for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        // {
        //     cell_population.GetNode(i)->ClearAppliedForce();
        //     //node_forces.push_back(zero_vector<double>(3));
        // }
        // // periodic_force->AddForceContribution(cell_population);
        // periodic_spring_force->AddForceContribution(cell_population);

        // for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        // {
        //     if( ! (cell_population.IsGhostNode(i)) )
        //     {
        //         PRINT_VARIABLE(i);
        //         PRINT_VECTOR(cell_population.GetNode(i)->rGetAppliedForce());
        //     }
        // }

        
    }     

};

#endif /*TEST3DBOXMODEL_HPP_*/

