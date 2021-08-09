#ifndef TestPeriodicBendingForce_FlatWntRadius_HPP_
#define TestPeriodicBendingForce_FlatWntRadius_HPP_

#include <cxxtest/TestSuite.h>
#include "Timer.hpp"

// Must be included before other cell_based headers
// #include "CellBasedSimulationArchiver.hpp"

// #include "TrianglesMeshReader.hpp"
#include "OffLatticeSimulation.hpp"
// #include "GeneralisedLinearSpringForce.hpp"
#include "DomMeshBasedCellPopulationWithGhostNodes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StromalCellMutationState.hpp"
// #include "CellsGenerator.hpp"

#include "AbstractTwoBodyInteractionForce.hpp"


// #include "RandomMotionForce.hpp"
#include "PeriodicCryptModelInteractionForceWithGhostNodes.hpp"
#include "PeriodicBendingForce3dHeightWithGhostNodes.hpp"
// #include "SloughingCellKiller3DWithGhostNodes.hpp"
#include "AnoikisCellKiller3DWithGhostNodes.hpp"
// #include "UniformCellKiller3dWithGhostNodes.hpp"
#include "PeriodicBoxBoundaryCondition3dGhosts.hpp"
#include "PeriodicStromalBoxBoundaryCondition3d.hpp"
// #include "TrianglesMeshWriter.hpp"
#include "Debug.hpp"
#include "DensityDependantCellKillerRadiual3DWithGhostNodesV2.hpp"

#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "CellIdWriter.hpp"
// #include "CellMutationStatesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellAncestor.hpp"
// #include "CellPopulationEpithelialWriter.hpp"
// #include "VoronoiDataWriter.hpp"
#include "NodeVelocityWriter.hpp"

#include "PetscSetupAndFinalize.hpp"
// #include "ForwardEulerNumericalMethod.hpp"
#include "PlanarDivisionRule.hpp"
#include "DriftPreventForce.hpp"

// #include "MeshModifier.hpp"
// #include "MeshRemeshModifier.hpp"
// #include "PeriodicRemeshCellsModifier.hpp"
#include "MeshRemeshCutoffModifier.hpp"

#include "DomSimpleWntCellCycleModel.hpp"

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
        EXIT_IF_PARALLEL;

    	// RandomNumberGenerator::Instance()->Reseed(1);

        std::vector<Node<3>*> nodes;

        std::string output_directory = "FlatWntRadius_20210804_03";

        unsigned width = 8;	   // x
        unsigned height = 10;      // y
        unsigned ghosts_bottom = 0;       // ghosts > depth
        unsigned ghosts_top = 2;       // ghosts > depth
        unsigned num_tissue_depth = 1;
        unsigned depth = num_tissue_depth + (ghosts_bottom + ghosts_top);        // z

        // Initialise the tissue in an equilibrum state
        double width_space = 1.0;
        double height_space = 1.0*sqrt(0.75);
        double ghost_sep = 1.0;
        double depth_space = 0.738431690356779*1.0; //Magic number for z-spaceing... 
        unsigned cells_per_layer = width*height;
        unsigned cell_iter = 0;

        double periodic_width = (double) (width+0.0)*width_space;
        double periodic_height = (double) (height+0.0)*height_space;

        double tissue_base = 5.0; //Hieght of tissue to prevent drift
        double tissue_middle = 0.0;

        bool isGhost = false;
        double x_coordinate, y_coordinate, z_coordinate;

        // Get the dimensions for the non-zero target curvature region
        // double centre_x = periodic_width*0.5;
        // double centre_y = periodic_height*0.5;
        double centre_x = 30;
        double centre_y = 30;

        double spring_strength = 20.0;
        double spring_cuttoff = 1.75;

        double radius =  2.0;//periodic_width+1.0;
        double target_curvature = -0.2; //maximum curvature is 0.2066 -> higher curvature means smaller sphere
        double beta_parameter = 1.0*spring_strength;
        double alpha_parameter = 1.2;

        double time_step = 0.001;
        double end_time = 120;
        double plot_step = 10;

        bool include_springs = true;
        bool include_bending = true;

        int is_transit[depth*height*width];
        int num_real_nodes = 0;

        std::vector<unsigned> ghost_node_indices, real_node_indices;
        for (unsigned k=0; k<=(depth); k++) //+1 puts a layer of ghosts on the bottom
        {
            isGhost = true;
            if(k < num_tissue_depth + 1)
            {
                isGhost = false;
            }
            // if(k == depth)
            // {
            //     isGhost = true;
            // }

            for (unsigned j=0; j<height; j++)
            {    
                for (unsigned i=0; i<width; i++)
                {
                    is_transit[cell_iter] = 0;
                        
                    c_vector<double, 3> node_i_new_location;

                    x_coordinate = (double) (i + 0.5*(j%2 + k%2))*width_space ;
                    y_coordinate = (double) j*height_space ;

                    // z_coordinate = (double) tissue_base + k*depth_space;
                    
                    if( k == depth)
                    {
                        z_coordinate = (double) tissue_base + (-1.0)*depth_space;
                    }
                    else
                    {
                        z_coordinate = (double) tissue_base + k*depth_space;
                    }

                    if( pow(x_coordinate - 0.5*periodic_width,2)+ pow(y_coordinate - 0.5*periodic_height ,2) <= pow(1.0,2) )
                    {
                        is_transit[cell_iter] = 1;
                    }
                    if(isGhost)
                    {
                        ghost_node_indices.push_back(cell_iter);
                    }
                    else
                    {
                        num_real_nodes++;
                        tissue_middle = tissue_middle + z_coordinate;
                        real_node_indices.push_back(cell_iter);
                    }
                    nodes.push_back(new Node<3>(cell_iter,  false,  x_coordinate, y_coordinate, z_coordinate));
                    
                    // PRINT_4_VARIABLES(cell_iter, x_coordinate,y_coordinate,z_coordinate);
                    // std::cout << z_coordinate << " ";
                    cell_iter++;

                }
            }
        }
        tissue_middle = tissue_middle/num_real_nodes;

        MutableMesh<3,3> mesh(nodes);

        std::vector<CellPtr> cells;

    	// First we sort the real node indices into increasing order (the last ones will correspond to the
    	// epithelial nodes)
    	sort(real_node_indices.begin(), real_node_indices.end());

        // Set up cells by iterating through the mesh nodes
        // unsigned num_nodes = mesh.GetNumAllNodes();
        unsigned num_epithelial_cells = (width)*(height);
        unsigned num_tissue_cells = (width)*(height)*(num_tissue_depth);
        unsigned num_ghosts_bottom = (width)*(height)*(ghosts_bottom);
        unsigned num_ghosts_top = (width)*(height)*(ghosts_top);

        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_stromal_state(new StromalCellMutationState);

        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        // Initialise Tissue cells (Stromal)
		for (unsigned i=0; i<real_node_indices.size()-num_epithelial_cells; i++)
		{
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();

            CellPtr p_differentiated_cell(new Cell(p_stromal_state, p_model));
            p_differentiated_cell->SetCellProliferativeType(p_differentiated_type);
            p_model->SetDimension(3);

            double birth_time = - 10.0;
            p_differentiated_cell->SetBirthTime(birth_time);

			cells.push_back(p_differentiated_cell);
        }


        // Initialise Epithelial cells
        for (unsigned i=real_node_indices.size()-num_epithelial_cells; i<real_node_indices.size(); i++)
        {

        	// FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            // p_model->SetMaxTransitGenerations(1);
            // p_model->SetSDuration(2); 
            // p_model->SetG2Duration(2); 
            // p_model->SetMDuration(2); 
            // p_model->SetDimension(3);
            
            DomSimpleWntCellCycleModel* p_model = new DomSimpleWntCellCycleModel();
            
            p_model->SetDimension(3);
            // p_model->SetMaxTransitGenerations(100);

            p_model->SetTransitCellG1Duration(4);
            p_model->SetSDuration(2);
            p_model->SetG2Duration(1);
            p_model->SetMDuration(1);

            // p_model->SetTransitCellG1Duration(11);
            // p_model->SetSDuration(8);
            // p_model->SetG2Duration(4);
            // p_model->SetMDuration(1);

            // p_model->SetTransitCellG1Duration(0.5);
            // p_model->SetSDuration(0.25);
            // p_model->SetG2Duration(0.25);
            // p_model->SetMDuration(0.5);
            
            CellPtr p_epithelial_cell(new Cell(p_state, p_model));
            // double birth_time = -8;
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*(p_model->GetTransitCellG1Duration() + p_model->GetSG2MDuration());

            p_epithelial_cell->SetCellProliferativeType(p_transit_type);
            // p_epithelial_cell->SetCellProliferativeType(p_differentiated_type);
			
            p_epithelial_cell->SetBirthTime(birth_time);
            p_epithelial_cell->SetApoptosisTime(1);
            p_epithelial_cell->InitialiseCellCycleModel();

            MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (i));
            p_epithelial_cell->SetAncestor(p_cell_ancestor);
            
			cells.push_back(p_epithelial_cell);
        }

        // std::cout<< "number of cells comp = " << real_node_indices.size() << "\n";
        // std::cout<< "number of cells all  = " << num_epithelial_cells+num_tissue_cells << "\n";
        // std::cout<< "number of ghosts     = " << ghost_node_indices.size() << "\n";

        DomMeshBasedCellPopulationWithGhostNodes<3> cell_population(mesh, cells, real_node_indices); //ghost_sep
        // assert(cell_population.GetNumRealCells() != 0);
        
        // Set the division rule for our population to be the random direction division rule
        boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_division_rule_to_set(new PlanarDivisionRule());
        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

        OffLatticeSimulation<3> simulator(cell_population);
        
        DomWntConcentration<3>::Instance()->SetType(DomRADIAL);
        DomWntConcentration<3>::Instance()->SetCellPopulation(cell_population);
        DomWntConcentration<3>::Instance()->SetCryptLength(20.0);
        DomWntConcentration<3>::Instance()->SetCryptCentreX(0.5*periodic_width);
        DomWntConcentration<3>::Instance()->SetCryptCentreY(0.5*periodic_height);
        DomWntConcentration<3>::Instance()->SetCryptRadius(1.5);
        DomWntConcentration<3>::Instance()->SetWntConcentrationParameter(2.0);

        // Pass an adaptive numerical method to the simulation
        // boost::shared_ptr<AbstractNumericalMethod<3,3> > p_method(new ForwardEulerNumericalMethod<3,3>());
        // p_method->SetUseAdaptiveTimestep(false);
        // simulator.SetNumericalMethod(p_method);

        // cell_population.InitialiseCells();

        // Make sure we have a Voronoi tessellation to begin with
        // cell_population.CreateVoronoiTessellation();
        
        cell_population.AddCellWriter<CellIdWriter>();
        // cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();
        // cell_population.AddPopulationWriter<CellPopulationEpithelialWriter>();
        cell_population.AddPopulationWriter<NodeVelocityWriter>();

        //cell_population.AddPopulationWriter<VoronoiDataWriter>(); // paraview is pretty pointless at viewing this, worth looking into
        
        // cell_population.WriteVtkResultsToFile(output_directory);
        //To fix paraview
        cell_population.SetWriteVtkAsPointsDom(true);
        //std::cout<<cell_population.GetWriteVtkAsPoints() << "\n";
        // cell_population.SetOutputMeshInVtkDom(false);



        std::map<Node<3>*, c_vector<double,3> > node_locations_before;
        for (std::list<CellPtr>::iterator cell_iter = cell_population.rGetCells().begin();
             cell_iter != cell_population.rGetCells().end();
             ++cell_iter)
        {
            Node<3>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
            node_locations_before[p_node] = p_node->rGetLocation();
            //PRINT_VECTOR(node_locations_before[p_node]);
        }

        MAKE_PTR_ARGS(PeriodicBoxBoundaryCondition3dGhosts, boundary_condition, (&cell_population));
        boundary_condition->SetCellPopulationWidth(periodic_width);
        boundary_condition->SetCellPopulationDepth(periodic_height);
        boundary_condition->SetMaxHeightForPinnedCells(0.0);       			      
        boundary_condition->ImposeBoundaryCondition(node_locations_before);
        simulator.AddCellPopulationBoundaryCondition(boundary_condition);

        MAKE_PTR_ARGS(PeriodicStromalBoxBoundaryCondition3d, stromal_boundary_condition, (&cell_population));
        stromal_boundary_condition->SetCellPopulationWidth(periodic_width);
        stromal_boundary_condition->SetCellPopulationDepth(periodic_height);
        stromal_boundary_condition->SetMaxHeightForPinnedCells(0.0);
        stromal_boundary_condition->ImposeBoundaryCondition(node_locations_before);
        simulator.AddCellPopulationBoundaryCondition(stromal_boundary_condition);

		// Create periodic spring force law
        MAKE_PTR(PeriodicCryptModelInteractionForceWithGhostNodes<3>, periodic_spring_force);
        periodic_spring_force->SetUseOneWaySprings(false); //turning this on makes the stromal cells act as ghosts..
        periodic_spring_force->SetCutOffLength(spring_cuttoff);
        //                     SetEpithelialStromalCellDependentSprings(ind , Ep-Ep, Str-Str, Ep-Str, apcTwoHitStromalMultiplier);
        periodic_spring_force->SetEpithelialStromalCellDependentSprings(true, 1.0,     1.0,     0.5,    1.0);
        periodic_spring_force->SetPeriodicDomainWidth(periodic_width);
        periodic_spring_force->SetPeriodicDomainDepth(periodic_height);
        periodic_spring_force->SetMeinekeSpringStiffness(spring_strength);
        if(include_springs)
        {
            simulator.AddForce(periodic_spring_force);
        }

		// Create periodic basement membrane force law
        MAKE_PTR(PeriodicBendingForce3dHeightWithGhostNodes, periodic_bending_force);
        periodic_bending_force->SetOutputDirectory(output_directory);
        periodic_bending_force->SetHeightDependantCurvatureParameter(1.2);
        periodic_bending_force->SetBasementMembraneParameter(beta_parameter);
        periodic_bending_force->SetExponentParameter(alpha_parameter);
        // periodic_bending_force->SetCircularNonZeroTargetCurvatureRegion(true, target_curvature, radius, 20, 20);
        periodic_bending_force->SetCircularNonZeroTargetCurvatureRegion(true, target_curvature, radius, centre_x, centre_y);
        periodic_bending_force->SetPeriodicDomainWidth(periodic_width);
        periodic_bending_force->SetPeriodicDomainDepth(periodic_height);
        if(include_bending)
        {
            simulator.AddForce(periodic_bending_force);
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////
        // THIS NEEDS TO BE LOOKED AT!!!!
        // MAKE_PTR(DriftPreventForce, p_drift_force);
        // p_drift_force->SetTissueMiddle(tissue_middle);
        // p_drift_force->SetMovementParameter(time_step);
        // simulator.AddForce(p_drift_force);
        ////////////////////////////////////////////////////////////////////////////////////////////////

        // Prevents getting stuck in a local minimums -> used to help break symmetry in cell anoikus
        // MAKE_PTR(RandomMotionForce<3>, p_random_force);
        // p_random_force->SetMovementParameter(0.001); //0.1 causes dissasociation, 0.001 is not enough
        // simulator.AddForce(p_random_force);
        
        double cut_off = 1.5;
        double density_threshold = 0.95;
        double domain_tol = 3.0;
        // Add anoikis cell killer
        MAKE_PTR_ARGS(AnoikisCellKiller3DWithGhostNodes, anoikis, (&cell_population, cut_off, periodic_width, periodic_height));
        anoikis->SetOutputDirectory(output_directory);
        simulator.AddCellKiller(anoikis);

        MAKE_PTR_ARGS(DensityDependantCellKillerRadiual3DWithGhostNodesV2, density, (&cell_population, domain_tol, density_threshold, periodic_width, periodic_height));
        density->SetOutputDirectory(output_directory);
        simulator.AddCellKiller(density);

        simulator.SetOutputDirectory(output_directory);	 

        // MAKE_PTR(MeshRemeshModifier<3>, p_modifier);
        // p_modifier->SetOutputDirectory(output_directory + "/results_from_time_0");
        // p_modifier->SetWidth(periodic_width);
        // p_modifier->SetDepth(periodic_height);
        // simulator.AddSimulationModifier(p_modifier);

        // MAKE_PTR(PeriodicRemeshCellsModifier<3>, p_per_modifier);
        // p_per_modifier->SetOutputDirectory(output_directory + "/results_from_time_0");
        // p_per_modifier->SetWidth(periodic_width);
        // p_per_modifier->SetDepth(periodic_height);
        // simulator.AddSimulationModifier(p_per_modifier);

        // Add random cell killer for death at the edges
        //                                                              ProbabilityOfDeathInAnHour,    MinXBoundary,                MaxXBoundary,     MinYBoundary,                    MaxYBoundary
        // MAKE_PTR_ARGS(UniformCellKiller3dWithGhostNodes, random_cell_death, (&cell_population, 1.0, 1.5*width_space,  periodic_width-1.5*width_space, 1.5*height_space,  periodic_height-1.5*height_space, num_epithelial_cells+num_tissue_cells));
        // // MAKE_PTR_ARGS(UniformCellKiller3dWithGhostNodes, random_cell_death, (&cell_population, 1.0, periodic_width + 1,  0 - 1, periodic_height + 1,  0 - 1 , num_epithelial_cells+num_tissue_cells ));
		// simulator.AddCellKiller(random_cell_death);              
    	
        // Create force law
	    // MAKE_PTR(GeneralisedLinearSpringForce<3>, linear_force);
	    // linear_force->SetCutOffLength(1.5);
        // simulator.AddForce(linear_force);

        MAKE_PTR(MeshRemeshCutoffModifier<3>, pc_modifier);
        pc_modifier->SetOutputDirectory(output_directory + "/results_from_time_0");
        pc_modifier->SetWidth(periodic_width);
        pc_modifier->SetDepth(periodic_height);
        pc_modifier->SetCutoff(spring_cuttoff+0.001);

        simulator.AddSimulationModifier(pc_modifier);

        simulator.SetSamplingTimestepMultiple(plot_step);			// Every hour
		simulator.SetEndTime(end_time);
        simulator.SetDt(time_step);

        // auto t1 = std::chrono::high_resolution_clock::now();
        Timer::Reset();
        simulator.Solve();
        Timer::Print("Time Ellapsed");
        // auto t2 = std::chrono::high_resolution_clock::now();
        // auto duration = pow(10.0,-6)*(std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count());

        // std::cout << "\nTime taken = " << duration << " seconds\n\n";
        
    }

};

#endif /*TEST3DBOXMODEL_HPP_*/

