#ifndef TEST3DBOXMODEL_HPP_
#define TEST3DBOXMODEL_HPP_

#include <cxxtest/TestSuite.h>
#include <chrono>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "TrianglesMeshReader.hpp"
#include "OffLatticeSimulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
//#include "FixedDurationGenerationBasedCellCycleModel.hpp"
//#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "DomMeshBasedCellPopulationWithGhostNodes.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellsGenerator.hpp"
#include "CryptModelInteractionForce.hpp"
#include "PeriodicCryptModelInteractionForceWithGhostNodes.hpp"
//#include "BasementMembraneForce3d.hpp"
#include "PeriodicBendingForce3dWithGhostNodes.hpp"
#include "CryptSimulation3dGhosts.hpp"
#include "SloughingCellKiller3DWithGhostNodes.hpp"
#include "AnoikisCellKiller3DWithGhostNodes.hpp"
#include "RandomRadialCellKiller3d.hpp"
#include "PeriodicBoxBoundaryCondition3d.hpp"
#include "PeriodicStromalBoxBoundaryCondition3d.hpp"
#include "TrianglesMeshWriter.hpp"
#include "Debug.hpp"
//#include "FixedBoundaries3d.hpp"

#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "CellIdWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "VoronoiDataWriter.hpp"

#include "PetscSetupAndFinalize.hpp"

// Tests have been copied directly from TestOffLatticeSimulation3d.hpp - most likely that bits will get changed
// and broken by me along the way

class Test3dBoxModel : public AbstractCellBasedTestSuite
{
private:

    MutableMesh<3,3>* Make3dMesh(unsigned width=3, unsigned height=3, unsigned depth=3)
    {
		/*   	   _ _ _ _ _
		 *        /        /|
		 *       /        / |
		 * 	    /_ _ _ _ /  | depth (z-direction)
		 * 	   |         |  |
		 *     |         |  |
		 *     |         |  /
		 *     |         | / height (y-direction)
		 * 	   |_ _ _ _ _|/
		 *        width
		 *    (x-direction)
		 *
		 *    width, depth and height correspond to number of elements in each direction
		 *
		 */		
		
        MutableMesh<3,3>* p_mesh = new MutableMesh<3,3>;
        p_mesh->ConstructCuboid(width-1, height-1, depth-1);
        TrianglesMeshWriter<3,3> mesh_writer("","3dSpringMesh");
        mesh_writer.WriteFilesUsingMesh(*p_mesh);

        return p_mesh;
    }

//    double mLastStartTime;
//    void setUp()
//    {
//        mLastStartTime = std::clock();
//        AbstractCellBasedTestSuite::setUp();
//    }
//    void tearDown()
//    {
//        double time = std::clock();
//        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
//        std::cout << "Elapsed time: " << elapsed_time << std::endl;
//        AbstractCellBasedTestSuite::tearDown();
//    }

public:
    
    /* A cube that consists of a block of stromal cells, and a single layer of epithelial cells.
     * Target curvature for the basement membrane is zero, and multiple division events occur.
     *
     */
    void TestPeriodicCubeWithGhosts() throw (Exception)
    {
    	RandomNumberGenerator::Instance()->Reseed(1);

        std::vector<Node<3>*> nodes;

        unsigned width = 10;	   // x
        unsigned height = 10;      // y
        unsigned ghosts_bottom = 0;       // ghosts > depth
        unsigned ghosts_top = 2;       // ghosts > depth
        unsigned num_tissue_depth = 6;
        unsigned num_tissue_beneith_crypt = 2; // If num_tissue_depth < 4 , use 0

        unsigned depth = ghosts_bottom + num_tissue_depth + 1.0 + ghosts_top ;        // +1 for epithelial cells -> which we always only have a monolayer

        // Initialise the tissue in an equilibrum state
        double width_space = 1.0;
        double height_space = 1.0;
        double depth_space = 1.0; //Magic number for z-spaceing... 
        unsigned cells_per_layer = width*height;
        unsigned cell_iter = 0;

        double periodic_width = (double) (width+0.0)*width_space;
        double periodic_height = (double) (height+0.0)*height_space;

        double off_set = 2.0;

        bool isGhost = false;
        double x_coordinate, y_coordinate, z_coordinate;

        // Get the dimensions for the non-zero target curvature region
        double centre_x = double (periodic_width)*0.5 - 0.5;
        double centre_y = double (periodic_height)*0.5;
        // double radius =  2.0;//periodic_width+1.0;
        double radius =  0;//periodic_width+1.0;
        double target_curvature = 0.15; //maximum curvature is 0.2066 -> higher curvature means smaller sphere
        double beta_parameter = 1.0;
        double alpha_parameter = 2.0;

        double epithelial_height = depth_space*num_tissue_depth + off_set;
        // double epithelial_radius = 4.25;
        // double ghost_radius = 3.4;
        double epithelial_radius = 2.5;
        double ghost_radius = 1.5;


        c_vector<unsigned,3> is_epithelial_cell;
        std::vector<c_vector<unsigned,3> > epithelial_cells_indicator;

        std::vector<unsigned> ghost_node_indices, real_node_indices;
        for (unsigned k=0; k<depth; k++)
        {
            // isGhost = false;
            // if(k < ghosts_bottom || k > ghosts_bottom + num_tissue_depth)
            // {
            //     isGhost = true;
            // }
            for (unsigned j=0; j<height; j++)
            {    
                for (unsigned i=0; i<width; i++)
                {
                    isGhost = false;
                    if(k < ghosts_bottom || k > ghosts_bottom + num_tissue_depth)
                    {
                        isGhost = true;
                    }

                    is_epithelial_cell[0] = cell_iter;
                    is_epithelial_cell[1] = 0;
                    is_epithelial_cell[2] = 0;

                    c_vector<double, 3> node_i_new_location;

                    //x_coordinate = (double) (i + 0.5*((j%2 + k%2)%2))*width_space;
                    x_coordinate = (double) i*width_space;
                    y_coordinate = (double) j*height_space;
                    z_coordinate = (double) off_set + k*depth_space;
                    
                    nodes.push_back(new Node<3>(cell_iter,  false,  x_coordinate, y_coordinate, z_coordinate));
                    if(isGhost == false)
                    {
                        // The bottom stromal layer
                        if (k == (ghosts_bottom + num_tissue_beneith_crypt) - 1.0)//(z_coordinate == off_set + depth_space*(ghosts_bottom + 1.0))
                        {
                            is_epithelial_cell[1] = 0;
                        }
                        // The base of the crypt
                        else if ( (k == (ghosts_bottom + num_tissue_beneith_crypt + 1.0)  ) && 
                                ( x_coordinate >= centre_x-ghost_radius &&  x_coordinate <= centre_x+ghost_radius ) &&
                                ( y_coordinate >= centre_y-ghost_radius &&  y_coordinate <= centre_y+ghost_radius ) )
                        {
                            is_epithelial_cell[1] = 1;
                            is_epithelial_cell[2] = 1;
                        }
                        // The crypt walls
                        else if ( (k > (ghosts_bottom + num_tissue_beneith_crypt + 1.0) ) && 
                                (k <= (ghosts_bottom + num_tissue_depth )) && 
                                ( x_coordinate >= centre_x-epithelial_radius &&  x_coordinate <= centre_x+epithelial_radius ) &&
                                ( y_coordinate >= centre_y-epithelial_radius &&  y_coordinate <= centre_y+epithelial_radius ) )
                        {
                            if (( x_coordinate >= centre_x-ghost_radius &&  x_coordinate <= centre_x+ghost_radius ) &&
                                ( y_coordinate >= centre_y-ghost_radius &&  y_coordinate <= centre_y+ghost_radius ) )
                            {
                                isGhost = true;
                            }
                            else
                            {
                                is_epithelial_cell[1] = 1;
                                is_epithelial_cell[2] = 1;
                            }
                            
                        }
                        // The epithelial layer above the stromal
                        else if (k == (ghosts_bottom + num_tissue_depth))
                        {
                            is_epithelial_cell[1] = 1;
                            is_epithelial_cell[2] = 0;
                        }
                        
                    }

                    if (isGhost == true)
                    {
                        ghost_node_indices.push_back(cell_iter);
                    }
                    else if(isGhost == false)
                    {
                        real_node_indices.push_back(cell_iter);
                        epithelial_cells_indicator.push_back(is_epithelial_cell);
                    }
                    

                    // Check if the cell is in the domain
                    

                    cell_iter++;
                }
            }
        }

        MutableMesh<3,3> mesh(nodes);
        // mesh.Translate(0.001, 0.001);


        std::vector<CellPtr> cells;

    	// First we sort the real node indices into increasing order (the last ones will correspond to the
    	// epithelial nodes)
    	sort(real_node_indices.begin(), real_node_indices.end());

        // Set up cells by iterating through the mesh nodes
        // unsigned num_nodes = mesh.GetNumAllNodes();
        // unsigned num_epithelial_cells = (width)*(height);
        // unsigned num_tissue_cells = (width)*(height)*(num_tissue_depth);
        // unsigned num_ghosts_bottom = (width)*(height)*(ghosts_bottom);
        // unsigned num_ghosts_top = (width)*(height)*(ghosts_top);

        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_stromal_state(new StromalCellMutationState);

        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        /*
        // Initialise Tissue cells (Stromal)
		for (unsigned i=0; i<real_node_indices.size()-num_epithelial_cells; i++)
		{
			//StochasticDurationGenerationBasedCellCycleModel* p_model = new StochasticDurationGenerationBasedCellCycleModel();
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();

            CellPtr p_differentiated_cell(new Cell(p_stromal_state, p_model));
            p_differentiated_cell->SetCellProliferativeType(p_differentiated_type);
            p_model->SetDimension(3);

			cells.push_back(p_differentiated_cell);
        }


        // Initialise Epithelial cells
        for (unsigned i=real_node_indices.size()-num_epithelial_cells; i<real_node_indices.size(); i++)
        {

        	FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_epithelial_cell(new Cell(p_state, p_model));
            p_epithelial_cell->SetCellProliferativeType(p_transit_type);
            p_model->SetDimension(3);
            
            //p_model->SetTransitCellG1Duration(5);
            //double birth_time = -RandomNumberGenerator::Instance()->ranf()*(p_model->GetTransitCellG1Duration() + p_model->GetSG2MDuration());
            
            p_model->SetMaxTransitGenerations(1);
            p_model->SetTransitCellG1Duration(100);       
            double birth_time = 0.0;
            //std::cout << "cell i = " << i << " birth_time = " << birth_time << "\n";
			
            p_epithelial_cell->SetBirthTime(birth_time);
            
            p_epithelial_cell->InitialiseCellCycleModel();
            
			cells.push_back(p_epithelial_cell);
        }
        */

        for(unsigned i=0; i<real_node_indices.size(); i++)
        {
            if (epithelial_cells_indicator[i][1] == 0)
            {
                //StochasticDurationGenerationBasedCellCycleModel* p_model = new StochasticDurationGenerationBasedCellCycleModel();
                FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();

                CellPtr p_differentiated_cell(new Cell(p_stromal_state, p_model));
                p_differentiated_cell->SetCellProliferativeType(p_differentiated_type);
                p_model->SetDimension(3);

                cells.push_back(p_differentiated_cell);
            }
            else if (epithelial_cells_indicator[i][1] == 1)
            {
                FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
                CellPtr p_epithelial_cell(new Cell(p_state, p_model));
                // if (epithelial_cells_indicator[i][2] == 1)
                // {
                //     p_epithelial_cell->SetCellProliferativeType(p_transit_type);
                //     // p_model->SetMaxTransitGenerations(1);
                //     p_model->SetTransitCellG1Duration(3);    
                // }
                // else if (epithelial_cells_indicator[i][2] == 0)
                // {
                //     p_epithelial_cell->SetCellProliferativeType(p_differentiated_type);
                // }
                // Set all epithelial as differentiated
                p_epithelial_cell->SetCellProliferativeType(p_differentiated_type);
                p_model->SetDimension(3);
                
                //p_model->SetTransitCellG1Duration(5);
                double birth_time = -RandomNumberGenerator::Instance()->ranf()*(p_model->GetTransitCellG1Duration() + p_model->GetSG2MDuration());
                   
                // double birth_time = 0.0;
                //std::cout << "cell i = " << i << " birth_time = " << birth_time << "\n";
                
                p_epithelial_cell->SetBirthTime(birth_time);
                
                p_epithelial_cell->InitialiseCellCycleModel();
                
                cells.push_back(p_epithelial_cell);

                // std::cout << p_epithelial_cell->GetCellProliferativeType()->GetColour() << "\n";
            }
        }

        DomMeshBasedCellPopulationWithGhostNodes<3> cell_population(mesh, cells, real_node_indices);
        //MeshBasedCellPopulationWithGhostNodes<3> cell_population(mesh, cells, real_node_indices);
//        MeshBasedCellPopulation<3> cell_population(mesh, cells);
        assert(cell_population.GetNumRealCells() != 0);


        //CryptSimulation3dGhosts simulator(cell_population, false, true);
        //CryptSimulation3d(rCellPopulation, bool deleteCellPopulationAndForceCollection, bool initialiseCells)
        OffLatticeSimulation<3> simulator(cell_population);
        std::string output_directory = "TestSimulation_MyBendingForce";
        simulator.SetOutputDirectory(output_directory);
//        cell_population.InitialiseCells();

        // Make sure we have a Voronoi tessellation to begin with
        //cell_population.CreateVoronoiTessellation();
        
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        // cell_population.AddCellWriter<CellAncestorWriter>();
        // cell_population.AddPopulationWriter<NodeVelocityWriter>();
        //cell_population.AddPopulationWriter<VoronoiDataWriter>(); // paraview is pretty pointless at viewing this, worth looking into
        
        //cell_population.WriteVtkResultsToFile(output_directory);
        //To fix paraview
        cell_population.SetWriteVtkAsPointsDom(true);
        //std::cout<<cell_population.GetWriteVtkAsPoints() << "\n";
        //PRINT_VARIABLE(cell_population.GetWriteVtkAsPoints());
        //cell_population.SetOutputMeshInVtk(true);



        std::map<Node<3>*, c_vector<double,3> > node_locations_before;
        for (std::list<CellPtr>::iterator cell_iter = cell_population.rGetCells().begin();
             cell_iter != cell_population.rGetCells().end();
             ++cell_iter)
        {
            Node<3>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
            node_locations_before[p_node] = p_node->rGetLocation();
            //PRINT_VECTOR(node_locations_before[p_node]);
        }

        

        MAKE_PTR_ARGS(PeriodicBoxBoundaryCondition3d, boundary_condition, (&cell_population));
        boundary_condition->SetCellPopulationWidth(periodic_width);
        boundary_condition->SetCellPopulationDepth(periodic_height);
        boundary_condition->SetMaxHeightForPinnedCells(0.0);       			      
        boundary_condition->ImposeBoundaryCondition(node_locations_before);
        simulator.AddCellPopulationBoundaryCondition(boundary_condition);

        MAKE_PTR_ARGS(PeriodicStromalBoxBoundaryCondition3d, stromal_boundary_condition, (&cell_population));
        stromal_boundary_condition->SetCellPopulationWidth(periodic_width);
        stromal_boundary_condition->SetCellPopulationDepth(periodic_height);
        stromal_boundary_condition->SetMaxHeightForPinnedCells(-5.0);
        stromal_boundary_condition->ImposeBoundaryCondition(node_locations_before);
        simulator.AddCellPopulationBoundaryCondition(stromal_boundary_condition);


		// Create periodic spring force law
        MAKE_PTR(PeriodicCryptModelInteractionForceWithGhostNodes<3>, periodic_spring_force);
     // PeriodicCryptModelInteractionForce<3> periodic_spring_force;  // Variable spring strengths
        periodic_spring_force->SetUseOneWaySprings(false);
        periodic_spring_force->SetCutOffLength(1.5);
        //                     SetEpithelialStromalCellDependentSprings(ind , Ep-Ep, Str-Str, Ep-Str, apcTwoHitStromalMultiplier);
        periodic_spring_force->SetEpithelialStromalCellDependentSprings(true, 1.0,     0.01,     1.0,    1.0);
        periodic_spring_force->SetPeriodicDomainWidth(periodic_width);
        periodic_spring_force->SetPeriodicDomainDepth(periodic_height);
        periodic_spring_force->SetMeinekeSpringStiffness(10.0);
        simulator.AddForce(periodic_spring_force);


	    // Create force law
	    //MAKE_PTR(GeneralisedLinearSpringForce<3>, linear_force);
	    //linear_force->SetCutOffLength(1.5);
        //simulator.AddForce(linear_force);

		// Create periodic basement membrane force law
        MAKE_PTR(PeriodicBendingForce3dWithGhostNodes, periodic_bending_force);
     // PeriodicBasementMembraneForce3d periodic_basement_membrane;
        periodic_bending_force->SetBasementMembraneParameter(beta_parameter);
        periodic_bending_force->SetExponentParameter(alpha_parameter);
        //periodic_bending_force->SetCircularNonZeroTargetCurvatureRegion(true, 0.15, radius, centre_x, centre_y);
        periodic_bending_force->SetCircularNonZeroTargetCurvatureRegion(true, target_curvature, radius, centre_x, centre_y);
        periodic_bending_force->SetPeriodicDomainWidth(periodic_width);
        periodic_bending_force->SetPeriodicDomainDepth(periodic_height);
        simulator.AddForce(periodic_bending_force);

 
        // Add cell sloughing -> dont need if using periodic conditions...
        // MAKE_PTR_ARGS(SloughingCellKiller3DWithGhostNodes, sloughing, (&cell_population, periodic_width, periodic_height));
        // simulator.AddCellKiller(sloughing);

        double cut_off = 2.5;
        // Add anoikis cell killer
        MAKE_PTR_ARGS(AnoikisCellKiller3DWithGhostNodes, anoikis, (&cell_population, cut_off, periodic_width, periodic_height));
        simulator.AddCellKiller(anoikis);

        //Test Forces
        
        // for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        // {
        //     cell_population.GetNode(i)->ClearAppliedForce();
        //     //node_forces.push_back(zero_vector<double>(3));
        // }

        // SimulationTime::Destroy();
        // SimulationTime::Instance()->SetStartTime(0.0);
        // SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(0.001, 1);

        // periodic_bending_force->AddForceContribution(cell_population);
        // //periodic_force->AddForceContribution(cell_population);

        // SimulationTime::Destroy();

        // for (unsigned i=0; i<real_node_indices.size(); i++)
        // {
        //     //std::cout<< "Cell i = " << i << "\n";
        //     PRINT_VECTOR(cell_population.GetNode(i)->rGetAppliedForce());
        // }
        
        
        // Run for a short time to allow it to deform		        
    	simulator.SetSamplingTimestepMultiple(10);			// Every hour
		simulator.SetEndTime(1);
        simulator.SetDt(0.001);

        auto t1 = std::chrono::high_resolution_clock::now();
        simulator.Solve();
        auto t2 = std::chrono::high_resolution_clock::now();
        auto duration = pow(10.0,-6)*(std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count());
        
        std::cout << "\nTime taken = " << duration << " seconds\n\n";
        
    }     
     
           
};

#endif /*TEST3DBOXMODEL_HPP_*/

