#ifndef TEST3DBOXMODEL_HPP_
#define TEST3DBOXMODEL_HPP_

#include <cxxtest/TestSuite.h>

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

public:

    /*
    * This checks a small region in the middle of the tissue to ensure the 
    * bending force is implemented on the patch of cells in an expected way.
    */
    void TestPeriodicCubeWithGhosts_SmallRegion() throw (Exception)
    {
    	RandomNumberGenerator::Instance()->Reseed(1);

        std::vector<Node<3>*> nodes;

        unsigned width = 10;	   // x
        unsigned height = 10;      // y
        unsigned ghosts_bottom = 0;       // ghosts > depth
        unsigned ghosts_top = 2;       // ghosts > depth
        unsigned num_tissue_depth = 1;
        unsigned depth = num_tissue_depth + (ghosts_bottom + ghosts_top) + 1;        // z

        // Initialise the tissue in an equilibrum state
        double width_space = 1.0;
        double height_space = sqrt(0.75);
        double depth_space = 0.738431690356779; //Magic number for z-spaceing... 
        unsigned cells_per_layer = width*height;
        unsigned cell_iter = 0;

        double periodic_width = (double) (width+0.0)*width_space + 0.00;
        double periodic_height = (double) (height+0.0)*height_space + 0.00;

        bool isGhost = false;
        double x_coordinate, y_coordinate, z_coordinate;

        std::vector<unsigned> ghost_node_indices, real_node_indices;
        for (unsigned k=0; k<depth; k++)
        {
            isGhost = false;
            if(k < ghosts_bottom || k > ghosts_bottom + num_tissue_depth)
            {
                isGhost = true;
            }
            for (unsigned j=0; j<height; j++)
            {    
                for (unsigned i=0; i<width; i++)
                {
                    c_vector<double, 3> node_i_new_location;

                    //x_coordinate = (double) (i + 0.5*((j%2 + k%2)%2))*width_space;
                    x_coordinate = (double) (i + 0.5*(j%2 + k%2))*width_space;
                    y_coordinate = (double) j*height_space;
                    z_coordinate = (double) 5.0 + k*depth_space;
                    
                    /*
                    if (cell_iter >= 100 && cell_iter < 200)
                    {
                        z_coordinate = (double) 5.0 + k*depth_space + 0.5*(2.0*RandomNumberGenerator::Instance()->ranf()-1.0);
                    }
                    
                    if (cell_iter == 154)
                    {
                        z_coordinate = (double) 5.0 + k*depth_space + 0.25;
                    }
                    */
                    nodes.push_back(new Node<3>(cell_iter,  false,  x_coordinate, y_coordinate, z_coordinate));
                    
                    if (isGhost == true)
                    {
                        ghost_node_indices.push_back(cell_iter);
                    }
                    else if(isGhost == false)
                    {
                        real_node_indices.push_back(cell_iter);
                    }

                    cell_iter++;
                }
            }
        }

        MutableMesh<3,3> mesh(nodes);
        //mesh.Translate(0.5, 0.5);

        std::vector<CellPtr> cells;

    	// First we sort the real node indices into increasing order (the last ones will correspond to the
    	// epithelial nodes)
    	sort(real_node_indices.begin(), real_node_indices.end());

        // Set up cells by iterating through the mesh nodes
        unsigned num_nodes = mesh.GetNumAllNodes();
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

			cells.push_back(p_differentiated_cell);
        }

        // Initialise Epithelial cells
        for (unsigned i=real_node_indices.size()-num_epithelial_cells; i<real_node_indices.size(); i++)
        {

        	FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_epithelial_cell(new Cell(p_state, p_model));
            p_epithelial_cell->SetCellProliferativeType(p_transit_type);
            p_model->SetDimension(3);
            
            p_model->SetMaxTransitGenerations(1);
            p_model->SetTransitCellG1Duration(1);
            double birth_time = 0.0;
			
            p_epithelial_cell->SetBirthTime(birth_time);
            
            p_epithelial_cell->InitialiseCellCycleModel();
            
			cells.push_back(p_epithelial_cell);
        }

        DomMeshBasedCellPopulationWithGhostNodes<3> cell_population(mesh, cells, real_node_indices);
        //MeshBasedCellPopulationWithGhostNodes<3> cell_population(mesh, cells, real_node_indices);
        //  MeshBasedCellPopulation<3> cell_population(mesh, cells);
        assert(cell_population.GetNumRealCells() != 0);


        //CryptSimulation3dGhosts simulator(cell_population, false, true);
        //CryptSimulation3d(rCellPopulation, bool deleteCellPopulationAndForceCollection, bool initialiseCells)
        OffLatticeSimulation<3> simulator(cell_population);
        std::string output_directory = "TestBendingForceCheck";
        simulator.SetOutputDirectory(output_directory);
        // ell_population.InitialiseCells();

        // Make sure we have a Voronoi tessellation to begin with
        //cell_population.CreateVoronoiTessellation();      
        
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();
        cell_population.AddPopulationWriter<NodeVelocityWriter>();
        //cell_population.AddPopulationWriter<VoronoiDataWriter>(); // paraview is pretty pointless at viewing this, worth looking into

        //cell_population.WriteVtkResultsToFile(output_directory);
        //To fix paraview
        cell_population.SetWriteVtkAsPointsDom(true);


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
        // MAKE_PTR(PeriodicCryptModelInteractionForceWithGhostNodes<3>, periodic_spring_force);
        // periodic_spring_force->SetUseOneWaySprings(false);
        // periodic_spring_force->SetCutOffLength(1.5);
        // //              SetEpithelialStromalCellDependentSprings(ind , Ep-Ep, Str-Str, Ep-Str, apcTwoHitStromalMultiplier);
        // periodic_spring_force->SetEpithelialStromalCellDependentSprings(true, 1.0,     1.0,     1.0,    1.0);
        // periodic_spring_force->SetPeriodicDomainWidth(periodic_width);
        // periodic_spring_force->SetPeriodicDomainDepth(periodic_height);
        // periodic_spring_force->SetMeinekeSpringStiffness(10.0);
        //simulator.AddForce(periodic_spring_force);


        // Get the dimensions for the non-zero target curvature region
        double centre_x = double (width)*0.5;
        double centre_y = double (height)*0.5;
        double radius =  1.75;
        double target_curvature = 0.15; //maximum curvature is 0.2066 -> higher curvature means smaller sphere
        double beta_parameter = 5.0;
        double alpha_parameter = 2.0;

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

 
        // Add cell sloughing
        MAKE_PTR_ARGS(SloughingCellKiller3DWithGhostNodes, sloughing, (&cell_population, periodic_width, periodic_height));
        simulator.AddCellKiller(sloughing);

        double cut_off = 10.0;
        // Add anoikis cell killer
        MAKE_PTR_ARGS(AnoikisCellKiller3DWithGhostNodes, anoikis, (&cell_population, cut_off));
        simulator.AddCellKiller(anoikis);


        //Test Forces
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
            //node_forces.push_back(zero_vector<double>(3));
        }

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(0.001, 1);

        periodic_bending_force->AddForceContribution(cell_population);

        SimulationTime::Destroy();

        /*
        * A check to see that the force model we have implemented is actually working.
        * This does require that:
        * target_curvature = 0.15
        * beta_parameter   = 5.0;
        * alpha_parameter  = 2.0;
        *
        * Forces are found using the Mathematica script " TestingBendingForce.nb "
        */
        for (AbstractCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
         cell_iter != cell_population.End();
         ++cell_iter)
        {

            unsigned real_node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
            c_vector<double, 3> real_node_location = cell_population.GetLocationOfCellCentre(*cell_iter);

            CellPtr p_cell_i_ext = cell_population.GetCellUsingLocationIndex(real_node_index);

            if ( p_cell_i_ext->GetMutationState()->IsType<WildTypeCellMutationState>() == true &&
            (pow(real_node_location[0] - centre_x,2) + pow(real_node_location[1] - centre_y,2) < radius*radius) )
            {
                TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[0], 0.0, 1e-4);
                TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[1], 0.0, 1e-4);
                TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[2], -3.7452, 1e-4);
            }
            else
            {
                TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[0], 0.0, 1e-4);
                TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[1], 0.0, 1e-4);
                TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[2], 0.0, 1e-4);
            }

        }
        
    }

    /*
    * This checks to make sure the periodic BCs in the bending force work as expected
    */
    void TestPeriodicCubeWithGhosts_WholeTissue() throw (Exception)
    {
    	RandomNumberGenerator::Instance()->Reseed(1);

        std::vector<Node<3>*> nodes;

        unsigned width = 10;	   // x
        unsigned height = 10;      // y
        unsigned ghosts_bottom = 0;       // ghosts > depth
        unsigned ghosts_top = 2;       // ghosts > depth
        unsigned num_tissue_depth = 1;
        unsigned depth = num_tissue_depth + (ghosts_bottom + ghosts_top) + 1;        // z

        // Initialise the tissue in an equilibrum state
        double width_space = 1.0;
        double height_space = sqrt(0.75);
        double depth_space = 0.738431690356779; //Magic number for z-spaceing... 
        unsigned cells_per_layer = width*height;
        unsigned cell_iter = 0;

        double periodic_width = (double) (width+0.0)*width_space + 0.00;
        double periodic_height = (double) (height+0.0)*height_space + 0.00;

        bool isGhost = false;
        double x_coordinate, y_coordinate, z_coordinate;

        std::vector<unsigned> ghost_node_indices, real_node_indices;
        for (unsigned k=0; k<depth; k++)
        {
            isGhost = false;
            if(k < ghosts_bottom || k > ghosts_bottom + num_tissue_depth)
            {
                isGhost = true;
            }
            for (unsigned j=0; j<height; j++)
            {    
                for (unsigned i=0; i<width; i++)
                {
                    c_vector<double, 3> node_i_new_location;

                    //x_coordinate = (double) (i + 0.5*((j%2 + k%2)%2))*width_space;
                    x_coordinate = (double) (i + 0.5*(j%2 + k%2))*width_space;
                    y_coordinate = (double) j*height_space;
                    z_coordinate = (double) 5.0 + k*depth_space;
                    
                    /*
                    if (cell_iter >= 100 && cell_iter < 200)
                    {
                        z_coordinate = (double) 5.0 + k*depth_space + 0.5*(2.0*RandomNumberGenerator::Instance()->ranf()-1.0);
                    }
                    
                    if (cell_iter == 154)
                    {
                        z_coordinate = (double) 5.0 + k*depth_space + 0.25;
                    }
                    */
                    nodes.push_back(new Node<3>(cell_iter,  false,  x_coordinate, y_coordinate, z_coordinate));
                    
                    if (isGhost == true)
                    {
                        ghost_node_indices.push_back(cell_iter);
                    }
                    else if(isGhost == false)
                    {
                        real_node_indices.push_back(cell_iter);
                    }

                    cell_iter++;
                }
            }
        }

        MutableMesh<3,3> mesh(nodes);
        //mesh.Translate(0.5, 0.5);


        std::vector<CellPtr> cells;

    	// First we sort the real node indices into increasing order (the last ones will correspond to the
    	// epithelial nodes)
    	sort(real_node_indices.begin(), real_node_indices.end());

        // Set up cells by iterating through the mesh nodes
        unsigned num_nodes = mesh.GetNumAllNodes();
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

			cells.push_back(p_differentiated_cell);
        }


        // Initialise Epithelial cells
        for (unsigned i=real_node_indices.size()-num_epithelial_cells; i<real_node_indices.size(); i++)
        {

        	FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_epithelial_cell(new Cell(p_state, p_model));
            p_epithelial_cell->SetCellProliferativeType(p_transit_type);
            p_model->SetDimension(3);
                        
            p_model->SetMaxTransitGenerations(1);
            p_model->SetTransitCellG1Duration(1);
            double birth_time = 0.0;
			
            p_epithelial_cell->SetBirthTime(birth_time);
            
            p_epithelial_cell->InitialiseCellCycleModel();
            
			cells.push_back(p_epithelial_cell);
        }

        DomMeshBasedCellPopulationWithGhostNodes<3> cell_population(mesh, cells, real_node_indices);
        //MeshBasedCellPopulationWithGhostNodes<3> cell_population(mesh, cells, real_node_indices);
        //  MeshBasedCellPopulation<3> cell_population(mesh, cells);
        assert(cell_population.GetNumRealCells() != 0);


        //CryptSimulation3dGhosts simulator(cell_population, false, true);
        //CryptSimulation3d(rCellPopulation, bool deleteCellPopulationAndForceCollection, bool initialiseCells)
        OffLatticeSimulation<3> simulator(cell_population);
        std::string output_directory = "TestBendingForceCheck";
        simulator.SetOutputDirectory(output_directory);

        // Make sure we have a Voronoi tessellation to begin with
        //cell_population.CreateVoronoiTessellation();      
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();
        cell_population.AddPopulationWriter<NodeVelocityWriter>();
        //cell_population.AddPopulationWriter<VoronoiDataWriter>(); // paraview is pretty pointless at viewing this, worth looking into
        cell_population.SetWriteVtkAsPointsDom(true);




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
        // MAKE_PTR(PeriodicCryptModelInteractionForceWithGhostNodes<3>, periodic_spring_force);
        // periodic_spring_force->SetUseOneWaySprings(false);
        // periodic_spring_force->SetCutOffLength(1.5);
        // //              SetEpithelialStromalCellDependentSprings(ind , Ep-Ep, Str-Str, Ep-Str, apcTwoHitStromalMultiplier);
        // periodic_spring_force->SetEpithelialStromalCellDependentSprings(true, 1.0,     1.0,     1.0,    1.0);
        // periodic_spring_force->SetPeriodicDomainWidth(periodic_width);
        // periodic_spring_force->SetPeriodicDomainDepth(periodic_height);
        // periodic_spring_force->SetMeinekeSpringStiffness(10.0);
        //simulator.AddForce(periodic_spring_force);

        // Get the dimensions for the non-zero target curvature region
        double centre_x = double (width)*0.5;
        double centre_y = double (height)*0.5;
        // Whole tissue
        double radius =  pow(periodic_width,2) + pow(periodic_height,2);
        double target_curvature = 0.15; //maximum curvature is 0.2066 -> higher curvature means smaller sphere
        double beta_parameter = 5.0;
        double alpha_parameter = 2.0;

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
 
        // Add cell sloughing
        MAKE_PTR_ARGS(SloughingCellKiller3DWithGhostNodes, sloughing, (&cell_population, periodic_width, periodic_height));
        simulator.AddCellKiller(sloughing);

        double cut_off = 10.0;
        // Add anoikis cell killer
        MAKE_PTR_ARGS(AnoikisCellKiller3DWithGhostNodes, anoikis, (&cell_population, cut_off));
        simulator.AddCellKiller(anoikis);


        //Test Forces
        
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
            //node_forces.push_back(zero_vector<double>(3));
        }

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(0.001, 1);

        periodic_bending_force->AddForceContribution(cell_population);

        SimulationTime::Destroy();

        /*
        * A check to see that the force model we have implemented is actually working.
        * This does require that:
        * target_curvature = 0.15
        * beta_parameter   = 5.0;
        * alpha_parameter  = 2.0;
        *
        * Forces are found using the Mathematica script " TestingBendingForce.nb "
        */
        for (AbstractCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
         cell_iter != cell_population.End();
         ++cell_iter)
        {

            unsigned real_node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
            c_vector<double, 3> real_node_location = cell_population.GetLocationOfCellCentre(*cell_iter);

            CellPtr p_cell_i_ext = cell_population.GetCellUsingLocationIndex(real_node_index);
            
            if ( p_cell_i_ext->GetMutationState()->IsType<WildTypeCellMutationState>() == true &&
            (pow(real_node_location[0] - centre_x,2) + pow(real_node_location[1] - centre_y,2) < radius*radius) )
            {
                TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[0], 0.0, 1e-4);
                TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[1], 0.0, 1e-4);
                TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[2], -3.7452, 1e-4);
            }
            else
            {
                TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[0], 0.0, 1e-4);
                TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[1], 0.0, 1e-4);
                TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[2], 0.0, 1e-4);
            }

        }

        
    }

    /*
    * This checks the force to flatten when a single cell has popped up
    */
    void TestPeriodicCubeWithGhosts_PoppedUpCell() throw (Exception)
    {
    	RandomNumberGenerator::Instance()->Reseed(1);

        std::vector<Node<3>*> nodes;

        unsigned width = 10;	   // x
        unsigned height = 10;      // y
        unsigned ghosts_bottom = 0;       // ghosts > depth
        unsigned ghosts_top = 2;       // ghosts > depth
        unsigned num_tissue_depth = 1;
        unsigned depth = num_tissue_depth + (ghosts_bottom + ghosts_top) + 1;        // z

        // Initialise the tissue in an equilibrum state
        double width_space = 1.0;
        double height_space = sqrt(0.75);
        double depth_space = 0.738431690356779; //Magic number for z-spaceing... 
        unsigned cells_per_layer = width*height;
        unsigned cell_iter = 0;

        double periodic_width = (double) (width+0.0)*width_space + 0.00;
        double periodic_height = (double) (height+0.0)*height_space + 0.00;

        bool isGhost = false;
        double x_coordinate, y_coordinate, z_coordinate;

        std::vector<unsigned> ghost_node_indices, real_node_indices;
        for (unsigned k=0; k<depth; k++)
        {
            isGhost = false;
            if(k < ghosts_bottom || k > ghosts_bottom + num_tissue_depth)
            {
                isGhost = true;
            }
            for (unsigned j=0; j<height; j++)
            {    
                for (unsigned i=0; i<width; i++)
                {
                    c_vector<double, 3> node_i_new_location;

                    //x_coordinate = (double) (i + 0.5*((j%2 + k%2)%2))*width_space;
                    x_coordinate = (double) (i + 0.5*(j%2 + k%2))*width_space;
                    y_coordinate = (double) j*height_space;
                    z_coordinate = (double) 5.0 + k*depth_space;
                    
                    /*
                    if (cell_iter >= 100 && cell_iter < 200)
                    {
                        z_coordinate = (double) 5.0 + k*depth_space + 0.5*(2.0*RandomNumberGenerator::Instance()->ranf()-1.0);
                    }
                    */
                    if (cell_iter == 154)
                    {
                        z_coordinate = (double) 5.0 + k*depth_space + 0.25;
                    }
                    
                    nodes.push_back(new Node<3>(cell_iter,  false,  x_coordinate, y_coordinate, z_coordinate));
                    
                    if (isGhost == true)
                    {
                        ghost_node_indices.push_back(cell_iter);
                    }
                    else if(isGhost == false)
                    {
                        real_node_indices.push_back(cell_iter);
                    }

                    cell_iter++;
                }
            }
        }

        MutableMesh<3,3> mesh(nodes);
        //mesh.Translate(0.5, 0.5);


        std::vector<CellPtr> cells;

    	// First we sort the real node indices into increasing order (the last ones will correspond to the
    	// epithelial nodes)
    	sort(real_node_indices.begin(), real_node_indices.end());

        // Set up cells by iterating through the mesh nodes
        unsigned num_nodes = mesh.GetNumAllNodes();
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

			cells.push_back(p_differentiated_cell);
        }


        // Initialise Epithelial cells
        for (unsigned i=real_node_indices.size()-num_epithelial_cells; i<real_node_indices.size(); i++)
        {

        	FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_epithelial_cell(new Cell(p_state, p_model));
            p_epithelial_cell->SetCellProliferativeType(p_transit_type);
            p_model->SetDimension(3);
                        
            p_model->SetMaxTransitGenerations(1);
            p_model->SetTransitCellG1Duration(1);
            double birth_time = 0.0;
			
            p_epithelial_cell->SetBirthTime(birth_time);
            
            p_epithelial_cell->InitialiseCellCycleModel();
            
			cells.push_back(p_epithelial_cell);
        }

        DomMeshBasedCellPopulationWithGhostNodes<3> cell_population(mesh, cells, real_node_indices);
        //MeshBasedCellPopulationWithGhostNodes<3> cell_population(mesh, cells, real_node_indices);
        // MeshBasedCellPopulation<3> cell_population(mesh, cells);
        assert(cell_population.GetNumRealCells() != 0);


        //CryptSimulation3dGhosts simulator(cell_population, false, true);
        //CryptSimulation3d(rCellPopulation, bool deleteCellPopulationAndForceCollection, bool initialiseCells)
        OffLatticeSimulation<3> simulator(cell_population);
        std::string output_directory = "TestBendingForceCheck";
        simulator.SetOutputDirectory(output_directory);

        // Make sure we have a Voronoi tessellation to begin with
        //cell_population.CreateVoronoiTessellation();      
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();
        cell_population.AddPopulationWriter<NodeVelocityWriter>();
        //cell_population.AddPopulationWriter<VoronoiDataWriter>(); // paraview is pretty pointless at viewing this, worth looking into
        cell_population.SetWriteVtkAsPointsDom(true);




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
        // MAKE_PTR(PeriodicCryptModelInteractionForceWithGhostNodes<3>, periodic_spring_force);
        // periodic_spring_force->SetUseOneWaySprings(false);
        // periodic_spring_force->SetCutOffLength(1.5);
        // //              SetEpithelialStromalCellDependentSprings(ind , Ep-Ep, Str-Str, Ep-Str, apcTwoHitStromalMultiplier);
        // periodic_spring_force->SetEpithelialStromalCellDependentSprings(true, 1.0,     1.0,     1.0,    1.0);
        // periodic_spring_force->SetPeriodicDomainWidth(periodic_width);
        // periodic_spring_force->SetPeriodicDomainDepth(periodic_height);
        // periodic_spring_force->SetMeinekeSpringStiffness(10.0);
        //simulator.AddForce(periodic_spring_force);

        // Get the dimensions for the non-zero target curvature region
        double centre_x = double (width)*0.5;
        double centre_y = double (height)*0.5;
        // Whole tissue
        double radius =  0.0;
        double target_curvature = 0.15; //maximum curvature is 0.2066 -> higher curvature means smaller sphere
        double beta_parameter = 5.0;
        double alpha_parameter = 2.0;

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
 
        // Add cell sloughing
        MAKE_PTR_ARGS(SloughingCellKiller3DWithGhostNodes, sloughing, (&cell_population, periodic_width, periodic_height));
        simulator.AddCellKiller(sloughing);

        double cut_off = 10.0;
        // Add anoikis cell killer
        MAKE_PTR_ARGS(AnoikisCellKiller3DWithGhostNodes, anoikis, (&cell_population, cut_off));
        simulator.AddCellKiller(anoikis);


        //Test Forces
        
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
            //node_forces.push_back(zero_vector<double>(3));
        }

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(0.001, 1);

        periodic_bending_force->AddForceContribution(cell_population);

        SimulationTime::Destroy();

        /*
        * A check to see that the force model we have implemented is actually working.
        * This does require that:
        * target_curvature = 0.15
        * beta_parameter   = 5.0;
        * alpha_parameter  = 2.0;
        *
        * Forces are found using the Mathematica script " TestingBendingForce.nb "
        */
        for (AbstractCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
         cell_iter != cell_population.End();
         ++cell_iter)
        {

            unsigned real_node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
            c_vector<double, 3> real_node_location = cell_population.GetLocationOfCellCentre(*cell_iter);

            CellPtr p_cell_i_ext = cell_population.GetCellUsingLocationIndex(real_node_index);
            
            if ( real_node_index == 154 )
            {
                TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[0], 0.0, 1e-4);
                TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[1], 0.0, 1e-4);
                TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[2], -2.6867, 1e-4);
            }
            // else
            // {
            //     TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[0], 0.0, 1e-4);
            //     TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[1], 0.0, 1e-4);
            //     TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[2], 0.0, 1e-4);
            // }

        }

        
    }


    /*
    * This checks the force to flatten when a single cell has popped Down
    */
    void TestPeriodicCubeWithGhosts_PoppedDownCell() throw (Exception)
    {
    	RandomNumberGenerator::Instance()->Reseed(1);

        std::vector<Node<3>*> nodes;

        unsigned width = 10;	   // x
        unsigned height = 10;      // y
        unsigned ghosts_bottom = 0;       // ghosts > depth
        unsigned ghosts_top = 2;       // ghosts > depth
        unsigned num_tissue_depth = 1;
        unsigned depth = num_tissue_depth + (ghosts_bottom + ghosts_top) + 1;        // z

        // Initialise the tissue in an equilibrum state
        double width_space = 1.0;
        double height_space = sqrt(0.75);
        double depth_space = 0.738431690356779; //Magic number for z-spaceing... 
        unsigned cells_per_layer = width*height;
        unsigned cell_iter = 0;

        double periodic_width = (double) (width+0.0)*width_space + 0.00;
        double periodic_height = (double) (height+0.0)*height_space + 0.00;

        bool isGhost = false;
        double x_coordinate, y_coordinate, z_coordinate;

        std::vector<unsigned> ghost_node_indices, real_node_indices;
        for (unsigned k=0; k<depth; k++)
        {
            isGhost = false;
            if(k < ghosts_bottom || k > ghosts_bottom + num_tissue_depth)
            {
                isGhost = true;
            }
            for (unsigned j=0; j<height; j++)
            {    
                for (unsigned i=0; i<width; i++)
                {
                    c_vector<double, 3> node_i_new_location;

                    //x_coordinate = (double) (i + 0.5*((j%2 + k%2)%2))*width_space;
                    x_coordinate = (double) (i + 0.5*(j%2 + k%2))*width_space;
                    y_coordinate = (double) j*height_space;
                    z_coordinate = (double) 5.0 + k*depth_space;
                    
                    /*
                    if (cell_iter >= 100 && cell_iter < 200)
                    {
                        z_coordinate = (double) 5.0 + k*depth_space + 0.5*(2.0*RandomNumberGenerator::Instance()->ranf()-1.0);
                    }
                    */
                    if (cell_iter == 154)
                    {
                        z_coordinate = (double) 5.0 + k*depth_space - 0.25;
                    }
                    
                    nodes.push_back(new Node<3>(cell_iter,  false,  x_coordinate, y_coordinate, z_coordinate));
                    
                    if (isGhost == true)
                    {
                        ghost_node_indices.push_back(cell_iter);
                    }
                    else if(isGhost == false)
                    {
                        real_node_indices.push_back(cell_iter);
                    }

                    cell_iter++;
                }
            }
        }

        MutableMesh<3,3> mesh(nodes);
        //mesh.Translate(0.5, 0.5);


        std::vector<CellPtr> cells;

    	// First we sort the real node indices into increasing order (the last ones will correspond to the
    	// epithelial nodes)
    	sort(real_node_indices.begin(), real_node_indices.end());

        // Set up cells by iterating through the mesh nodes
        unsigned num_nodes = mesh.GetNumAllNodes();
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

			cells.push_back(p_differentiated_cell);
        }


        // Initialise Epithelial cells
        for (unsigned i=real_node_indices.size()-num_epithelial_cells; i<real_node_indices.size(); i++)
        {

        	FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_epithelial_cell(new Cell(p_state, p_model));
            p_epithelial_cell->SetCellProliferativeType(p_transit_type);
            p_model->SetDimension(3);
                        
            p_model->SetMaxTransitGenerations(1);
            p_model->SetTransitCellG1Duration(1);
            double birth_time = 0.0;
			
            p_epithelial_cell->SetBirthTime(birth_time);
            
            p_epithelial_cell->InitialiseCellCycleModel();
            
			cells.push_back(p_epithelial_cell);
        }

        DomMeshBasedCellPopulationWithGhostNodes<3> cell_population(mesh, cells, real_node_indices);
        //MeshBasedCellPopulationWithGhostNodes<3> cell_population(mesh, cells, real_node_indices);
        // MeshBasedCellPopulation<3> cell_population(mesh, cells);
        assert(cell_population.GetNumRealCells() != 0);


        //CryptSimulation3dGhosts simulator(cell_population, false, true);
        //CryptSimulation3d(rCellPopulation, bool deleteCellPopulationAndForceCollection, bool initialiseCells)
        OffLatticeSimulation<3> simulator(cell_population);
        std::string output_directory = "TestBendingForceCheck";
        simulator.SetOutputDirectory(output_directory);

        // Make sure we have a Voronoi tessellation to begin with
        //cell_population.CreateVoronoiTessellation();      
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();
        cell_population.AddPopulationWriter<NodeVelocityWriter>();
        //cell_population.AddPopulationWriter<VoronoiDataWriter>(); // paraview is pretty pointless at viewing this, worth looking into
        cell_population.SetWriteVtkAsPointsDom(true);




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
        // MAKE_PTR(PeriodicCryptModelInteractionForceWithGhostNodes<3>, periodic_spring_force);
        // periodic_spring_force->SetUseOneWaySprings(false);
        // periodic_spring_force->SetCutOffLength(1.5);
        // //              SetEpithelialStromalCellDependentSprings(ind , Ep-Ep, Str-Str, Ep-Str, apcTwoHitStromalMultiplier);
        // periodic_spring_force->SetEpithelialStromalCellDependentSprings(true, 1.0,     1.0,     1.0,    1.0);
        // periodic_spring_force->SetPeriodicDomainWidth(periodic_width);
        // periodic_spring_force->SetPeriodicDomainDepth(periodic_height);
        // periodic_spring_force->SetMeinekeSpringStiffness(10.0);
        //simulator.AddForce(periodic_spring_force);

        // Get the dimensions for the non-zero target curvature region
        double centre_x = double (width)*0.5;
        double centre_y = double (height)*0.5;
        // Whole tissue
        double radius =  0.0;
        double target_curvature = 0.15; //maximum curvature is 0.2066 -> higher curvature means smaller sphere
        double beta_parameter = 5.0;
        double alpha_parameter = 2.0;

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
 
        // Add cell sloughing
        MAKE_PTR_ARGS(SloughingCellKiller3DWithGhostNodes, sloughing, (&cell_population, periodic_width, periodic_height));
        simulator.AddCellKiller(sloughing);

        double cut_off = 10.0;
        // Add anoikis cell killer
        MAKE_PTR_ARGS(AnoikisCellKiller3DWithGhostNodes, anoikis, (&cell_population, cut_off));
        simulator.AddCellKiller(anoikis);


        //Test Forces
        
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
            //node_forces.push_back(zero_vector<double>(3));
        }

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(0.001, 1);

        periodic_bending_force->AddForceContribution(cell_population);

        SimulationTime::Destroy();

        /*
        * A check to see that the force model we have implemented is actually working.
        * This does require that:
        * target_curvature = 0.15
        * beta_parameter   = 5.0;
        * alpha_parameter  = 2.0;
        *
        * Forces are found using the Mathematica script " TestingBendingForce.nb "
        */
        for (AbstractCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
         cell_iter != cell_population.End();
         ++cell_iter)
        {

            unsigned real_node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
            c_vector<double, 3> real_node_location = cell_population.GetLocationOfCellCentre(*cell_iter);

            CellPtr p_cell_i_ext = cell_population.GetCellUsingLocationIndex(real_node_index);
            
            if ( real_node_index == 154 )
            {
                TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[0], 0.0, 1e-4);
                TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[1], 0.0, 1e-4);
                TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[2], 2.6867, 1e-4);
            }
            // else
            // {
            //     TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[0], 0.0, 1e-4);
            //     TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[1], 0.0, 1e-4);
            //     TS_ASSERT_DELTA(cell_population.GetNode(real_node_index)->rGetAppliedForce()[2], 0.0, 1e-4);
            // }

        }

        
    }
     
           
};

#endif /*TEST3DBOXMODEL_HPP_*/

