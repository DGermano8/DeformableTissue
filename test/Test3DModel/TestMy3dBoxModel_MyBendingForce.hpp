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
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellsGenerator.hpp"
#include "CryptModelInteractionForce.hpp"
#include "PeriodicCryptModelInteractionForce.hpp"
//#include "BasementMembraneForce3d.hpp"
#include "PeriodicBendingForce3d.hpp"
#include "CryptSimulation3d.hpp"
#include "SloughingCellKiller3D.hpp"
#include "AnoikisCellKiller3D.hpp"
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
    void TestSimpleCubeBending() throw (Exception)
    {
    	RandomNumberGenerator::Instance()->Reseed(1);

        unsigned width = 10;	   // x
        unsigned height = 10;      // y
        unsigned depth = 3;        // z

        // Get the dimensions for the non-zero target curvature region
        double centre_x = double (width)*0.5;
        double centre_y = double (height)*0.5;
        double radius =  1.75;
        double target_curvature = 0.2066;
        double beta_parameter = 5.0;
        double alpha_parameter = 1.5;

//        double periodic_width = (double) width + 0.0;
//        double periodic_height = (double) height + 0.0;
        
        MutableMesh<3,3>* p_mesh = Make3dMesh(width, height, depth);
        p_mesh->Translate(0.5, 0.5);
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_nodes = p_mesh->GetNumAllNodes();
        
        unsigned num_epithelial_cells = (width)*(height);

        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellMutationState> p_stromal_state(new StromalCellMutationState);

        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        std::vector<CellPtr> cells;

        // Initialise Tissue cells (Stromal)
		for (unsigned i=0; i<num_nodes-num_epithelial_cells; i++)
		{
			//StochasticDurationGenerationBasedCellCycleModel* p_model = new StochasticDurationGenerationBasedCellCycleModel();
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();

            CellPtr p_differentiated_cell(new Cell(p_stromal_state, p_model));
            p_differentiated_cell->SetCellProliferativeType(p_differentiated_type);
            p_model->SetDimension(3);

			cells.push_back(p_differentiated_cell);
        }

        // Initialise Epithelial cells
        for (unsigned i=num_nodes-num_epithelial_cells; i<num_nodes; i++)
        {
        	FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_epithelial_cell(new Cell(p_state, p_model));
            p_epithelial_cell->SetCellProliferativeType(p_transit_type);
            p_model->SetDimension(3);
            
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*(p_model->GetTransitCellG1Duration() + p_model->GetSG2MDuration());
            
            //p_model->SetMaxTransitGenerations(100);
            //p_model->SetTransitCellG1Duration(100);
            //double birth_time = 0.0;
            //std::cout << "cell i = " << i << " birth_time = " << birth_time << "\n";
			
            p_epithelial_cell->SetBirthTime(birth_time);
            
            p_epithelial_cell->InitialiseCellCycleModel();
            
			cells.push_back(p_epithelial_cell);
        }

        MeshBasedCellPopulation<3> cell_population(*p_mesh, cells);
        assert(cell_population.GetNumRealCells() != 0);

        // Initialise the tissue in an equilibrum state
        double width_space = 1.0;
        double height_space = sqrt(0.75);
        double depth_space = 0.738431690356779; //Magic number for z-spaceing... 
        unsigned cells_per_layer = width*height;
        unsigned cell_iter = 0;

        double periodic_width = (double) (width+0.0)*width_space;
        double periodic_height = (double) (height+0.0)*height_space;
        for (unsigned k=0; k<depth; k++)
        {
            for (unsigned j=0; j<height; j++)
            {    
                for (unsigned i=0; i<width; i++)
                {
                    c_vector<double, 3> node_i_new_location;

                    node_i_new_location[0] = (double) (i + 0.5*((j%2 + k%2)%2))*width_space;
                    node_i_new_location[1] = (double) j*height_space;
                    node_i_new_location[2] = (double) 5.0 + k*depth_space;

                    ChastePoint<3> node_i_new_point(node_i_new_location);
                    Node<3>* p_node_i = p_mesh->GetNode(cell_iter);
                    p_node_i->SetPoint(node_i_new_point);

                    cell_iter++;
                }
            }
        }

        CryptSimulation3d simulator(cell_population, false, true);

        simulator.SetOutputDirectory("TestSimulation_MyBendingForce");
//        cell_population.InitialiseCells();

        // Make sure we have a Voronoi tessellation to begin with
        //cell_population.CreateVoronoiTessellation();      
        
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();
        cell_population.AddPopulationWriter<NodeVelocityWriter>();
//        cell_population.AddPopulationWriter<VoronoiDataWriter>(); // paraview is pretty pointless at viewing this, worth looking into

        //To fix paraview
        cell_population.SetWriteVtkAsPoints(true);
        //cell_population.SetOutputMeshInVtk(true);

        // CryptSimulation3d(rCellPopulation, bool deleteCellPopulationAndForceCollection, bool initialiseCells)

        std::map<Node<3>*, c_vector<double,3> > node_locations_before;
        for (std::list<CellPtr>::iterator cell_iter = cell_population.rGetCells().begin();
             cell_iter != cell_population.rGetCells().end();
             ++cell_iter)
        {
            Node<3>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
            node_locations_before[p_node] = p_node->rGetLocation();
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
        MAKE_PTR(PeriodicCryptModelInteractionForce<3>, periodic_force);
     // PeriodicCryptModelInteractionForce<3> periodic_force;  // Variable spring strengths
        periodic_force->SetUseOneWaySprings(false);
        periodic_force->SetCutOffLength(1.5);
        //              SetEpithelialStromalCellDependentSprings(ind , Ep-Ep, Str-Str, Ep-Str, apcTwoHitStromalMultiplier);
        periodic_force->SetEpithelialStromalCellDependentSprings(true, 1.0,   2.0,     1.0,    1.0);
        periodic_force->SetPeriodicDomainWidth(periodic_width);
        periodic_force->SetPeriodicDomainDepth(periodic_height);
        periodic_force->SetMeinekeSpringStiffness(10.0);
        simulator.AddForce(periodic_force);
/*        
		// Create periodic basement membrane force law
        MAKE_PTR(PeriodicBendingForce3d, periodic_bending_force);
     // PeriodicBasementMembraneForce3d periodic_basement_membrane;
        periodic_bending_force->SetBasementMembraneParameter(beta_parameter);
        periodic_bending_force->SetExponentParameter(alpha_parameter);
        //periodic_bending_force->SetCircularNonZeroTargetCurvatureRegion(true, 0.15, radius, centre_x, centre_y);
        periodic_bending_force->SetCircularNonZeroTargetCurvatureRegion(true, target_curvature, radius, centre_x, centre_y);
        periodic_bending_force->SetPeriodicDomainWidth(periodic_width);
        periodic_bending_force->SetPeriodicDomainDepth(periodic_height);
        simulator.AddForce(periodic_bending_force);
 */
        // Add cell sloughing
        MAKE_PTR_ARGS(SloughingCellKiller3D, sloughing, (&cell_population, periodic_width, periodic_height));
     // SloughingCellKiller3D sloughing(&cell_population, periodic_width, periodic_width);
        simulator.AddCellKiller(sloughing);

                // Add anoikis cell killer
        MAKE_PTR_ARGS(AnoikisCellKiller3D, anoikis, (&cell_population));
     // AnoikisCellKiller3D anoikis(&cell_population);
        simulator.AddCellKiller(anoikis);


        // Run for a short time to allow it to deform		        
    	simulator.SetSamplingTimestepMultiple(1);			// Every hour
		simulator.SetEndTime(0.002);							// Firstly run to a steady state (i.e. til the thing is deformed as much as it will be)
        simulator.SetDt(0.001);	
        simulator.Solve();



        //Test Forces
        /*
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
            //node_forces.push_back(zero_vector<double>(3));
        }

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(0.01, 1);

        periodic_basement_membrane->AddForceContribution(cell_population);
        //periodic_force->AddForceContribution(cell_population);

        SimulationTime::Destroy();

        for (unsigned i=0; i<num_nodes; i++)
        {
            std::cout<< "Cell i = " << i << "\n";
            PRINT_VECTOR(cell_population.GetNode(i)->rGetAppliedForce());
        }
        */
    }     
     
           
};

#endif /*TEST3DBOXMODEL_HPP_*/

