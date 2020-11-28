#ifndef PERIODICBENDINGFORCE3D_HPP_
#define PERIODICBENDINGFORCE3D_HPP_

#include "GeneralisedLinearSpringForce.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "AbstractVanLeeuwen2009WntSwatCellCycleModel.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "Debug.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "StromalCellMutationState.hpp"
#include <cmath>
#include <list>
#include <fstream>

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/*
 * Mechanics class that defines how the basement membrane force should act in the 3D Crypt Model.
 *
 */

class PeriodicBendingForce3d : public AbstractForce<3>
{
    friend class TestPeriodicForces3d;

private :

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<AbstractForce<3> >(*this);
        archive & mBasementMembraneParameter;
        archive & mExponentParameter;
        archive & mPeriodicDomainWidth;
        archive & mPeriodicDomainDepth;
        archive & mSetNonZeroTargetCurvatureRegion;
        archive & mNonZeroTargetCurvature;
        archive & mRadiusOfNonZeroTargetCurvatureRegion;
        archive & mXCoordinateOfCentreOfNonZeroTargetCurvatureRegion;
        archive & mYCoordinateOfCentreOfNonZeroTargetCurvatureRegion;
    }


protected :
 
    /** Parameter that multiplies the curvature to give the basal lamina force */
    double mBasementMembraneParameter;

    double mExponentParameter;

    /** Width of the periodic domain. */
    double mPeriodicDomainWidth;

    /** Depth of the periodic domain. */
    double mPeriodicDomainDepth;

    /** An extended mesh, used to implement periodicity. */
    MutableMesh<3,3>* mpExtendedMesh;

    /** A map from node indices in mpExtendedMesh to node indices in the cell population. */
    std::map<unsigned, unsigned> mExtendedMeshNodeIndexMap;

    /** Whether or not you will prescribe a non-zero target curvature region */
    bool mSetNonZeroTargetCurvatureRegion;

    /** Target curvature for the prescribed circular region corresponding to the crypt base */
    double mNonZeroTargetCurvature;

    /** Radius of circular non-zero target curvature region */
    double mRadiusOfNonZeroTargetCurvatureRegion;

    /** x-coordinate of the centre of the circular non-zero target curvature region */
    double mXCoordinateOfCentreOfNonZeroTargetCurvatureRegion;

    /** y-coordinate of the centre of the circular non-zero target curvature region */
    double mYCoordinateOfCentreOfNonZeroTargetCurvatureRegion;


public :

    /**
     * Constructor.
     */
    PeriodicBendingForce3d();

    /**
     * Destructor.
     */
    ~PeriodicBendingForce3d();

    /* Set method for Basement Membrane Parameter
     */
    void SetBasementMembraneParameter(double basementMembraneParameter);

    /* Get method for Basement Membrane Parameter
     */
    double GetBasementMembraneParameter();


    void SetExponentParameter(double exponentParameter);

    double GetExponentParameter();

    /* Setting the boundaries for the non-zero target curvature region
     *
     * @param setNonZeroTargetCurvatureRegion whether you want a non-zero target curvature region (defaults to false in constructor)
     * @param nonZeroTargetCurvature the value of the non-zero target curvature
     * @param radiusOfNonZeroTargetCurvatureRegion the radius of the non-zero target curvature region
     * @param xCoordinateOfCentreOfNonZeroTargetCurvatureRegion the x-coordinate of the centre of the non-zero target curvature region
     * @param yCoordinateOfCentreOfNonZeroTargetCurvatureRegion the y-coordinate of the centre of the non-zero target curvature region
     */
    void SetCircularNonZeroTargetCurvatureRegion(bool setNonZeroTargetCurvatureRegion, double nonZeroTargetCurvature,
												double radiusOfNonZeroTargetCurvatureRegion, double xCoordinateOfCentreOfNonZeroTargetCurvatureRegion,
												double yCoordinateOfCentreOfNonZeroTargetCurvatureRegion);


    /* Get method for whether or not you want a non-zero target curvature region */
    bool GetWhetherToSetNonZeroTargetCurvatureRegion();

    /* Get method for the value of the non zero target curvature */
    double GetNonZeroTargetCurvature();

    /* Get the radius of the non-zero target curvature region */
    double GetRadiusOfNonZeroTargetCurvatureRegion();

    /* Get the x-coordinate of the centre of the non-zero target curvature region */
    double GetXCoordinateOfCentreOfNonZeroTargetCurvatureRegion();

    /* Get the y-coordinate of the centre of the non-zero target curvature region */
    double GetYCoordinateOfCentreOfNonZeroTargetCurvatureRegion();

    /* Removing duplicated entries of a vector */
    void RemoveDuplicates1D(std::vector<unsigned>& rVectorWithDuplicates);

    /* Returns a vector of all epithelial-tissue node pairings
     */
    std::vector<c_vector<unsigned, 2> > GetEpithelialTissuePairs(AbstractCellPopulation<3>& rCellPopulation);
    
    ////
    std::vector<c_vector<unsigned, 10> > GetEpithelialNeighbours(std::vector<c_vector<unsigned, 3> > rEpithelialMeshVector, int number_of_cells);

    std::vector<c_vector<unsigned, 3> > GetEpithelialMesh(AbstractCellPopulation<3>& rCellPopulation);

    std::vector<c_vector<double, 3> > FitPlaneAndFindImage(AbstractCellPopulation<3>& rCellPopulation, std::vector<unsigned> second_order_neighs, unsigned cell_i);

    c_vector<double, 3> GetForceDueToDiscreteCurvature(AbstractCellPopulation<3>& rCellPopulation, std::vector<c_vector<unsigned, 3> > rEpithelialMeshVector, std::vector<unsigned> first_order_neighs,  std::vector<unsigned> second_order_neighs, std::vector<c_vector<double, 3> > image_location, unsigned cell_i, unsigned num_cells);

    ////

//    double GetCurvatureDueToSpringMidpoints(AbstractCellPopulation<3>& rCellPopulation, unsigned epithelialNodeIndex, unsigned tissueNodeIndex);

//Doms add
    double GetCurvatureDueToSpringMidpoints(AbstractCellPopulation<3>& rCellPopulation, unsigned epithelialNodeIndex, unsigned tissueNodeIndex);

    /* Returns the area of the triangle ABC given two vectors: AB and AC
     */
    double GetAreaOfTriangle(c_vector<double,3> vectorOneToTwo, c_vector<double,3> vectorOneToThree);

    /* Returns a boolean for whether the element contains ghost nodes
     */
    bool DoesElementContainGhostNodes(AbstractCellPopulation<3>& rCellPopulation, unsigned elementIndex);

    /* Returns a boolean for whether this element has a long edge
     */
    bool DoesElementContainLongEdge(AbstractCellPopulation<3>& rCellPopulation, unsigned elementIndex, double maxEdgeLength);

    /**
     * Overridden AddForceContribution method.
     *
     * @param rForces reference to vector of forces on nodes
     * @param rCellPopulation reference to the tissue
     */
    void AddForceContribution(AbstractCellPopulation<3>& rCellPopulation);
    //void AddForceContribution(std::vector<c_vector<double, 3> >& rForces,
    //                          AbstractCellPopulation<3>& rCellPopulation);
    
    /**
     * Returns the initial width which we use to define the periodic boundaries.
     */
    double GetPeriodicDomainWidth();

    /**
     * Returns the initial width which we use to define the periodic boundaries.
     */
    void SetPeriodicDomainWidth(double periodicDomainWidth);

    /**
     * Returns the initial depth (y-direction) which we use to define the periodic boundaries (used in 3d only).
     */
    double GetPeriodicDomainDepth();

    /**
     * Returns the initial depth (y-direction) which we use to define the periodic boundaries (used in 3d only).
     */
    void SetPeriodicDomainDepth(double periodicDomainDepth);


    /**
     * Outputs force Parameters to file
	 *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(PeriodicBendingForce3d)

#endif /* PERIODICBENDINGFORCE3D_HPP_ */
