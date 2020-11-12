#ifndef BendingForce2DGhosts_HPP_
#define BendingForce2DGhosts_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractOffLatticeCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * A force class to model random cell movement.
 */
template<unsigned DIM>
class BendingForce2DGhosts : public AbstractForce<DIM>
{
private :

    double mBendingCoefficient;
    double mExponentParameter;
    double mDomainWidth;
    bool   mSetNonZeroTargetCurvatureRegion;
    double mNonZeroTargetCurvature;
    double mRadiusOfNonZeroTargetCurvatureRegion;
    double mXCoordinateOfCentreOfNonZeroTargetCurvatureRegion;
    



    /**
     * Archiving.
     */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mBendingCoefficient;
        archive & mExponentParameter;
        archive & mDomainWidth;
        archive & mSetNonZeroTargetCurvatureRegion;
        archive & mNonZeroTargetCurvature;
        archive & mRadiusOfNonZeroTargetCurvatureRegion;
        archive & mXCoordinateOfCentreOfNonZeroTargetCurvatureRegion;
    }

public :

    /**
     * Constructor.
     */
    BendingForce2DGhosts();

    /**
     * Destructor.
     */
    ~BendingForce2DGhosts();

    void SetBendingCoefficient(double bendingCoefficient);

    double GetBendingCoefficient();

    void SetExponentParameter(double exponentParameter);

    double GetExponentParameter();

    void SetDomainWidth(double domainWidth);

    double GetDomainWidth();


    void SetCircularNonZeroTargetCurvatureRegion(bool setNonZeroTargetCurvatureRegion, double nonZeroTargetCurvature,
		double radiusOfNonZeroTargetCurvatureRegion, double xCoordinateOfCentreOfNonZeroTargetCurvatureRegion);
    
    c_vector<double, DIM> GetForceDueToDiscreteCurvature(AbstractCellPopulation<DIM>&  rCellPopulation, unsigned cell_i, std::vector<unsigned> first_neighs, std::vector<unsigned> second_neighs, std::vector<c_vector<double, DIM> > imag_loc);

    std::vector<c_vector<double, DIM> > FitPlaneAndFindImage(AbstractCellPopulation<DIM>&  rCellPopulation, std::vector<unsigned> second_neighs, c_vector<double, DIM> cell_i_loc);


    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the tissue
     *
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BendingForce2DGhosts)

#endif /*BendingForce2DGhosts_HPP_*/
