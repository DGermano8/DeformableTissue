/*
*/

#ifndef HistoryDepMutableMesh_HPP_
#define HistoryDepMutableMesh_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>

#include "NodeMap.hpp"
#include "TetrahedralMesh.hpp"
#include "MutableMesh.hpp"
/**
 * A concrete mutable mesh class.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class HistoryDepMutableMesh : public MutableMesh<ELEMENT_DIM, SPACE_DIM>
{
    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Save the mesh, along with attributes if nodes have attributes.
     * Saving of attributes is covered in TestNodesOnlyMesh.
     * @param archive the archive to save to.
     * @param version the version number.
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        archive & boost::serialization::base_object<MutableMesh<ELEMENT_DIM, SPACE_DIM> >(*this);

        // // Assume that the first node is indicative of the rest.
        // bool does_have_attributes = this->mNodes[0]->HasNodeAttributes();

        // archive & does_have_attributes;

        // if (does_have_attributes)
        // {
        //     for (unsigned i=0; i<this->mNodes.size(); i++)
        //     {
        //         double radius;
        //         radius = this->mNodes[i]->GetRadius();
        //         archive & radius;

        //         bool is_particle;
        //         is_particle = this->mNodes[i]->IsParticle();
        //         archive & is_particle;
        //     }
        // }
    }

    


    /**
     * Load the mesh, along with attributes if saved nodes have attributes.
     * Loading of attributes is covered in TestNodesOnlyMesh.
     * @param archive the archive to save to.
     * @param version the version number.
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<MutableMesh<ELEMENT_DIM, SPACE_DIM> >(*this);

        bool does_have_attributes;

        archive & does_have_attributes;

        if (does_have_attributes)
        {
            for (unsigned i=0; i<this->mNodes.size(); i++)
            {
                double radius;
                archive & radius;
                this->mNodes[i]->SetRadius(radius);

                bool is_particle;
                archive & is_particle;

                if (is_particle)
                {
                    this->mNodes[i]->SetIsParticle(true);
                }
            }
        }

        // If ELEMENT_DIM==SPACE_DIM do a remesh after archiving has finished to get right number of boundary nodes etc.
        // NOTE - Subclasses must archive their member variables BEFORE calling this method.
        if (ELEMENT_DIM == SPACE_DIM)
        {
            NodeMap map(this->GetNumNodes());
            this->ReMesh(map);
            assert(map.IsIdentityMap()); // Otherwise the mesh will get VERY confused.
        }
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()



public:

    /**
     * Constructor.
     */
    HistoryDepMutableMesh();

    /**
     * Constructor which takes in a vector of nodes.
     *
     * @param nodes  a vector of nodes
     */
    HistoryDepMutableMesh(std::vector<Node<SPACE_DIM> *> nodes);

    void DeleteMesh();
    void AssignNewMesh(HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>* New_Mesh);
    void AddANewNodeBehindBoundary();
    void CreateNewMesh(HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>* New_Mesh, std::map<unsigned, c_vector<double, SPACE_DIM> > mInitalPositionOfRemeshedNodes);
    

    /**
     * Destructor.
     */
    virtual ~HistoryDepMutableMesh();
};


#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(HistoryDepMutableMesh)

#endif /*HistoryDepMutableMesh_HPP_*/

