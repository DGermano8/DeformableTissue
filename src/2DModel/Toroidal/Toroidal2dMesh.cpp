/*

Copyright (c) 2005-2020, University of Oxford.
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

#include <map>

/* These lines are very useful for debugging (visualize with 'showme').
#include "TrianglesMeshWriter.hpp"
TrianglesMeshWriter<2,2> mesh_writer("Toroidal2dMeshDebug", "mesh", false);
mesh_writer.WriteFilesUsingMesh(*this);
*/
#include "Toroidal2dMesh.hpp"
#include "Exception.hpp"
#include "Debug.hpp"

Toroidal2dMesh::Toroidal2dMesh(double width, double depth)
  : MutableMesh<2,2>(),
    mWidth(width),
    mDepth(depth)
{
    assert(width > 0.0);
    assert(depth > 0.0);
}

Toroidal2dMesh::~Toroidal2dMesh()
{
}

Toroidal2dMesh::Toroidal2dMesh(double width, double depth, std::vector<Node<2>* > nodes)
  : MutableMesh<2,2>(),
    mWidth(width),
    mDepth(depth)
{
    assert(width > 0.0);
    assert(depth > 0.0);

//    mMismatchedBoundaryElements = false;
    for (unsigned index=0; index<nodes.size(); index++)
    {
        Node<2>* p_temp_node = nodes[index];
        double x = p_temp_node->rGetLocation()[0];
        UNUSED_OPT(x); // Fix optimised build
        assert( 0 <= x && x < width);
        mNodes.push_back(p_temp_node);
    }

    NodeMap node_map(nodes.size());
    ReMesh(node_map);
}

//bool Toroidal2dMesh::GetInstanceOfMismatchedBoundaryNodes()
//{
//    return mMismatchedBoundaryElements;
//}


void Toroidal2dMesh::CreateMirrorNodes()
{
    double x_half_way = 0.5*mWidth;
    double y_half_way = 0.5*mDepth;

    mLeftOriginals.clear();
    mLeftImages.clear();
    mImageToLeftOriginalNodeMap.clear();
    mLeftPeriodicBoundaryElementIndices.clear();
    
    mRightOriginals.clear();
    mRightImages.clear();
    mImageToRightOriginalNodeMap.clear();
    mRightPeriodicBoundaryElementIndices.clear();

    mTopOriginals.clear();
    mTopImages.clear();
    mImageToTopOriginalNodeMap.clear();
    mTopPeriodicBoundaryElementIndices.clear();

    mBottomOriginals.clear();
    mBottomImages.clear();
    mImageToBottomOriginalNodeMap.clear();
    mBottomPeriodicBoundaryElementIndices.clear();


    for (AbstractMesh<2,2>::NodeIterator node_iter = GetNodeIteratorBegin();
         node_iter != GetNodeIteratorEnd();
         ++node_iter)
    {
        c_vector<double, 2> location;
        location = node_iter->rGetLocation();
        unsigned this_node_index = node_iter->GetIndex();
        double this_node_y_location = location[1];

        // Check the mesh currently conforms to the dimensions given
        assert(0.0 <= this_node_y_location);
        assert(this_node_y_location <= mDepth);

        // Put the nodes which are to be mirrored in the relevant vectors
        if (this_node_y_location < y_half_way)
        {
            mBottomOriginals.push_back(this_node_index);
        }
        else
        {
            mTopOriginals.push_back(this_node_index);
        }
    }

    // For each Bottom original node, create an image node and record its new index
    for (unsigned i=0; i<mBottomOriginals.size(); i++)
    {
        c_vector<double, 2> location;
        location = mNodes[mBottomOriginals[i]]->rGetLocation();
        location[1] = location[1] + mDepth;

        unsigned new_node_index = MutableMesh<2,2>::AddNode(new Node<2>(0, location));
        mBottomImages.push_back(new_node_index);
        mImageToBottomOriginalNodeMap[new_node_index] = mBottomOriginals[i];
    }

    // For each Top original node, create an image node and record its new index
    for (unsigned i=0; i<mTopOriginals.size(); i++)
    {
        c_vector<double, 2> location;
        location = mNodes[mTopOriginals[i]]->rGetLocation();
        location[1] = location[1] - mDepth;

        unsigned new_node_index = MutableMesh<2,2>::AddNode(new Node<2>(0, location));
        mTopImages.push_back(new_node_index);
        mImageToTopOriginalNodeMap[new_node_index] = mTopOriginals[i];
    }

    assert(mTopOriginals.size() == mTopImages.size());
    assert(mBottomOriginals.size() == mBottomImages.size());
    assert(mImageToBottomOriginalNodeMap.size() == mBottomOriginals.size());
    assert(mImageToTopOriginalNodeMap.size() == mTopOriginals.size());


    for (AbstractMesh<2,2>::NodeIterator node_iter = GetNodeIteratorBegin();
         node_iter != GetNodeIteratorEnd();
         ++node_iter)
    {
        c_vector<double, 2> location;
        location = node_iter->rGetLocation();
        unsigned this_node_index = node_iter->GetIndex();
        double this_node_x_location = location[0];

        // Check the mesh currently conforms to the dimensions given
        assert(0.0 <= this_node_x_location);
        assert(this_node_x_location<= mWidth);

        // Put the nodes which are to be mirrored in the relevant vectors
        if (this_node_x_location < x_half_way)
        {
            mLeftOriginals.push_back(this_node_index);
        }
        else
        {
            mRightOriginals.push_back(this_node_index);
        }
    }

    // For each left original node, create an image node and record its new index
    for (unsigned i=0; i<mLeftOriginals.size(); i++)
    {
        c_vector<double, 2> location;
        location = mNodes[mLeftOriginals[i]]->rGetLocation();
        location[0] = location[0] + mWidth;

        unsigned new_node_index = MutableMesh<2,2>::AddNode(new Node<2>(0, location));
        mLeftImages.push_back(new_node_index);
        mImageToLeftOriginalNodeMap[new_node_index] = mLeftOriginals[i];
    }

    // For each right original node, create an image node and record its new index
    for (unsigned i=0; i<mRightOriginals.size(); i++)
    {
        c_vector<double, 2> location;
        location = mNodes[mRightOriginals[i]]->rGetLocation();
        location[0] = location[0] - mWidth;

        unsigned new_node_index = MutableMesh<2,2>::AddNode(new Node<2>(0, location));
        mRightImages.push_back(new_node_index);
        mImageToRightOriginalNodeMap[new_node_index] = mRightOriginals[i];
    }

    assert(mRightOriginals.size() == mRightImages.size());
    assert(mLeftOriginals.size() == mLeftImages.size());
    assert(mImageToLeftOriginalNodeMap.size() == mLeftOriginals.size());
    assert(mImageToRightOriginalNodeMap.size() == mRightOriginals.size());

}

void Toroidal2dMesh::ReMesh(NodeMap& rMap)
{
    unsigned old_num_all_nodes = GetNumAllNodes();

    rMap.Resize(old_num_all_nodes);
    rMap.ResetToIdentity();

    // Flag the deleted nodes as deleted in the map
    for (unsigned i=0; i<old_num_all_nodes; i++)
    {
        if (mNodes[i]->IsDeleted())
        {
            rMap.SetDeleted(i);
        }
    }

  
    // Create mirrored nodes for the normal remesher to work with
    CreateMirrorNodes();

    /*
     * The mesh now has messed up boundary elements, but this
     * doesn't matter as the ReMesh() below doesn't read them in
     * and reconstructs the boundary elements.
     *
     * Call ReMesh() on the parent class. Note that the mesh now has lots
     * of extra nodes which will be deleted, hence the name 'big_map'.
     */
    NodeMap big_map(GetNumAllNodes());
    MutableMesh<2,2>::ReMesh(big_map);

    /*
     * If the big_map isn't the identity map, the little map ('map') needs to be
     * altered accordingly before being passed to the user. Not sure how this all
     * works, so deal with this bridge when we get to it.
     */
    assert(big_map.IsIdentityMap());

    // Re-index the vectors according to the big nodemap, and set up the maps.
    mImageToLeftOriginalNodeMap.clear();
    mImageToRightOriginalNodeMap.clear();

    assert(mLeftOriginals.size() == mLeftImages.size());
    assert(mRightOriginals.size() == mRightImages.size());

    for (unsigned i=0; i<mLeftOriginals.size(); i++)
    {
        mLeftOriginals[i] = big_map.GetNewIndex(mLeftOriginals[i]);
        mLeftImages[i] = big_map.GetNewIndex(mLeftImages[i]);
        mImageToLeftOriginalNodeMap[mLeftImages[i]] = mLeftOriginals[i];
    }

    for (unsigned i=0; i<mRightOriginals.size(); i++)
    {
        mRightOriginals[i] = big_map.GetNewIndex(mRightOriginals[i]);
        mRightImages[i] = big_map.GetNewIndex(mRightImages[i]);
        mImageToRightOriginalNodeMap[mRightImages[i]] = mRightOriginals[i];
    }

    /*
     * Check that elements crossing the periodic boundary have been meshed
     * in the same way at each side.
     */
    CorrectCylindricalNonPeriodicMesh();

    /*
     * Take the double-sized mesh, with its new boundary elements, and
     * remove the relevant nodes, elements and boundary elements to leave
     * a proper periodic mesh.
     */
    ReconstructCylindricalMesh();


    // /*
    //  * Create a random boundary element between two nodes of the first
    //  * element if it is not deleted. This is a temporary measure to get
    //  * around re-index crashing when there are no boundary elements.
    //  */
    // unsigned num_elements = GetNumAllElements();
    // bool boundary_element_made = false;
    // unsigned elem_index = 0;

    // while (elem_index<num_elements && !boundary_element_made)
    // {
    //     Element<2,2>* p_element = GetElement(elem_index);
    //     if (!p_element->IsDeleted())
    //     {
    //         boundary_element_made = true;
    //         std::vector<Node<2>*> nodes;
    //         nodes.push_back(p_element->GetNode(0));
    //         nodes.push_back(p_element->GetNode(1));
    //         BoundaryElement<1,2>* p_boundary_element = new BoundaryElement<1,2>(0, nodes);
    //         p_boundary_element->RegisterWithNodes();
    //         mBoundaryElements.push_back(p_boundary_element);
    //         this->mBoundaryElementWeightedDirections.push_back(zero_vector<double>(2));
    //         this->mBoundaryElementJacobianDeterminants.push_back(0.0);
    //     }
    //     elem_index++;
    // }

    // Now call ReIndex() to remove the temporary nodes which are marked as deleted
    NodeMap reindex_map(GetNumAllNodes());
    ReIndex(reindex_map);
    assert(!reindex_map.IsIdentityMap());  // maybe don't need this

    /*
     * Go through the reindex map and use it to populate the original NodeMap
     * (the one that is returned to the user).
     */
    for (unsigned i=0; i<rMap.GetSize(); i++) // only going up to be size of map, not size of reindex_map
    {
        if (reindex_map.IsDeleted(i))
        {
            /*
             * i < num_original_nodes and node is deleted, this should correspond to
             * a node that was labelled as before the remeshing, so should have already
             * been set as deleted in the map above.
             */
            assert(rMap.IsDeleted(i));
        }
        else
        {
            rMap.SetNewIndex(i, reindex_map.GetNewIndex(i));
        }
    }

    // We can now clear the index vectors and maps; they are only used for remeshing
    mLeftOriginals.clear();
    mLeftImages.clear();
    mImageToLeftOriginalNodeMap.clear();
    mRightOriginals.clear();
    mRightImages.clear();
    mImageToRightOriginalNodeMap.clear();
    mLeftPeriodicBoundaryElementIndices.clear();
    mRightPeriodicBoundaryElementIndices.clear();

    // At this point we have an extended cylindrical mesh now we need to stitch the top and bottom together
    
    // Re-index the vectors according to the big nodemap, and set up the maps.
    mImageToBottomOriginalNodeMap.clear();
    mImageToTopOriginalNodeMap.clear();

    assert(mBottomOriginals.size() == mBottomImages.size());
    assert(mTopOriginals.size() == mTopImages.size());

    for (unsigned i=0; i<mBottomOriginals.size(); i++)
    {
        mBottomOriginals[i] = big_map.GetNewIndex(mBottomOriginals[i]);
        mBottomImages[i] = big_map.GetNewIndex(mBottomImages[i]);
        mImageToBottomOriginalNodeMap[mBottomImages[i]] = mBottomOriginals[i];
    }

    for (unsigned i=0; i<mTopOriginals.size(); i++)
    {
        mTopOriginals[i] = big_map.GetNewIndex(mTopOriginals[i]);
        mTopImages[i] = big_map.GetNewIndex(mTopImages[i]);
        mImageToTopOriginalNodeMap[mTopImages[i]] = mTopOriginals[i];
    }

    /*
     * Check that elements crossing the periodic boundary have been meshed
     * in the same way at each side.
     */
    CorrectToroidalNonPeriodicMesh();

    /*
     * Take the double-sized mesh, with its new boundary elements, and
     * remove the relevant nodes, elements and boundary elements to leave
     * a proper periodic mesh.
     */
    ReconstructToroidalMesh();


    /*
     * Create a random boundary element between two nodes of the first
     * element if it is not deleted. This is a temporary measure to get
     * around re-index crashing when there are no boundary elements.
     */
    unsigned num_elements = GetNumAllElements();
    bool boundary_element_made = false;
    unsigned elem_index = 0;

    while (elem_index<num_elements && !boundary_element_made)
    {
        Element<2,2>* p_element = GetElement(elem_index);
        if (!p_element->IsDeleted())
        {
            boundary_element_made = true;
            std::vector<Node<2>*> nodes;
            nodes.push_back(p_element->GetNode(0));
            nodes.push_back(p_element->GetNode(1));
            BoundaryElement<1,2>* p_boundary_element = new BoundaryElement<1,2>(0, nodes);
            p_boundary_element->RegisterWithNodes();
            mBoundaryElements.push_back(p_boundary_element);
            this->mBoundaryElementWeightedDirections.push_back(zero_vector<double>(2));
            this->mBoundaryElementJacobianDeterminants.push_back(0.0);
        }
        elem_index++;
    }

    // Now call ReIndex() to remove the temporary nodes which are marked as deleted
    NodeMap reindex_map_2(GetNumAllNodes());
    ReIndex(reindex_map_2);
    assert(!reindex_map_2.IsIdentityMap());  // maybe don't need this

    /*
     * Go through the reindex map and use it to populate the original NodeMap
     * (the one that is returned to the user).
     */
    for (unsigned i=0; i<rMap.GetSize(); i++) // only going up to be size of map, not size of reindex_map
    {
        if (reindex_map.IsDeleted(i))
        {
            /*
             * i < num_original_nodes and node is deleted, this should correspond to
             * a node that was labelled as before the remeshing, so should have already
             * been set as deleted in the map above.
             */
            assert(rMap.IsDeleted(i));
        }
        else
        {
            rMap.SetNewIndex(i, reindex_map.GetNewIndex(i));
        }
    }

    // We can now clear the index vectors and maps; they are only used for remeshing
    mBottomOriginals.clear();
    mBottomImages.clear();
    mImageToBottomOriginalNodeMap.clear();
    mTopOriginals.clear();
    mTopImages.clear();
    mImageToTopOriginalNodeMap.clear();
    mBottomPeriodicBoundaryElementIndices.clear();
    mTopPeriodicBoundaryElementIndices.clear();

}

void Toroidal2dMesh::ReconstructCylindricalMesh()
{
    /*
     * Figure out which elements have real nodes and image nodes in them
     * and replace image nodes with corresponding real ones.
     */
    for (MutableMesh<2,2>::ElementIterator elem_iter = GetElementIteratorBegin();
         elem_iter != GetElementIteratorEnd();
         ++elem_iter)
    {
        // Left images are on the right of the mesh
        unsigned number_of_left_image_nodes = 0;
        unsigned number_of_right_image_nodes = 0;

        for (unsigned i=0; i<3; i++)
        {
            unsigned this_node_index = elem_iter->GetNodeGlobalIndex(i);

            if (mImageToLeftOriginalNodeMap.find(this_node_index) != mImageToLeftOriginalNodeMap.end())
            {
                number_of_left_image_nodes++;
            }
            else if (mImageToRightOriginalNodeMap.find(this_node_index) != mImageToRightOriginalNodeMap.end())
            {
                number_of_right_image_nodes++;
            }
        }

        // Delete all the elements on the left hand side (images of right)...
        if (number_of_right_image_nodes >= 1)
        {
            elem_iter->MarkAsDeleted();
            mDeletedElementIndices.push_back(elem_iter->GetIndex());
        }

        // Delete only purely imaginary elements on the right (images of left nodes)
        if (number_of_left_image_nodes == 3)
        {
            elem_iter->MarkAsDeleted();
            mDeletedElementIndices.push_back(elem_iter->GetIndex());
        }

        /*
         * If some are images then replace them with the real nodes. There
         * can be elements with either two image nodes on the right (and one
         * real) or one image node on the right (and two real).
         */
        if (number_of_left_image_nodes == 1 || number_of_left_image_nodes == 2)
        {
            for (unsigned i=0; i<3; i++)
            {
                unsigned this_node_index = elem_iter->GetNodeGlobalIndex(i);
                std::map<unsigned, unsigned>::iterator it = mImageToLeftOriginalNodeMap.find(this_node_index);
                if (it != mImageToLeftOriginalNodeMap.end())
                {
                    elem_iter->ReplaceNode(mNodes[this_node_index], mNodes[it->second]);
                }
            }
        }
    }

    /*
     * Figure out which boundary elements have real nodes and image nodes in them
     * and replace image nodes with corresponding real ones.
     */
    for (unsigned elem_index = 0; elem_index<GetNumAllBoundaryElements(); elem_index++)
    {
        BoundaryElement<1,2>* p_boundary_element = GetBoundaryElement(elem_index);
        if (!p_boundary_element->IsDeleted())
        {
            unsigned number_of_image_nodes = 0;
            for (unsigned i=0; i<2; i++)
            {
                unsigned this_node_index = p_boundary_element->GetNodeGlobalIndex(i);

                if (mImageToLeftOriginalNodeMap.find(this_node_index) != mImageToLeftOriginalNodeMap.end())
                {
                    number_of_image_nodes++;
                }
                else if (mImageToRightOriginalNodeMap.find(this_node_index) != mImageToRightOriginalNodeMap.end())
                {
                    number_of_image_nodes++;
                }
            }

            if (number_of_image_nodes == 2)
            {
                p_boundary_element->MarkAsDeleted();
                mDeletedBoundaryElementIndices.push_back(p_boundary_element->GetIndex());
            }

            /*
             * To avoid having two copies of the boundary elements on the periodic
             * boundaries we only deal with the elements on the left image and
             * delete the ones on the right image.
             */
            if (number_of_image_nodes == 1)
            {
                for (unsigned i=0; i<2; i++)
                {
                    unsigned this_node_index = p_boundary_element->GetNodeGlobalIndex(i);
                    std::map<unsigned, unsigned>::iterator it = mImageToLeftOriginalNodeMap.find(this_node_index);
                    if (it != mImageToLeftOriginalNodeMap.end())
                    {
                        p_boundary_element->ReplaceNode(mNodes[this_node_index], mNodes[it->second]);
                    }
                    else
                    {
                        it = mImageToRightOriginalNodeMap.find(this_node_index);
                        if (it != mImageToRightOriginalNodeMap.end())
                        {
                            p_boundary_element->MarkAsDeleted();
                            mDeletedBoundaryElementIndices.push_back(p_boundary_element->GetIndex());
                        }
                    }
                }
            }
        }
    }

    // Delete all image nodes unless they have already gone 
    for (unsigned i=0; i<mLeftImages.size(); i++)
    {
        mNodes[mLeftImages[i]]->MarkAsDeleted();
        mDeletedNodeIndices.push_back(mLeftImages[i]);
    }

    for (unsigned i=0; i<mRightImages.size(); i++)
    {
        mNodes[mRightImages[i]]->MarkAsDeleted();
        mDeletedNodeIndices.push_back(mRightImages[i]);
    }
}

void Toroidal2dMesh::ReconstructToroidalMesh()
{
    /*
     * Figure out which elements have real nodes and image nodes in them
     * and replace image nodes with corresponding real ones.
     */
    for (MutableMesh<2,2>::ElementIterator elem_iter = GetElementIteratorBegin();
         elem_iter != GetElementIteratorEnd();
         ++elem_iter)
    {
        // Bottom images are on the top of the mesh
        unsigned number_of_bottom_image_nodes = 0;
        unsigned number_of_top_image_nodes = 0;

        for (unsigned i=0; i<3; i++)
        {
            unsigned this_node_index = elem_iter->GetNodeGlobalIndex(i);

            if (mImageToBottomOriginalNodeMap.find(this_node_index) != mImageToBottomOriginalNodeMap.end())
            {
                number_of_bottom_image_nodes++;
            }
            else if (mImageToTopOriginalNodeMap.find(this_node_index) != mImageToTopOriginalNodeMap.end())
            {
                number_of_top_image_nodes++;
            }
        }

        // Delete all the elements on the bottom hand side (images of top)...
        if (number_of_top_image_nodes >= 1)
        {
            elem_iter->MarkAsDeleted();
            mDeletedElementIndices.push_back(elem_iter->GetIndex());
        }

        // Delete only purely imaginary elements on the top (images of bottom nodes)
        if (number_of_bottom_image_nodes == 3)
        {
            elem_iter->MarkAsDeleted();
            mDeletedElementIndices.push_back(elem_iter->GetIndex());
        }

        /*
         * If some are images then replace them with the real nodes. There
         * can be elements with either two image nodes on the top (and one
         * real) or one image node on the top (and two real).
         */
        if (number_of_bottom_image_nodes == 1 || number_of_bottom_image_nodes == 2)
        {
            for (unsigned i=0; i<3; i++)
            {
                unsigned this_node_index = elem_iter->GetNodeGlobalIndex(i);
                std::map<unsigned, unsigned>::iterator it = mImageToBottomOriginalNodeMap.find(this_node_index);
                if (it != mImageToBottomOriginalNodeMap.end())
                {
                    elem_iter->ReplaceNode(mNodes[this_node_index], mNodes[it->second]);
                }
            }
        }
    }

    /*
     * Figure out which boundary elements have real nodes and image nodes in them
     * and replace image nodes with corresponding real ones.
     */
    for (unsigned elem_index = 0; elem_index<GetNumAllBoundaryElements(); elem_index++)
    {
        BoundaryElement<1,2>* p_boundary_element = GetBoundaryElement(elem_index);
        if (!p_boundary_element->IsDeleted())
        {
            unsigned number_of_image_nodes = 0;
            for (unsigned i=0; i<2; i++)
            {
                unsigned this_node_index = p_boundary_element->GetNodeGlobalIndex(i);
                if (mImageToBottomOriginalNodeMap.find(this_node_index) != mImageToBottomOriginalNodeMap.end())
                {
                    number_of_image_nodes++;
                }
                else if (mImageToTopOriginalNodeMap.find(this_node_index) != mImageToTopOriginalNodeMap.end())
                {
                    number_of_image_nodes++;
                }
            }

            if (number_of_image_nodes == 2)
            {
                p_boundary_element->MarkAsDeleted();
                mDeletedBoundaryElementIndices.push_back(p_boundary_element->GetIndex());
            }

            /*
             * To avoid having two copies of the boundary elements on the periodic
             * boundaries we only deal with the elements on the bottom image and
             * delete the ones on the top image.
             */
            if (number_of_image_nodes == 1)
            {
                for (unsigned i=0; i<2; i++)
                {
                    unsigned this_node_index = p_boundary_element->GetNodeGlobalIndex(i);
                    std::map<unsigned, unsigned>::iterator it = mImageToBottomOriginalNodeMap.find(this_node_index);
                    if (it != mImageToBottomOriginalNodeMap.end())
                    {
                        p_boundary_element->ReplaceNode(mNodes[this_node_index], mNodes[it->second]);
                    }
                    else
                    {
                        it = mImageToTopOriginalNodeMap.find(this_node_index);
                        if (it != mImageToTopOriginalNodeMap.end())
                        {
                            p_boundary_element->MarkAsDeleted();
                            mDeletedBoundaryElementIndices.push_back(p_boundary_element->GetIndex());
                        }
                    }
                }
            }
        }
    }

    // Delete all image nodes unless they have already gone 
    for (unsigned i=0; i<mBottomImages.size(); i++)
    {
        mNodes[mBottomImages[i]]->MarkAsDeleted();
        mDeletedNodeIndices.push_back(mBottomImages[i]);
    }

    for (unsigned i=0; i<mTopImages.size(); i++)
    {
        mNodes[mTopImages[i]]->MarkAsDeleted();
        mDeletedNodeIndices.push_back(mTopImages[i]);
    }
}

c_vector<double, 2> Toroidal2dMesh::GetVectorFromAtoB(const c_vector<double, 2>& rLocation1, const c_vector<double, 2>& rLocation2)
{
    // assert(0); //TODO Fix for both periodicity

    assert(mWidth > 0.0);
    assert(mDepth > 0.0);

    c_vector<double, 2> vector = rLocation2 - rLocation1;
    vector[0] = fmod(vector[0], mWidth);
    vector[1] = fmod(vector[1], mDepth);

    /*
     * Handle the toroidal condition here: if the points are more
     * than halfway around the torus apart, measure the other way.
     */
    if (vector[0] > 0.5*mWidth)
    {
        vector[0] -= mWidth;
    }
    else if (vector[0] < -0.5*mWidth)
    {
        vector[0] += mWidth;
    }

    if (vector[1] > 0.5*mDepth)
    {
        vector[1] -= mDepth;
    }
    else if (vector[1] < -0.5*mDepth)
    {
        vector[1] += mDepth;
    }
    return vector;
}

void Toroidal2dMesh::SetNode(unsigned index, ChastePoint<2> point, bool concreteMove)
{
    assert(0); //TODO Fix for both periodicity

    // Perform a periodic movement if necessary
    if (point.rGetLocation()[0] >= mWidth)
    {
        // Move point to the left
        point.SetCoordinate(0, point.rGetLocation()[0] - mWidth);
    }
    else if (point.rGetLocation()[0] < 0.0)
    {
        // Move point to the right
        point.SetCoordinate(0, point.rGetLocation()[0] + mWidth);
    }

    if (point.rGetLocation()[1] >= mDepth)
    {
        // Move point to the left
        point.SetCoordinate(1, point.rGetLocation()[1] - mDepth);
    }
    else if (point.rGetLocation()[1] < 0.0)
    {
        // Move point to the right
        point.SetCoordinate(1, point.rGetLocation()[1] + mDepth);
    }

    // Update the node's location
    MutableMesh<2,2>::SetNode(index, point, concreteMove);
}

double Toroidal2dMesh::GetWidth(const unsigned& rDimension) const
{
    double width = 0.0;
    assert(rDimension==0 || rDimension==1);
    if (rDimension==0)
    {
        width = mWidth;
    }
    else
    {
        width = mDepth;
    }
    return width;
}

unsigned Toroidal2dMesh::AddNode(Node<2>* pNewNode)
{
    unsigned node_index = MutableMesh<2,2>::AddNode(pNewNode);

    // If necessary move it to be back on the cylinder
    ChastePoint<2> new_node_point = pNewNode->GetPoint();
    SetNode(node_index, new_node_point, false);

    return node_index;
}

void Toroidal2dMesh::CorrectCylindricalNonPeriodicMesh()
{
    GenerateVectorsOfElementsStraddlingCylindricalPeriodicBoundaries();

    /*
     * Copy the member variables into new vectors, which we modify
     * by knocking out elements which pair up on each side.
     */
    std::set<unsigned> temp_left_hand_side_elements = mLeftPeriodicBoundaryElementIndices;
    std::set<unsigned> temp_right_hand_side_elements = mRightPeriodicBoundaryElementIndices;

//    if ((mLeftPeriodicBoundaryElementIndices.size()!=mRightPeriodicBoundaryElementIndices.size())
//            || (temp_left_hand_side_elements.size() <= 2)
//            || (temp_right_hand_side_elements.size() <= 2) )
//    {
//        mMismatchedBoundaryElements = true;
//    }
    assert(mLeftPeriodicBoundaryElementIndices.size()==mRightPeriodicBoundaryElementIndices.size());

    // Go through all of the elements on the left periodic boundary
    for (std::set<unsigned>::iterator left_iter = mLeftPeriodicBoundaryElementIndices.begin();
         left_iter != mLeftPeriodicBoundaryElementIndices.end();
         ++left_iter)
    {
        unsigned elem_index = *left_iter;

        Element<2,2>* p_element = GetElement(elem_index);

        /*
         * Make lists of the nodes which the elements on the left contain and
         * the nodes which should be in a corresponding element on the right.
         */
        c_vector<unsigned,3> original_element_node_indices;
        c_vector<unsigned,3> corresponding_element_node_indices;
        for (unsigned i=0; i<3; i++)
        {
            original_element_node_indices[i] = p_element->GetNodeGlobalIndex(i);
            corresponding_element_node_indices[i] = GetCorrespondingCylindricalNodeIndex(original_element_node_indices[i]);
        }

        // Search the right hand side elements for the corresponding element
        for (std::set<unsigned>::iterator right_iter = mRightPeriodicBoundaryElementIndices.begin();
             right_iter != mRightPeriodicBoundaryElementIndices.end();
             ++right_iter)
        {
            unsigned corresponding_elem_index = *right_iter;

            Element<2,2>* p_corresponding_element = GetElement(corresponding_elem_index);

            bool is_corresponding_node = true;

            for (unsigned i=0; i<3; i++)
            {
                if ((corresponding_element_node_indices[i] != p_corresponding_element->GetNodeGlobalIndex(0)) &&
                    (corresponding_element_node_indices[i] != p_corresponding_element->GetNodeGlobalIndex(1)) &&
                    (corresponding_element_node_indices[i] != p_corresponding_element->GetNodeGlobalIndex(2)) )
                {
                    is_corresponding_node = false;
                    break;
                }
            }

            if (is_corresponding_node)
            {
                // Remove original and corresponding element from sets
                temp_left_hand_side_elements.erase(elem_index);
                temp_right_hand_side_elements.erase(corresponding_elem_index);
            }
        }
    }

    /*
     * If either of these ever throw you have more than one situation where the mesher has an option
     * of how to mesh. If it does ever throw you need to be cleverer and match up the
     * elements into as many pairs as possible on the left hand and right hand sides.
     */
 //    assert(temp_left_hand_side_elements.size() <= 2);
 //    assert(temp_right_hand_side_elements.size() <= 2);

    /*
     * Now we just have to use the first pair of elements and copy their info over to the other side.
     * First we need to get hold of both elements on either side.
     */
    if (temp_left_hand_side_elements.empty() || temp_right_hand_side_elements.empty())
    {
        assert(temp_right_hand_side_elements.empty());
        assert(temp_left_hand_side_elements.empty());
    }
    else
    {
        if (temp_right_hand_side_elements.size() == 2 && temp_left_hand_side_elements.size() == 2)
        {
            if (temp_right_hand_side_elements.size() == 2)
            {
                // Use the right hand side meshing and map to left
                UseTheseElementsToDecideCylindricalMeshing(temp_right_hand_side_elements);
            }
            else
            {
                /*
                 * If you get here there are more than two mixed up elements on the periodic edge.
                 * We need to knock the pair out and then rerun this function. This shouldn't be
                 * too hard to do but is as yet unnecessary.
                 */
                NEVER_REACHED;
            }
        }
    }
}

void Toroidal2dMesh::CorrectToroidalNonPeriodicMesh()
{
    GenerateVectorsOfElementsStraddlingToroidalPeriodicBoundaries();

    /*
     * Copy the member variables into new vectors, which we modify
     * by knocking out elements which pair up on each side.
     */
    std::set<unsigned> temp_bottom_hand_side_elements = mBottomPeriodicBoundaryElementIndices;
    std::set<unsigned> temp_top_hand_side_elements = mTopPeriodicBoundaryElementIndices;

//    if ((mBottomPeriodicBoundaryElementIndices.size()!=mTopPeriodicBoundaryElementIndices.size())
//            || (temp_bottom_hand_side_elements.size() <= 2)
//            || (temp_top_hand_side_elements.size() <= 2) )
//    {
//        mMismatchedBoundaryElements = true;
//    }
    assert(mBottomPeriodicBoundaryElementIndices.size()==mTopPeriodicBoundaryElementIndices.size());

    // Go through all of the elements on the bottom periodic boundary
    for (std::set<unsigned>::iterator bottom_iter = mBottomPeriodicBoundaryElementIndices.begin();
         bottom_iter != mBottomPeriodicBoundaryElementIndices.end();
         ++bottom_iter)
    {
        unsigned elem_index = *bottom_iter;

        Element<2,2>* p_element = GetElement(elem_index);

        /*
         * Make lists of the nodes which the elements on the bottom contain and
         * the nodes which should be in a corresponding element on the top.
         */
        c_vector<unsigned,3> original_element_node_indices;
        c_vector<unsigned,3> corresponding_element_node_indices;
        for (unsigned i=0; i<3; i++)
        {
            original_element_node_indices[i] = p_element->GetNodeGlobalIndex(i);
            corresponding_element_node_indices[i] = GetCorrespondingToroidalNodeIndex(original_element_node_indices[i]);
        }

        // Search the top hand side elements for the corresponding element
        for (std::set<unsigned>::iterator top_iter = mTopPeriodicBoundaryElementIndices.begin();
             top_iter != mTopPeriodicBoundaryElementIndices.end();
             ++top_iter)
        {
            unsigned corresponding_elem_index = *top_iter;

            Element<2,2>* p_corresponding_element = GetElement(corresponding_elem_index);

            bool is_corresponding_node = true;

            for (unsigned i=0; i<3; i++)
            {
                if ((corresponding_element_node_indices[i] != p_corresponding_element->GetNodeGlobalIndex(0)) &&
                    (corresponding_element_node_indices[i] != p_corresponding_element->GetNodeGlobalIndex(1)) &&
                    (corresponding_element_node_indices[i] != p_corresponding_element->GetNodeGlobalIndex(2)) )
                {
                    is_corresponding_node = false;
                    break;
                }
            }

            if (is_corresponding_node)
            {
                // Remove original and corresponding element from sets
                temp_bottom_hand_side_elements.erase(elem_index);
                temp_top_hand_side_elements.erase(corresponding_elem_index);
            }
        }
    }

    /*
     * If either of these ever throw you have more than one situation where the mesher has an option
     * of how to mesh. If it does ever throw you need to be cleverer and match up the
     * elements into as many pairs as possible on the left hand and right hand sides.
     */
 //    assert(temp_bottom_hand_side_elements.size() <= 2);
 //    assert(temp_top_hand_side_elements.size() <= 2);

    /*
     * Now we just have to use the first pair of elements and copy their info over to the other side.
     * First we need to get hold of both elements on either side.
     */
    if (temp_bottom_hand_side_elements.empty() || temp_top_hand_side_elements.empty())
    {
        assert(temp_top_hand_side_elements.empty());
        assert(temp_bottom_hand_side_elements.empty());
    }
    else
    {
        if (temp_top_hand_side_elements.size() == 2 && temp_bottom_hand_side_elements.size() == 2)
        {
            if (temp_top_hand_side_elements.size() == 2)
            {
                // Use the top hand side meshing and map to bottom
                UseTheseElementsToDecideToroidalMeshing(temp_top_hand_side_elements);
            }
            else
            {
                /*
                 * If you get here there are more than two mixed up elements on the periodic edge.
                 * We need to knock the pair out and then rerun this function. This shouldn't be
                 * too hard to do but is as yet unnecessary.
                 */
                NEVER_REACHED;
            }
        }
    }
}

void Toroidal2dMesh::UseTheseElementsToDecideCylindricalMeshing(std::set<unsigned>& rMainSideElements)
{
    assert(rMainSideElements.size() == 2);

    // We find the four nodes surrounding the dodgy meshing, on each side.
    std::set<unsigned> main_four_nodes;
    for (std::set<unsigned>::iterator elem_iter = rMainSideElements.begin();
         elem_iter != rMainSideElements.end();
         ++elem_iter)
    {
        unsigned elem_index = *elem_iter;
        Element<2,2>* p_element = GetElement(elem_index);
        for (unsigned i=0; i<3; i++)
        {
            unsigned index = p_element->GetNodeGlobalIndex(i);
            main_four_nodes.insert(index);
        }
    }
    assert(main_four_nodes.size() == 4);

    std::set<unsigned> other_four_nodes;
    for (std::set<unsigned>::iterator iter = main_four_nodes.begin();
         iter != main_four_nodes.end();
         ++iter)
    {
        other_four_nodes.insert(GetCorrespondingCylindricalNodeIndex(*iter));
    }
    assert(other_four_nodes.size() == 4);

    /*
     * Find the elements surrounded by the nodes on the right
     * and change them to match the elements on the left.
     */
    std::vector<unsigned> corresponding_elements;

    // Loop over all elements
    for (MutableMesh<2,2>::ElementIterator elem_iter = GetElementIteratorBegin();
         elem_iter != GetElementIteratorEnd();
         ++elem_iter)
    {
        // Loop over the nodes of the element
        if (!(other_four_nodes.find(elem_iter->GetNodeGlobalIndex(0))==other_four_nodes.end()) &&
            !(other_four_nodes.find(elem_iter->GetNodeGlobalIndex(1))==other_four_nodes.end()) &&
            !(other_four_nodes.find(elem_iter->GetNodeGlobalIndex(2))==other_four_nodes.end()) )
        {
            corresponding_elements.push_back(elem_iter->GetIndex());
            elem_iter->MarkAsDeleted();
            mDeletedElementIndices.push_back(elem_iter->GetIndex());
        }
    }
    assert(corresponding_elements.size() == 2);

    // Now corresponding_elements contains the two elements which are going to be replaced by rMainSideElements
    unsigned num_elements = GetNumAllElements();
    for (std::set<unsigned>::iterator iter = rMainSideElements.begin();
         iter != rMainSideElements.end();
         ++iter)
    {
        Element<2,2>* p_main_element = GetElement(*iter);
        std::vector<Node<2>*> nodes;

        // Put corresponding nodes into a std::vector
        for (unsigned i=0; i<3; i++)
        {
            unsigned main_node = p_main_element->GetNodeGlobalIndex(i);
            nodes.push_back(this->GetNode(GetCorrespondingCylindricalNodeIndex(main_node)));
        }

        // Make a new element
        Element<2,2>* p_new_element = new Element<2,2>(num_elements, nodes);
        this->mElements.push_back(p_new_element);
        this->mElementJacobians.push_back(zero_matrix<double>(2,2));
        this->mElementInverseJacobians.push_back(zero_matrix<double>(2,2));
        this->mElementJacobianDeterminants.push_back(0.0);
        num_elements++;
    }

    // Reindex to get rid of extra elements indices
    NodeMap map(GetNumAllNodes());
    this->ReIndex(map);
}

void Toroidal2dMesh::UseTheseElementsToDecideToroidalMeshing(std::set<unsigned>& rMainSideElements)
{
    assert(rMainSideElements.size() == 2);

    // We find the four nodes surrounding the dodgy meshing, on each side.
    std::set<unsigned> main_four_nodes;
    for (std::set<unsigned>::iterator elem_iter = rMainSideElements.begin();
         elem_iter != rMainSideElements.end();
         ++elem_iter)
    {
        unsigned elem_index = *elem_iter;
        Element<2,2>* p_element = GetElement(elem_index);
        for (unsigned i=0; i<3; i++)
        {
            unsigned index = p_element->GetNodeGlobalIndex(i);
            main_four_nodes.insert(index);
        }
    }
    assert(main_four_nodes.size() == 4);

    std::set<unsigned> other_four_nodes;
    for (std::set<unsigned>::iterator iter = main_four_nodes.begin();
         iter != main_four_nodes.end();
         ++iter)
    {
        other_four_nodes.insert(GetCorrespondingToroidalNodeIndex(*iter));
    }
    assert(other_four_nodes.size() == 4);

    /*
     * Find the elements surrounded by the nodes on the right
     * and change them to match the elements on the left.
     */
    std::vector<unsigned> corresponding_elements;

    // Loop over all elements
    for (MutableMesh<2,2>::ElementIterator elem_iter = GetElementIteratorBegin();
         elem_iter != GetElementIteratorEnd();
         ++elem_iter)
    {
        // Loop over the nodes of the element
        if (!(other_four_nodes.find(elem_iter->GetNodeGlobalIndex(0))==other_four_nodes.end()) &&
            !(other_four_nodes.find(elem_iter->GetNodeGlobalIndex(1))==other_four_nodes.end()) &&
            !(other_four_nodes.find(elem_iter->GetNodeGlobalIndex(2))==other_four_nodes.end()) )
        {
            corresponding_elements.push_back(elem_iter->GetIndex());
            elem_iter->MarkAsDeleted();
            mDeletedElementIndices.push_back(elem_iter->GetIndex());
        }
    }
    assert(corresponding_elements.size() == 2);

    // Now corresponding_elements contains the two elements which are going to be replaced by rMainSideElements
    unsigned num_elements = GetNumAllElements();
    for (std::set<unsigned>::iterator iter = rMainSideElements.begin();
         iter != rMainSideElements.end();
         ++iter)
    {
        Element<2,2>* p_main_element = GetElement(*iter);
        std::vector<Node<2>*> nodes;

        // Put corresponding nodes into a std::vector
        for (unsigned i=0; i<3; i++)
        {
            unsigned main_node = p_main_element->GetNodeGlobalIndex(i);
            nodes.push_back(this->GetNode(GetCorrespondingToroidalNodeIndex(main_node)));
        }

        // Make a new element
        Element<2,2>* p_new_element = new Element<2,2>(num_elements, nodes);
        this->mElements.push_back(p_new_element);
        this->mElementJacobians.push_back(zero_matrix<double>(2,2));
        this->mElementInverseJacobians.push_back(zero_matrix<double>(2,2));
        this->mElementJacobianDeterminants.push_back(0.0);
        num_elements++;
    }

    // Reindex to get rid of extra elements indices
    NodeMap map(GetNumAllNodes());
    this->ReIndex(map);
}

void Toroidal2dMesh::GenerateVectorsOfElementsStraddlingCylindricalPeriodicBoundaries()
{
    mLeftPeriodicBoundaryElementIndices.clear();
    mRightPeriodicBoundaryElementIndices.clear();

    unsigned incidences_of_zero_left_image_nodes = 0;
    unsigned incidences_of_zero_right_image_nodes = 0;

    for (MutableMesh<2,2>::ElementIterator elem_iter = GetElementIteratorBegin();
         elem_iter != GetElementIteratorEnd();
         ++elem_iter)
    {
        // Left images are on the right of the mesh
        unsigned number_of_left_image_nodes = 0;
        unsigned number_of_right_image_nodes = 0;

        for (unsigned i=0; i<3; i++)
        {
            unsigned this_node_index = elem_iter->GetNodeGlobalIndex(i);

            if (mImageToLeftOriginalNodeMap.find(this_node_index) != mImageToLeftOriginalNodeMap.end())
            {
                number_of_left_image_nodes++;
            }
            else if (mImageToRightOriginalNodeMap.find(this_node_index) != mImageToRightOriginalNodeMap.end())
            {
                number_of_right_image_nodes++;
            }
        }

        if ((number_of_left_image_nodes == 0) && (number_of_right_image_nodes == 1 || number_of_right_image_nodes == 2) )
        {
            incidences_of_zero_left_image_nodes++;
        }
        if ((number_of_right_image_nodes == 0) && (number_of_left_image_nodes == 1 || number_of_left_image_nodes == 2) )
        {
            incidences_of_zero_right_image_nodes++;
        }

        /* SJ - Have checked the following:
         * - It is never the case that number_of_left_image_nodes + number_of_right_image_nodes > 3
         * - There are 1264 incidences of zero left image nodes, and 1252 incidences of zero right image nodes
         * - There are 450 incidences of zero left image nodes and non-zero right image nodes; 438 incidences of zero right image nodes and non-zero left
         *   image nodes
         * - There are 40 incidences of 0 left image nodes, and 1/2 right image nodes; 39 incidences of 0 right image nodes and 1/2 left image nodes
         * As this corresponds to 40 left-hand boundary elements and 39 right-hand boundary elements, then it is this extra occurrence of a 0 left
         * image node with 1/2 right image nodes that is responsible for the extra element.
         */

        // Elements on the left hand side (images of right nodes)
        if (number_of_right_image_nodes == 1 || number_of_right_image_nodes == 2)
        {
            mLeftPeriodicBoundaryElementIndices.insert(elem_iter->GetIndex());
        }

        // Elements on the right (images of left nodes)
        if (number_of_left_image_nodes == 1 || number_of_left_image_nodes == 2)
        {
            mRightPeriodicBoundaryElementIndices.insert(elem_iter->GetIndex());
        }
    }

//    if (mLeftPeriodicBoundaryElementIndices.size() != mRightPeriodicBoundaryElementIndices.size())
//    {
//        mMismatchedBoundaryElements = true;
//        // In here - if you hit this case, we want to stop the test and move on, so we work with a stopping event
//
//    }

    // Every boundary element on the left must have a corresponding element on the right
    assert(mLeftPeriodicBoundaryElementIndices.size() == mRightPeriodicBoundaryElementIndices.size());
}

void Toroidal2dMesh::GenerateVectorsOfElementsStraddlingToroidalPeriodicBoundaries()
{
    mBottomPeriodicBoundaryElementIndices.clear();
    mTopPeriodicBoundaryElementIndices.clear();

    unsigned incidences_of_zero_bottom_image_nodes = 0;
    unsigned incidences_of_zero_top_image_nodes = 0;

    for (MutableMesh<2,2>::ElementIterator elem_iter = GetElementIteratorBegin();
         elem_iter != GetElementIteratorEnd();
         ++elem_iter)
    {
        // Left images are on the right of the mesh
        unsigned number_of_bottom_image_nodes = 0;
        unsigned number_of_top_image_nodes = 0;

        for (unsigned i=0; i<3; i++)
        {
            unsigned this_node_index = elem_iter->GetNodeGlobalIndex(i);

            if (mImageToBottomOriginalNodeMap.find(this_node_index) != mImageToBottomOriginalNodeMap.end())
            {
                number_of_bottom_image_nodes++;
            }
            else if (mImageToTopOriginalNodeMap.find(this_node_index) != mImageToTopOriginalNodeMap.end())
            {
                number_of_top_image_nodes++;
            }
        }

        if ((number_of_bottom_image_nodes == 0) && (number_of_top_image_nodes == 1 || number_of_top_image_nodes == 2) )
        {
            incidences_of_zero_bottom_image_nodes++;
        }
        if ((number_of_top_image_nodes == 0) && (number_of_bottom_image_nodes == 1 || number_of_bottom_image_nodes == 2) )
        {
            incidences_of_zero_top_image_nodes++;
        }
        
        // Elements on the bottom hand side (images of top nodes)
        if (number_of_top_image_nodes == 1 || number_of_top_image_nodes == 2)
        {
            mBottomPeriodicBoundaryElementIndices.insert(elem_iter->GetIndex());
        }

        // Elements on the right (images of bottom nodes)
        if (number_of_bottom_image_nodes == 1 || number_of_bottom_image_nodes == 2)
        {
            mTopPeriodicBoundaryElementIndices.insert(elem_iter->GetIndex());
        }
    }

//    if (mBottomPeriodicBoundaryElementIndices.size() != mTopPeriodicBoundaryElementIndices.size())
//    {
//        mMismatchedBoundaryElements = true;
//        // In here - if you hit this case, we want to stop the test and move on, so we work with a stopping event
//
//    }

    // Every boundary element on the left must have a corresponding element on the right
    assert(mBottomPeriodicBoundaryElementIndices.size() == mTopPeriodicBoundaryElementIndices.size());
}

unsigned Toroidal2dMesh::GetCorrespondingCylindricalNodeIndex(unsigned nodeIndex)
{
    unsigned corresponding_node_index = UINT_MAX;

    // If nodeIndex is a member of mRightOriginals, then find the corresponding node index in mRightImages
    std::vector<unsigned>::iterator right_orig_iter = std::find(mRightOriginals.begin(), mRightOriginals.end(), nodeIndex);
    if (right_orig_iter != mRightOriginals.end())
    {
        corresponding_node_index = mRightImages[right_orig_iter - mRightOriginals.begin()];
    }
    else
    {
        // If nodeIndex is a member of mRightImages, then find the corresponding node index in mRightOriginals
        std::vector<unsigned>::iterator right_im_iter = std::find(mRightImages.begin(), mRightImages.end(), nodeIndex);
        if (right_im_iter != mRightImages.end())
        {
            corresponding_node_index = mRightOriginals[right_im_iter - mRightImages.begin()];
        }
        else
        {
            // If nodeIndex is a member of mLeftOriginals, then find the corresponding node index in mLeftImages
            std::vector<unsigned>::iterator left_orig_iter = std::find(mLeftOriginals.begin(), mLeftOriginals.end(), nodeIndex);
            if (left_orig_iter != mLeftOriginals.end())
            {
                corresponding_node_index = mLeftImages[left_orig_iter - mLeftOriginals.begin()];
            }
            else
            {
                // If nodeIndex is a member of mLeftImages, then find the corresponding node index in mLeftOriginals
                std::vector<unsigned>::iterator left_im_iter = std::find(mLeftImages.begin(), mLeftImages.end(), nodeIndex);
                if (left_im_iter != mLeftImages.end())
                {
                    corresponding_node_index = mLeftOriginals[left_im_iter - mLeftImages.begin()];
                }
            }
        }
    }

    // We must have found the corresponding node index
    assert(corresponding_node_index != UINT_MAX);
    return corresponding_node_index;
}


unsigned Toroidal2dMesh::GetCorrespondingToroidalNodeIndex(unsigned nodeIndex)
{
    unsigned corresponding_node_index = UINT_MAX;

    // If nodeIndex is a member of mTopOriginals, then find the corresponding node index in mTopImages
    std::vector<unsigned>::iterator top_orig_iter = std::find(mTopOriginals.begin(), mTopOriginals.end(), nodeIndex);
    if (top_orig_iter != mTopOriginals.end())
    {
        corresponding_node_index = mTopImages[top_orig_iter - mTopOriginals.begin()];
    }
    else
    {
        // If nodeIndex is a member of mTopImages, then find the corresponding node index in mTopOriginals
        std::vector<unsigned>::iterator top_im_iter = std::find(mTopImages.begin(), mTopImages.end(), nodeIndex);
        if (top_im_iter != mTopImages.end())
        {
            corresponding_node_index = mTopOriginals[top_im_iter - mTopImages.begin()];
        }
        else
        {
            // If nodeIndex is a member of mBottomOriginals, then find the corresponding node index in mBottomImages
            std::vector<unsigned>::iterator bottom_orig_iter = std::find(mBottomOriginals.begin(), mBottomOriginals.end(), nodeIndex);
            if (bottom_orig_iter != mBottomOriginals.end())
            {
                corresponding_node_index = mBottomImages[bottom_orig_iter - mBottomOriginals.begin()];
            }
            else
            {
                // If nodeIndex is a member of mBottomImages, then find the corresponding node index in mBottomOriginals
                std::vector<unsigned>::iterator bottom_im_iter = std::find(mBottomImages.begin(), mBottomImages.end(), nodeIndex);
                if (bottom_im_iter != mBottomImages.end())
                {
                    corresponding_node_index = mBottomOriginals[bottom_im_iter - mBottomImages.begin()];
                }
            }
        }
    }

    // We must have found the corresponding node index
    assert(corresponding_node_index != UINT_MAX);
    return corresponding_node_index;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Toroidal2dMesh)
