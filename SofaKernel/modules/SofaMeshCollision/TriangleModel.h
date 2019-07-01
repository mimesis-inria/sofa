/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2019 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_COMPONENT_COLLISION_TRIANGLEMODEL_H
#define SOFA_COMPONENT_COLLISION_TRIANGLEMODEL_H
#include "config.h"

#include <sofa/core/CollisionModel.h>
#include <SofaMeshCollision/LocalMinDistanceFilter.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <SofaBaseTopology/TopologyData.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/defaulttype/VecTypes.h>
#include <SofaMeshCollision/PointModel.h>
#include <map>

namespace sofa
{

namespace component
{

namespace collision
{

template<class DataTypes>
class TriangleCollisionModel;

class TriangleLocalMinDistanceFilter;

template<class TDataTypes>
class TTriangle : public core::TCollisionElementIterator< TriangleCollisionModel<TDataTypes> >
{
public:
    typedef TDataTypes DataTypes;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef TriangleCollisionModel<DataTypes> ParentModel;
	typedef typename DataTypes::Real Real;

    TTriangle(ParentModel* model, int index);
    TTriangle() {}
    explicit TTriangle(const core::CollisionElementIterator& i);
	TTriangle(ParentModel* model, int index, helper::ReadAccessor<typename DataTypes::VecCoord>& /*x*/);

    const Coord& p1() const;
    const Coord& p2() const;
    const Coord& p3() const;

    const Coord& p(int i)const;

    int p1Index() const;
    int p2Index() const;
    int p3Index() const;

    const Coord& p1Free() const;
    const Coord& p2Free() const;
    const Coord& p3Free() const;

    const Coord& operator[](int i) const;

    const Deriv& v1() const;
    const Deriv& v2() const;
    const Deriv& v3() const;
    const Deriv& v(int i) const;


    const Deriv& n() const;
    Deriv& n();

    /// Return true if the element stores a free position vector
    bool hasFreePosition() const;

    int flags() const;

	TTriangle& shape() { return *this; }
    const TTriangle& shape() const { return *this; }

    Coord interpX(defaulttype::Vec<2,Real> bary) const
    {
		return (p1()*(1-bary[0]-bary[1])) + (p2()*bary[0]) + (p3()*bary[1]);
	}
};

/**
 * This class will create collision elements based on a triangle and/or quad mesh.
 * It uses directly the information of the topology and the dof to compute the triangle normals, BB and BoundingTree.
 * The class \sa TTriangle is used to access specific triangle of this collision Model.
 */
template<class TDataTypes>
class SOFA_MESH_COLLISION_API TriangleCollisionModel : public core::CollisionModel
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(TriangleCollisionModel, TDataTypes), core::CollisionModel);

    typedef TDataTypes DataTypes;
    typedef DataTypes InDataTypes;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef TTriangle<DataTypes> Element;
    friend class TTriangle<DataTypes>;

    enum TriangleFlag
    {
        FLAG_P1  = 1<<0, ///< Point 1  is attached to this triangle
        FLAG_P2  = 1<<1, ///< Point 2  is attached to this triangle
        FLAG_P3  = 1<<2, ///< Point 3  is attached to this triangle
        FLAG_E23 = 1<<3, ///< Edge 2-3 is attached to this triangle
        FLAG_E31 = 1<<4, ///< Edge 3-1 is attached to this triangle
        FLAG_E12 = 1<<5, ///< Edge 1-2 is attached to this triangle
        FLAG_BE23 = 1<<6, ///< Edge 2-3 is attached to this triangle and is a boundary
        FLAG_BE31 = 1<<7, ///< Edge 3-1 is attached to this triangle and is a boundary
        FLAG_BE12 = 1<<8, ///< Edge 1-2 is attached to this triangle and is a boundary
        FLAG_POINTS  = FLAG_P1|FLAG_P2|FLAG_P3,
        FLAG_EDGES   = FLAG_E12|FLAG_E23|FLAG_E31,
        FLAG_BEDGES  = FLAG_BE12|FLAG_BE23|FLAG_BE31,
    };

	enum { NBARY = 2 };

    Data<bool> d_bothSide; ///< to activate collision on both side of the triangle model
    Data<bool> d_computeNormals; ///< set to false to disable computation of triangles normal

protected:
    core::behavior::MechanicalState<DataTypes>* m_mstate; ///< Pointer to the corresponding MechanicalState
    sofa::core::topology::BaseMeshTopology* m_topology; ///< Pointer to the corresponding Topology

    VecDeriv m_normals; ///< Vector of normal direction per triangle.

    /** Pointer to the triangle array of this collision model.
     * Will point directly to the topology triangle buffer if only triangles are present. If topology is using/mixing quads and triangles,
     * This pointer will target \sa m_internalTriangles
     * @brief m_triangles
     */
    const sofa::core::topology::BaseMeshTopology::SeqTriangles* m_triangles;

    sofa::core::topology::BaseMeshTopology::SeqTriangles m_internalTriangles; ///< Internal Buffer of triangles to combine quads splitted and other triangles.

    bool m_needsUpdate; ///< parameter storing the info boundingTree has to be recomputed.
    int m_topologyRevision; ///< internal revision number to check if topology has changed.

    PointModel* m_pointModels;

    TriangleLocalMinDistanceFilter *m_lmdFilter;

protected:

    TriangleCollisionModel();

    virtual void updateFromTopology();
    virtual void updateFlags(int ntri=-1);
    virtual void updateNormals();

public:
    void init() override;

    // -- CollisionModel interface

    void resize(int size) override;

    void computeBoundingTree(int maxDepth=0) override;

    void computeContinuousBoundingTree(double dt, int maxDepth=0) override;

    void draw(const core::visual::VisualParams*,int index) override;

    void draw(const core::visual::VisualParams* vparams) override;

    bool canCollideWithElement(int index, CollisionModel* model2, int index2) override;

    core::behavior::MechanicalState<DataTypes>* getMechanicalState() { return m_mstate; }
    const core::behavior::MechanicalState<DataTypes>* getMechanicalState() const { return m_mstate; }

    const VecCoord& getX() const { return(getMechanicalState()->read(core::ConstVecCoordId::position())->getValue()); }
    const sofa::core::topology::BaseMeshTopology::SeqTriangles& getTriangles() const { return *m_triangles; }
    const VecDeriv& getNormals() const { return m_normals; }
    int getTriangleFlags(sofa::core::topology::BaseMeshTopology::TriangleID i);

    TriangleLocalMinDistanceFilter *getFilter() const;

    void setFilter(TriangleLocalMinDistanceFilter * /*lmdFilter*/);

    Deriv velocity(int index)const;


    /// Pre-construction check method called by ObjectFactory.
    /// Check that DataTypes matches the MechanicalState.
    template<class T>
    static bool canCreate(T*& obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
    {
        if (dynamic_cast<core::behavior::MechanicalState<DataTypes>*>(context->getMechanicalState()) == NULL)
            return false;
        return BaseObject::canCreate(obj, context, arg);
    }

    virtual std::string getTemplateName() const override
    {
        return templateName(this);
    }

    static std::string templateName(const TriangleCollisionModel<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

    void computeBBox(const core::ExecParams* params, bool onlyVisible=false) override;
};

template <class TDataTypes> using TTriangleModel [[deprecated("The TTriangleModel is now deprecated please use TriangleCollisionModel instead.")]] = TriangleCollisionModel<TDataTypes>;

template<class DataTypes>
inline TTriangle<DataTypes>::TTriangle(ParentModel* model, int index)
    : core::TCollisionElementIterator<ParentModel>(model, index)
{}

template<class DataTypes>
inline TTriangle<DataTypes>::TTriangle(const core::CollisionElementIterator& i)
    : core::TCollisionElementIterator<ParentModel>(static_cast<ParentModel*>(i.getCollisionModel()), i.getIndex())
{}

template<class DataTypes>
inline TTriangle<DataTypes>::TTriangle(ParentModel* model, int index, helper::ReadAccessor<typename DataTypes::VecCoord>& x)
    : TTriangle(model, index)
{
    SOFA_UNUSED(x);
}

template<class DataTypes>
inline const typename DataTypes::Coord& TTriangle<DataTypes>::p1() const { return this->model->m_mstate->read(core::ConstVecCoordId::position())->getValue()[(*(this->model->m_triangles))[this->index][0]]; }
template<class DataTypes>
inline const typename DataTypes::Coord& TTriangle<DataTypes>::p2() const { return this->model->m_mstate->read(core::ConstVecCoordId::position())->getValue()[(*(this->model->m_triangles))[this->index][1]]; }
template<class DataTypes>
inline const typename DataTypes::Coord& TTriangle<DataTypes>::p3() const { return this->model->m_mstate->read(core::ConstVecCoordId::position())->getValue()[(*(this->model->m_triangles))[this->index][2]]; }
template<class DataTypes>
inline const typename DataTypes::Coord& TTriangle<DataTypes>::p(int i) const {
    return this->model->m_mstate->read(core::ConstVecCoordId::position())->getValue()[(*(this->model->m_triangles))[this->index][i]];
}
template<class DataTypes>
inline const typename DataTypes::Coord& TTriangle<DataTypes>::operator[](int i) const {
    return this->model->m_mstate->read(core::ConstVecCoordId::position())->getValue()[(*(this->model->m_triangles))[this->index][i]];
}

template<class DataTypes>
inline const typename DataTypes::Coord& TTriangle<DataTypes>::p1Free() const { return (*this->model->m_mstate->read(sofa::core::ConstVecCoordId::freePosition())->getValue())[(*(this->model->m_triangles))[this->index][0]]; }
template<class DataTypes>
inline const typename DataTypes::Coord& TTriangle<DataTypes>::p2Free() const { return (*this->model->m_mstate->read(sofa::core::ConstVecCoordId::freePosition())->getValue())[(*(this->model->m_triangles))[this->index][1]]; }
template<class DataTypes>
inline const typename DataTypes::Coord& TTriangle<DataTypes>::p3Free() const { return (*this->model->m_mstate->read(sofa::core::ConstVecCoordId::freePosition())->getValue())[(*(this->model->m_triangles))[this->index][2]]; }

template<class DataTypes>
inline int TTriangle<DataTypes>::p1Index() const { return (*(this->model->m_triangles))[this->index][0]; }
template<class DataTypes>
inline int TTriangle<DataTypes>::p2Index() const { return (*(this->model->m_triangles))[this->index][1]; }
template<class DataTypes>
inline int TTriangle<DataTypes>::p3Index() const { return (*(this->model->m_triangles))[this->index][2]; }

template<class DataTypes>
inline const typename DataTypes::Deriv& TTriangle<DataTypes>::v1() const { return (this->model->m_mstate->read(core::ConstVecDerivId::velocity())->getValue())[(*(this->model->m_triangles))[this->index][0]]; }
template<class DataTypes>
inline const typename DataTypes::Deriv& TTriangle<DataTypes>::v2() const { return this->model->m_mstate->read(core::ConstVecDerivId::velocity())->getValue()[(*(this->model->m_triangles))[this->index][1]]; }
template<class DataTypes>
inline const typename DataTypes::Deriv& TTriangle<DataTypes>::v3() const { return this->model->m_mstate->read(core::ConstVecDerivId::velocity())->getValue()[(*(this->model->m_triangles))[this->index][2]]; }
template<class DataTypes>
inline const typename DataTypes::Deriv& TTriangle<DataTypes>::v(int i) const { return this->model->m_mstate->read(core::ConstVecDerivId::velocity())->getValue()[(*(this->model->m_triangles))[this->index][i]]; }

template<class DataTypes>
inline const typename DataTypes::Deriv& TTriangle<DataTypes>::n() const { return this->model->m_normals[this->index]; }
template<class DataTypes>
inline       typename DataTypes::Deriv& TTriangle<DataTypes>::n()       { return this->model->m_normals[this->index]; }

template<class DataTypes>
inline int TTriangle<DataTypes>::flags() const { return this->model->getTriangleFlags(this->index); }

template<class DataTypes>
inline bool TTriangle<DataTypes>::hasFreePosition() const { return this->model->m_mstate->read(core::ConstVecCoordId::freePosition())->isSet(); }

template<class DataTypes>
inline typename DataTypes::Deriv TriangleCollisionModel<DataTypes>::velocity(int index) const { return (m_mstate->read(core::ConstVecDerivId::velocity())->getValue()[(*(m_triangles))[index][0]] + m_mstate->read(core::ConstVecDerivId::velocity())->getValue()[(*(m_triangles))[index][1]] +
                                                                                                m_mstate->read(core::ConstVecDerivId::velocity())->getValue()[(*(m_triangles))[index][2]])/((Real)(3.0)); }
typedef TriangleCollisionModel<sofa::defaulttype::Vec3Types> TriangleModel;
typedef TTriangle<sofa::defaulttype::Vec3Types> Triangle;

#if  !defined(SOFA_COMPONENT_COLLISION_TRIANGLEMODEL_CPP)
extern template class SOFA_MESH_COLLISION_API TriangleCollisionModel<defaulttype::Vec3Types>;
#endif

} // namespace collision

} // namespace component

} // namespace sofa

#endif
