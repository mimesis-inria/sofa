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
#ifndef SOFA_COMPONENT_PROJECTIVECONSTRAINTSET_FIXEDCONSTRAINT_INL
#define SOFA_COMPONENT_PROJECTIVECONSTRAINTSET_FIXEDCONSTRAINT_INL

#include <sofa/core/topology/BaseMeshTopology.h>
#include <SofaBoundaryCondition/FixedConstraint.h>
#include <SofaBaseLinearSolver/SparseMatrix.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/Simulation.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <iostream>
#include <SofaBaseTopology/TopologySubsetData.inl>

#include <sofa/core/objectmodel/BaseObject.h>
using sofa::core::objectmodel::ComponentState;


namespace sofa
{

namespace component
{

namespace projectiveconstraintset
{

// Define TestNewPointFunction
template< class DataTypes>
bool FixedConstraint<DataTypes>::FCPointHandler::applyTestCreateFunction(unsigned int, const sofa::helper::vector<unsigned int> &, const sofa::helper::vector<double> &)
{
    if (fc)
    {
        return false;
    }
    else
    {
        return false;
    }
}

// Define RemovalFunction
template< class DataTypes>
void FixedConstraint<DataTypes>::FCPointHandler::applyDestroyFunction(unsigned int pointIndex, value_type &)
{
    if (fc)
    {
        fc->removeConstraint((unsigned int) pointIndex);
    }
}

template <class DataTypes>
FixedConstraint<DataTypes>::FixedConstraint()
    : core::behavior::ProjectiveConstraintSet<DataTypes>(NULL)
    , d_indices( initData(&d_indices,"indices","Indices of the fixed points") )
    , d_fixAll( initData(&d_fixAll,false,"fixAll","filter all the DOF to implement a fixed object") )
    , d_showObject(initData(&d_showObject,true,"showObject","draw or not the fixed constraints"))
    , d_drawSize( initData(&d_drawSize,(SReal)0.0,"drawSize","0 -> point based rendering, >0 -> radius of spheres") )
    , d_projectVelocity( initData(&d_projectVelocity,false,"activate_projectVelocity","activate project velocity to set velocity") )
    , data(new FixedConstraintInternalData<DataTypes>())
{
    // default to indice 0
    d_indices.beginEdit()->push_back(0);
    d_indices.endEdit();

    pointHandler = new FCPointHandler(this, &d_indices);
}


template <class DataTypes>
FixedConstraint<DataTypes>::~FixedConstraint()
{
    if (pointHandler)
        delete pointHandler;

    delete data;
}

template <class DataTypes>
void FixedConstraint<DataTypes>::clearConstraints()
{
    d_indices.beginEdit()->clear();
    d_indices.endEdit();
}

template <class DataTypes>
void FixedConstraint<DataTypes>::addConstraint(unsigned int index)
{
    d_indices.beginEdit()->push_back(index);
    d_indices.endEdit();
}

template <class DataTypes>
void FixedConstraint<DataTypes>::removeConstraint(unsigned int index)
{
    removeValue(*d_indices.beginEdit(),index);
    d_indices.endEdit();
}

// -- Constraint interface


template <class DataTypes>
void FixedConstraint<DataTypes>::init()
{
    this->m_componentstate = ComponentState::Invalid;
    this->core::behavior::ProjectiveConstraintSet<DataTypes>::init();

    if (!this->mstate.get())
    {
        msg_warning() << "Missing mstate, cannot initialize the component.";
        return;
    }

    topology = this->getContext()->getMeshTopology();
    if (!topology)
        msg_warning() << "Can not find the topology, won't be able to handle topological changes";

    // Initialize topological functions
    d_indices.createTopologicalEngine(topology, pointHandler);
    d_indices.registerTopologicalData();

    this->checkIndices();
    this->m_componentstate = ComponentState::Valid;
}

template <class DataTypes>
void  FixedConstraint<DataTypes>::reinit()
{
    this->checkIndices();
}

template <class DataTypes>
void  FixedConstraint<DataTypes>::checkIndices()
{
    // Check value of given indices
    unsigned int maxIndex=this->mstate->getSize();

    const SetIndexArray & indices = d_indices.getValue();
    SetIndexArray invalidIndices;
    for (unsigned int i=0; i<indices.size(); ++i)
    {
        const unsigned int index=indices[i];
        if (index >= maxIndex)
        {
            msg_warning() << "Index " << index << " not valid, should be [0,"<< maxIndex <<"]. Constraint will be removed.";
            invalidIndices.push_back(index);
        }
    }

    // if invalid indices, sort them and remove in decreasing order as removeConstraint perform a swap and pop_back.
    if (!invalidIndices.empty())
    {
        std::sort( invalidIndices.begin(), invalidIndices.end(), std::greater<unsigned int>() );
        int max = invalidIndices.size()-1;
        for (int i=max; i>= 0; i--)
        {
            removeConstraint(invalidIndices[i]);
        }
    }
}

template <class DataTypes>
void FixedConstraint<DataTypes>::projectMatrix( sofa::defaulttype::BaseMatrix* M, unsigned offset )
{
    static const unsigned blockSize = DataTypes::deriv_total_size;

    if( d_fixAll.getValue() )
    {
        unsigned size = this->mstate->getSize();
        for( unsigned i=0; i<size; i++ )
        {
            M->clearRowsCols( offset + i * blockSize, offset + (i+1) * (blockSize) );
        }
    }
    else
    {
        // clears the rows and columns associated with fixed particles
        for(SetIndexArray::const_iterator it= d_indices.getValue().begin(), iend=d_indices.getValue().end(); it!=iend; it++ )
        {
            M->clearRowsCols( offset + (*it) * blockSize, offset + (*it+1) * (blockSize) );
        }
    }
}


template <class DataTypes>
void FixedConstraint<DataTypes>::projectResponse(const core::MechanicalParams* mparams, DataVecDeriv& resData)
{
    helper::WriteAccessor<DataVecDeriv> res ( mparams, resData );
    const SetIndexArray & indices = d_indices.getValue(mparams);

    if( d_fixAll.getValue() )
    {
        // fix everything
        typename VecDeriv::iterator it;
        for( it = res.begin(); it != res.end(); ++it )
        {
            *it = Deriv();
        }
    }
    else
    {
        for (SetIndexArray::const_iterator it = indices.begin(); it != indices.end(); ++it)
        {
            res[*it] = Deriv();
        }
    }
}

template <class DataTypes>
void FixedConstraint<DataTypes>::projectJacobianMatrix(const core::MechanicalParams* mparams, DataMatrixDeriv& cData)
{
    helper::WriteAccessor<DataMatrixDeriv> c ( mparams, cData );
    const SetIndexArray & indices = d_indices.getValue(mparams);

    MatrixDerivRowIterator rowIt = c->begin();
    MatrixDerivRowIterator rowItEnd = c->end();

    if( d_fixAll.getValue() )
    {
        // fix everything
        while (rowIt != rowItEnd)
        {
            rowIt.row().clear();
            ++rowIt;
        }
    }
    else
    {
        while (rowIt != rowItEnd)
        {
            for (SetIndexArray::const_iterator it = indices.begin(); it != indices.end(); ++it)
            {
                rowIt.row().erase(*it);
            }
            ++rowIt;
        }
    }
}

// projectVelocity applies the same changes on velocity vector as projectResponse on position vector :
// Each fixed point received a null velocity vector.
// When a new fixed point is added while its velocity vector is already null, projectVelocity is not usefull.
// But when a new fixed point is added while its velocity vector is not null, it's necessary to fix it to null. If not, the fixed point is going to drift.
template <class DataTypes>
void FixedConstraint<DataTypes>::projectVelocity(const core::MechanicalParams* mparams, DataVecDeriv& vData)
{
    if(!d_projectVelocity.getValue()) return;
    const SetIndexArray & indices = this->d_indices.getValue();
    helper::WriteAccessor<DataVecDeriv> res ( mparams, vData );

    if( d_fixAll.getValue() )    // fix everyting
    {
        for( unsigned i=0; i<res.size(); i++ )
            res[i] = Deriv();
    }
    else
    {
        for (SetIndexArray::const_iterator it = indices.begin(); it != indices.end(); ++it)
        {
            res[*it] = Deriv();
        }
    }
}


template <class DataTypes>
void FixedConstraint<DataTypes>::projectPosition(const core::MechanicalParams* /*mparams*/, DataVecCoord& /*xData*/)
{

}

// Matrix Integration interface
template <class DataTypes>
void FixedConstraint<DataTypes>::applyConstraint(const core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix)
{
    core::behavior::MultiMatrixAccessor::MatrixRef r = matrix->getMatrix(this->mstate.get(mparams));
    if(r)
    {
        const unsigned int N = Deriv::size();

        if( d_fixAll.getValue() )
        {
            unsigned size = this->mstate->getSize();
            for(unsigned int i=0; i<size; i++)
            {
                // Reset Fixed Row and Col
                for (unsigned int c=0; c<N; ++c)
                    r.matrix->clearRowCol(r.offset + N * i + c);
                // Set Fixed Vertex
                for (unsigned int c=0; c<N; ++c)
                    r.matrix->set(r.offset + N * i + c, r.offset + N * i + c, 1.0);
            }
        }
        else
        {
            const SetIndexArray & indices = d_indices.getValue();

            for (SetIndexArray::const_iterator it = indices.begin(); it != indices.end(); ++it)
            {
                // Reset Fixed Row and Col
                for (unsigned int c=0; c<N; ++c)
                    r.matrix->clearRowCol(r.offset + N * (*it) + c);
                // Set Fixed Vertex
                for (unsigned int c=0; c<N; ++c)
                    r.matrix->set(r.offset + N * (*it) + c, r.offset + N * (*it) + c, 1.0);
            }
        }
    }
}

template <class DataTypes>
void FixedConstraint<DataTypes>::applyConstraint(const core::MechanicalParams* mparams, defaulttype::BaseVector* vect, const sofa::core::behavior::MultiMatrixAccessor* matrix)
{
    int o = matrix->getGlobalOffset(this->mstate.get(mparams));
    if (o >= 0)
    {
        unsigned int offset = (unsigned int)o;
        const unsigned int N = Deriv::size();

        if( d_fixAll.getValue() )
        {
            for( int i=0; i<vect->size(); i++ )
            {
                for (unsigned int c=0; c<N; ++c)
                    vect->clear(offset + N * i + c);
            }
        }
        else
        {
            const SetIndexArray & indices = d_indices.getValue();
            for (SetIndexArray::const_iterator it = indices.begin(); it != indices.end(); ++it)
            {
                for (unsigned int c=0; c<N; ++c)
                    vect->clear(offset + N * (*it) + c);
            }
        }
    }
}

template <class DataTypes>
void FixedConstraint<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (this->m_componentstate!=ComponentState::Valid) return;
    if (!vparams->displayFlags().getShowBehaviorModels()) return;
    if (!d_showObject.getValue()) return;
    if (!this->isActive()) return;

    vparams->drawTool()->saveLastState();

    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();
    const SetIndexArray & indices = d_indices.getValue();

    if( d_drawSize.getValue() == 0) // old classical drawing by points
    {
        std::vector< sofa::defaulttype::Vector3 > points;
        sofa::defaulttype::Vector3 point;

        if( d_fixAll.getValue() )
            for (unsigned i=0; i<x.size(); i++ )
            {
                point = DataTypes::getCPos(x[i]);
                points.push_back(point);
            }
        else
        {
            for (SetIndexArray::const_iterator it = indices.begin(); it != indices.end(); ++it)
            {
                point = DataTypes::getCPos(x[*it]);
                points.push_back(point);
            }
        }
        vparams->drawTool()->drawPoints(points, 10, sofa::defaulttype::Vec<4,float>(1,0.5,0.5,1));
    }
    else // new drawing by spheres
    {
        vparams->drawTool()->setLightingEnabled(true);

        std::vector< sofa::defaulttype::Vector3 > points;
        sofa::defaulttype::Vector3 point;
        if( d_fixAll.getValue()==true )
            for (unsigned i=0; i<x.size(); i++ )
            {
                point = DataTypes::getCPos(x[i]);
                points.push_back(point);
            }
        else
        {
            for (SetIndexArray::const_iterator it = indices.begin(); it != indices.end(); ++it)
            {
                point = DataTypes::getCPos(x[*it]);
                points.push_back(point);
            }
        }
        vparams->drawTool()->drawSpheres(points, (float)d_drawSize.getValue(), sofa::defaulttype::Vec<4,float>(1.0f,0.35f,0.35f,1.0f));
    }

    vparams->drawTool()->restoreLastState();
}

} // namespace constraint

} // namespace component

} // namespace sofa

#endif


