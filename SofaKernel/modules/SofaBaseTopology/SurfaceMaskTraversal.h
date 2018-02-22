#ifndef SOFA_COMPONENT_TOPOLOGY_SURFACEMASKTRAVERSAL_H
#define SOFA_COMPONENT_TOPOLOGY_SURFACEMASKTRAVERSAL_H

#include <sofa/core/objectmodel/BaseObject.h>
#include <SofaBaseTopology/SurfaceTopologyContainer.h>
#include <sofa/core/behavior/MechanicalState.h>

namespace sofa
{

namespace component
{

namespace topology
{

template<class TDataTypes>
class SurfaceMaskTraversal : public virtual core::objectmodel::BaseObject
{

public:
    SOFA_CLASS(SOFA_TEMPLATE(SurfaceMaskTraversal, TDataTypes), core::objectmodel::BaseObject);

    using SurfaceTopology = sofa::component::topology::SurfaceTopologyContainer;
    using Face = SurfaceTopology::Face;
    typedef unsigned int Index;
    typedef typename TDataTypes::VecCoord VecCoord;

protected:
    SurfaceMaskTraversal()
        :BaseObject()
    , mstate(initLink("mstate", "MechanicalState used by this ForceField"))
    {}
	virtual ~SurfaceMaskTraversal()
	{}

public:    
    SurfaceTopology* topology_;
    VecCoord x_;

    SingleLink<SurfaceMaskTraversal<TDataTypes>,sofa::core::behavior::MechanicalState<TDataTypes>,BaseLink::FLAG_STRONGLINK> mstate;

    virtual void init()
    {
        this->getContext()->get(topology_);
        x_ = mstate->read(core::ConstVecCoordId::position())->getValue();
    }


    // Operator() called by the thread
    bool operator() (cgogn::Dart face)
    {
        const auto& dofs = topology_->get_dofs(Face(face));

        Index a = dofs[0];
        Index b = dofs[1];
        Index c = dofs[2];

        return x_[a][2] <= 12. && x_[b][2] <= 12. && x_[c][2] <= 12. ;
    }
};

//template<class TDataTypes>
//class FixedConstraintMask : public SurfaceMaskTraversal
//{

//public:
//    SOFA_CLASS(SOFA_TEMPLATE(FixedConstraintMask, TDataTypes), SOFA_TEMPLATE(SurfaceMaskTraversal, TDataTypes));

//protected:
//	FixedConstraintMask():BaseObject() {}
//	virtual ~FixedConstraintMask()
//	{}

////public:
////	bool select(cgogn::Dart d) override;
//};


} //end namespace topology

} //end namespace component

} //end namespace sofa

#endif // SOFA_COMPONENT_TOPOLOGY_SURFACEMASKTRAVERSAL_H
