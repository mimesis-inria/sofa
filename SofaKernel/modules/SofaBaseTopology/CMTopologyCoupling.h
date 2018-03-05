#ifndef SOFA_COMPONENT_TOPOLOGY_CMTOPOLOGYCOUPLING_H
#define SOFA_COMPONENT_TOPOLOGY_CMTOPOLOGYCOUPLING_H

#include <sofa/core/objectmodel/BaseObject.h>
#include <SofaBaseTopology/SurfaceTopologyContainer.h>


namespace sofa
{

namespace component
{

namespace topology
{


class CMTopologyCoupling : public virtual core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(CMTopologyCoupling, BaseObject);

    using SurfaceTopology = sofa::component::topology::SurfaceTopologyContainer;
    using Face = SurfaceTopology::Face;

protected:
    FixedConstraintMask()
        :BaseObject()
    {}

public:
    virtual void init() override
    {
        this->getContext()->get(topology_);
//        this->getContext()->get(grah_);

    }

   void computeCoupling()
   {
       topology_->foreach_cell([&](Face f)
       {
//            graph_->foreach_cell([&](Node n)
//            {
                //compute m_f the center of f = 1/card(f) * (sum (f_i))

                //compute the cost(m_f, n) = ||m_f, n ||

                //save n in f if costs are minimal up to now
//            });
       });
   }

private:
    SurfaceTopology* topology_;
    //GraphTopology* graph_;
};

} //end namespace topology

} //end namespace component

} //end namespace sofa

#endif // SOFA_COMPONENT_TOPOLOGY_SURFACEMASKTRAVERSAL_H
