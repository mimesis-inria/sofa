/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
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
#pragma once

#include <sofa/core/behavior/BaseConstraintSet.h>
#include <sofa/core/fwd.h>

namespace sofa::core::behavior
{

/**
 *  \brief Component computing constraints within a simulated body.
 *
 *  This class defines the abstract API common to all constraints.
 *  A BaseConstraint computes constraints applied to one or more simulated body
 *  given its current position and velocity.
 *
 *  Constraints can be internal to a given body (attached to one MechanicalState,
 *  see the LagrangianConstraint class), or link several bodies together (such as contacts,
 *  see the InteractionConstraint class).
 *
 */
class SOFA_CORE_API BaseLagrangianConstraint : public BaseConstraintSet
{
public:
    SOFA_ABSTRACT_CLASS(BaseLagrangianConstraint, BaseConstraintSet);
    SOFA_BASE_CAST_IMPLEMENTATION(BaseLagrangianConstraint)

protected:
    BaseLagrangianConstraint() {}
    ~BaseLagrangianConstraint() override {}

private:
    BaseLagrangianConstraint(const BaseLagrangianConstraint& n) = delete ;
    BaseLagrangianConstraint& operator=(const BaseLagrangianConstraint& n) = delete ;

public:
    /// Get the ID of the group containing this constraint. This ID is used to specify which constraints are solved by which solver, by specifying in each solver which groups of constraints it should handle.
    int getGroup() const;

    /// Set the ID of the group containing this constraint. This ID is used to specify which constraints are solved by which solver, by specifying in each solver which groups of constraints it should handle.
    void setGroup(int g);

    typedef long long PersistentID;
    typedef type::vector<PersistentID> VecPersistentID;

    class ConstraintBlockInfo
    {
    public:
        BaseLagrangianConstraint* parent;
        int const0; ///< index of first constraint
        int nbLines; ///< how many dofs (i.e. lines in the matrix) are used by each constraint
        int nbGroups; ///< how many groups of constraints are active
        bool hasId; ///< true if this constraint has persistent ID information
        int offsetId; ///< index of first constraint group info in vector of persistent ids and coordinates
        ConstraintBlockInfo() : parent(nullptr), const0(0), nbLines(1), nbGroups(0), hasId(false), offsetId(0)
        {}
    };
    typedef type::vector<ConstraintBlockInfo> VecConstraintBlockInfo;

    /// Get information for each constraint: pointer to parent BaseConstraint, unique persistent ID, 3D position
    /// \param cParams defines the state vectors to use for positions and velocities. Also defines the order of the constraint (POS, VEL, ACC) and resolution parameters (smoothness, ...)
    virtual void getConstraintInfo(const ConstraintParams* cParams, VecConstraintBlockInfo& blocks, VecPersistentID& ids);


    //DEPRECATED(v25.06, v25.12)
    typedef sofa::type::vector<type::Vec<3,double>> VecConstCoord;
    typedef sofa::type::vector<type::Vec<3,double>> VecConstDeriv;
    typedef sofa::type::vector<double> VecConstArea;
    SOFA_ATTRIBUTE_DEPRECATED__DELETED_ARGUMENTS()
    virtual void getConstraintInfo(const core::ConstraintParams* cParams, VecConstraintBlockInfo& blocks, VecPersistentID& ids, VecConstCoord& positions, VecConstDeriv& directions, VecConstArea& areas) final;

    /// Add the corresponding ConstraintResolution using the offset parameter
    /// \param cParams defines the state vectors to use for positions and velocities. Also defines the order of the constraint (POS, VEL, ACC) and resolution parameters (smoothness, ...)
    /// \param resTab is the result vector that contains the constraint resolution algorithms
    virtual void getConstraintResolution(const ConstraintParams* cParams, std::vector<ConstraintResolution*> &resTab, unsigned int &offset);

    virtual void getConstraintResolution(std::vector<ConstraintResolution*> &resTab, unsigned int &offset);

    type::vector<std::string> getIdentifiers()
    {
        type::vector<std::string> ids = getBaseConstraintIdentifiers();
        ids.push_back("Base");
        return ids;
    }

    virtual type::vector<std::string> getBaseConstraintIdentifiers() = 0;


    /// Store the constraint lambda at the constraint dofs at the given VecDerivId location. 
    /// res = J^t * lambda. 
    /// J is the sparse matrix containing the constraint jacobian that was used to build the constraint matrix ( see BaseConstraintSet::buildConstraintMatrix ).
    /// \param cParams stores the id of the state vectors used during the constraint solving step. Mostly it helps retrieving the MatrixDerivId where
    ///        the constraint jacobian J is stored.
    /// \param res is the state vector Id where to store the result.
    /// \param lambda is the vector of scalar constraint impulses. The direction are stored in the MatrixDerivId stored in the cParams.
    virtual void storeLambda(const ConstraintParams* cParams, MultiVecDerivId res, const sofa::linearalgebra::BaseVector* lambda) = 0;
};

} // namespace sofa::core::behavior
