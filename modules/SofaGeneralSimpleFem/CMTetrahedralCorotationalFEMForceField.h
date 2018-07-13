/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2016 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_COMPONENT_FORCEFIELD_CMTETRAHEDRALCOROTATIONALFEMFORCEFIELD_H
#define SOFA_COMPONENT_FORCEFIELD_CMTETRAHEDRALCOROTATIONALFEMFORCEFIELD_H
#include "config.h"

#include <sofa/core/behavior/ForceField.h>
#include <SofaBaseTopology/VolumeTopologyContainer.h>
#include <SofaBaseTopology/CMTopologyData.inl>
#include <sofa/helper/vector.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/helper/map.h>

// corotational tetrahedron from
// @InProceedings{NPF05,
//   author       = "Nesme, Matthieu and Payan, Yohan and Faure, Fran\c{c}ois",
//   title        = "Efficient, Physically Plausible Finite Elements",
//   booktitle    = "Eurographics (short papers)",
//   month        = "august",
//   year         = "2005",
//   editor       = "J. Dingliana and F. Ganovelli",
//   keywords     = "animation, physical model, elasticity, finite elements",
//   url          = "http://www-evasion.imag.fr/Publications/2005/NPF05"
// }


namespace sofa
{

namespace component
{

namespace cm_forcefield
{


/** Compute Finite Element forces based on tetrahedral elements.
 */
template<class DataTypes>
class CMTetrahedralCorotationalFEMForceField : public core::behavior::ForceField<DataTypes>
{
public:
	SOFA_CLASS(SOFA_TEMPLATE(CMTetrahedralCorotationalFEMForceField, DataTypes), SOFA_TEMPLATE(core::behavior::ForceField, DataTypes));

	typedef unsigned int Index;
	typedef typename DataTypes::VecCoord VecCoord;
	typedef typename DataTypes::VecDeriv VecDeriv;
	typedef typename DataTypes::VecReal VecReal;
	typedef VecCoord Vector;
	typedef typename DataTypes::Coord Coord;
	typedef typename DataTypes::Deriv Deriv;
	typedef typename Coord::value_type Real;
	using VolumeTopology = sofa::component::topology::VolumeTopologyContainer;
	using BaseVolume = core::topology::MapTopology::Volume;
	using Vertex = VolumeTopology::Vertex;
	using Volume = VolumeTopology::Volume;

	typedef core::objectmodel::Data<VecDeriv>    DataVecDeriv;
	typedef core::objectmodel::Data<VecCoord>    DataVecCoord;

	enum
	{
		SMALL = 0, // Symbol of small displacements tetrahedron solver
		LARGE = 1, // Symbol of large displacements tetrahedron solver
		POLAR = 2,  // Symbol of polar displacements tetrahedron solver
		PLARGE = 3 // Symbol of large displacements tetrahedron solver with parallel addForce
	};

protected:

	/// @name Per element (tetrahedron) data
	/// @{

	/// Displacement vector (deformation of the 4 corners of a tetrahedron
	typedef defaulttype::VecNoInit<12, Real> Displacement;

	/// Material stiffness matrix of a tetrahedron
	typedef defaulttype::Mat<6, 6, Real> MaterialStiffness;

	/// Strain-displacement matrix
	typedef defaulttype::Mat<12, 6, Real> StrainDisplacementTransposed;

	/// Rigid transformation (rotation) matrix
	typedef defaulttype::MatNoInit<3, 3, Real> Transformation;

	/// Stiffness matrix ( = RJKJtRt  with K the Material stiffness matrix, J the strain-displacement matrix, and R the transformation matrix if any )
	typedef defaulttype::Mat<12, 12, Real> StiffnessMatrix;

	/// @}

	/// The information stored for each tetrahedron
	class TetrahedronInformation
	{
	public:
		/// the indices of the DOF (vertices) of the tetrahedron
		helper::fixed_array<Index, 4> dofs;

		/// for all methods
		MaterialStiffness materialMatrix;
		StrainDisplacementTransposed strainDisplacementTransposedMatrix;

		/// for the large displacement method
		helper::fixed_array<Coord,4> rotatedInitialElements;
		Transformation rotation;

		/// for the polar method
		Transformation initialTransformation;

		TetrahedronInformation()
		{
		}

		/// Output stream
		inline friend std::ostream& operator<< ( std::ostream& os, const TetrahedronInformation& /*tri*/ )
		{
			return os;
		}

		/// Input stream
		inline friend std::istream& operator>> ( std::istream& in, TetrahedronInformation& /*tri*/ )
		{
			return in;
		}
	};
	/// container that stotes all requires information for each tetrahedron
	cm_topology::TetrahedronData<TetrahedronInformation> _volumeAttribute;

	/// @name Full system matrix assembly support
	/// @{

	typedef std::pair<int,Real> Col_Value;
	typedef helper::vector< Col_Value > CompressedValue;
	typedef helper::vector< CompressedValue > CompressedMatrix;

	CompressedMatrix _stiffnesses;
	/// @}

	SReal m_potentialEnergy;

	VolumeTopology* _topology;
public:
	class TetrahedronHandler : public cm_topology::TopologyDataHandler<core::topology::MapTopology::Volume, TetrahedronInformation>
	{
	public :
		typedef typename CMTetrahedralCorotationalFEMForceField<DataTypes>::TetrahedronInformation TetrahedronInformation;
		TetrahedronHandler(CMTetrahedralCorotationalFEMForceField<DataTypes>* forcefield,
						   cm_topology::TetrahedronData<TetrahedronInformation>* data)
			: cm_topology::TopologyDataHandler<core::topology::MapTopology::Volume, TetrahedronInformation>(data)
			,ff(forcefield)
		{

		}

		void applyCreateFunction(TetrahedronInformation &t, core::topology::MapTopology::Volume,
				const sofa::helper::vector<core::topology::MapTopology::Volume> &,
				const sofa::helper::vector<double> &);

	protected:
		CMTetrahedralCorotationalFEMForceField<DataTypes>* ff;

	};
public:
	int method;
	Data<std::string> f_method; ///< the computation method of the displacements
	Data<Real> _poissonRatio;
	Data<Real> _youngModulus;
	Data<VecReal> _localStiffnessFactor;
	Data<bool> _updateStiffnessMatrix;
	Data<bool> _assembling;
	Data<bool> f_drawing;
	Data<bool> _displayWholeVolume;
	Data<defaulttype::Vec4f> drawColor1;
	Data<defaulttype::Vec4f> drawColor2;
	Data<defaulttype::Vec4f> drawColor3;
	Data<defaulttype::Vec4f> drawColor4;
	Data<std::map < std::string, sofa::helper::vector<double> > > _volumeGraph;
protected:
	CMTetrahedralCorotationalFEMForceField();
	TetrahedronHandler* tetrahedronHandler;
public:

	void setPoissonRatio(Real val) { this->_poissonRatio.setValue(val); }

	void setYoungModulus(Real val) { this->_youngModulus.setValue(val); }

	void setMethod(int val) { method = val; }

	void setUpdateStiffnessMatrix(bool val) { this->_updateStiffnessMatrix.setValue(val); }

	void setComputeGlobalMatrix(bool val) { this->_assembling.setValue(val); }

	virtual void init();
	virtual void reinit();

	virtual void addForce(const core::MechanicalParams* mparams, DataVecDeriv& d_f, const DataVecCoord& d_x, const DataVecDeriv& d_v);
	virtual void addDForce(const core::MechanicalParams* mparams, DataVecDeriv& d_df, const DataVecDeriv& d_dx);
	virtual SReal getPotentialEnergy(const core::MechanicalParams* /*mparams*/, const DataVecCoord&  /* x */) const
	{
		serr << "Get potentialEnergy not implemented" << sendl;
		return 0.0;
	}

	virtual void addKToMatrix(sofa::defaulttype::BaseMatrix *m, SReal kFactor, unsigned int &offset);

	// Getting the rotation of the vertex by averaing the rotation of neighboring elements
	void getRotation(Transformation& R, Vertex vertex);
	void getRotations() {}
	void getElementRotation(Transformation& R, unsigned int elementIdx);

	// Getting the stiffness matrix of index i
	void getElementStiffnessMatrix(Real* stiffness, Volume w);
	void getElementStiffnessMatrix(Real* stiffness, const core::topology::BaseMeshTopology::Tetrahedron& t);

	void draw(const core::visual::VisualParams* vparams);

protected:

	void computeStrainDisplacement( StrainDisplacementTransposed &J, Coord a, Coord b, Coord c, Coord d );
	Real peudo_determinant_for_coef ( const defaulttype::Mat<2, 3, Real>&  M );

	void computeStiffnessMatrix( StiffnessMatrix& S,StiffnessMatrix& SR,const MaterialStiffness &K, const StrainDisplacementTransposed &J, const Transformation& Rot );

	void computeMaterialStiffness(TetrahedronInformation& info);

	/// overloaded by classes with non-uniform stiffness
	virtual void computeMaterialStiffness(MaterialStiffness& materialMatrix, Index&a, Index&b, Index&c, Index&d, SReal localStiffnessFactor=1);

	void computeForce( Displacement &F, const Displacement &Depl, const MaterialStiffness &K, const StrainDisplacementTransposed &J );
	void computeForce( Displacement &F, const Displacement &Depl, const MaterialStiffness &K, const StrainDisplacementTransposed &J, SReal fact );

	////////////// small displacements method
	void initSmall(const VecCoord& X0, TetrahedronInformation& info);
	void accumulateForceSmall(Vector& f, const Vector & p, TetrahedronInformation& info);
	void applyStiffnessSmall(Vector& f, const Vector& x, const TetrahedronInformation& info, SReal fact=1.0 );

	////////////// large displacements method
	void initLarge(const VecCoord& X0, TetrahedronInformation& info);
	void computeRotationLarge( Transformation &r, const Vector &p, const Index &a, const Index &b, const Index &c);
	void accumulateForceLarge(Vector& f, const Vector & p, TetrahedronInformation& info);
	void applyStiffnessLarge(Vector& f, const Vector& x, const TetrahedronInformation& info, SReal fact=1.0 );

	////////////// polar decomposition method
	void initPolar(const VecCoord& X0, TetrahedronInformation& info);
	void accumulateForcePolar(Vector& f, const Vector& p, TetrahedronInformation& info);
	void applyStiffnessPolar(Vector& f, const Vector& x, const TetrahedronInformation& info, SReal fact=1.0 );

	void printStiffnessMatrix(Volume idTetra);

public:
	virtual std::string getClassName() const { return "CMTetrahedralCorotationalFEMForceField"; }

	static std::string className(const CMTetrahedralCorotationalFEMForceField* /* ptr = nullptr */)
	{
	  return "CMTetrahedralCorotationalFEMForceField";
	}
};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_FORCEFIELD_CMTETRAHEDRALCOROTATIONALFEMFORCEFIELD_CPP)

#ifndef SOFA_FLOAT
extern template class SOFA_GENERAL_SIMPLE_FEM_API CMTetrahedralCorotationalFEMForceField<sofa::defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_GENERAL_SIMPLE_FEM_API CMTetrahedralCorotationalFEMForceField<sofa::defaulttype::Vec3fTypes>;
#endif

#endif // defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_FORCEFIELD_CMTETRAHEDRALCOROTATIONALFEMFORCEFIELD_CPP)


} // namespace cm_forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_CMTETRAHEDRALCOROTATIONALFEMFORCEFIELD_H
