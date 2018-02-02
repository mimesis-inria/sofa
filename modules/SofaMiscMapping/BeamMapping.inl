/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2017 INRIA, USTL, UJF, CNRS, MGH                    *
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
#ifndef SOFA_COMPONENT_MAPPING_BEAMMAPPING_INL
#define SOFA_COMPONENT_MAPPING_BEAMMAPPING_INL

#include <SofaMiscMapping/BeamMapping.h>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/helper/io/MassSpringLoader.h>
#include <sofa/helper/io/SphereLoader.h>
#include <sofa/helper/io/Mesh.h>
#include <sofa/helper/gl/template.h>

#include <sofa/simulation/Simulation.h>

#include <string>

namespace sofa
{

namespace component
{

namespace mapping
{


template <class TIn, class TOut>
void BeamMapping<TIn, TOut>::computeMatrix(MatElement& matrixA, Real ksi, Real eta, Real zeta, Real l)
{
    matrixA[0][0] = 1 - ksi;                                  matrixA[1][0] = 0;                                        matrixA[2][0] = 0;
    matrixA[0][1] = 6 * (ksi - ksi*ksi) * eta;                matrixA[1][1] = 1 - 3 * ksi*ksi + 2*ksi*ksi*ksi;          matrixA[2][1] = 0;
    matrixA[0][2] = 6 * (ksi - ksi*ksi) * zeta;               matrixA[1][2] = 0;                                        matrixA[2][2] = 1 - 3 * ksi*ksi + 2*ksi*ksi*ksi;
    matrixA[0][3] = 0;                                        matrixA[1][3] = -(1-ksi) * l * zeta;                      matrixA[2][3] = -(1-ksi) * l * eta;
    matrixA[0][4] = (1 - 4*ksi + 3*ksi*ksi) * l * zeta;       matrixA[1][4] = 0;                                        matrixA[2][4] = (-ksi + 2*ksi*ksi-ksi*ksi*ksi)*l;
    matrixA[0][5] = (-1 + 4 * ksi - 3 * ksi*ksi) * l * eta;   matrixA[1][5] = (ksi - 2 * ksi*ksi + ksi*ksi*ksi) * l;    matrixA[2][5] = 0;
    matrixA[0][6] = ksi;                                      matrixA[1][6] = 0;                                        matrixA[2][6] = 0;
    matrixA[0][7] = 6 * (-ksi + ksi*ksi) * eta;               matrixA[1][7] = 3*ksi*ksi - 2 * ksi*ksi*ksi;              matrixA[2][7] = 0;
    matrixA[0][8] = 6 * (-ksi + ksi*ksi) * zeta;              matrixA[1][8] = 0;                                        matrixA[2][8] = 3*ksi*ksi - 2 * ksi*ksi*ksi;
    matrixA[0][9] = 0;                                        matrixA[1][9] = -l * ksi * zeta;                          matrixA[2][9] = -l * ksi * eta;
    matrixA[0][10] = (-2 * ksi + 3 * ksi*ksi) * l * zeta;     matrixA[1][10] = 0;                                       matrixA[2][10] = (ksi*ksi - ksi*ksi*ksi) * l;
    matrixA[0][11] = (2 * ksi - 3 * ksi*ksi) * l * eta;       matrixA[1][11] = (-ksi*ksi + ksi*ksi*ksi) * l;            matrixA[2][11] = 0;
}

template <class TIn, class TOut>
void BeamMapping<TIn, TOut>::init()
{
    if (this->bary_points.empty() && this->toModel!=NULL)
    {
        const typename In::VecCoord& xfrom = this->fromModel->read(core::ConstVecCoordId::position())->getValue();
        beamLength.resize(xfrom.size() - 1);

        Real totalBeamLength(0);
        for (unsigned int i=0; i<xfrom.size()-1; i++)
        {
            beamLength[i] = (Real)((xfrom[i]-xfrom[i+1]).norm());
            totalBeamLength += beamLength[i];
        }

        const VecCoord& x = this->toModel->read(core::ConstVecCoordId::position())->getValue();
        bary_points.resize(x.size());

        for(unsigned int i = 0 ; i < x.size() ; i++)
        {
            Coord p = x[i] / totalBeamLength; //(ksi, eta, zeta)
            bary_points[i] = p;
        }
    }



    Inherit::init();
}

template <class TIn, class TOut>
void BeamMapping<TIn, TOut>::apply(const core::MechanicalParams * /*mparams*/, Data< typename Out::VecCoord >& _out, const Data< typename In::VecCoord >& _in)
{
    helper::WriteAccessor< Data< typename Out::VecCoord > > out = _out;
    helper::ReadAccessor< Data< typename In::VecCoord > > in = _in;


    for(unsigned int i = 0 ; i < out.size() ; i++)
    {
        computeMatrix(matrixA, bary_points[i][0], bary_points[i][1], bary_points[i][2], beamLength[i]);
        defaulttype::Vec<12, Real> U(in[i][0], in[i][1], in[i][2], in[i][3], in[i][4], in[i][5],
                                    in[i+1][0], in[i+1][1], in[i+1][2], in[i+1][3], in[i+1][4], in[i+1][5]);
        out[i] = matrixA * U;
    }

}

template <class TIn, class TOut>
void BeamMapping<TIn, TOut>::applyJ(const core::MechanicalParams * /*mparams*/, Data< typename Out::VecDeriv >& _out, const Data< typename In::VecDeriv >& _in)
{
    helper::WriteAccessor< Data< typename Out::VecDeriv > > out = _out;
    helper::ReadAccessor< Data< typename In::VecDeriv > > in = _in;

}

template <class TIn, class TOut>
void BeamMapping<TIn, TOut>::applyJT(const core::MechanicalParams * /*mparams*/, Data< typename In::VecDeriv >& _out, const Data< typename Out::VecDeriv >& _in)
{
    helper::WriteAccessor< Data< typename In::VecDeriv > > out = _out;
    helper::ReadAccessor< Data< typename Out::VecDeriv > > in = _in;

}


template <class TIn, class TOut>
void BeamMapping<TIn, TOut>::applyJT(const core::ConstraintParams * /*cparams*/, Data< typename In::MatrixDeriv >& _out, const Data< typename Out::MatrixDeriv >& _in)
{
    typename In::MatrixDeriv* out = _out.beginEdit();
    const typename Out::MatrixDeriv& in = _in.getValue();


}


template <class TIn, class TOut>
void BeamMapping<TIn, TOut>::draw(const core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowMappings()) return;
    std::vector< sofa::defaulttype::Vector3 > points;
//    sofa::defaulttype::Vector3 point;

    const typename In::VecCoord& xfrom = this->fromModel->read(core::ConstVecCoordId::position())->getValue();

    helper::WriteOnlyAccessor<Data<VecCoord>> out(*this->toModel->write(core::VecCoordId::position())); //->getValue();

    for (unsigned int i=0; i<out.size(); i++)
    {
        computeMatrix(matrixA, bary_points[i][0], bary_points[i][1], bary_points[i][2], beamLength[i]);
        defaulttype::Vec<12, Real> U(xfrom[i][0], xfrom[i][1], xfrom[i][2], xfrom[i][3], xfrom[i][4], xfrom[i][5],
                                    xfrom[i+1][0], xfrom[i+1][1], xfrom[i+1][2], xfrom[i+1][3], xfrom[i+1][4], xfrom[i+1][5]);
        out[i] = matrixA * U;
        points.push_back(out[i]);
    }

    vparams->drawTool()->drawPoints(points, 7, sofa::defaulttype::Vec<4,float>(1,1,0,1));
}


template <class TIn, class TOut>
const sofa::defaulttype::BaseMatrix* BeamMapping<TIn, TOut>::getJ()
{

    const unsigned int  inStateSize = this->fromModel->getSize();
//    const unsigned int outStateSize = points.size();

}


} // namespace mapping

} // namespace component

} // namespace sofa

#endif
