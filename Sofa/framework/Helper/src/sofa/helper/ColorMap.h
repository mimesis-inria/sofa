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
#ifndef SOFA_HELPER_COLORMAP_H
#define SOFA_HELPER_COLORMAP_H

#include <sofa/helper/config.h>

#include <sofa/type/vector.h>
#include <sofa/type/RGBAColor.h>
#include <sofa/helper/rmath.h>
#include <sofa/type/Vec.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
//#include <sofa/helper/OptionsGroup.h>


namespace sofa::helper
{
    
class SOFA_HELPER_API ColorMap 
{
public:
    typedef sofa::type::vector<type::RGBAColor> VecColor;
    
    ColorMap(unsigned int paletteSize = 256, const std::string& colorScheme = "HSV");
    virtual ~ColorMap();

    template<class Real>
    class evaluator
    {
    public:
        evaluator()
            : map(nullptr), vmin(0), vmax(0), vscale(0)
        {}

        evaluator(const ColorMap* map, Real vmin, Real vmax)
            : map(map), vmin(vmin), vmax(vmax),
              // An exact vmax == vmin test is fragile: two values meant to
              // represent "no variation" (e.g. a field that is uniformly
              // zero) routinely differ by more than a couple of ULPs once
              // they've gone through independent floating-point paths (e.g.
              // mapping between two meshes that represent the same geometry
              // but aren't bit-identical -- observed noise on one real case
              // was ~4 ULPs, i.e. a tight few-epsilon tolerance isn't safe
              // margin), leaving a near-zero but nonzero vmax-vmin. That
              // produces a huge vscale, which amplifies that noise across
              // the whole palette instead of the flat color a genuinely-
              // uniform field should get. Tolerate a relative (scaled by
              // the larger endpoint magnitude, floored at 1 so near-zero
              // ranges aren't over-sensitive) difference of up to 1e4
              // epsilons as "equal" instead.
              vscale((std::abs(vmax - vmin) <= Real(1e4) * std::numeric_limits<Real>::epsilon()
                          * std::max({Real(1), std::abs(vmax), std::abs(vmin)}))
                         ? (Real)0
                         : (map->entries.size()-1)/(vmax-vmin)) {}

        auto operator()(Real r) const
        {
            Real e = (r-vmin)*vscale;
            if (e<0) return map->entries.front();

            unsigned int i = (unsigned int)(e);
            if (i>=map->entries.size()-1) return map->entries.back();

            const auto& c1 = map->entries[i];
            const auto& c2 = map->entries[i+1];
            return c1+(c2-c1)*(e-i);
        }
    protected:
        const ColorMap* map;
        Real vmin;
        Real vmax;
        Real vscale;
    };

    unsigned int m_paletteSize;
    std::string m_colorScheme;
    
    VecColor entries;
    
    void init();
    void reinit();
    
    unsigned int getPaletteSize() const { return m_paletteSize;  }
    void setPaletteSize(unsigned int paletteSize) { m_paletteSize = paletteSize; }

    const std::string& getColorScheme() const { return m_colorScheme;  }
    void setColorScheme(const std::string& colorScheme) { m_colorScheme = colorScheme; }

    unsigned int getNbColors() const { return (unsigned int) entries.size(); }
    type::RGBAColor getColor(unsigned int i) {
        if (i < entries.size()) return entries[i];
        return type::RGBAColor(0.f, 0.f, 0.f, 0.f);
    }

    static ColorMap* getDefault();

    template<class Real>
    evaluator<Real> getEvaluator(Real vmin, Real vmax)
    {
        if (!entries.empty()) {
            return evaluator<Real>(this, vmin, vmax);
        } else {
            return evaluator<Real>(getDefault(), vmin, vmax);
        }
    }

    inline friend std::ostream& operator << (std::ostream& out, const ColorMap& m )
    {
        out << m.entries;
        return out;
    }

    inline friend std::istream& operator >> (std::istream& in, ColorMap& m )
    {
        in >> m.entries;
        return in;
    }
};


} // namespace sofa::helper


#endif
