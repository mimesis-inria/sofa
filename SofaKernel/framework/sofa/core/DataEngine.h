/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2018 INRIA, USTL, UJF, CNRS, MGH                    *
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
#ifndef SOFA_CORE_DATAENGINE_H
#define SOFA_CORE_DATAENGINE_H

#include <sofa/core/core.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/DataTracker.h>

namespace sofa
{

namespace core
{

/**
 *  \brief from a set of Data inputs computes a set of Data outputs
 *
 * Implementation good rules:
 *
 * void init()
 * {
 *    addInput // indicate all inputs
 *    addOutput // indicate all outputs
 * }
 *
 * // optional (called each time a data is modified in the gui)
 * // it is not always desired
 * void reinit()
 * {
 *    update();
 * }
 *
 * void doUpdate() override
 * {
 *    access your inputs, set your outputs...
 * }
 *
 */
class SOFA_CORE_API DataEngine : public core::DataTrackerDDGNode, public virtual core::objectmodel::BaseObject
{
public:
    SOFA_ABSTRACT_CLASS2(DataEngine, core::objectmodel::BaseObject, core::DataTrackerDDGNode);
    SOFA_BASE_CAST_IMPLEMENTATION(DataEngine)
protected:
    /// Constructor
    DataEngine();

    /// Destructor. Do nothing
    virtual ~DataEngine();

private:
	DataEngine(const DataEngine& n) ;
	DataEngine& operator=(const DataEngine& n) ;

    /// Called in update(), back-propagates the data update
    /// in the data dependency graph
    void updateAllInputs();

protected:
    /// Where you put your engine's impl
    virtual void doUpdate() = 0;

public:
    /// Updates your inputs and calls cleanDirty() for you.
    /// User implementation moved to doUpdate()
    virtual void update() final;

    /// Add a new input to this engine
    /// Automatically adds the input fields to the datatracker
    void addInput(sofa::core::objectmodel::BaseData* data);

    /// Add a new output to this engine
    void addOutput(objectmodel::BaseData* n);

    // The methods below must be redefined because of the
    // double inheritance from Base and DDGNode

    /// @name Class reflection system
    /// @{

    template<class T>
    static std::string typeName(const T* ptr= NULL)
    {
        return core::objectmodel::BaseObject::typeName(ptr);
    }

    /// Helper method to get the class name of a type derived from this class
    ///
    /// This method should be used as follow :
    /// \code  T* ptr = NULL; std::string type = T::className(ptr); \endcode
    /// This way derived classes can redefine the className method
    template<class T>
    static std::string className(const T* ptr= NULL)
    {
        return core::objectmodel::BaseObject::className(ptr);
    }

    /// Helper method to get the namespace name of a type derived from this class
    ///
    /// This method should be used as follow :
    /// \code  T* ptr = NULL; std::string type = T::namespaceName(ptr); \endcode
    /// This way derived classes can redefine the namespaceName method
    template<class T>
    static std::string namespaceName(const T* ptr= NULL)
    {
        return core::objectmodel::BaseObject::namespaceName(ptr);
    }

    /// Helper method to get the template name of a type derived from this class
    ///
    /// This method should be used as follow :
    /// \code  T* ptr = NULL; std::string type = T::templateName(ptr); \endcode
    /// This way derived classes can redefine the templateName method
    template<class T>
    static std::string templateName(const T* ptr= NULL)
    {
        return core::objectmodel::BaseObject::templateName(ptr);
    }

    /// Helper method to get the shortname of a type derived from this class.
    /// The default implementation return the class name.
    ///
    /// This method should be used as follow :
    /// \code  T* ptr = NULL; std::string type = T::shortName(ptr); \endcode
    /// This way derived classes can redefine the shortName method
    template< class T>
    static std::string shortName( const T* ptr = NULL, core::objectmodel::BaseObjectDescription* desc = NULL )
    {
        return core::objectmodel::BaseObject::shortName(ptr,desc);
    }

    template<class T>
    static void dynamicCast(T*& ptr, Base* b)
    {
        core::objectmodel::BaseObject::dynamicCast(ptr, b);
    }

    /// @}

    /// This method is needed by DDGNode
    const std::string& getName() const override
    {
        return objectmodel::BaseObject::getName();
    }

    /// This method is needed by DDGNode
    objectmodel::Base* getOwner() const override
    {
        return const_cast<DataEngine*>(this);
    }

    /// This method is needed by DDGNode
    objectmodel::BaseData* getData() const override
    {
        return NULL;
    }

    /// Add a link.
    void addLink(objectmodel::BaseLink* l)
    {
        objectmodel::BaseObject::addLink(l);
    }
};

} // namespace core

} // namespace sofa

#endif
