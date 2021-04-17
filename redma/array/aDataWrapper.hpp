// Reduced Modeling of Arteries (RedMA)
// Copyright (C) 2019  Luca Pegolotti
//
// RedMA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RedMA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef aDATAWRAPPER_HPP
#define aDATAWRAPPER_HPP

#include <redma/RedMA.hpp>
#include <redma/array/Datatypes.hpp>
#include <redma/utils/Exception.hpp>
#include <redma/array/TypesUtils.hpp>

#include <memory>

namespace RedMA
{

/*! \brief Abstract class for a data wrapper.
 *
 * The data is stored as shared pointer.
 */
class aDataWrapper
{
public:

    /// Virtual distructor.
    virtual ~aDataWrapper() {}

    /*! \brief Getter for the internal type.
     *
     * \return The internal Datatype.
     */
    virtual Datatype type() const = 0;

    /*! \brief Getter for the data.
     *
     * \return Shared pointer to the data.
     */
    virtual shp<void> data() const = 0;

    /*! \brief Setter for the data.
     *
     * \param data Shared pointer to the data.
     */
    virtual void setData(shp<void> data) = 0;

    /*! \brief Returns true if the data container is zero.
     *
     * The functions returns true even if the data is not set.
     * \return True if the data is zero or not set.
     */
    virtual bool isZero() const = 0;

    /*! \brief Clones the wrapper
     *
     * \return Raw pointer to a copy of the wrapper.
     */
    virtual aDataWrapper* clone() const = 0;

    /*! \brief Shallow copy.
     *
     * The shared pointer to the data of the argument is copied.
     *
     * \param other Shared pointer to another aDataWrapper.
     */
    virtual void shallowCopy(shp<aDataWrapper> other) = 0;

    /*! \brief Deep copy.
     *
     * The data of the argument is copied.
     *
     * \param other Shared pointer to another aDataWrapper.
     */
    virtual void deepCopy(shp<aDataWrapper> other) = 0;

    /*! \brief Save content of the wrapper to file.
     *
     * \param namefile Name of the output file.
     */
    virtual void dump(std::string namefile) const = 0;
};

}

#endif // aDATAWRAPPER_HPP
