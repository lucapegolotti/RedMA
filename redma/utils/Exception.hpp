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

#ifndef EXCEPTION_HPP
#define EXCEPTION_HPP

#include <exception>
#include <string>
#include <redma/utils/PrintLog.hpp>
namespace RedMA
{

class Exception: public std::exception
{
public:
    explicit Exception(const char* message) :
      M_msg(message)
    {
        printMsg();
    }

    explicit Exception(const std::string& message) :
      M_msg(message)
    {
        printMsg();
    }

    virtual ~Exception() throw ()
    {
    }

    virtual const char* what() const throw ()
    {
       return M_msg.c_str();
    }

    void printMsg()
    {
        std::string msg = "[Exception] ";
        msg += M_msg;
        msg += "\n";
        printlog(RED, msg);
    }

protected:
    std::string M_msg;
};

}  // namespace RedMA

#endif  // EXCEPTION_HPP
