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

#ifndef DATACONTAINER_HPP
#define DATACONTAINER_HPP

#include <redma/RedMA.hpp>
#include <lifev/core/filter/GetPot.hpp>

namespace RedMA
{

class DataContainer
{
public:
    DataContainer();

    void setDatafile(const std::string& datafile);

    void setInflow(const std::function<double(double)>& inflow);

    void setVerbose(bool verbose);

    inline GetPot getDatafile() const {return *M_datafile;}

    inline std::function<double(double)> getInflow() const {return M_inflow;}

    inline bool getVerbose() const {return M_verbose;}

    std::string operator()(std::string location, std::string defValue) const;

    int operator()(std::string location, int defValue) const;

    double operator()(std::string location, double defValue) const;

protected:
    SHP(GetPot)                     M_datafile;
    std::function<double(double)>   M_inflow;
    bool                            M_verbose;
};

}

#endif // DATACONTAINER_HPP
