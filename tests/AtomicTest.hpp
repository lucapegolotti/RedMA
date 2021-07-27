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

#include <redma/utils/Exception.hpp>
#include <redma/utils/PrintLog.hpp>

#include <functional>

#define SUCCESS 0
#define FAILURE 1
#define TZERO   1e-14

/// Class testing a specific feature.
class AtomicTest
{
public:
    /*! \brief Constructor.
     *
     * \param name The name of the test.
     * \param test The test; the return code determines failure or success.
     */
    AtomicTest(std::string name,
               std::function<int(void)> test);

    /*! \brief Run the test.
     *
     * \return A code; by default, 0 if success and 1 if failure.
     */
    int run();

private:
    std::string                     M_testName;
    std::function<bool(void)>       M_test;

};
