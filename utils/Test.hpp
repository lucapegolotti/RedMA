// Reduced Modeling of Arteries (ReMA)
// Copyright (C) 2019  Luca Pegolotti
//
// ReMA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ReMA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef TEST_HPP
#define TEST_HPP

#include <iostream>
#include <string>
#include <vector>

#include "PrintLog.hpp"

namespace ReMA
{

class Test
{
public:
  Test(std::string testName);

  void addSubTest(void (*subTest)(Test&));

  void assertTrue(bool statement);

  void run();

private:
  Test() {};

  std::vector<void (*)(Test&)> M_subTests;
  std::string M_testName;

  unsigned int M_nTests;
  unsigned int M_successes;
};

}  // namespace ReMA

#endif  // TEST_HPP
