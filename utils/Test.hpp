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
