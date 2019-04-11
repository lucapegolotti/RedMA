#include <iostream>
#include <string>
#include <Test.hpp>
#include <Tube.hpp>

using namespace ReMA;

void subTest1(Test& test)
{
  Tube tube;
  try
  {
      tube.setParameterValue("abab", 2.0);
      test.assertTrue(false);
  }
  catch (std::exception& e)
  {
      test.assertTrue(true);
  }
}

int main()
{
  Test test("TubeTest");

  test.addSubTest(*subTest1);

  test.run();

  return 0;
}
