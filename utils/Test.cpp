#include "test.hpp"

Test::Test(std::string testName) :
  M_subTests(),
  M_testName(testName),
  M_nTests(0),
  M_successes(0)
{

}

void Test::addSubTest(void (*subTest)(Test&))
{
  M_subTests.push_back(subTest);
}

void Test::assert(bool statement)
{
  M_successes = statement? M_successes + 1 : M_successes;
  if (!statement)
    printlog(RED, "ERROR: assertion failed\n");
  M_nTests++;
}


void Test::run()
{

  std::string msg = "\nRunning test " + M_testName + "\n";
  printlog(MAGENTA, msg);
  printlog(WHITE, "-----------------------\n");

  int count = 0;
  for (std::vector<void (*)(Test&)>::iterator it = M_subTests.begin();
       it < M_subTests.end(); it++)
  {
    count++;
    printlog(BLUE, "\tTest ");
    printlog(BLUE, count);
    printlog(BLUE, " ... \n");

    (*(*it))(*this);
  }

  if (M_successes == M_nTests)
  {
    printlog(GREEN, "Successful asserts: ");
    printlog(GREEN, M_successes);
    printlog(GREEN, "/");
    printlog(GREEN, M_nTests);
    printlog(GREEN, "\n");
  }
  else
  {
    printlog(RED, "Successful asserts: ");
    printlog(RED, M_successes);
    printlog(RED, "/");
    printlog(RED, M_nTests);
    printlog(RED, "\n");
  }
}
