#include "Test.hpp"

namespace RedMA
{

Test::
Test(std::string testName, std::shared_ptr<Epetra_Comm> comm) :
  M_subTests(),
  M_testName(testName),
  M_nTests(0),
  M_successes(0),
  M_comm(comm)
{
}

void
Test::
addSubTest(void (*subTest)(Test&))
{
    M_subTests.push_back(subTest);
}

void
Test::
assertTrue(bool statement)
{
    M_successes = statement? M_successes + 1 : M_successes;
    if (!statement)
        printlog(RED, "ERROR: assertion failed\n");
    M_nTests++;
}

void
Test::
run()
{

    std::string msg = "\nRunning test " + M_testName + "\n";
    printlog(MAGENTA, msg);
    printlog(WHITE, "-----------------------\n");

    int count = 0;
    for (std::vector<void (*)(Test&)>::iterator it = M_subTests.begin();
       it < M_subTests.end(); it++)
    {
        count++;
        std::string msgBlue = std::string("\tTest ") + std::to_string(count)
                              + " ... \n";
        printlog(BLUE, msgBlue);

        (*(*it))(*this);
    }

    if (M_successes == M_nTests)
    {
        std::string msgGreen = std::string("Successful asserts: ") +
                               std::to_string(M_successes) + "/" +
                               std::to_string(M_nTests) + "\n";
        printlog(GREEN, msgGreen);
    }
    else
    {
        std::string msgRed = std::string("Successful asserts: ") +
                             std::to_string(M_successes) + "/" +
                             std::to_string(M_nTests) + "\n";
        printlog(RED, msgRed);
    }
}

std::shared_ptr<Epetra_Comm>&
Test::
getComm()
{
    return M_comm;
}

}  // namespace RedMA
