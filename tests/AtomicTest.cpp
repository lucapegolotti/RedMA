#include "AtomicTest.hpp"

AtomicTest::
AtomicTest(std::string name, std::function<int(void)> test) :
  M_testName(name),
  M_test(test)
{
}

int
AtomicTest::
run()
{
    std::string msg = "Running atomic test ";
    msg += M_testName;
    msg += "\n";

    printlog(RedMA::GREEN, msg, true, false);

    int status = M_test();

    if (status != SUCCESS)
    {
        msg = "Atomic test ";
        msg += M_testName;
        msg += " has failed.\n";
        printlog(RedMA::RED, msg, true, false);
    }
    else
    {
        printlog(RedMA::GREEN, "SUCCESS\n", true, false);
    }

    return status;
}
