#include "driver.h"

int main(int narg, char **arg)
{
  Driver *driver = new Driver(narg, arg);
  delete driver;

return 0;
}
