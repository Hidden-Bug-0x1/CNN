#include <stdlib.h>

u_long pow(u_long base, u_long exp)
{
    if (exp == 0)
        return 1;
    else if (exp == 1)
        return base;

    u_long out = base;
    while (exp--)
    {
        out *= base;
    }

    return out;
}