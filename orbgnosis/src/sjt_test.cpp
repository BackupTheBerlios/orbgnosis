#include "SJT.h"
#include <iostream>
using namespace std;

int main (void)
{
    int numTargets;
    cout << "This is a test of the SJT class.\n";
    cout << "Enter number of targets [1-" << SJT_MAX << "]: ";
    cin >> numTargets;

    if (numTargets < 1)
    {
        cout << "Error: bad input.";
        return 1;
    }

    SJT blah(numTargets);
    blah.print();

    return EXIT_SUCCESS;
}
