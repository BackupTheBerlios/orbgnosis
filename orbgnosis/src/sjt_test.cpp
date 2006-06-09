#include "SJT.h"
#include <iostream>
using namespace std;

int main (void)
{
    int numTargets;
    cout << "This is a test of the SJT class.\n";
    cout << "Enter number of targets [1-" << SJT_MAX << "]: ";
    cin >> numTargets;

    SJT blah(numTargets);
    blah.print();

    return EXIT_SUCCESS;
}
