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
    //blah.print();

    int* tNow       = new int[blah.getCols()];
    int* tLast      = new int[blah.getCols()];
    int* tMoveCount = new int[blah.getCols()];

    for (int i = 0; i < blah.getCols(); i++)
    {
        tLast[i] = i;
        tMoveCount[i] = 0;
    }


    for (int i = 0; i < blah.getRows(); i++)
    {
        for (int j = 0; j < blah.getCols(); j++)
        {
            tNow[blah.getElement(i, j)] = j; // find position of each number
        }

        // compare with last position to see who moved
        for (int j = 0; j < blah.getCols(); j++)
        {
            if (tNow[j] != tLast[j])
                tMoveCount[j] += 1;
        }

        for (int j = 0; j < blah.getCols(); j++)
            tLast[j] = tNow[j];
    }

    cout << "Status report for SJT, n = " << blah.getCols() << "\n";
    for (int i = 0; i < blah.getCols(); i++)
        cout << i+1 << ", " << tMoveCount[i] << "\n";

    return EXIT_SUCCESS;
}
