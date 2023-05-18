#include <iostream>
#include "cl_matrix.h"


int main()
{

    Matrix A(5,5);
    int k = 0;
    for (int i= 0; i < 5; ++i)
        for (int j = 0 ; j< 5; ++j)
        {
            A.set_elem(i,j,k);
            ++k;
        }
    std::cout << A;

    !A;
    std::cout << A;

}

