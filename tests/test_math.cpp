#include <iostream>
#include <my_math.hpp>
#include <my_math.inl>

using namespace matrix;

int main(void)
{
    Matrix<2, 2> m1(2, 3, 1, -2);
    printf("Matrix 1:\nDeterminant: %.2f\n\n", m1.det2());
    m1.print();
    m1.inv2().print();
    (m1*m1.inv2()).print();

    printf("\n");

    Matrix<3, 3> m2(2, 5, 3, -1, 2, 5, 2, 1, 3);
    printf("Matrix 2:\nDeterminant: %.2f\n\n", m2.det3());
    m2.print();
    m2.inv3().print();
    (m2*m2.inv3()).print();

    printf("\n");

    Matrix<4, 4> m3(1, 0, 4, -6, 2, 5, 0, 3, -1, 2, 3, 5, 2, 1, -2, 3);
    printf("Matrix 3:\nDeterminant: %.2f\n\n", m3.det4());
    m3.print();
    m3.inv4().print();
    (m3*m3.inv4()).print();

    Matrix<4, 4> m4(true);
    m4.print();

    return 0;
}
