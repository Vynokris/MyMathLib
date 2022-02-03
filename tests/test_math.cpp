#include <iostream>
#include <my_math.hpp>
#include <my_math.inl>

using namespace matrix;
using namespace render3D;

int main(void)
{
    Matrix<2, 2> m1(2, 3, 1, -2);
    m1.print();
    printf("%.2f\n", m1.det2());
    m1.inv2().print();

    Matrix<3, 3> m2(2, 5, 3, -1, 2, 5, 2, 1, 3);
    m2.print();
    printf("%.2f\n", m2.det3());
    m2.inv3().print();

    Matrix<4, 4> m3(1, 0, 4, -6, 2, 5, 0, 3, -1, 2, 3, 5, 2, 1, -2, 3);
    m3.print();
    printf("%.2f\n", m3.det4());
    m3.inv4().print();

    return 0;
}