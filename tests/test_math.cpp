#include <iostream>
#include <my_math.hpp>
#include <my_math.inl>

using namespace arithmetic;
using namespace render3D;

int main(void)
{
    Matrix<4, 4> m0 = getTranslationMatrix({0, 0, 0});
    Matrix<4, 4> m1 = getXRotationMatrix(PI);
    Matrix<4, 4> m2 = getYRotationMatrix(0.f);
    Matrix<4, 4> m3 = getZRotationMatrix(0.f);
    Matrix<4, 4> m4 = getScaleMatrix({0, 0, 0});
    Matrix<4, 4> m = getTransformMatrix({1, 1, 1}, {PI, PI, PI}, {1, 1, 1});
    m0.print();
    m1.print();
    m2.print();
    m3.print();
    m4.print();
    m.print();
    return 0;
}