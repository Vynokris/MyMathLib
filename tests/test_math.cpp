#include <iostream>
#include <my_math.hpp>
#include <my_math.inl>

using namespace arithmetic;

int main(void)
{
    Matrix<4, 4> m0( 
        { 1, 0, 4, 6 },
        { 2, 5, 0, 3 },
        { -1, 2, 3, 5 },
        { 2, 1, -2, 3 }
    );
    m0.print();
    Vector4 v0(1, 1, 1, 1);
    std::cout << "v0 * m0: " << (v0 * m0).x << ", " << (v0 * m0).y << ", " << (v0 * m0).z << ", " << (v0 * m0).w << "\n";

    return 0;
}