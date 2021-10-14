#include <my_math.h>
#include <stdio.h>

int main()
{
    printf("\n----- TESTING MATH FUNCTIONS -----\n");
    printf("abs(-1) = %d\n", (int)absolute(-1));
    printf("2^2 = %d\n", (int)sqpow(2));
    printf("2^4 = %d\n", (int)power(2, 4));
    printf("2^-2 = %f\n", power(2, -2));
    printf("sqrt(81) = %d\n", (int)sqroot(8));
    printf("ceil(0.2) = %d\n", ceilVal(0.2));
    printf("floor(0.8) = %d\n", floorVal(0.8));
    printf("round(0.2) = %d\n", roundVal(0.2));
    printf("round(0.8) = %d\n", roundVal(0.8));
    printf("copySign(-10, -5) = %d\n", (int)copySign(-10, -5));
    printf("max(5, 10) = %d\n", (int)max(5, 10));
    printf("min(5, 10) = %d\n", (int)min(5, 10));

    printf("\n----- TESTING MATH UTILITY FUNCTIONS -----\n");
    printf("degToRad(180) = %f\n", degToRad(180));
    printf("radToDeg(PI) = %f\n", radToDeg(PI));
    printf("clamp(11, 0, 10) = %d\n", (int)clamp(11, 0, 10));
    printf("clamp(-1, 0, 10) = %d\n", (int)clamp(-1, 0, 10));
    printf("clampUnder(11, 10) = %d\n", (int)clampUnder(11, 10));
    printf("clampAbove(-1, 0) = %d\n", (int)clampAbove(-1, 0));

    printf("\n----- TESTING VECTOR MATH FUNCTIONS -----\n");
    printf("vectorZero = { %d, %d }\n", (int)vectorZero().x, (int)vectorZero().y);
    Vector vec_add = vectorAdd((Vector){1, 2}, (Vector){3, 4});
    printf("vectorAdd({1, 2}, {3, 4}) = { %d, %d }\n", (int)vec_add.x, (int)vec_add.y);
    Vector vec_add_val = vectorAddVal((Vector){1, 2}, 3);
    printf("vectorAddVal({1, 2}, 3) = { %d, %d }\n", (int)vec_add_val.x, (int)vec_add_val.y);
    Vector vec_sub = vectorSubstract((Vector){1, 2}, (Vector){3, 4});
    printf("vectorSubstract({1, 2}, {3, 4}) = { %d, %d }\n", (int)vec_sub.x, (int)vec_sub.y);
    Vector vec_mult = vectorMultiply((Vector){1, 2}, (Vector){3, 4});
    printf("vectorMultiply({1, 2}, {3, 4}) = { %d, %d }\n", (int)vec_mult.x, (int)vec_mult.y);
    Vector vec_mult_val = vectorMultiplyVal((Vector){1, 2}, 3);
    printf("vectorMultiplyVal({1, 2}, 3) = { %d, %d }\n", (int)vec_mult_val.x, (int)vec_mult_val.y);
    Vector vec_div = vectorDivide((Vector){1, 2}, (Vector){3, 4});
    printf("vectorDivide({1, 2}, {3, 4}) = { %f, %f }\n", vec_div.x, vec_div.y);
    printf("vectorLength({1, 1}) = %f\n", vectorLength((Vector){1, 1}));
    Vector vec_neg = vectorNegate((Vector){1, -2});
    printf("vectorNegate({1, -2}) = { %d, %d }\n", (int)vec_neg.x, (int)vec_neg.y);
    Vector vec_norm = vectorNormalize((Vector){1, 2});
    printf("vectorNormalize({1, 2}) = { %f, %f } of length = %d\n", vec_norm.x, vec_norm.y, roundVal(vectorLength(vec_norm)));

    return 0;
}