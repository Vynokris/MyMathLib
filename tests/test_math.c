#include <my_math.h>
#include <stdio.h>

int main()
{
    printf("\n----- TESTING MATH UTILITY FUNCTIONS -----\n");
    printf("sqpow(2)           = %d\n", (int)sqpow(2));
    printf("degToRad(180)      = %f\n", degToRad(180));
    printf("radToDeg(PI)       = %d\n", roundInt(radToDeg(PI)));
    printf("clamp(11, 0, 10)   = %d\n", (int)clamp(11, 0, 10));
    printf("clamp(-1, 0, 10)   = %d\n", (int)clamp(-1, 0, 10));
    printf("clampUnder(11, 10) = %d\n", (int)clampUnder(11, 10));
    printf("clampAbove(-1, 0)  = %d\n", (int)clampAbove(-1, 0));

    printf("\n----- TESTING Vector2 MATH FUNCTIONS -----\n");
    printf("Vector2Zero                        = { %d, %d }\n", (int)Vector2Zero().x, (int)Vector2Zero().y);
    Vector2 vec_from_rot = Vector2FromAngle(PI, 1);
    printf("Vector2FromRotation(pi, 1)         = { %d, %d }\n", roundInt(vec_from_rot.x), roundInt(vec_from_rot.y));
    printf("Vector2GetAngle(vec_from_rot)  = %f\n", Vector2GetAngle(vec_from_rot));
    Vector2 vec_add = Vector2Add((Vector2){1, 2}, (Vector2){3, 4});
    printf("Vector2Add         ({1, 2}, {3, 4}) = { %d, %d }\n", (int)vec_add.x, (int)vec_add.y);
    Vector2 vec_add_val = Vector2AddVal((Vector2){1, 2}, 3);
    printf("Vector2AddVal      ({1, 2}, 3)      = { %d, %d }\n", (int)vec_add_val.x, (int)vec_add_val.y);
    Vector2 vec_sub = Vector2Substract((Vector2){1, 2}, (Vector2){3, 4});
    printf("Vector2Substract   ({1, 2}, {3, 4}) = { %d, %d }\n", (int)vec_sub.x, (int)vec_sub.y);
    Vector2 vec_mult = Vector2Multiply((Vector2){1, 2}, (Vector2){3, 4});
    printf("Vector2Multiply    ({1, 2}, {3, 4}) = { %d, %d }\n", (int)vec_mult.x, (int)vec_mult.y);
    Vector2 vec_mult_val = Vector2MultiplyVal((Vector2){1, 2}, 3);
    printf("Vector2MultiplyVal ({1, 2}, 3)      = { %d, %d }\n", (int)vec_mult_val.x, (int)vec_mult_val.y);
    Vector2 vec_div = Vector2Divide((Vector2){1, 2}, (Vector2){3, 4});
    printf("Vector2Divide      ({1, 2}, {3, 4}) = { %f, %f }\n", vec_div.x, vec_div.y);
    printf("Vector2Length      ({1, 1})         = %f\n", Vector2Length((Vector2){1, 1}));
    Vector2 vec_neg = Vector2Negate((Vector2){1, -2});
    printf("Vector2Negate      ({1, -2})        = { %d, %d }\n", (int)vec_neg.x, (int)vec_neg.y);
    Vector2 vec_norm = Vector2Normalize((Vector2){1, 2});
    printf("Vector2Normalize   ({1, 2})         = { %f, %f } of length = %d\n", vec_norm.x, vec_norm.y, roundInt(Vector2Length(vec_norm)));
    printf("Vector2DotProduct  ({0, 1}, {1, 0}) = %d\n", roundInt(Vector2DotProduct((Vector2){0, 1}, (Vector2){1, 0})));
    printf("Vector2CrossProduct({0, 1}, {0, 1}) = %d\n", roundInt(Vector2CrossProduct((Vector2){0, 1}, (Vector2){0, 1})));
    printf("Vector2Angle       ({0, 1}, {1, 0}) = %d\n", roundInt(radToDeg(Vector2Angle((Vector2){0, 1}, (Vector2){1, 0}))));
    Vector2 vec_rotate = Vector2Rotate((Vector2){0, 1}, PI/2);
    printf("Vector2Rotate      ({0, 1}, pi/2)   = { %d, %d }\n", roundInt(vec_rotate.x), roundInt(vec_rotate.y));

    printf("\n");

    return 0;
}