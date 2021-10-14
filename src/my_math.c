#include <my_math.h>
#include <stdio.h>


// ---------- MATH FUNCTIONS ---------- //

double absolute(double val)
{
    if (val < 0) val = -val;
    return val;
}

double sqpow(double val)
{
    return val * val;
}

double power(double val, int power)
{
    if (power == 0) return 1;

    else if (power > 0) {
        for (int i = 1; i < power; i++) {
            val *= val;
        }
        return val;
    }

    else {
        double output = 1;
        for (int i = 0; i > power; i--) {
            output /= val;
        }
        return output;
    }
}

double sqroot(double val)
{
    assert(val >= 0);

    double root = val;
    double precision = 0.000001f;

    while (absolute(val - sqpow(root)) > precision) {
        root = (root + val / root) / 2;
    }

    return root;
}

int ceilVal(double val)
{
    return (int)val + 1;
}

int floorVal(double val)
{
    return (int)val;
}

int roundVal(double val)
{
    return (val + 0.5f > ceilVal(val) ? ceilVal(val) : floorVal(val));
}

double copySign(double val, int sign)
{
    return absolute(val) * sign / absolute(sign);
}

double max(double val1, double val2)
{
    return (val1 > val2 ? val1 : val2);
}

double min(double val1, double val2)
{
    return (val1 < val2 ? val1 : val2);
}


// ---------- MATH UTILITY FUNCTIONS ---------- //

double degToRad(double deg)
{
    return deg * (PI/180.0f);
}

double radToDeg(double rad)
{
    return rad * (180.0f/PI);
}

double clamp(double val, double min, double max)
{
    assert(min <= max);

    if (val < min) {
        val = min;
    }
    if (val > max) {
        val = max;
    }

    return val;
}

double clampUnder(double val, double max)
{
    if (val > max) {
        val = max;
    }

    return val;
}

double clampAbove(double val, double min)
{
    if (val < min) {
        val = min;
    }

    return val;
}

// TODO: figure out how these work.
double Lerp(double val, double start, double end);
double Normalize(double val, double start, double end);
double Remap(double val, double inputStart, double inputEnd, double outputStart, double outputEnd);


// ---------- VECTOR MATH FUNCTIONS ---------- //

Vector vectorZero()
{
    return (Vector){ 0, 0 };
}

Vector vectorAdd(Vector v1, Vector v2)
{
    return (Vector){ v1.x + v2.x, v1.y + v2.y };
}

Vector vectorAddVal(Vector v, double val)
{
    return (Vector){ v.x + val, v.y + val };
}

Vector vectorSubstract(Vector v1, Vector v2)
{
    return (Vector){ v1.x - v2.x, v1.y - v2.y };
}

Vector vectorMultiply(Vector v1, Vector v2)
{
    return (Vector){ v1.x * v2.x, v1.y * v2.y };
}

Vector vectorMultiplyVal(Vector v, double val)
{
    return (Vector){ v.x * val, v.y * val };
}

Vector vectorDivide(Vector v1, Vector v2)
{
    return (Vector){ v1.x / v2.x, v1.y / v2.y };
}

double vectorLength(Vector v)
{
    return sqroot(sqpow(v.x) + sqpow(v.y));
}

Vector vectorNegate(Vector v)
{
    return (Vector){ -v.x, -v.y };
}

Vector vectorNormalize(Vector v)
{
    return (Vector){ v.x / vectorLength(v), v.y / vectorLength(v) };
}
