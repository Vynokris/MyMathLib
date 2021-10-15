#pragma once
#include <assert.h>
#include <math.h>
#include <stdlib.h>


// ---------- DEFINITIONS --------- //

#define PI 3.14159265358979323846f



// ---------- STRUCTURES ---------- //

// Point structure that holds values for x and y.
typedef struct Point {
    double x, y;
} Point;

// Vector structure that holds values for x and y.
typedef struct Vector2 {
    double x, y;
} Vector2;

// Sgment structure that holds values for the starting point and the end point.
typedef struct Segment {
    Point a, b;
} Segment;

// Triangle structure that holds values for 3 points.
typedef struct Triangle {
    Point a, b, c;
} Triangle;

// Rectangle structure that holds values for the origin point, width and height.
typedef struct Rectangle {
    Point origin;
    double width, height;
} Rectangle;

// Polygon structure that holds values for the origin point, the radius and the number of sides.
typedef struct Polygon {
    Point origin;
    double radius;
    int sides;
} Polygon;

// Circle structure that holds values for the origin point and radius.
typedef struct Circle {
    Point origin;
    double radius;
} Circle;



// ---------- MATH FUNCTIONS ---------- //

// Rounds the given value to the nearest integer.
int roundInt(double val)
{
    return (int)round(val);
}

// Rounds down the given value.
int floorInt(double val)
{
    return (int)floor(val);
}

// Rounds up the given value.
int ceilInt(double val)
{
    return (int)ceil(val);
}

// Returns the square of the given value.
double sqpow(double val)
{
    return val * val;
}

// Returns 1 if val is positive or null, -1 if it is negative.
int signOf(double val)
{
    if (val == 0) return 1;
    return val / abs(val);
}

// Converts degrees to radians.
double degToRad(double deg)
{
    return deg * (PI/180.0f);
}

// Converts radians to degrees.
double radToDeg(double rad)
{
    return rad * (180.0f/PI);
}

// Clamps a value between two values.
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

// Clamps a value to be under a value.
double clampUnder(double val, double max)
{
    if (val > max) {
        val = max;
    }

    return val;
}

// Clamps a value to be above a value.
double clampAbove(double val, double min)
{
    if (val < min) {
        val = min;
    }

    return val;
}

// TODO: figure out how these work.
// Linear interpolation between two values.
double lerp(double val, double start, double end);
// Normalize value in a range.
double normalize(double val, double start, double end);
// Remap value from a range to another range.
double remap(double val, double inputStart, double inputEnd, double outputStart, double outputEnd);



// ---------- VECTOR MATH FUNCTIONS ---------- //

// Returns a vector of values { 0, 0 }.
Vector2 Vector2Zero()
{
    return (Vector2){ 0, 0 };
}

// Creates a vector given a rotation and a length.
Vector2 Vector2FromAngle(double rad, double length)
{
    return (Vector2){ cos(rad) * length, sin(rad) * length };
}

// Returns the angle of the given vector.
double Vector2GetAngle(Vector2 v)
{
    return acos(v.x);
}

// Adds two vectors together.
Vector2 Vector2Add(Vector2 v1, Vector2 v2)
{
    return (Vector2){ v1.x + v2.x, v1.y + v2.y };
}

// Adds a value to a vector.
Vector2 Vector2AddVal(Vector2 v, double val)
{
    return (Vector2){ v.x + val, v.y + val };
}

// Substract two vectors together.
Vector2 Vector2Substract(Vector2 v1, Vector2 v2)
{
    return (Vector2){ v1.x - v2.x, v1.y - v2.y };
}

// Multiplies two vectors together.
Vector2 Vector2Multiply(Vector2 v1, Vector2 v2)
{
    return (Vector2){ v1.x * v2.x, v1.y * v2.y };
}

// Multiplies a vector by a value.
Vector2 Vector2MultiplyVal(Vector2 v, double val)
{
    return (Vector2){ v.x * val, v.y * val };
}

// Divides a vector by another.
Vector2 Vector2Divide(Vector2 v1, Vector2 v2)
{
    return (Vector2){ v1.x / v2.x, v1.y / v2.y };
}

// Returns the length of the given vector.
double Vector2Length(Vector2 v)
{
    return sqrt(sqpow(v.x) + sqpow(v.y));
}

// Negates the values of the given vector.
Vector2 Vector2Negate(Vector2 v)
{
    return (Vector2){ -v.x, -v.y };
}

// Normalizes the given vector.
Vector2 Vector2Normalize(Vector2 v)
{
    return (Vector2){ v.x / Vector2Length(v), v.y / Vector2Length(v) };
}

// Returns the dot product of two vectors.
double Vector2DotProduct(Vector2 v1, Vector2 v2)
{
    return (v1.x * v2.x) + (v1.y * v2.y);
}

// Returns the cross product of two vectors.
double Vector2CrossProduct(Vector2 v1, Vector2 v2)
{
    return (v1.x * v2.y) - (v1.y * v2.x);
}

// Returns the angle between two vectors.
double Vector2Angle(Vector2 v1, Vector2 v2)
{
    double v1_angle = Vector2GetAngle(v1);
    double v2_angle = Vector2GetAngle(v2);
    return (v1_angle >= v2_angle ? (v1_angle - v2_angle) : (v2_angle - v1_angle));
}

// Rotates the given vector by the given angle (in radians).
Vector2 Vector2Rotate(Vector2 v, double angle)
{
    double v_length = Vector2Length(v);
    double v_angle = Vector2GetAngle(v);
    return Vector2FromAngle(v_angle + angle, v_length);
}
