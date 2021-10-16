#pragma once
#ifndef MY_MATH_H
#define MY_MATH_H

#include <assert.h>
#include <raylib.h>
#include <math.h>
#include <stdlib.h>

// ---------- DEFINITIONS --------- //

#define PI 3.14159265358979323846f
#define TO_RAY_V2(v) (*(*Vector2) & v)
#define TO_RAY_REC(r) (*(*Rectangle) & r)

// ---------- STRUCTURES ---------- //

// Vector structure that holds values for x and y.
typedef struct MyVector2 {
    double x, y;
} MyVector2;

// Sgment structure that holds values for the starting point and the end point.
typedef struct Segment {
    MyVector2 a, b;
} Segment;

// Triangle structure that holds values for 3 points.
typedef struct Triangle {
    MyVector2 a, b, c;
} Triangle;

// Rectangle structure that holds values for the origin point, width and height.
typedef struct MyRectangle {
    MyVector2 origin;
    double width, height;
} MyRectangle;

// Polygon structure that holds values for the origin point, the radius and the number of sides.
typedef struct Polygon {
    MyVector2 origin;
    double radius;
    int sides;
} Polygon;

// Circle structure that holds values for the origin point and radius.
typedef struct Circle {
    MyVector2 origin;
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
    if (val == 0)
        return 1;
    return val / abs((int)val);
}

// Converts degrees to radians.
double degToRad(double deg)
{
    return deg * (PI / 180.0f);
}

// Converts radians to degrees.
double radToDeg(double rad)
{
    return rad * (180.0f / PI);
}

// Clamps a value between two values.
double clamp(double val, double min, double max)
{
    assert(min <= max);

    if (val < min)
    {
        val = min;
    }
    if (val > max)
    {
        val = max;
    }

    return val;
}

// Clamps a value to be under a value.
double clampUnder(double val, double max)
{
    if (val > max)
    {
        val = max;
    }

    return val;
}

// Clamps a value to be above a value.
double clampAbove(double val, double min)
{
    if (val < min)
    {
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

// Returns a vector of the given values.
MyVector2 Vector2Create(double x, double y)
{
    return (MyVector2){x, y};
}

// Returns a vector of values { 0, 0 }.
MyVector2 Vector2Zero()
{
    return (MyVector2){0, 0};
}

// Adds two vectors together.
MyVector2 Vector2Add(MyVector2 v1, MyVector2 v2)
{
    return (MyVector2){v1.x + v2.x, v1.y + v2.y};
}

// Adds a value to a vector.
MyVector2 Vector2AddVal(MyVector2 v, double val)
{
    return (MyVector2){v.x + val, v.y + val};
}

// Substract two vectors together.
MyVector2 Vector2Substract(MyVector2 v1, MyVector2 v2)
{
    return (MyVector2){v1.x - v2.x, v1.y - v2.y};
}

// Multiplies two vectors together.
MyVector2 Vector2Multiply(MyVector2 v1, MyVector2 v2)
{
    return (MyVector2){v1.x * v2.x, v1.y * v2.y};
}

// Multiplies a vector by a value.
MyVector2 Vector2MultiplyVal(MyVector2 v, double val)
{
    return (MyVector2){v.x * val, v.y * val};
}

// Divides a vector by another.
MyVector2 Vector2Divide(MyVector2 v1, MyVector2 v2)
{
    return (MyVector2){v1.x / v2.x, v1.y / v2.y};
}

// Returns the length of the given vector.
double Vector2Length(MyVector2 v)
{
    return sqrt(sqpow(v.x) + sqpow(v.y));
}

// Normalizes the given vector.
MyVector2 Vector2Normalize(MyVector2 v)
{
    return (MyVector2){v.x / Vector2Length(v), v.y / Vector2Length(v)};
}

// Creates a vector given a rotation and a length.
MyVector2 Vector2FromAngle(double rad, double length)
{
    return (MyVector2){cos(rad) * length, sin(rad) * length};
}

// Returns the angle of the given vector.
double Vector2GetAngle(MyVector2 v)
{
    return acos(Vector2Normalize(v).x);
}

// Resizes the given vector to the given length.
MyVector2 Vector2SetLength(MyVector2 v, double length)
{
    return Vector2FromAngle(Vector2GetAngle(v), length);
}

// Negates the values of the given vector.
MyVector2 Vector2Negate(MyVector2 v)
{
    return (MyVector2){-v.x, -v.y};
}

// Returns the dot product of two vectors.
double Vector2DotProduct(MyVector2 v1, MyVector2 v2)
{
    return (v1.x * v2.x) + (v1.y * v2.y);
}

// Returns the cross product of two vectors.
double Vector2CrossProduct(MyVector2 v1, MyVector2 v2)
{
    return (v1.x * v2.y) - (v1.y * v2.x);
}

// Returns the angle between two vectors.
double Vector2Angle(MyVector2 v1, MyVector2 v2)
{
    double v1_angle = Vector2GetAngle(v1);
    double v2_angle = Vector2GetAngle(v2);
    return (v1_angle >= v2_angle ? (v1_angle - v2_angle) : (v2_angle - v1_angle));
}

// Rotates the given vector by the given angle (in radians).
MyVector2 Vector2Rotate(MyVector2 v, double angle)
{
    double v_length = Vector2Length(v);
    double v_angle = Vector2GetAngle(v);
    return Vector2FromAngle(v_angle + angle, v_length);
}

#endif
