#pragma once
#include <assert.h>


// ---------- DEFINITIONS --------- //

#define PI 3.14159265358979323846f



// ---------- STRUCTURES ---------- //

// Point structure that holds values for x and y.
typedef struct Point {
    double x;
    double y;
} Point;

// Vector structure that holds values for x and y.
typedef struct Vector {
    double x;
    double y;
} Vector;

// Sgment structure that holds values for the starting point and the end point.
typedef struct Segment {
    Point start;
    Point end;
} Segment;

// Triangle structure that holds values for 3 points.
typedef struct Triangle {
    Point p1;
    Point p2;
    Point p3;
} Triangle;

// Rectangle structure that holds values for the origin point, width and height.
typedef struct Rectangle {
    Point orig;
    double width;
    double height;
} Rectangle;

// Polygon structure that holds values for the origin point, the radius and the number of sides.
typedef struct Polygon {
    Point orig;
    double radius;
    int sides;
} Polygon;

// Circle structure that holds values for the origin point and radius.
typedef struct Circle {
    Point orig;
    double radius;
} Circle;



// ---------- MATH FUNCTIONS ---------- //

// Returns the absolute value of the given value.
double absolute(double value);

// Returns the square of the given value.
double sqpow(double val);

// Returns the given value risen to the given power.
double power(double val, int power);

// Returns the square root of a value.
double sqroot(double val);

// Rounds up the given value.
int ceilVal(double val);

// Rounds down the given value.
int floorVal(double val);

// Rounds the given value.
int roundVal(double val);

// Returns the given value with the given sign.
double copySign(double val, int sign);

// Returns the highest of the two given values.
double max(double val1, double val2);

// Returns the lowest of the two given values.
double min(double val1, double val2);



// ---------- MATH UTILITY FUNCTIONS ---------- //

// Converts degrees to radians.
double degToRad(double deg);

// Converts radians to degrees.
double radToDeg(double rad);

// Clamps a value between two values.
double clamp(double val, double min, double max);

// Clamps a value to be under a value.
double clampUnder(double val, double max);

// Clamps a value to be above a value.
double clampAbove(double val, double min);

// Linear interpolation between two values.
double lerp(double val, double start, double end);

// Normalize value in a range.
double normalize(double val, double start, double end);

// Remap value from a range to another range.
double remap(double val, double inputStart, double inputEnd, double outputStart, double outputEnd);



// ---------- VECTOR MATH FUNCTIONS ---------- //

// Returns a vector of values { 0, 0 }.
Vector vectorZero();

// Adds two vectors together.
Vector vectorAdd(Vector v1, Vector v2);

// Adds a value to a vector.
Vector vectorAddVal(Vector v1, double val);

// Substract two vectors together.
Vector vectorSubstract(Vector v1, Vector v2);

// Multiplies two vectors together.
Vector vectorMultiply(Vector v1, Vector v2);

// Multiplies a vector by a value.
Vector vectorMultiplyVal(Vector v1, double val);

// Divides a vector by another.
Vector vectorDivide(Vector v1, Vector v2);

// Returns the length of the given vector.
double vectorLength(Vector v);

// Negates the values of the given vector.
Vector vectorNegate(Vector v);

// Normalizes the given vector.
Vector vectorNormalize(Vector v);
