#pragma once
#ifndef MY_MATH_H
#define MY_MATH_H

#include <assert.h>
#include <raylib.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

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

// ---------- CONVERSIONS --------- //

// Converts a my_math 2D vector to a raylib 2D vector.
static inline Vector2 toRayVec(MyVector2 v)
{
    return (Vector2){v.x, v.y};
}

// Converts a my_math rectangle to a raylib rectangle.
static inline Rectangle toRayRec(MyRectangle r)
{
    return (Rectangle){r.origin.x, r.origin.y, r.width, r.height};
}

// ---------- OBJECT CREATION ---------- //

// Creates a 2D vector from two values.
static inline MyVector2 Vector2Create(double x, double y)
{
    return (MyVector2){x, y};
}

// Creates a 2D vector from one point to another.
static inline MyVector2 Vector2FromPoints(MyVector2 p1, MyVector2 p2)
{
    return (MyVector2){p2.x - p1.x, p2.y - p1.y};
}

// Creates a 2D vector given an angle and a length.
static inline MyVector2 Vector2FromAngle(double rad, double length)
{
    return (MyVector2){cos(rad) * length, sin(rad) * length};
}

// Creates a segment from one point to another.
static inline Segment SegmentCreate(MyVector2 a, MyVector2 b)
{
    return (Segment){a, b};
}

// Creates a triangle given three points.
static inline Triangle TriangleCreate(MyVector2 a, MyVector2 b, MyVector2 c)
{
    return (Triangle){a, b, c};
}

// Creates a rectangle given an origin point, a width and a height.
static inline MyRectangle RectangleCreate(MyVector2 origin, double width, double height)
{
    return (MyRectangle){origin, width, height};
}

// Create a polygon given an origin point, a radius and a number of sides.
static inline Polygon PolygonCreate(MyVector2 origin, double radius, double sides)
{
    return (Polygon){origin, radius, sides};
}

// Create a circle given an origin point and a radius.
static inline Circle CircleCreate(MyVector2 origin, double radius)
{
    return (Circle){origin, radius};
}

// ---------- MATH FUNCTIONS ---------- //

// Rounds the given value to the nearest int.
static inline int roundInt(double val)
{
    return (int)round(val);
}

// Rounds down the given value.
static inline int floorInt(double val)
{
    return (int)floor(val);
}

// Rounds up the given value.
static inline int ceilInt(double val)
{
    return (int)ceil(val);
}

// Returns the sqare power of the given value.
static inline double sqpow(double val)
{
    return val * val;
}

// Returns 1 if the given value is positive or null, and -1 if it is negative.
static inline int signOf(double val)
{
    if (val == 0)
        return 1;
    return val / abs((int)val);
}

// Converts the given angle from degrees to radians.
static inline double degToRad(double deg)
{
    return deg * (PI / 180.0f);
}

// Converts the given angle from radians to degrees.
static inline double radToDeg(double rad)
{
    return rad * (180.0f / PI);
}

// Clamps the given value to be superior or equal to the minimum value and inferior or equal to the maximum value.
static inline double clamp(double val, double min, double max)
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

// Clamps the given value to be inferior or equal to the maximum value.
static inline double clampUnder(double val, double max)
{
    if (val > max)
    {
        val = max;
    }

    return val;
}

// Clamps the given value to be superior or equal to the minimum value.
static inline double clampAbove(double val, double min)
{
    if (val < min)
    {
        val = min;
    }

    return val;
}

// Remaps the given value from one range to another.
static inline double remap(double val, double inputStart, double inputEnd, double outputStart, double outputEnd)
{
    // Source: https://stackoverflow.com/a/3451607/13858872
    // Find how far you are into the first range, scale that distance by the ratio of sizes of the ranges, and that's how far you should be into the second range.
    return outputStart + (val - inputStart) * (outputEnd - outputStart) / (inputEnd - inputStart);
}

// ---------- VECTOR MATH FUNCTIONS ---------- //

// Returns a vector of coordinates { 0, 0 }.
static inline MyVector2 Vector2Zero()
{
    return (MyVector2){0, 0};
}

// Adds two vectors together.
static inline MyVector2 Vector2Add(MyVector2 v1, MyVector2 v2)
{
    return (MyVector2){v1.x + v2.x, v1.y + v2.y};
}

// Adds a value to a vector.
static inline MyVector2 Vector2AddVal(MyVector2 v, double val)
{
    return (MyVector2){v.x + val, v.y + val};
}

// Substracts a vector by another.
static inline MyVector2 Vector2Substract(MyVector2 v1, MyVector2 v2)
{
    return (MyVector2){v1.x - v2.x, v1.y - v2.y};
}

// Substracts a value from a vector.
static inline MyVector2 Vector2SubstractVal(MyVector2 v, double val)
{
    return (MyVector2){v.x - val, v.y - val};
}

// Multiplies two vectors together.
static inline MyVector2 Vector2Multiply(MyVector2 v1, MyVector2 v2)
{
    return (MyVector2){v1.x * v2.x, v1.y * v2.y};
}

// Multiplies a vector by a value.
static inline MyVector2 Vector2MultiplyVal(MyVector2 v, double val)
{
    return (MyVector2){v.x * val, v.y * val};
}

// Divides a vector by another.
static inline MyVector2 Vector2Divide(MyVector2 v1, MyVector2 v2)
{
    return (MyVector2){v1.x / v2.x, v1.y / v2.y};
}

// Divides a vector by a value.
static inline MyVector2 Vector2DivideVal(MyVector2 v, double val)
{
    return (MyVector2){v.x / val, v.y / val};
}

// Returns the length of the given vetor.
static inline double Vector2Length(MyVector2 v)
{
    return sqrt(sqpow(v.x) + sqpow(v.y));
}

// Normalizes the given vector so that its length is 1.
static inline MyVector2 Vector2Normalize(MyVector2 v)
{
    return (MyVector2){v.x / Vector2Length(v), v.y / Vector2Length(v)};
}

// Returns the angle (in radians) of the given vector.
static inline double Vector2GetAngle(MyVector2 v)
{
    return acos(Vector2Normalize(v).x);
}

// Modifies the length of the given vector to correspond to the given value.
static inline MyVector2 Vector2SetLength(MyVector2 v, double length)
{
    return Vector2FromAngle(Vector2GetAngle(v), length);
}

// Negates both of the coordinates of the given vector.
static inline MyVector2 Vector2Negate(MyVector2 v)
{
    return (MyVector2){-v.x, -v.y};
}

// Returns the dot product of the given vectors.
static inline double Vector2DotProduct(MyVector2 v1, MyVector2 v2)
{
    return (v1.x * v2.x) + (v1.y * v2.y);
}

// Returns the cross product of the given vectors.
static inline double Vector2CrossProduct(MyVector2 v1, MyVector2 v2)
{
    return (v1.x * v2.y) - (v1.y * v2.x);
}

// Returns the angle (in radians) between two vectors.
static inline double Vector2Angle(MyVector2 v1, MyVector2 v2)
{
    // TODO: test this.
    double v1_angle = Vector2GetAngle(v1);
    double v2_angle = Vector2GetAngle(v2);
    return (v1_angle >= v2_angle ? (v1_angle - v2_angle) : (v2_angle - v1_angle));
}

// Rotates the given vector by the given angle.
static inline MyVector2 Vector2Rotate(MyVector2 v, double angle)
{
    double v_length = Vector2Length(v);
    double v_angle = Vector2GetAngle(v);
    return Vector2FromAngle(v_angle + angle, v_length);
}

// ---------- POLYGON MATH FUNCTIONS ---------- //

// Returns the side of the given polygon that corresponds to the given index.
static inline Segment PolygonGetSide(Polygon poly, int index)
{
    // TODO: Test this.
    assert(index < poly.sides);

    double corner_angle = 360 / poly.sides;
    double angle_offset = index * corner_angle;

    MyVector2 poly_point_a = Vector2FromAngle(angle_offset, poly.radius);
    MyVector2 poly_point_b = Vector2FromAngle(angle_offset + corner_angle, poly.radius);

    return SegmentCreate(poly_point_a, poly_point_b);
}

// ---------- DISTANCE FUNCTIONS ----------- //

// Returns the distance between two points.
static inline double distancePoints(MyVector2 p1, MyVector2 p2)
{
    return Vector2Length(Vector2FromPoints(p1, p2));
}

// ---------- COLLISIONS ---------- //

// Returns true if the given point is colliding with the given segment.
static inline bool collisionSegmentPoint(Segment s, MyVector2 p)
{
    MyVector2 vec_ap = Vector2FromPoints(s.a, p);
    MyVector2 vec_bp = Vector2FromPoints(s.b, p);
    return (roundInt(Vector2CrossProduct(vec_ap, vec_bp)) == 0 ? true : false);
}

// Returns true if the given point is colliding with the sides of the given triangle.
static inline bool collisionTrianglePoint(Triangle t, MyVector2 p)
{
    Segment ab = SegmentCreate(t.a, t.b);
    Segment bc = SegmentCreate(t.b, t.c);
    Segment ca = SegmentCreate(t.c, t.a);

    if (collisionSegmentPoint(ab, p) || collisionSegmentPoint(bc, p) || collisionSegmentPoint(ca, p))
    {
        return true;
    }
    return false;
}

// Returns true if the given point is colliding with the given rectangle.
static inline bool collisionRectanglePoint(MyRectangle r, MyVector2 p)
{
    if (p.x >= r.origin.x && p.x <= r.origin.x + r.width &&
        p.y >= r.origin.y && p.y <= r.origin.y + r.height)
    {
        return true;
    }
    return false;
}

// Returns true if the given point is colliding with the sides of the given polygon.
static inline bool collisionPolygonPoint(Polygon poly, MyVector2 p)
{
    Segment poly_seg;
    for (int i = 0; i < poly.sides; i++)
    {
        poly_seg = PolygonGetSide(poly, i);
        if (collisionSegmentPoint(poly_seg, p))
        {
            return true;
        }
    }
    return false;
}

// Returns true if the given point is colliding with the given circle.
static inline bool collisionCirclePoint(Circle c, MyVector2 p)
{
    return (distancePoints(c.origin, p) <= c.radius ? true : false);
}

// Returns true if the given circles are in collision.
static inline bool collisionCircles(Circle c1, Circle c2)
{
    return (distancePoints(c1.origin, c2.origin) <= c1.radius + c2.radius ? true : false);
}

#endif