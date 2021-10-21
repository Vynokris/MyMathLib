#pragma once
#ifndef MY_MATH_H
#define MY_MATH_H

#include "raylib.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>


// ---------- DEFINES ---------- //

static bool __debug_shapes = false;
static bool __debug_bounding_boxes = false;
static bool __debug_axes = false;
static bool __debug_projections = false;
static bool __debug_failed_projections = false;


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
    double rotation;
    int sides;
} Polygon;

// Circle structure that holds values for the origin point and radius.
typedef struct Circle {
    MyVector2 origin;
    double radius;
} Circle;

// Union that can contain any shape.
typedef union Shape {
    MyVector2 vector;
    Segment segment;
    Triangle triangle;
    MyRectangle rectangle;
    Polygon polygon;
    Circle circle;
} Shape;

// Shape types enum.
typedef enum ShapeTypes {
    VECTOR2,
    SEGMENT,
    TRIANGLE,
    RECTANGLE,
    POLYGON,
    CIRCLE,
} ShapeTypes;

typedef struct ShapeInfo {
    ShapeTypes type;
    Shape data;
} ShapeInfo;


// ---------- CONVERSIONS --------- //

// Converts a my_math 2D vector to a raylib 2D vector.
static inline Vector2
toRayVec(MyVector2 v)
{
    return (Vector2){v.x, v.y};
}

// Converts a my_math rectangle to a raylib rectangle.
static inline Rectangle toRayRec(MyRectangle r)
{
    return (Rectangle){r.origin.x, r.origin.y, r.width, r.height};
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


// ---------- OBJECT CREATION ---------- //

// Creates a 2D vector from two values.
static inline MyVector2 Vector2Create(double x, double y)
{
    return (MyVector2){x, y};
}

// Creates a 2D vector from one point to another.
static inline MyVector2 Vector2FromPoints(MyVector2 p1, MyVector2 p2)
{
    return Vector2Create(p2.x - p1.x, p2.y - p1.y);
}

// Creates a 2D vector given an angle and a length.
static inline MyVector2 Vector2FromAngle(double rad, double length)
{
    return Vector2Create(cos(rad) * length, sin(rad) * length);
}

// Creates a 2D vector from a segement.
static inline MyVector2 Vector2FromSegment(Segment s)
{
    return Vector2FromPoints(s.a, s.b);
}

// Creates a segment from one point to another.
static inline Segment SegmentCreate(MyVector2 a, MyVector2 b)
{
    return (Segment){a, b};
}

// Creates a segment given an origin point and a vector.
static inline Segment SegmentFromVector2(MyVector2 origin, MyVector2 v)
{
    return SegmentCreate(origin, Vector2Create(origin.x + v.x, origin.y + v.y));
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
static inline Polygon PolygonCreate(MyVector2 origin, double radius, double rotation, int sides)
{
    return (Polygon){origin, radius, rotation, sides};
}

// Create a circle given an origin point and a radius.
static inline Circle CircleCreate(MyVector2 origin, double radius)
{
    return (Circle){origin, radius};
}

// ---------- VECTOR MATH FUNCTIONS ---------- //

// Returns a vector of coordinates { 0, 0 }.
static inline MyVector2 Vector2Zero()
{
    return Vector2Create(0, 0);
}

// Adds two vectors together.
static inline MyVector2 Vector2Add(MyVector2 v1, MyVector2 v2)
{
    return Vector2Create(v1.x + v2.x, v1.y + v2.y);
}

// Adds a value to a vector.
static inline MyVector2 Vector2AddVal(MyVector2 v, double val)
{
    return Vector2Create(v.x + val, v.y + val);
}

// Substracts a vector by another.
static inline MyVector2 Vector2Substract(MyVector2 v1, MyVector2 v2)
{
    return Vector2Create(v1.x - v2.x, v1.y - v2.y);
}

// Substracts a value from a vector.
static inline MyVector2 Vector2SubstractVal(MyVector2 v, double val)
{
    return Vector2Create(v.x - val, v.y - val);
}

// Multiplies two vectors together.
static inline MyVector2 Vector2Multiply(MyVector2 v1, MyVector2 v2)
{
    return Vector2Create(v1.x * v2.x, v1.y * v2.y);
}

// Multiplies a vector by a value.
static inline MyVector2 Vector2MultiplyVal(MyVector2 v, double val)
{
    return Vector2Create(v.x * val, v.y * val);
}

// Divides a vector by another.
static inline MyVector2 Vector2Divide(MyVector2 v1, MyVector2 v2)
{
    return Vector2Create(v1.x / v2.x, v1.y / v2.y);
}

// Divides a vector by a value.
static inline MyVector2 Vector2DivideVal(MyVector2 v, double val)
{
    return Vector2Create(v.x / val, v.y / val);
}

// Returns the length of the given vector.
static inline double Vector2Length(MyVector2 v)
{
    return sqrt(sqpow(v.x) + sqpow(v.y));
}

// Returns the middle of the given vector
static inline MyVector2 Vector2Middle(MyVector2 v)
{
    return Vector2Create(v.x / 2, v.y / 2);
}

// Normalizes the given vector so that its length is 1.
static inline MyVector2 Vector2Normalize(MyVector2 v)
{
    return Vector2Create(v.x / Vector2Length(v), v.y / Vector2Length(v));
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
    return Vector2Create(-v.x, -v.y);
}

// Copies the signs from the source vector to the destination vector.
static inline MyVector2 Vector2Copysign(MyVector2 dest, MyVector2 source)
{
    return Vector2Create(copysign(dest.x, source.x), copysign(dest.y, source.y));
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

// Returns the normal of a given vector.
static inline MyVector2 Vector2Normal(MyVector2 vector)
{
    return Vector2Create(-vector.y, vector.x);
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


// ---------- MISC. FUNCTIONS ----------- //

// Returns the distance between two points.
static inline double distancePoints(MyVector2 p1, MyVector2 p2)
{
    return Vector2Length(Vector2FromPoints(p1, p2));
}

// Returns the side of the given polygon that corresponds to the given index.
static inline Segment PolygonGetSide(Polygon poly, int index)
{
    assert(index < poly.sides);

    double corner_angle = degToRad(360 / poly.sides);
    double angle_offset = PI/2 + (index * corner_angle);

    MyVector2 poly_point_a = Vector2Add(Vector2FromAngle(angle_offset + poly.rotation, poly.radius), poly.origin);
    MyVector2 poly_point_b = Vector2Add(Vector2FromAngle(angle_offset + corner_angle + poly.rotation, poly.radius), poly.origin);

    return SegmentCreate(poly_point_a, poly_point_b);
}


// ---------- CENTER OF MASS ---------- //

// Returns the center of mass of a given segment.
static inline MyVector2 SegmentCenterOfMass(Segment segment)
{
    return Vector2Create((segment.a.x + segment.b.x) / 2, (segment.a.y + segment.b.y) / 2);
}

// Returns the center of mass of a given triangle.
static inline MyVector2 TriangleCenterOfMass(Triangle triangle)
{
    return Vector2Create((triangle.a.x + triangle.b.x + triangle.c.x) / 3, (triangle.a.y + triangle.b.y + triangle.c.y) / 3);
}

// Returns the center of mass of a given rectangle.
static inline MyVector2 RectangleCenterOfMass(MyRectangle rectangle)
{
    return Vector2Create(rectangle.origin.x + rectangle.width / 2, rectangle.origin.y + rectangle.height / 2);
}

// Returns the center of mass of a given polygon.
static inline MyVector2 PolygonCenterOfMass(Polygon poly)
{
    return poly.origin;
}

// Returns the center of mass of a given circle.
static inline MyVector2 CircleCenterOfMass(Circle circle)
{
    return circle.origin;
}

// Returns the center of mass of the given shape (returns a vector of coordinates 1,000,000 if the shape type isn't supported).
static inline MyVector2 ShapeCenterOfMass(ShapeInfo shape)
{
    switch (shape.type)
    {
    case SEGMENT:
        return SegmentCenterOfMass(shape.data.segment);
    case TRIANGLE:
        return TriangleCenterOfMass(shape.data.triangle);
    case RECTANGLE:
        return RectangleCenterOfMass(shape.data.rectangle);
    case POLYGON:
        return PolygonCenterOfMass(shape.data.polygon);
    case CIRCLE:
        return CircleCenterOfMass(shape.data.circle);
    default:
        return Vector2Create(1000000, 1000000);
    }
}


// ---------- DRAWING FUNCTIONS ---------- //

// Draws a point in a raylib window.
static inline void DrawPoint(MyVector2 p, Color color)
{
    DrawCircle(p.x, p.y, 2, color);
}

// Draws a vector at a certain origin point of a raylib window.
static inline void DrawVector2(MyVector2 v, MyVector2 origin, Color color)
{
    DrawLine(origin.x, origin.y, origin.x + v.x, origin.y + v.y, color);
    DrawPoly(toRayVec(Vector2Add(origin, v)), 3, 4, radToDeg(Vector2GetAngle(v) - PI/2), color);
}

// Draws a segment in a raylib window.
static inline void DrawSegment(Segment s, Color color)
{
    DrawLine(s.a.x, s.a.y, s.b.x, s.b.y, color);
    DrawCircle(s.a.x, s.a.y, 2, color);
    DrawCircle(s.b.x, s.b.y, 2, color);
}

// Draws a triangle in a raylib window.
static inline void DrawMyTriangle(Triangle t, Color color)
{
    DrawTriangleLines(toRayVec(t.a), toRayVec(t.b), toRayVec(t.c), color);
}

// Draws a rectangle in a raylib window.
static inline void DrawMyRectangle(MyRectangle r, Color color)
{
    DrawRectangleLines(r.origin.x, r.origin.y, r.width, r.height, color);
}

// Draws a polygon in a raylib window.
static inline void DrawMyPolygon(Polygon poly, Color color)
{
    DrawPolyLines(toRayVec(poly.origin), poly.sides, poly.radius, radToDeg(poly.rotation), color);
}

// Draws a circle in a raylib window.
static inline void DrawMyCircle(Circle c, Color color)
{
    DrawCircleLines(c.origin.x, c.origin.y, c.radius, color);
}

// Draws any shape in a raylib window.
static inline void DrawShape(ShapeInfo shape, MyVector2 origin, Color color)
{
    switch (shape.type)
    {
    case VECTOR2:
        DrawVector2(shape.data.vector, origin, color);
        break;
    case SEGMENT:
        DrawSegment(shape.data.segment, color);
        break;
    case TRIANGLE:
        DrawMyTriangle(shape.data.triangle, color);
        break;
    case RECTANGLE:
        DrawMyRectangle(shape.data.rectangle, color);
        break;
    case POLYGON:
        DrawMyPolygon(shape.data.polygon, color);
        break;
    case CIRCLE:
        DrawMyCircle(shape.data.circle, color);
        break;
    default:
        break;
    }
}


// ---------- SHAPE VERTICES AND SIDES ---------- //

// Returns the number of sides of a given shape (returns 2 for rectangles).
static inline int getSidesNum(ShapeInfo shape)
{
    switch (shape.type)
    {
    case SEGMENT:
        return 1;
    case TRIANGLE:
        return 3;
    case RECTANGLE:
        return 2; // There are only two axes to check for collision in a rectangle.
    case POLYGON:
        return shape.data.polygon.sides;
    case CIRCLE:
        return 1;
    default:
        return 0;
    }
}

// Returns the sides of a given shape as segments.
static inline Segment *getSides(ShapeInfo shape)
{
    Segment *sides;

    switch (shape.type)
    {
    case SEGMENT:
        sides = malloc(sizeof(Segment));
        sides[0] = shape.data.segment;
        break;

    case TRIANGLE:
        sides = malloc(sizeof(Segment) * 3);
        sides[0] = SegmentCreate(shape.data.triangle.a, shape.data.triangle.b);
        sides[1] = SegmentCreate(shape.data.triangle.b, shape.data.triangle.c);
        sides[2] = SegmentCreate(shape.data.triangle.c, shape.data.triangle.a);
        break;

    case RECTANGLE: // There are only two axes to check for collision in a rectangle.
        sides = malloc(sizeof(Segment) * 2);
        sides[0] = SegmentCreate(shape.data.rectangle.origin,
                                 Vector2Create(shape.data.rectangle.origin.x + shape.data.rectangle.width, shape.data.rectangle.origin.y));
        sides[1] = SegmentCreate(shape.data.rectangle.origin,
                                 Vector2Create(shape.data.rectangle.origin.x, shape.data.rectangle.origin.y + shape.data.rectangle.height));
        break;

    case POLYGON:
        sides = malloc(sizeof(Segment) * shape.data.polygon.sides);
        for (int i = 0; i < shape.data.polygon.sides; i++)
        {
            sides[i] = PolygonGetSide(shape.data.polygon, i);
        }
        break;

    default:
        break;
    }

    return sides;
}

// Returns the number of vertices of a given shape.
static inline int getVerticesNum(ShapeInfo shape)
{
    int vertices_num;

    switch (shape.type)
    {
    case SEGMENT:
        return 2;
    case TRIANGLE:
        return 3;
    case RECTANGLE:
        return 4;
    case POLYGON:
        return shape.data.polygon.sides;
    default:
        return 0;
    }

    return vertices_num;
}

// Returns the vertices of a given shape as vectors.
static inline MyVector2 *getVertices(ShapeInfo shape)
{
    MyVector2 *vertices;

    switch (shape.type)
    {
    case SEGMENT:
        vertices = malloc(sizeof(MyVector2) * 2);
        vertices[0] = shape.data.segment.a;
        vertices[1] = shape.data.segment.b;
        break;
    case TRIANGLE:
        vertices = malloc(sizeof(MyVector2) * 3);
        vertices[0] = shape.data.triangle.a;
        vertices[1] = shape.data.triangle.b;
        vertices[2] = shape.data.triangle.c;
        break;
    case RECTANGLE:
        vertices = malloc(sizeof(MyVector2) * 4);
        vertices[0] = shape.data.rectangle.origin;
        vertices[1] = Vector2Create(shape.data.rectangle.origin.x + shape.data.rectangle.width, shape.data.rectangle.origin.y);
        vertices[2] = Vector2Create(shape.data.rectangle.origin.x + shape.data.rectangle.width, shape.data.rectangle.origin.y + shape.data.rectangle.height);
        vertices[3] = Vector2Create(shape.data.rectangle.origin.x, shape.data.rectangle.origin.y + shape.data.rectangle.height);
        break;
    case POLYGON:
        vertices = malloc(sizeof(MyVector2) * shape.data.polygon.sides);
        for (int i = 0; i < shape.data.polygon.sides; i++)
        {
            vertices[i] = PolygonGetSide(shape.data.polygon, i).a;
        }
    default:
        break;
    }

    return vertices;
}


// ---------- COLLISIONS ---------- //

// Returns the smallest rectangle that contanins the given shape.
static inline MyRectangle getBoundingBox(ShapeInfo shape)
{
    // If the shape is a circle.
    if (shape.type == CIRCLE)
    {
        //! Debug render.
        if (__debug_bounding_boxes) {
            DrawMyRectangle(RectangleCreate(Vector2SubstractVal(shape.data.circle.origin, shape.data.circle.radius), shape.data.circle.radius * 2, shape.data.circle.radius * 2), GRAY);
        }

        return RectangleCreate(Vector2SubstractVal(shape.data.circle.origin, shape.data.circle.radius), shape.data.circle.radius * 2, shape.data.circle.radius * 2);
    }

    // Get the shape's vertices information.
    int vertices_num = getVerticesNum(shape);
    MyVector2 *vertices = getVertices(shape);

    // Create the min and max values for x and y.
    double xmin = vertices[0].x;
    double xmax = vertices[0].x;
    double ymin = vertices[0].y;
    double ymax = vertices[0].y;

    // Loop though the vertices and find the min and max values for x and y.
    for (int i = 0; i < vertices_num; i++)
    {
        if (vertices[i].x <= xmin)
            xmin = vertices[i].x;
        if (vertices[i].x >= xmax)
            xmax = vertices[i].x;
        if (vertices[i].y <= ymin)
            ymin = vertices[i].y;
        if (vertices[i].y >= ymax)
            ymax = vertices[i].y;
    }

    free(vertices);

    //! Debug render.
    if (__debug_bounding_boxes) {
        DrawMyRectangle(RectangleCreate(Vector2Create(xmin, ymin), xmax - xmin, ymax - ymin), GRAY);
    }

    return RectangleCreate(Vector2Create(xmin, ymin), xmax - xmin, ymax - ymin);
}

// Returns an axis that passes through the center of the given circle and the center of the given shape.
static inline Segment *getAxisCircleShape(Circle circle, ShapeInfo shape)
{
    Segment* axis = malloc(sizeof(Segment));
    axis[0] = SegmentFromVector2(circle.origin,
                                 Vector2Normalize(Vector2FromPoints(circle.origin, ShapeCenterOfMass(shape))));
        
    //! Debug render.
    if (__debug_axes) {
        DrawVector2(Vector2MultiplyVal(Vector2FromSegment(axis[0]), 100), axis[0].a, BLUE);
    }

    return axis;
}

// Get all the axes of two given shapes.
static inline Segment *getAxes(ShapeInfo shape1, ShapeInfo shape2)
{
    if (shape1.type == CIRCLE) {
        return getAxisCircleShape(shape1.data.circle, shape2);
    }
    if (shape2.type == CIRCLE) {
        return getAxisCircleShape(shape2.data.circle, shape1);
    }

    // Get all the sides of the given shape.
    Segment* sides1 = getSides(shape1);
    Segment* sides2 = getSides(shape2);
    Segment* axes = malloc(sizeof(Segment) * (getSidesNum(shape1) + getSidesNum(shape2)));

    // Loop over them and get their normals to get the axes.
    for (int i = 0; i < getSidesNum(shape1) + getSidesNum(shape2); i++)
    {
        if (i < getSidesNum(shape1)) {
            axes[i] = SegmentFromVector2(Vector2DivideVal(Vector2Add(sides1[i].a, sides1[i].b), 2),
                                         Vector2Normal(Vector2Normalize(Vector2FromSegment(sides1[i]))));
        }
        else {
            axes[i] = SegmentFromVector2(Vector2DivideVal(Vector2Add(sides2[i-getSidesNum(shape1)].a, sides2[i-getSidesNum(shape1)].b), 2),
                                         Vector2Normal(Vector2Normalize(Vector2FromSegment(sides2[i-getSidesNum(shape1)]))));
        }
        
        //! Debug render.
        if (__debug_axes) {
            DrawVector2(Vector2MultiplyVal(Vector2FromSegment(axes[i]), 100), axes[i].a, BLUE);
        }
    }
    free(sides1);
    free(sides2);
    return axes;
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

// Checks for collision between two rectangles.
static inline bool collisionAABB(MyRectangle rec1, MyRectangle rec2)
{
    if (rec1.origin.x + rec1.width >= rec2.origin.x &&
        rec1.origin.x <= rec2.origin.x + rec2.width &&
        rec1.origin.y + rec1.height >= rec2.origin.y &&
        rec1.origin.y <= rec2.origin.y + rec2.height) 
    {
        return true;
    }
    else {
        return false;
    }
}

// Project a shape onto a given axis.
static inline Segment projectShapeOnAxis(Segment axis, ShapeInfo shape)
{
    // Get the axis' vector.
    MyVector2 axis_vec = Vector2FromSegment(axis);

    // Handle circles.
    if (shape.type == CIRCLE)
    {
        // Project the circle's origin onto the axis.
        MyVector2 origin_projection = Vector2Add(axis.a, Vector2MultiplyVal(axis_vec, Vector2DotProduct(Vector2FromPoints(axis.a, shape.data.circle.origin), axis_vec)));

        // Create a segment of the circle's projection.
        Segment circle_projection = SegmentCreate(Vector2Substract(origin_projection, Vector2MultiplyVal(axis_vec, shape.data.circle.radius)),
                                                  Vector2Add      (origin_projection, Vector2MultiplyVal(axis_vec, shape.data.circle.radius)));

        //! Debug render.
        if (__debug_projections) {
            DrawSegment(circle_projection, ORANGE);
        }
        
        return circle_projection;
    }

    //* https://fr.wikipedia.org/wiki/Projection_orthogonale#Projet%C3%A9_orthogonal_sur_une_droite,_distance

    // Get all the vertices of the shape.
    int vertices_num = getVerticesNum(shape);
    MyVector2 *vertices = getVertices(shape);

    MyVector2 projected_points[vertices_num];

    // Loop over the vertices of the shape and get their projections onto the axis.
    for (int i = 0; i < vertices_num; i++)
    {
        projected_points[i] = Vector2Add(axis.a, Vector2MultiplyVal(axis_vec, Vector2DotProduct(Vector2FromPoints(axis.a, vertices[i]), axis_vec)));
    }
    free(vertices);

    // Find the closest and farthest points from the axis origin.
    MyVector2 min_point = projected_points[0];
    MyVector2 max_point = min_point;

    for (int i = 0; i < vertices_num; i++)
    {
        if (Vector2Copysign(projected_points[i], axis_vec).x > Vector2Copysign(max_point, axis_vec).x ||
            Vector2Copysign(projected_points[i], axis_vec).y > Vector2Copysign(max_point, axis_vec).y)
        {
            max_point = projected_points[i];
        }

        if (Vector2Copysign(projected_points[i], axis_vec).x < Vector2Copysign(min_point, axis_vec).x ||
            Vector2Copysign(projected_points[i], axis_vec).y < Vector2Copysign(min_point, axis_vec).y)
        {
            min_point = projected_points[i];
        }
    }

    MyVector2 axis_orig_to_min_point = Vector2FromPoints(axis.a, min_point);

    //! Debug render.
    if (__debug_projections) {
        DrawSegment(SegmentFromVector2(Vector2Add(axis.a, axis_orig_to_min_point), Vector2FromPoints(min_point, max_point)), ORANGE);
    }

    return SegmentFromVector2(Vector2Add(axis.a, axis_orig_to_min_point),
                            Vector2FromPoints(min_point, max_point));
}

// Returns true if the given point is colliding with the given segment.
static inline bool collisionSegmentPoint(Segment segment, MyVector2 point)
{
    if (roundInt(Vector2CrossProduct(Vector2FromSegment(segment), Vector2FromPoints(segment.a, point))) == 0)
    {
        if ((point.x >= segment.a.x && point.x <= segment.b.x) || (point.y >= segment.a.y && point.y <= segment.b.y) ||
            (point.x <= segment.a.x && point.x >= segment.b.x) || (point.y <= segment.a.y && point.y >= segment.b.y))
        {
            return true;
        }
    }
    return false;
}

// Returns true if the given projections are colliding each others
static inline bool collisionProjections(Segment projection1, Segment projection2)
{
    if (collisionSegmentPoint(projection1, projection2.a) ||
        collisionSegmentPoint(projection1, projection2.b) ||
        collisionSegmentPoint(projection2, projection1.a) ||
        collisionSegmentPoint(projection2, projection1.b))
    {
        return true;
    }
    return false;
}

// Checks for collision between two given shapes.
static inline bool collisionSAT(ShapeInfo shape1, ShapeInfo shape2)
{
    // If both the shapes are, circles, don't use SAT.
    if (shape1.type == CIRCLE && shape2.type == CIRCLE)
        collisionCircles(shape1.data.circle, shape2.data.circle);

    // Check for collisions on the shapes' bounding boxes to not have to check if they are not in collision.
    if (collisionAABB(getBoundingBox(shape1), getBoundingBox(shape2)))
    {
        //! Debug render.
        if (__debug_shapes) {
            DrawShape(shape1, Vector2Zero(), GREEN); DrawShape(shape2, Vector2Zero(), GREEN);
        }

        // Get the number of sides of both shapes.
        int sides = getSidesNum(shape1) + getSidesNum(shape2);
        if (shape1.type == CIRCLE || shape2.type == CIRCLE) sides = 1;

        // Get the axes for both shapes.
        Segment *axes = getAxes(shape1, shape2);

        // Loop over all of the axes.
        for (int i = 0; i < sides; i++)
        {
            // Project both shapes onto the axis.
            Segment projection1 = projectShapeOnAxis(axes[i], shape1);
            Segment projection2 = projectShapeOnAxis(axes[i], shape2);

            // If the projections don't overlap, the shapes are not in collision.
            if (!collisionProjections(projection1, projection2))
            {
                //! Debug render.
                if (__debug_failed_projections) {
                    DrawSegment(projection1, PINK); DrawSegment(projection2, PINK);
                }

                free(axes);
                return false;
            }
        }
        free(axes);
        return true;
    }
    //! Debug render.
    if (__debug_shapes) {
        DrawShape(shape1, Vector2Zero(), GREEN); DrawShape(shape2, Vector2Zero(), GREEN);
    }

    return false;
}

#endif