#pragma once

#include "raylib.h"
#include <cmath>
#include <cassert>

using namespace std;


// ----------------- MY MATH FUNCTIONS -------------------- //

namespace MyMathLib
{

// ----------------------- DEFINES ----------------------- //

static bool __debug_shapes = false;
static bool __debug_bounding_boxes = false;
static bool __debug_axes = false;
static bool __debug_projections = false;
static bool __debug_failed_projections = false;
static bool __debug_points = false;

// --------------------- ARITHMECTIC ---------------------- //

namespace arithmetic
{
    // Rounds the given value to the nearest int.
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

    // Returns the sqare power of the given value.
    double sqpow(double val)
    {
        return val * val;
    }

    // Returns 1 if the given value is positive or null, and -1 if it is negative.
    int signOf(double val)
    {
        if (val == 0) return 1;
        return val / abs((int)val);
    }

    // Converts the given angle from degrees to radians.
    double degToRad(double deg)
    {
        return deg * (PI / 180.0f);
    }

    // Converts the given angle from radians to degrees.
    double radToDeg(double rad)
    {
        return rad * (180.0f / PI);
    }

    // Clamps the given value to be superior or equal to the minimum value and inferior or equal to the maximum value.
    double clamp(double val, double min, double max)
    {
        assert (min <= max);

        if (val < min) val = min;
        if (val > max) val = max;

        return val;
    }

    // Clamps the given value to be inferior or equal to the maximum value.
    double clampUnder(double val, double max)
    {
        if (val > max) val = max;

        return val;
    }

    // Clamps the given value to be superior or equal to the minimum value.
    double clampAbove(double val, double min)
    {
        if (val < min) val = min;

        return val;
    }

    // Remaps the given value from one range to another.
    double remap(double val, double inputStart, double inputEnd, double outputStart, double outputEnd)
    {
        // Source: https://stackoverflow.com/a/3451607/13858872
        // Find how far you are into the first range, scale that distance by the ratio of sizes of the ranges, and that's how far you should be into the second range.
        return outputStart + (val - inputStart) * (outputEnd - outputStart) / (inputEnd - inputStart);
    }
}

// ----------------------- GEOMETRY ----------------------- //

namespace geometry
{
    // Vector structure that holds values for x and y.
    class MyVector2 {
        private:
            MyVector2* self = this;
        
        public :
            // Attributes
            double x, y;

            // Builder
            MyVector2(double new_x, new_y)
            {
                x = new_x;
                y = new_y;
            };

        // Adds another vector to this vector.
        MyVector2 add(MyVector2 v)
        {
            x += v.x;
            y += v.y;
        }

        // Adds a value to a vector.
        MyVector2 addVal(double val)
        {
            x += val;
            y += val;
        }

        // Substracts a vector by another.
        MyVector2 substract(MyVector2 v)
        {
            x -= v.x;
            y -= v.y;
        }

        // Substracts a value from a vector.
        MyVector2 substractVal(double val)
        {
            x -= val;
            y -= val;
        }

        // Multiplies two vectors together.
        MyVector2 multiply(MyVector2 v)
        {
            x *= v.x;
            y *= v.y;
        }

        // Multiplies a vector by a value.
        MyVector2 multiplyVal(double val)
        {
            x *= val;
            y *= val;
        }

        // Divides a vector by another.
        MyVector2 divide(MyVector2 v)
        {
            x /= v.x;
            y /= v.y;
        }

        // Divides a vector by a value.
        MyVector2 divideVal(double val)
        {
            x /= val;
            y /= val;
        }

        // Returns the length of the given vector.
        double getLength()
        {
            return sqrt(sqpow(x) + sqpow(y));
        }

        // Returns the middle of the given vector.
        MyVector2 getMiddle()
        {
            return MyVector2(x / 2, y / 2);
        }

        // Normalizes the given vector so that its length is 1.
        void normalize()
        {
            x = x / self->getLength();
            y = y / self->getLength();
        }

        // Normalizes the given vector so that its length is 1.
        MyVector2 getNormalized()
        {
            return MyVector2(x / self->getLength(), y / self->getLength());
        }

        // Returns the angle (in radians) of the given vector.
        double getAngle()
        {
            return copysign(acos(self->getNormalized().x), asin(self->getNormalized().y));
        }

        // Modifies the length of the given vector to correspond to the given value.
        void setLength(double length)
        {
            *(self) = Vector2FromAngle(self->getAngle(), length);
        }

        // Negates both of the coordinates of the given vector.
        void negate()
        {
            *(self) = MyVector2(-x, -y);
        }

        // Copies the signs from the source vector to the destination vector.
        void copysign(MyVector2 source)
        {
            *(self) = MyVector2(copysign(x, source.x), copysign(y, source.y));
        }

        // Returns the dot product of the given vectors.
        double dotProduct(MyVector2 v)
        {
            return (x * v.x) + (y * v.y);
        }

        // Returns the cross product of the given vectors.
        double crossProduct(MyVector2 v)
        {
            return (x * v.y) - (y * v.x);
        }

        // Returns the normal of a given vector.
        MyVector2 getNormal()
        {
            return MyVector2(-y, x);
        }

        // Returns the angle (in radians) between two vectors.
        double angleWith(MyVector2 v)
        {
            // TODO: test this.
            double self_angle = self->getAngle();
            double v_angle = v.getAngle();
            return (self_angle >= v_angle ? (self_angle - v_angle) : (v_angle - self_angle));
        }

        // Rotates the given vector by the given angle.
        void rotate(double angle)
        {
            double self_length = self->getLength();
            double self_angle = self->getAngle();
            *(this) = Vector2FromAngle(self_angle + angle, self_length);
        }

        // Draws a vector at a certain origin point of a raylib window.
        static inline void draw(MyVector2 origin, Color color)
        {
            DrawLine(origin.x, origin.y, origin.x + x, origin.y + y, color);
            DrawPoly(toRayVec(origin.add(self))), 3, 4, radToDeg(self.getAngle() - PI/2), color);
        }

        // Draws a point in a raylib window.
        static inline void drawAsPoint(Color color)
        {
            DrawCircle(x, y, 2, color);
        }
    };

    // Segment structure that holds values for the starting point and the end point.
    class MySegment {
        private:
            // Self reference
            Segment *self = this;
        
        public :
            // Attributes
            MyVector2 a, b;

            // Builder
            Segment(MyVector2 new_a, MyVector2 new_b)
            {
                a = new_a;
                b = new_b;
            };

        // Returns the center of mass of a given segment.
        MyVector2 getCenterOfMass()
        {
            return MyVector2((a.x + b.x) / 2, (a.y +b.y) / 2);
        }

        // Returns the vertex of the given segment that corresponds to the given index.
        MyVector2 getVertex(int index)
        {
            assert (index < 2);

            switch (index)
            {
            case 0:
                return a;
            case 1:
                return b;
            default:
                // This should never happen thanks to the assert.
                return MyVector2(100000, 1000000);
            }
        }

        // Draws a segment in a raylib window.
        void draw(Color color)
        {
            DrawLine(a.x, a.y, b.x, b.y, color);
            DrawCircle(a.x, a.y, 2, color);
            DrawCircle(b.x, b.y, 2, color);
        }
    };

    // Triangle structure that holds values for 3 points.
    class MyTriangle {
        private:
            // Self reference
            MyTriangle *self = this;
        
        public:
            // Attributes
            MyVector2 a, b, c;

            // Builder
            MyTriangle(MyVector2 new_a, MyVector2 new_b, MyVector2 new_c)
            {
                a = new_a;
                b = new_b;
                c = new_c;

            }

        // Returns the center of mass of a given triangle.
            MyVector2 getCenterOfMass()
            {

            }

            }Bd
TrialMyVecor2a, MyVector2 b, MyVector2 c)/ Returns the side of the given tiangle that corresponds to the given index.
     MySeg(n{        {
            assert (index < 3);

            switch (index)
            {
            case 0:
                return MySegment(a, b);
            case 1:
                return MySegment(b, c);
            case 2:
                return MySegment(c, a);
            default:
                // This will never happen thanks to the assert.
                return MySegment(Vector2Zero(), Vector2Zero());
            }
        }

        // Returns the vertex of the given triangle that corresponds to the given index.
        MyVector2 getVertex(int index)
        {
            assert (index < 3);

            switch (index)
            {
            case 0:
                return a;
            case 1:
                return b;
            case 2:
                return c;
            default:
                // This should never happen thanks to the assert.
                return MyVector2(1000000, 1000000);
            }
        }

        // Draws a triangle in a raylib window.
        void draw(Color color)
        {
            DrawTriangleLines(toRayVec(a), toRayVec(b), toRayVec(c), color);
        }
    };

    // Rectangle structure that holds values for the origin point, width and height.
    class MyRectangle {
        private:
            // Self referance
            MyRectangle *self = this;
        
        public:
            // Attributes
            MyVector2     origin;
            double width, height;

            // Builder
           eight)
            {
                origin = new_origin;
                width = new_width;
                height = new_height;
            };

        // Returns the center of mass of a given rectangle.
        MyVector2 getCenterOfMass()
        {
            return MyVector2(origin.x + width / 2, origin.y + height / 2);
        }

        // Returns the side of the given rectangle that corresponds to the given index.
        MySegment getSide(int index)
        {
            assert (index < 4);

            switch (index)
            {
            case 0:
                return MySegment(MyVector2(origin.x + width, origin.y), origin);
            case 1:
                return MySegment(origin, MyVector2(origin.x, origin.y + height));
            case 2:
                return MySegment(MyVector2(origin.x, origin.y + height), MyVector2(origin.x + width, origin.y + height));
            case 3:
                return MySegment(MyVector2(origin.x + width, origin.y + height), MyVector2(origin.x + width, origin.y));
            default:
                // This should never happen thanks to the assert.
                return MySegment(Vector2Zero(), Vector2Zero());
            }
        }

        // Returns the vertex of the given rectangle that corresponds to the given index.
        MyVector2 getVertex(int index)
        {
            assert (index < 4);

            switch (index)
            {
            case 0:
                return MyVector2(origin.x + width, origin.y);
            case 1:
                return origin;
            case 2:
                return MyVector2(origin.x, origin.y + height);
            case 3:
                return MyVector2(origin.x + width, origin.y + height);
            default:
                // This should never happen thanks to the assert.
                return MyVector2(1000000, 1000000);
            }
        }

        // Draws a rectangle in a raylib window.
        void draw(Color color)
        {
            DrawRectangleLines(origin.x, origin.y, width, height, color);
        }
    };

    // Polygon structure that holds values for the origin point, the radius and the number of sides.
    class MyPolygon {
        private:
            // Self reference
            MyPolygon *self = this;
        
        public:
            MyVector2 origin;
            double radius, rotation;
            int sides;

        // Returns the center of mass of a given polygon.
        MyVector2 getCenterOfMass()
        {
            return origin;
        }

        // Returns the side of the given polygon that corresponds to the given index.
        MySegment getSide(int index)
        {
            assert(index < sides);

            double corner_angle = degToRad(360 / sides);
            double angle_offset = PI/2 + (index * corner_angle);

            MyVector2 poly_point_a = Vector2Add(Vector2FromAngle(angle_offset + rotation, radius), origin);
            MyVector2 poly_point_b = Vector2Add(Vector2FromAngle(angle_offset + corner_angle + rotation, radius), origin);

            return SegmentCreate(poly_point_a, poly_point_b);
        }

        // Returns the vertex of the given polygon that corresponds to the given index.
        MyVector2 getVertex(int index)
        {
            assert (index < sides);
            return getSide(self, index).a;
        }

        // Draws a polygon in a raylib window.
        void draw(Color color)
        {
            DrawPolyLines(toRayVec(origin), sides, radius, radToDeg(rotation), color);
        }
    };

    // Circle structure that holds values for the origin point and radius.
    class MyCircle {
        private:
            // Self reference
            MyCircle *self = this;
        
        public :
            MyVector2 origin;
            double radius;

        // Returns the center of mass of a given circle.
        MyVector2 getCenterOfMass()
        {
            return origin;
        }

        // Draws a circle in a raylib window.
        void draw(Color color)
        {
            DrawCircleLines(origin.x, origin.y, radius, color);
        }
    };

    // Union that can contain any shape.
    union Shape {
        MyVector2 vector;
        MySegment segment;
        MyTriangle triangle;
        MyRectangle rectangle;
        MyPolygon polygon;
        MyCircle circle;
    };

    // Shape types enum.
    enum class ShapeTypes {
        VECTOR2,
        SEGMENT,
        TRIANGLE,
        RECTANGLE,
        POLYGON,
        CIRCLE,
    };

    // Structure for shape info, holds shape type and data.
    class ShapeInfo {
        public :
            ShapeTypes type;
            Shape data;

        // Returns the center of mass of the given shape.
        MyVector2 getCenterOfMass()
        {
            assert(type != VECTOR2)

            switch (type)
            {
            case SEGMENT:
                return data.segment.getCenterOfMass();
            case TRIANGLE:
                return data.triangle.getCenterOfMass();
            case RECTANGLE:
                return data.rectangle.getCenterOfMass();
            case POLYGON:
                return data.polygon.getCenterOfMass();
            case CIRCLE:
                return data.circle.getCenterOfMass();
            default:
                break;
            }

            // This should never happen thanks to the assert.
            // It is there to shut the wrning up.
            return Vector2Zero();
        }

        // Returns the number of sides of a given shape.
        int getSidesNum()
        {
            switch (type)
            {
            case SEGMENT:
                return 1;
            case TRIANGLE:
                return 3;
            case RECTANGLE:
                return 4;
            case POLYGON:
                return data.polygon.sides;
            case CIRCLE:
                return 1;
            default:
                return 0;
            }
        }

        // Returns the side of the given shape that corresponds to the given index.
        // Returns a (0, 0) segment if the shape type is not supported (circle and vector).
        MySegment getSide(int index)
        {
            switch (type)
            {
            case SEGMENT:
                assert (index < 1);
                return data.segment;
            case TRIANGLE:
                return data.triangle.getSide(index);
            case RECTANGLE:
                return data.rectangle.getSide(index);
            case POLYGON:
                return data.polygon.getSide(index);
            default:
                return SegmentCreate(Vector2Zero(), Vector2Zero());
            }
        }

        // Returns the number of vertices of a given shape.
        int getVerticesNum()
        {
            switch (type)
            {
            case SEGMENT:
                return 2;
            case TRIANGLE:
                return 3;
            case RECTANGLE:
                return 4;
            case POLYGON:
                return data.polygon.sides;
            default:
                return 0;
            }
        }

        // Returns the vertex of the given shape that corresponds to the given index.
        MyVector2 ShapeGetVertex(int index)
        {
            switch (type)
            {
            case SEGMENT:
                return data.segment.getVertex(index);
            case TRIANGLE:
                return data.triangle.getVertex(index);
            case RECTANGLE:
                return data.rectangle.getVertex(index);
            case POLYGON:
                return data.polygon.getVertex(index);
            default:
                return MyVector2(1000000, 1000000);
            }
        }

        // Draws any shape in a raylib window.
        void draw(Color color, MyVector2 origin = Vector2Zero())
        {
            switch (type)
            {
            case VECTOR2:
                data.vector.draw(origin, color);
                break;
            case SEGMENT:
                data.segment.draw(color);
                break;
            case TRIANGLE:
                data.triangle.draw(color);
                break;
            case RECTANGLE:
                data.rectangle.draw(color);
                break;
            case POLYGON:
                data.polygon.draw(color);
                break;
            case CIRCLE:
                data.circle.draw(color);
                break;
            default:
                break;
            }
        }
    };

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

    // Returns a vector of coordinates { 0, 0 }.
    static inline MyVector2 Vector2Zero()
    {
        return MyVector2(0, 0);
    }

    // Creates a 2D vector from one point to another.
    static inline MyVector2 Vector2FromPoints(MyVector2 p1, MyVector2 p2)
    {
        return MyVector2(p2.x - p1.x, p2.y - p1.y);
    }

    // Creates a 2D vector given an angle and a length.
    static inline MyVector2 Vector2FromAngle(double rad, double length)
    {
        return MyVector2(cos(rad) * length, sin(rad) * length);
    }

    // Creates a 2D vector from a segement.
    static inline MyVector2 Vector2FromSegment(MySegment s)
    {
        return Vector2FromPoints(s.a, s.b);
    }

    // Creates a segment from one point to another.
    static inline MySegment SegmentCreate(MyVector2 a, MyVector2 b)
    {
        return (MySegment){a, b};
    }

    // Creates a segment given an origin point and a vector.
    static inline MySegment SegmentFromVector2(MyVector2 origin, MyVector2 v)
    {
        return SegmentCreate(origin, MyVector2(origin.x + v.x, origin.y + v.y));
    }

    // Creates a triangle given three points.
    static inline MyTriangle TriangleCreate(MyVector2 a, MyVector2 b, MyVector2 c)
    {
        return (MyTriangle){a, b, c};
    }

    // Creates a rectangle given an origin point, a width and a height.
    static inline MyRectangle RectangleCreate(MyVector2 origin, double width, double height)
    {
        return (MyRectangle){origin, width, height};
    }

    // Create a polygon given an origin point, a radius and a number of sides.
    static inline MyPolygon PolygonCreate(MyVector2 origin, double radius, double rotation, int sides)
    {
        return (MyPolygon){origin, radius, rotation, sides};
    }

    // Create a circle given an origin point and a radius.
    static inline MyCircle CircleCreate(MyVector2 origin, double radius)
    {
        return (MyCircle){origin, radius};
    }

    // ---------- MISC. FUNCTIONS ----------- //

    // Returns the distance between two points.
    static inline double distancePoints(MyVector2 p1, MyVector2 p2)
    {
        return Vector2FromPoints(p1, p2).getLength();
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

        // Create the min and max values for x and y.
        MyVector2 vertex = ShapeGetVertex(shape, 0);
        double xmin = vertex.x;
        double xmax = vertex.x;
        double ymin = vertex.y;
        double ymax = vertex.y;

        // Loop though the vertices and find the min and max values for x and y.
        for (int i = 1; i < vertices_num; i++)
        {
            vertex = ShapeGetVertex(shape, i);
            if (vertex.x <= xmin)
                xmin = vertex.x;
            if (vertex.x >= xmax)
                xmax = vertex.x;
            if (vertex.y <= ymin)
                ymin = vertex.y;
            if (vertex.y >= ymax)
                ymax = vertex.y;
        }

        // Create the shape's bounding box.
        MyRectangle bounding_box = RectangleCreate(MyVector2(xmin, ymin), xmax - xmin, ymax - ymin);

        //! Debug render.
        if (__debug_bounding_boxes) {
            DrawMyRectangle(bounding_box, GRAY);
        }

        return bounding_box;
    }

    // Returns an axis that passes through the center of the given circle and the center of the given shape.
    static inline MySegment CircleGetAxis(MyCircle circle, ShapeInfo shape)
    {
        // Make a segment that starts at the center of the circle, goes in the direction of the center of the shape and is of length 1.
        return SegmentFromVector2(circle.origin,
                                Vector2Normalize(Vector2FromPoints(circle.origin, ShapeCenterOfMass(shape))));
    }

    // Returns the axis of the given shapes that corresponds to the given index.
    static inline MySegment ShapesGetAxis(ShapeInfo shape1, ShapeInfo shape2, int index)
    {
        assert (index < getSidesNum(shape1) + getSidesNum(shape2));

        MySegment side;
        MySegment axis;

        // If the given index refers to an axis of the first shape...
        if (index < getSidesNum(shape1))
        {
            // If the first shape is not a circle, get the side pointed to by the index and calculate its normal.
            if (shape1.type != CIRCLE) {
                side = ShapeGetSide(shape1, index);
                axis = SegmentFromVector2(Vector2DivideVal(Vector2Add(side.a, side.b), 2),
                                            Vector2Normal(Vector2Normalize(Vector2FromSegment(side))));
            }
            // If the first shape is a circle, get its axis.
            else
                axis = CircleGetAxis(shape1.data.circle, shape2);
        }
        // If the given index refers to an axis of the second shape...
        else
        {
            // If the second shape is not a circle, get the side pointed to by the index and calculate its normal.
            if (shape2.type != CIRCLE) {
                side = ShapeGetSide(shape2, index - getSidesNum(shape1));
                axis = SegmentFromVector2(Vector2DivideVal(Vector2Add(side.a, side.b), 2),
                                            Vector2Normal(Vector2Normalize(Vector2FromSegment(side))));
            }
            // If the second shape is a circle, get its axis.
            else
                axis = CircleGetAxis(shape2.data.circle, shape1);
        }

        //! Debug render.
        if (__debug_axes) {
            DrawMyVector2(Vector2MultiplyVal(Vector2FromSegment(axis), 100), axis.a, BLUE);
        }

        return axis;
    }

    // Returns true if the given point is colliding with the given circle.
    static inline bool collisionCirclePoint(MyCircle c, MyVector2 p)
    {
        return (distancePoints(c.origin, p) <= c.radius ? true : false);
    }

    // Returns true if the given circles are in collision.
    static inline bool collisionCircles(MyCircle c1, MyCircle c2)
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
    template <typename T>
    static inline MySegment projectShapeOnAxis(MySegment axis, T shape)
    {
        // Get the axis' vector.
        MyVector2 axis_vec = Vector2FromSegment(axis);

        // Handle circles.
        if (T == MyCircle)
        {
            // Project the circle's origin onto the axis.
            MyVector2 origin_projection = axis.a.add(axis_vec.muliplyVal(Vector2DotProduct(Vector2FromPoints(axis.a, shape.data.circle.origin), axis_vec)));

            // Create a segment of the circle's projection.
            MySegment circle_projection = SegmentCreate(Vector2Substract(origin_projection, Vector2MultiplyVal(axis_vec, shape.data.circle.radius)),
                                                    Vector2Add      (origin_projection, Vector2MultiplyVal(axis_vec, shape.data.circle.radius)));

            //! Debug render.
            if (__debug_points) {
                DrawMyPoint(shape.data.circle.origin, WHITE);
                DrawMyPoint(Vector2Add      (origin_projection, Vector2MultiplyVal(axis_vec, shape.data.circle.radius)), SKYBLUE);
                DrawMyPoint(Vector2Substract(origin_projection, Vector2MultiplyVal(axis_vec, shape.data.circle.radius)), BLUE);
            }

            //! Debug render.
            if (__debug_projections) {
                DrawMySegment(circle_projection, ORANGE);
            }
            
            return circle_projection;
        }

        //* https://fr.wikipedia.org/wiki/Projection_orthogonale#Projet%C3%A9_orthogonal_sur_une_droite,_distance

        // Get all the vertices of the shape.
        int vertices_num = getVerticesNum(shape);
        MyVector2 vertex;
        MyVector2 projected_points[vertices_num];

        // Loop over the vertices of the shape and get their projections onto the axis.
        for (int i = 0; i < vertices_num; i++)
        {
            vertex = ShapeGetVertex(shape, i);
            projected_points[i] = Vector2Add(axis.a, Vector2MultiplyVal(axis_vec, Vector2DotProduct(Vector2FromPoints(axis.a, vertex), axis_vec)));

            //! Debug render.
            if (__debug_points) {
                DrawCircle(projected_points[i].x, projected_points[i].y, 2, WHITE);
            }
        }

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

        //! Debug render.
        if (__debug_points) {
            DrawCircle(min_point.x, min_point.y, 2, SKYBLUE);
            DrawCircle(max_point.x, max_point.y, 2, BLUE);
        }

        MyVector2 axis_orig_to_min_point = Vector2FromPoints(axis.a, min_point);
        MySegment projection = SegmentFromVector2(Vector2Add(axis.a, axis_orig_to_min_point), 
                                                Vector2FromPoints(min_point, max_point));

        //! Debug render.
        if (__debug_projections) {
            DrawMySegment(projection, ORANGE);
        }

        return projection;
    }

    // Returns true if the given point is colliding with the given segment.
    static inline bool collisionSegmentPoint(MySegment segment, MyVector2 point)
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
    static inline bool collisionProjections(MySegment projection1, MySegment projection2)
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
    template<typename T1, typename T2>
    static inline bool collisionSAT(T1 shape1, T2 shape2)
    {
        // If both shapes are circles, don't use SAT.
        if (T1 == MyCircle && T2 == MyCircle)
            return collisionCircles(shape1, shape2);

        // If both shapes are rectangles, don't use SAT.
        else if (T1 == MyRectangle && T2 == MyRectangle)
            return collisionAABB(shape1, shape2);

        // Check for collisions on the shapes' bounding boxes to not have to check if they are not in collision.
        else if (collisionAABB(getBoundingBox(shape1), getBoundingBox(shape2)))
        {
            // Get the number of sides of both shapes.
            int sides = shape1.getSidesNum() + shape2.getSidesNum();

            // Loop over all of the axes.
            for (int i = 0; i < sides; i++)
            {
                // Project both shapes onto the axis.
                MySegment projection1 = projectShapeOnAxis(ShapesGetAxis(shape1, shape2, i), shape1);
                MySegment projection2 = projectShapeOnAxis(ShapesGetAxis(shape1, shape2, i), shape2);

                // If the projections don't overlap, the shapes are not in collision.
                if (!collisionProjections(projection1, projection2))
                {
                    //! Debug render.
                    if (__debug_failed_projections) {
                        projection1.draw(PINK); 
                        projection2.draw(PINK);
                    }
                    return false;
                }
            }
            return true;
        }

        //! Debug render.
        if (__debug_shapes) {
            shape1.draw(GREEN); 
            shape2.draw(GREEN);
        }

        return false;
    }

}

}