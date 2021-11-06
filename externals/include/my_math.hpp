/*************************************************************
 * MyMathLib, a custom matmematic library written by
 * Rémi Serra and Alexandre Perché, students at ISART Digital.
 *************************************************************/

// ---------------- INCLUDES / NAMESPACE ------------------ //

#pragma once
#include "raylib.h"
#include <cmath>
#include <cassert>
#include <stdint.h>
#include <typeinfo>

using namespace std;

namespace MyMathLib
{
// --------------------- ARITHMECTIC ---------------------- //

namespace arithmetic
{
    // Rounds the given value to the nearest int.
    int roundInt(double val)    { return (int)round(val); }

    // Rounds down the given value.
    int floorInt(double val)    { return (int)floor(val); }

    // Rounds up the given value.
    int ceilInt(double val)     { return (int)ceil(val); }

    // Returns the sqare power of the given value.
    double sqpow(double val)    { return val * val; }

    // Returns 1 if the given value is positive or null, and
    // -1 if it is negative.
    int signOf(double val)      { if (val == 0) return 1; return val / abs((int)val); }

    // Converts the given angle from degrees to radians.
    double degToRad(double deg) { return deg * (PI / 180.0f); }

    // Converts the given angle from radians to degrees.
    double radToDeg(double rad) { return rad * (180.0f / PI); }

    // Clamps the given value to be superior or equal to the 
    // minimum value and inferior or equal to the maximum value.
    double clamp(double val, double min, double max)
    {
        assert (min <= max); 
        if (val < min) val = min;
        if (val > max) val = max;
        return val;
    }

    // Clamps the given value to be inferior or equal to the maximum value.
    double clampUnder(double val, double max) { if (val > max) val = max; return val; }

    // Clamps the given value to be superior or equal to the minimum value.
    double clampAbove(double val, double min) { if (val < min) val = min; return val; }

    // Remaps the given value from one range to another.
    double remap(double val, double inputStart, double inputEnd, double outputStart, double outputEnd)
    {
        return outputStart + (val - inputStart) * (outputEnd - outputStart) / (inputEnd - inputStart);
    }

    // Return a random integer according to a seed and a state
    int32_t getNextRandomNumber(int32_t state)
    {
        uint64_t result = (((int64_t)1103515245) * state + 12345) % 2147483648;
        return result;
    }

    // Return a random value according to a range
    template<typename T>
    T RandomValue(T start, T end)
    {
        assert(start > end);
        T r = (T)GetRandomNumberFromRange(0);
        for (int i = 0; r >= start && r <= end; i++) r = (T)GetRandomNumberFromRange(i);
        return r;
    }
}

// ----------------------- GEOMETRY ----------------------- //

namespace geometry
{
    // ------------- DEFINES ------------- //

    static bool __debug_shapes             = false;
    static bool __debug_bounding_boxes     = false;
    static bool __debug_axes               = false;
    static bool __debug_projections        = false;
    static bool __debug_failed_projections = false;
    static bool __debug_points             = false;
    

    // ------------- FORWARD DECLARATIONS ------------- //

    class MyVector2;
    class MySegment;
    class MyTriangle;
    class MyRectangle;
    class MyPolygon;
    class MyCircle;

    // Makes sure that given object is a shape.
    template <typename T>
    bool isShape(T object, bool notVector2 = false)
    {
        return ((notVector2 ? true : typeid(T) == typeid(MyVector2))  ||
                                     typeid(T) == typeid(MySegment)   ||
                                     typeid(T) == typeid(MyTriangle)  ||
                                     typeid(T) == typeid(MyRectangle) ||
                                     typeid(T) == typeid(MyPolygon)   ||
                                     typeid(T) == typeid(MyCircle));
    }
    

    // ------------- CLASS CREATION FUNCTIONS ------------- //

    // Returns a vector of coordinates { 0, 0 }.
    MyVector2 Vector2Zero()                                     { return MyVector2(0, 0); }

    // Creates a 2D vector from one point to another.
    MyVector2 Vector2FromPoints(MyVector2 p1, MyVector2 p2)     { return MyVector2(p2.x - p1.x, p2.y - p1.y); }

    // Creates a 2D vector given an angle and a length.
    MyVector2 Vector2FromAngle(double rad, double length)       { return MyVector2(cos(rad) * length, sin(rad) * length); }

    // Creates a 2D vector from a segement.
    MyVector2 Vector2FromSegment(MySegment s)                   { return Vector2FromPoints(s.a, s.b); }

    // Creates a segment given an origin point and a vector.
    MySegment SegmentFromVector2(MyVector2 origin, MyVector2 v) { return MySegment(origin, MyVector2(origin.x + v.x, origin.y + v.y)); }

    // ---------------------- CLASSES --------------------- //

    // Vector class that holds values for x and y (2 dimensions).
    class MyVector2 {
        private:
            MyVector2* self = this;
        
        public :
            // Attributes.
            double x, y;

            // Builder.
            MyVector2(double new_x, double new_y) : x(new_x), y(new_y){};

        // ---------- VECTOR2 METHODS ---------- //

        // Vector2 addition.
        template <typename T>
        MyVector2 operator+(T val)
        {
            if (typeid(T) == typeid(MyVector2) || 
                typeid(T) == typeid(int)       || 
                typeid(T) == typeid(float)     || 
                typeid(T) == typeid(double));
            {
                if (typeid(T) == typeid(MyVector2)) return MyVector2(x + val.x, y + val.y);
                else                                return MyVector2(x + val,   y + val);
            }
        }

        // Vector2 substraction.
        template <typename T>
        MyVector2 operator-(T val)
        {
            if (typeid(T) == typeid(MyVector2) || 
                typeid(T) == typeid(int)       || 
                typeid(T) == typeid(float)     || 
                typeid(T) == typeid(double));
            {
                if (typeid(T) == typeid(MyVector2))  return MyVector2(x - val.x, y - val.y);
                else                                 return MyVector2(x - val,   y - val);
            }
        }

        // Vector2 multiplication.
        template <typename T>
        MyVector2 operator*(T val)
        {
            if (typeid(T) == typeid(MyVector2) || 
                typeid(T) == typeid(int)       || 
                typeid(T) == typeid(float)     || 
                typeid(T) == typeid(double));
            {
                if (typeid(T) == typeid(MyVector2)) return MyVector2(x * val.x, y * val.y);
                else                                return MyVector2(x * val,   y * val);
            }
        }

        // Vector2 division.
        template <typename T>
        MyVector2 operator/(T val)
        {
            if (typeid(T) == typeid(MyVector2) || 
                typeid(T) == typeid(int)       || 
                typeid(T) == typeid(float)     || 
                typeid(T) == typeid(double));
            {
                if (typeid(T) == typeid(MyVector2)) return MyVector2(x / val.x, y / val.y);
                else                                return MyVector2(x / val,   y / val);
            }
        }

        // Vector2 dot product.
        double operator&(MyVector2 val) { return (x * val.x) + (y * val.y); }

        // Vector2 cross product.
        double operator^(MyVector2 val) { return (x * val.y) - (y * val.x); }

        // Returns the length of the given vector.
        double getLength()              { return sqrt(arithmetic::sqpow(x) + arithmetic::sqpow(y)); }

        // Returns the angle (in radians) of the given vector.
        double getAngle()               { return std::copysign(std::acos(self->getNormalized().x), std::asin(self->getNormalized().y)); }

        // Returns the middle of the given vector.
        MyVector2 getMiddle()           { return MyVector2(x / 2, y / 2); }

        // Normalizes the given vector so that its length is 1.
        void normalize()                { x /= self->getLength(); y /= self->getLength(); }

        // Normalizes the given vector so that its length is 1.
        MyVector2 getNormalized()       { return MyVector2(x / self->getLength(), y / self->getLength()); }

        // Modifies the length of the given vector to correspond to the given value.
        void setLength(double length)   { *(self) = Vector2FromAngle(self->getAngle(), length); }

        // Negates both of the coordinates of the given vector.
        void negate()                   { *(self) = MyVector2(-x, -y); }

        // Copies the signs from the source vector to the destination vector.
        void copysign(MyVector2 source) { *(self) = MyVector2(std::copysign(x, source.x), std::copysign(y, source.y)); }

        // Copies the signs from the source vector to the destination vector.
        MyVector2 getCopiedSign(MyVector2 source) { return MyVector2(std::copysign(x, source.x), std::copysign(y, source.y)); }

        // Returns the normal of a given vector.
        MyVector2 getNormal()           { return MyVector2(-y, x); }

        // Returns the angle (in radians) between two vectors.
        double angleWith(MyVector2 v)
        {
            double self_angle = self->getAngle();
            double v_angle    = v.getAngle();
            return (self_angle >= v_angle ? (self_angle - v_angle) : (v_angle - self_angle));
        }

        // Rotates the given vector by the given angle.
        void rotate(double angle)
        {
            double self_length = self->getLength();
            double self_angle  = self->getAngle();
            *(self) = Vector2FromAngle(self_angle + angle, self_length);
        }

        // Draws a vector at a certain origin point of a raylib window.
        void draw(MyVector2 origin, Color color)
        {
            DrawLine(origin.x, origin.y, origin.x + x, origin.y + y, color);
            MyVector2 tmp = origin + *self;
            DrawPoly(tmp.toRayVec(), 3, 4, arithmetic::radToDeg(self->getAngle() - PI / 2), color);
        }

        // Draws a point in a raylib window.
        void drawAsPoint(Color color) { DrawCircle(x, y, 2, color); }

        // Converts a my_math 2D vector to a raylib 2D vector.
        Vector2 toRayVec() { return (Vector2){x, y}; }
    };

    // Segment structure that holds values for the starting point and the end point.
    class MySegment {
        private:
            // Self reference.
            MySegment *self = this;
        
        public :
            // Attributes.
            MyVector2 a, b;

            // Builder.
            MySegment(MyVector2 new_a, MyVector2 new_b) : a(new_a), b(new_b){};

        // ---------- SEGMENT METHODS ---------- //

        // Returns the center of mass of the segment.
        MyVector2 getCenterOfMass()
        {
            return MyVector2((a.x + b.x) / 2, (a.y +b.y) / 2);
        }

        // Returns the number of sides of the segment.
        int getSidesNum() { return 1; }

        // Returns the side of the segment that corresponds to the given index.
        MySegment getSide(int index)
        {
            assert (0 <= index && index < 1);
            return *self;
        }

        // Returns the number of vertices of the segment.
        int getVerticesNum() { return 2; }

        // Returns the vertex of the segment that corresponds to the given index.
        MyVector2 getVertex(int index)
        {
            assert (0 <= index && index < 2);

            switch (index)
            {
                case 0:  return a;
                case 1:  return b;
                default: break;
            }
            return Vector2Zero();
        }

        // Draws the segment in a raylib window.
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
            // Self reference.
            MyTriangle *self = this;
        
        public:
            // Attributes.
            MyVector2 a, b, c;

            // Builder.
            MyTriangle(MyVector2 new_a, MyVector2 new_b, MyVector2 new_c) : a(new_a), b(new_b), c(new_c){};

        // ---------- TRIANGLE METHODS ---------- //
        
        // Returns the center of mass of the triangle.
        MyVector2 getCenterOfMass() { return MyVector2((a.x + b.x + c.x) / 3, (a.y + b.y + c.y) / 3); }

        // Returns the number of sides of the triangle.
        int getSidesNum() { return 3; }

        // Returns the side of the triangle that corresponds to the given index.
        MySegment getSide(int index) 
        {
            assert (0<= index && index < 3);

            switch (index)
            {
                case 0:   return MySegment(a, b);
                case 1:   return MySegment(b, c);
                case 2:   return MySegment(c, a);
                default:  break;
            }

            return MySegment(Vector2Zero(), Vector2Zero());
        }

        // Returns the number of vertices of the triangle.
        int getVerticesNum() { return 3; }

        // Returns the vertex of the triangle that corresponds to the given index.
        MyVector2 getVertex(int index)
        {
            assert (0 <= index && index < 3);

            switch (index)
            {
            case 0:  return a;
            case 1:  return b;
            case 2:  return c;
            default: break;
            }

            return Vector2Zero();
        }

        // Draws the triangle in a raylib window.
        void draw(Color color) { DrawTriangleLines(a.toRayVec(), b.toRayVec(), c.toRayVec(), color); }
    };

    // Rectangle structure that holds values for the origin point, width and height.
    class MyRectangle {
        private:
            // Self referance.
            MyRectangle *self = this;
        
        public:            
            // Attributes.
            MyVector2 origin;
            double width, height;

            // Builder.
            MyRectangle(MyVector2 new_origin, double new_width, double new_height) : origin(new_origin), width(new_width), height(new_height){};

        // Returns the center of mass of the rectangle.
        MyVector2 getCenterOfMass() { return MyVector2(origin.x + width / 2, origin.y + height / 2); };

        // Returns the number of sides of the rectangle.
        int getSidesNum() { return 4; }

        // Returns the side of the rectangle that corresponds to the given index.
        MySegment getSide(int index)
        {
            assert (0 <= index && index < 4);

            switch (index)
            {
                case 0:  return MySegment(MyVector2(origin.x + width, origin.y), origin);
                case 1:  return MySegment(origin, MyVector2(origin.x, origin.y + height));
                case 2:  return MySegment(MyVector2(origin.x, origin.y + height), MyVector2(origin.x + width, origin.y + height));
                case 3:  return MySegment(MyVector2(origin.x + width, origin.y + height), MyVector2(origin.x + width, origin.y));
                default: break;
            }
            return MySegment(Vector2Zero(), Vector2Zero());
        }

        // Returns the number of vertices of the rectangle.
        int getVerticesNum() { return 4; }

        // Returns the vertex of the rectangle that corresponds to the given index.
        MyVector2 getVertex(int index)
        {
            assert (0 <= index && index < 4);

            switch (index)
            {
                case 0:  return MyVector2(origin.x + width, origin.y);
                case 1:  return origin;
                case 2:  return MyVector2(origin.x, origin.y + height);
                case 3:  return MyVector2(origin.x + width, origin.y + height);
                default: break;
            }

            return Vector2Zero();
        }

        // Draws the rectangle in a raylib window.
        void draw(Color color) { DrawRectangleLines(origin.x, origin.y, width, height, color); }

        // Converts the rectangle to a raylib rectangle.
        Rectangle toRayRec() { return (Rectangle){origin.x, origin.y, width, height}; }
    };

    // Polygon structure that holds values for the origin point, the radius and the number of sides.
    class MyPolygon {
        private:
            // Self reference.
            MyPolygon *self = this;
        
        public:
            // Attributes.
            MyVector2 origin;
            double    radius, rotation;
            int       sides;

            // Builder.
            MyPolygon(MyVector2 new_origin, double new_radius, double new_rotation, int new_sides) : origin(new_origin), radius(new_radius), rotation(new_rotation), sides(new_sides){};

        // ---------- POLYGON METHODS ---------- //

        // Returns the center of mass of the polygon.
        MyVector2 getCenterOfMass() { return origin; }

        // Returns the number of sides of the polygon.
        int getSidesNum() { return sides; }

        // Returns the side of the polygon that corresponds to the given index.
        MySegment getSide(int index)
        {
            assert(0 <= index && index < sides);

            double corner_angle    = arithmetic::degToRad(360 / sides);
            double angle_offset    = PI / 2 + (index * corner_angle);
            MyVector2 poly_point_a = origin + Vector2FromAngle(angle_offset + rotation, radius);
            MyVector2 poly_point_b = origin + Vector2FromAngle(angle_offset + corner_angle + rotation, radius);

            return MySegment(poly_point_a, poly_point_b);
        }

        // Returns the number of vertices of the polygon.
        int getVerticesNum() { return sides; }

        // Returns the vertex of the polygon that corresponds to the given index.
        MyVector2 getVertex(int index)
        {
            assert (0 <= index && index < sides);
            return self->getSide(index).a;
        }

        // Draws the polygon in a raylib window.
        void draw(Color color) { DrawPolyLines(origin.toRayVec(), sides, radius, arithmetic::radToDeg(rotation), color); }
    };

    // Circle structure that holds values for the origin point and radius.
    class MyCircle {
        private:
            // Self reference.
            MyCircle *self = this;
        
        public :
            // Attributes.
            MyVector2 origin;
            double    radius;

            // Builder.
            MyCircle(MyVector2 new_origin, double new_radius) : origin(new_origin), radius(new_radius){};

        // ---------- CIRCLE METHODS ---------- //

        // Returns the center of mass of the circle.
        MyVector2 getCenterOfMass() { return origin; }

        // Returns the number of sides of the circle.
        int getSidesNum() { return 1; }

        // Returns the number of vertices of the circle.
        int getVerticesNum() { return 0; }

        // Draws the circle in a raylib window.
        void draw(Color color) { DrawCircleLines(origin.x, origin.y, radius, color); }
    };

    // ----------- SHAPES UNION RELATED FUNCTIONS --------- //

    union Shape {
        MyVector2   vector;
        MySegment   segment;
        MyTriangle  triangle;
        MyRectangle rectangle;
        MyPolygon   polygon;
        MyCircle    circle;
    };

    // Shape types enum.
    enum class ShapeTypes {VECTOR2, SEGMENT, TRIANGLE, RECTANGLE, POLYGON, CIRCLE};

    // Structure for shape info, holds shape type and data.
    class ShapeInfo {
        public :
            ShapeTypes type;
            Shape data;

        // Returns the center of mass of the given shape.
        MyVector2 getCenterOfMass()
        {
            assert (type != ShapeTypes::VECTOR2);

            switch (type)
            {
                case ShapeTypes::SEGMENT:   return data.segment.getCenterOfMass();
                case ShapeTypes::TRIANGLE:  return data.triangle.getCenterOfMass();
                case ShapeTypes::RECTANGLE: return data.rectangle.getCenterOfMass();
                case ShapeTypes::POLYGON:   return data.polygon.getCenterOfMass();
                case ShapeTypes::CIRCLE:    return data.circle.getCenterOfMass();
                default:                    break;
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
                case ShapeTypes::SEGMENT:   return 1;
                case ShapeTypes::TRIANGLE:  return 3;
                case ShapeTypes::RECTANGLE: return 4;
                case ShapeTypes::POLYGON:   return data.polygon.sides;
                case ShapeTypes::CIRCLE:    return 1;
                default:                    return 0;
            }
        }

        // Returns the side of the given shape that corresponds to the given index.
        // Returns a (0, 0) segment if the shape type is not supported (circle and vector).
        MySegment getSide(int index)
        {
            switch (type)
            {
                case ShapeTypes::SEGMENT:   assert (index < 1); return data.segment;
                case ShapeTypes::TRIANGLE:  return data.triangle.getSide(index);
                case ShapeTypes::RECTANGLE: return data.rectangle.getSide(index);
                case ShapeTypes::POLYGON:   return data.polygon.getSide(index);
                default:                    return MySegment(Vector2Zero(), Vector2Zero());
            }
        }

        // Returns the number of vertices of a given shape.
        int getVerticesNum()
        {
            switch (type)
            {
                case ShapeTypes::SEGMENT:   return 2;
                case ShapeTypes::TRIANGLE:  return 3;
                case ShapeTypes::RECTANGLE: return 4;
                case ShapeTypes::POLYGON:   return data.polygon.sides;
                default:                    return 0;
            }
        }

        // Returns the vertex of the given shape that corresponds to the given index.
        MyVector2 ShapeGetVertex(int index)
        {
            switch (type)
            {
                case ShapeTypes::SEGMENT:   return data.segment.getVertex(index);
                case ShapeTypes::TRIANGLE:  return data.triangle.getVertex(index);
                case ShapeTypes::RECTANGLE: return data.rectangle.getVertex(index);
                case ShapeTypes::POLYGON:   return data.polygon.getVertex(index);
                default:                    return MyVector2(1000000, 1000000);
            }
        }

        // Draws any shape in a raylib window.
        void draw(Color color, MyVector2 origin = Vector2Zero())
        {
            switch (type)
            {
                case ShapeTypes::VECTOR2:   data.vector.draw(origin, color); break;
                case ShapeTypes::SEGMENT:   data.segment.draw(color); break;
                case ShapeTypes::TRIANGLE:  data.triangle.draw(color); break;
                case ShapeTypes::RECTANGLE: data.rectangle.draw(color); break;
                case ShapeTypes::POLYGON:   data.polygon.draw(color); break;
                case ShapeTypes::CIRCLE:    data.circle.draw(color); break;
                default:                    break;
            }
        }
    };

    // ----------------- MISC. FUNCTIONS ------------------ //

    // Returns the distance between two points.
    double distancePoints(MyVector2 p1, MyVector2 p2) { return Vector2FromPoints(p1, p2).getLength(); }

    // ------------------- COLLISIONS -------------------- //

    // Returns the smallest rectangle that contanins the given shape.
    template <typename T> MyRectangle getBoundingBox(T shape)
    {
        if (isShape(shape, true))
        {
            // If the shape is a circle.
            if (typeid(T) == typeid(MyCircle))
            {
                MyRectangle bounding_box = MyRectangle(shape.data.circle.origin - shape.data.circle.radius, 
                                                    shape.data.circle.radius * 2, 
                                                    shape.data.circle.radius * 2);

                //! Debug render.
                if (__debug_bounding_boxes) bounding_box.draw(GRAY);

                return bounding_box;
            }

            // Get the shape's vertices information.
            int vertices_num = shape.getVerticesNum();

            // Create the min and max values for x and y.
            MyVector2 vertex = shape.getVertex(0);
            double xmin      = vertex.x;
            double xmax      = vertex.x;
            double ymin      = vertex.y;
            double ymax      = vertex.y;

            // Loop though the vertices and find the min and max values for x and y.
            for (int i = 1; i < vertices_num; i++)
            {
                vertex = shape.getVertex(i);
                if (vertex.x <= xmin) xmin = vertex.x;
                if (vertex.x >= xmax) xmax = vertex.x;
                if (vertex.y <= ymin) ymin = vertex.y;
                if (vertex.y >= ymax) ymax = vertex.y;
            }

            // Create the shape's bounding box.
            MyRectangle bounding_box = MyRectangle(MyVector2(xmin, ymin), xmax - xmin, ymax - ymin);

            //! Debug render.
            if (__debug_bounding_boxes) bounding_box.draw(GRAY);

            return bounding_box;
        }
    }

    // Returns an axis that passes through the center of the given circle and the center of the given shape.
    template <typename T>
    MySegment CircleGetAxis(MyCircle circle, T shape)
    {
        if (isShape(shape, true))
        {
            // Make a segment that starts at the center of the circle, goes in the direction of the center of the shape and is of length 1.
            return SegmentFromVector2(circle.origin, Vector2FromPoints(circle.origin, shape.getCenterOfMass()).normalize());
        }
    }

    // Returns the axis of the given shapes that corresponds to the given index.
    template <typename T1, typename T2>
    MySegment ShapesGetAxis(T1 shape1, T2 shape2, int index)
    {
        if (isShape(shape1, true) && isShape(shape2, true))
        {
            assert (0 <= index && index < shape1.getSidesNum() + shape2.getSidesNum());

            MySegment side;
            MySegment axis;

            // If the given index refers to an axis of the first shape...
            if (index < shape1.getSidesNum())
            {
                // If the first shape is not a circle, get the side pointed to by the index and calculate its normal.
                if (T1 != MyCircle) {
                    side = shape1.ShapeGetSide(index);
                    axis = SegmentFromVector2(side.getCenterOfMass()),
                                            Vector2FromSegment(side).normalize().getNormal());
                }
                // If the first shape is a circle, get its axis.
                else
                    axis = CircleGetAxis(shape1, shape2);
            }
            // If the given index refers to an axis of the second shape...
            else
            {
                // If the second shape is not a circle, get the side pointed to by the index and calculate its normal.
                if (T2 != MyCircle) {
                    side = shape2.ShapeGetSide(index - getSidesNum(shape1);
                    axis = SegmentFromVector2(side.getCenterOfMass()),
                                                Vector2Normal(Vector2Normalize(Vector2FromSegment(side)));
                }
                // If the second shape is a circle, get its axis.
                else
                    axis = CircleGetAxis(shape2.data.circle, shape1);
            }

            //! Debug render.
            if (__debug_axes) DrawMyVector2(Vector2MultiplyVal(Vector2FromSegment(axis), 100), axis.a, BLUE);

            return axis;
        }
    }

    // Returns true if the given point is colliding with the given circle.
    bool collisionCirclePoint(MyCircle c, MyVector2 p) { return (distancePoints(c.origin, p) <= c.radius ? true : false); }

    // Returns true if the given circles are in collision.
    bool collisionCircles(MyCircle c1, MyCircle c2)    { return (distancePoints(c1.origin, c2.origin) <= c1.radius + c2.radius ? true : false); }

    // Checks for collision between two rectangles.
    bool collisionAABB(MyRectangle rec1, MyRectangle rec2)
    {
        return (rec1.origin.x + rec1.width  >= rec2.origin.x              &&
                rec1.origin.x               <= rec2.origin.x + rec2.width &&
                rec1.origin.y + rec1.height >= rec2.origin.y              &&
                rec1.origin.y               <= rec2.origin.y + rec2.height);
    }

    // Project a shape onto a given axis.
    template <typename T>
    MySegment projectShapeOnAxis(MySegment axis, T shape)
    {
        if (isShape(shape, true))
        {
            // Get the axis' vector.
            MyVector2 axis_vec = Vector2FromSegment(axis);

            // Handle circles.
            if (typeid(T) == typeid(MyCircle))
            {
                // Project the circle's origin onto the axis.
                MyVector2 origin_projection = axis.a + axis_vec * (Vector2FromPoints(axis.a, shape.origin) & axis_vec);

                // Create a segment of the circle's projection.
                MySegment circle_projection = MySegment(origin_projection - axis_vec * shape.radius,
                                                        origin_projection + axis_vec * shape.radius);

                //! Debug render.
                if (__debug_points)
                {
                    shape.origin.drawAsPoint(WHITE);
                    circle_projection.a.drawAsPoint(SKYBLUE);
                    circle_projection.b.drawAsPoint(BLUE);
                }
                if (__debug_projections) circle_projection.draw(ORANGE);
                
                return circle_projection;
            }

            //* https://fr.wikipedia.org/wiki/Projection_orthogonale#Projet%C3%A9_orthogonal_sur_une_droite,_distance

            // Get all the vertices of the shape.
            int vertices_num = shape.getVerticesNum();
            MyVector2 vertex;
            MyVector2 projected_points[vertices_num];

            // Loop over the vertices of the shape and get their projections onto the axis.
            for (int i = 0; i < vertices_num; i++)
            {
                vertex = shape.getVertex(i);
                projected_points[i] = axis.a + axis_vec * (Vector2FromPoints(axis.a, vertex) & axis_vec);

                //! Debug render.
                if (__debug_points) projected_points[i].drawAsPoint(WHITE);
            }

            // Find the closest and farthest points from the axis origin.
            MyVector2 min_point = projected_points[0];
            MyVector2 max_point = min_point;

            for (int i = 0; i < vertices_num; i++)
            {
                if (projected_points[i].getCopiedSign(axis_vec).x > max_point.getCopiedSign(axis_vec).x ||
                    projected_points[i].getCopiedSign(axis_vec).y > max_point.getCopiedSign(axis_vec).y)
                {
                    max_point = projected_points[i];
                }

                if (projected_points[i].getCopiedSign(axis_vec).x < min_point.getCopiedSign(axis_vec).x ||
                    projected_points[i].getCopiedSign(axis_vec).y < min_point.getCopiedSign(axis_vec).y)
                {
                    min_point = projected_points[i];
                }
            }

            MyVector2 axis_orig_to_min_point = Vector2FromPoints(axis.a, min_point);
            MySegment projection = SegmentFromVector2(axis.a + axis_orig_to_min_point, Vector2FromPoints(min_point, max_point));

            //! Debug render.
            if (__debug_points)
            {
                min_point.drawAsPoint(SKYBLUE);
                max_point.drawAsPoint(BLUE);
            }
            if (__debug_projections) projection.draw(ORANGE);

            return projection;
        }
    }

    // Returns true if the given point is colliding with the given segment.
    bool collisionSegmentPoint(MySegment segment, MyVector2 point)
    {
        if (arithmetic::roundInt(Vector2FromSegment(segment) ^ Vector2FromPoints(segment.a, point)) == 0)
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
    bool collisionProjections(MySegment projection1, MySegment projection2)
    {
        return (collisionSegmentPoint(projection1, projection2.a) ||
                collisionSegmentPoint(projection1, projection2.b) ||
                collisionSegmentPoint(projection2, projection1.a) ||
                collisionSegmentPoint(projection2, projection1.b));
    }

    // Checks for collision between two given shapes.
    template<typename T1, typename T2>
    bool collisionSAT(T1 shape1, T2 shape2)
    {
        if (isShape(shape1, true) && isShape(shape2, true))
        {
            // If both shapes are circles, don't use SAT.
            if (typeid(T1) == typeid(MyCircle) && 
                typeid(T2) == typeid(MyCircle))
                return collisionCircles(shape1, shape2);

            // If both shapes are rectangles, don't use SAT.
            else if (typeid(T1) == typeid(MyRectangle) && 
                    typeid(T2) == typeid(MyRectangle))
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
                        if (__debug_failed_projections)
                        {
                            projection1.draw(PINK); 
                            projection2.draw(PINK);
                        }
                        return false;
                    }
                }
                return true;
            }

            //! Debug render.
            if (__debug_shapes)
            {
                shape1.draw(GREEN); 
                shape2.draw(GREEN);
            }
        }

        return false;
    }
}

}