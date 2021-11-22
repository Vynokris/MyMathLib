/*************************************************************
 * MyMathLib, a custom mathematics library written by
 * Rémi Serra and Alexandre Perché, students at ISART Digital.
 *************************************************************/


// ------------------- INCLUDES / NAMESPACE -------------------- //

#pragma once
#include "raylib.h"
#include <cmath>
#include <cassert>
#include <stdint.h>
#include <typeinfo>
#include <cstdarg>
#include <vector>
#include <iostream> //! DEBUG

using namespace std;

namespace MyMathLib
{
    // --------------------- ARITHMECTIC ---------------------- //
    namespace arithmetic
    {
        // --------- ARITHMETIC FUNCTIONS ---------- //
        
        // Fast inverse square root from Quake III.
        float Q_rsqrt(float number);

        // Rounds the given value to the nearest int.
        int roundInt(double val);

        // Rounds down the given value.
        int floorInt(double val);

        // Rounds up the given value.
        int ceilInt(double val);

        // Returns the sqare power of the given value.
        double sqpow(double val);

        // Returns 1 if the given value is positive or null, and -1 if it is negative.
        int signOf(double val);

        // Converts the given angle from degrees to radians.
        double degToRad(double deg);

        // Converts the given angle from radians to degrees.
        double radToDeg(double rad);

        // Clamps the given value to be superior or equal to the minimum value and inferior or equal to the maximum value.
        double clamp(double val, double min, double max);

        // Clamps the given value to be inferior or equal to the maximum value.
        double clampUnder(double val, double max);

        // Clamps the given value to be superior or equal to the minimum value.
        double clampAbove(double val, double min);
        // Remaps the given value from one range to another.
        double remap(double val, double inputStart, double inputEnd, double outputStart, double outputEnd);

        // Returns true if the given number is a power of 2.
        bool is_power_of_two(int val);

        // Returns the closest power of 2 that is inferior or equal to val.
        int get_power_of_two_under(int val);

        // Returns the closest power of 2 that is superior or equal to val.
        int get_power_of_two_above(int val);

        // Return a random integer according to a seed and a state
        int32_t getNextRandomNumber(int32_t state);

        // Return a random value according to a range
        int getRandomValue(int start, int end);

        // ---------------- MATRIX ---------------- //

        // Matrix class 
        // (Web docs: https://www.bestprog.net/en/2019/08/23/c-an-example-of-creating-a-template-class-matrix-dynamic-memory-allocation/)
        template<typename T, int R, int C>
        class MyMatrix
        {
            private:
                T   matrix[R][C];

            public:
                // ----- Constructors & Destructor ----- //
                MyMatrix() { assert(R >= 2 && C >= 2); }
                MyMatrix(const MyMatrix<T,R,C>& _matrix) { assert(R >= 2 && C >= 2); for (int i = 0; i < R; i++) for (int j = 0; j < C; j++) matrix[i][j] = _matrix[i][j]; }
                ~MyMatrix() {}

                // ----- Methods ----- //

                // Getters.
                int getRows()                    { return R;            }
                int getColumns()                 { return C;            }  
                T   getMatrixValue(int i, int j) { return matrix[i][j]; }

                // Setters.
                void setMatrixValues(std::initializer_list<T> inputs)
                {
                    vector<T> list = inputs;
                    assert(list.size() == R * C);
                    int tmp = 0;

                    for (int i = 0; i < R; i++)
                        for (int j = 0; j < C; j++)
                        {
                            matrix[i][j] = list[tmp];
                            tmp++;
                        }
                }

                // Arithmetic.
                bool isSquare() { return R == C; }

                bool isIdentity()
                {
                    for (int i = 0; i < R; i++)
                        for (int j = 0; j < C; j++)
                            if ((i != j && matrix[i][j] != 0) || (i == j && matrix[i][j] != 1))
                                return false;
                    return true;
                }

                // ----- Operators ----- //

                // Matrix bracket operators.
                      T* operator[](int index)       { return matrix[index]; }
                const T* operator[](int index) const { return matrix[index]; }

                // Matrix copy.
                MyMatrix<T,R,C> operator=(const MyMatrix<T,R,C>& _matrix) const
                {
                    // Matrix content copy
                    for (int i = 0; i < R; i++) 
                        for (int j = 0; j < C; j++)
                            matrix[i][j] = _matrix[i][j];
                    
                    return *this;
                }

                // Matrix addition.
                MyMatrix<T,R,C> operator+(const T& val) const
                {
                    MyMatrix<T,R,C> tmp;
                    for (int i = 0; i < R; i++)
                        for (int j = 0; j < C; j++)
                            tmp[i][j] = matrix[i][j] + val;
                    return tmp;
                }

                template<int _R, int _C>
                MyMatrix<T,R,C> operator+(const MyMatrix<T,_R,_C>& _matrix) const
                {
                    assert(R == _R && C == _C); // Matrix must have the same dimension
                    MyMatrix<T,_R,_C> tmp;
                    for(int i = 0; i < R; i++)
                        for(int j = 0; j < C; j++)
                            tmp[i][j] = matrix[i][j] + _matrix[i][j];
                    return tmp;
                }

                // Matrix substraction.
                MyMatrix<T,R,C> operator-(const T& val) const
                {
                    MyMatrix<T,R,C> tmp;
                    for (int i = 0; i < R; i++)
                        for (int j = 0; j < C; j++)
                            tmp[i][j] = matrix[i][j] - val;
                    return tmp;
                }

                template<int _R, int _C>
                MyMatrix<T,R,C> operator-(const MyMatrix<T,_R,_C>& _matrix) const
                {
                    assert(R == _R && C == _C);
                    MyMatrix<T,_R,_C> tmp;
                    for (int i = 0; i < R; i++)
                        for (int j = 0; j < C; j++)
                            tmp[i][j] = matrix[i][j] - _matrix[i][j];
                    return tmp;
                }


                // Matrix multiplication.
                MyMatrix<T,R,C> operator*(const T& val) const
                {
                    MyMatrix<T,R,C> tmp;
                    for (int i = 0; i < R; i++)
                        for (int j = 0; j < C; j++)
                            tmp[i][j] = matrix[i][j] * val;
                    return tmp;
                }

                template<int _R, int _C>
                MyMatrix<T, (R > _R ? R : _R),(C > _C ? C : _C)> operator*(const MyMatrix<T,_R,_C>& _matrix) const
                {
                    assert(C == _R); // Size condition to calculate
    
                    MyMatrix<T, (R > _R ? R : _R),(C > _C ? C : _C)> result;
                    for (int i = 0; i < R; i++)
                        for (int j = 0; j < _C; j++)
                        {
                            result[i][j] = 0;
                            for (int k = 0; k < C; k++) {
                                result[i][j] = matrix[i][k] + _matrix[k][j];
                            }
                        }
                    return result;
                }

                // Matrix division by a scalar.
                MyMatrix<T,R,C> operator/(const T& val) const
                {
                    MyMatrix<T,R,C> tmp;
                    for (int i = 0; i < R; i++)
                        for (int j = 0; j < C; j++)
                            tmp[i][j] = matrix[i][j] / val;
                    return matrix;
                }

                // Matrix addition assignement.
                void operator+=(const T &val)
                {
                    for (int i = 0; i < R; i++)
                        for (int j = 0; j < C; j++)
                            matrix[i][j] += val;
                }

                template<int _R, int _C>
                void operator+=(const MyMatrix<T,_R,_C>& _matrix)
                {
                    assert(R == _R && C == _C);
                    for (int i = 0; i < R; i++)
                        for (int j = 0; j < C; j++)
                            matrix[i][j] += _matrix[i][j];
                }

                // Matrix substraction assignement.
                void operator-=(const T &val)
                {
                    for (int i = 0; i < R; i++)
                        for (int j = 0; j < C; j++)
                            matrix[i][j] -= val;
                }

                template<int _R, int _C>
                void operator-=(const MyMatrix<T,_R,_C>& _matrix)
                {
                    assert(R == _R && C == _C);
                    for (int i = 0; i < R; i++)
                        for (int j = 0; j < C; j++)
                            matrix[i][j] -= _matrix[i][j];
                }

                // Matrix multiplication assignement.
                void operator*=(const T &val)
                {
                    for (int i = 0; i < R; i++)
                        for (int j = 0; j < C; j++)
                            matrix[i][j] *= val;
                }

                template<int _R, int _C>
                void operator*=(const MyMatrix<T,_R,_C>& _matrix)
                {
                    assert(C == _R); // Size condition to calculate

                    MyMatrix<T, (R > _R ? R : _R),(C > _C ? C : _C)> result;
                    for (int i = 0; i < R; i++)
                        for (int j = 0; j < C; j++)
                        {
                            result[i][j] = 0;
                            for (int k = 0; k < C; k++) result[i][j] = matrix[i][k] + _matrix[k][j];
                        }
                }

                // ----- Others Methods ----- //
                void print()
                {
                    // Print data
                    std::cout << "MyMatrix (" << &matrix << ") | R: " << R << " C: " << C << std::endl;
                    
                    // Print content
                    for (int i = 0; i < R; i++)
                    {
                        for (int j = 0; j < C; j++) std::cout << (j == 0 ? "| " : "") << matrix[i][j] << " | ";
                        std::cout << std::endl;
                    }
                }

        };
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

        // Makes sure that the given object is a shape.
        template <typename T>
        bool isShape(T object, bool notVector2 = false)
        {
            return ((notVector2 ? false : typeid(T) == typeid(MyVector2))         ||
                                          typeid(T) == typeid(MySegment)          ||
                                          typeid(T) == typeid(MyTriangle)         ||
                                          typeid(T) == typeid(MyRectangle)        ||
                                          typeid(T) == typeid(MyPolygon)          ||
                                          typeid(T) == typeid(MyCircle)           ||
                    (notVector2 ? false : typeid(T) == typeid(const MyVector2&))  ||
                                          typeid(T) == typeid(const MySegment&)   ||
                                          typeid(T) == typeid(const MyTriangle&)  ||
                                          typeid(T) == typeid(const MyRectangle&) ||
                                          typeid(T) == typeid(const MyPolygon&)   ||
                                          typeid(T) == typeid(const MyCircle&));
        }

        // Makes sure that the given object is a number.
        template <typename T>
        bool isNumeral(T object)
        {
            return (typeid(T) == typeid(int)         ||
                    typeid(T) == typeid(float)       ||
                    typeid(T) == typeid(double)      ||
                    typeid(T) == typeid(const int)   ||
                    typeid(T) == typeid(const float) ||
                    typeid(T) == typeid(const double));
        }

        enum class ShapeTypes {
            SEGMENT,
            TRIANGLE,
            RECTANGLE,
            POLYGON,
            CIRCLE,
        };

        // ---------------------- CLASSES --------------------- //

        // Vector class that holds values for x and y (2 dimensions).
        class MyVector2
        {
            public :
                // Attributes.
                double x, y;

                // Constructors.
                MyVector2();                                                          // Null vector.
                MyVector2(const double new_x, const double new_y);                    // Vector with 2 coordinates.
                MyVector2(const MyVector2& p1, const MyVector2& p2);                  // Vector from 2 points.
                MyVector2(const MySegment& seg);                                      // Vector from semgent.
                MyVector2(const double rad, const double length, const bool isAngle); // Vector from angle (useless bool).

                // Destroyer.
                ~MyVector2();

                // ---------- VECTOR2 OPERATOR ---------- //

                // Copy constructor.
                void operator=(const MyVector2& other);

                // Vector2 addition.
                template <typename T> MyVector2 operator+(const T& val) const;

                // Vector2 substraction.
                template <typename T> MyVector2 operator-(const T& val) const;

                // Vector2 multiplication.
                template <typename T> MyVector2 operator*(const T& val) const;

                // Vector2 division.
                template <typename T> MyVector2 operator/(const T &val) const;

                // Vector2 addition assignement.
                template <typename T> void operator+=(const T &val);

                // Vector2 substraction assignement.
                template <typename T> void operator-=(const T &val);

                // Vector2 multiplication assignement.
                template <typename T> void operator*=(const T &val);

                // Vector2 division assignement.
                template <typename T> void operator/=(const T &val);

                // Vector2 dot product.
                double operator&(const MyVector2& val) const;

                // Vector2 cross product.
                double operator^(const MyVector2& val) const;

                // ---------- VECTOR2 METHODS ---------- //
                
                // Returns the length of the given vector.
                double getLength();

                // Returns the angle (in radians) of the given vector.
                double getAngle();

                // Returns the middle of the given vector.
                MyVector2 getMiddle();

                // Normalizes the given vector so that its length is 1.
                void normalize();

                // Normalizes the given vector so that its length is 1.
                MyVector2 getNormalized();

                // Modifies the length of the given vector to correspond to the given value.
                void setLength(double length);

                // Negates both of the coordinates of the given vector.
                void negate();

                // Copies the signs from the source vector to the destination vector.
                void copysign(MyVector2 source);

                // Copies the signs from the source vector to the destination vector.
                MyVector2 getCopiedSign(MyVector2 source);

                // Returns the normal of a given vector.
                MyVector2 getNormal();

                // Interprets the vector as a point and returns the distance to another point.
                double getDistanceFromPoint(const MyVector2 &p);

                // Rotates the given vector by the given angle.
                void rotate(double angle);

                // Returns the angle (in radians) between two vectors.
                double getAngleWithVector2(MyVector2 v);

                // Draws a vector at a certain origin point of a raylib window.
                void draw(MyVector2 origin, Color color);

                // Draws a point in a raylib window.
                void drawAsPoint(Color color);

                // Converts a my_math 2D vector to a raylib 2D vector.
                Vector2 toRayVec();
        };

        // Segment structure that holds values for the starting point and the end point.
        class MySegment
        {
            public :
                // Attributes.
                MyVector2  a, b;
                ShapeTypes type = ShapeTypes::SEGMENT;

                // Constructors.
                MySegment();                                                                 // Null segment.
                MySegment(const MyVector2& new_a,  const MyVector2& new_b);                  // Segment from points.
                MySegment(const MyVector2& origin, const MyVector2& vec, const bool vector); // Segment from point and vector.

                // Destroyer.
                ~MySegment();

                // ---------- SEGMENT METHODS ---------- //

                // Returns the center of mass of the segment.
                MyVector2 getCenterOfMass();

                // Returns the number of sides of the segment.
                int getSidesNum();

                // Returns the side of the segment that corresponds to the given index.
                MySegment getSide(int index);

                // Returns the number of vertices of the segment.
                int getVerticesNum();

                // Returns the vertex of the segment that corresponds to the given index.
                MyVector2 getVertex(int index);

                // Draws the segment in a raylib window.
                void draw(Color color);
        };

        // Triangle structure that holds values for 3 points.
        class MyTriangle 
        {
            public:
                // Attributes.
                MyVector2  a, b, c;
                ShapeTypes type = ShapeTypes::TRIANGLE;

                // Constructor.
                MyTriangle();                                                  // Null triangle.
                MyTriangle(MyVector2 new_a, MyVector2 new_b, MyVector2 new_c); // Triangle from points.

                // Destroyer.
                ~MyTriangle();

                // ---------- TRIANGLE METHODS ---------- //
                
                // Returns the center of mass of the triangle.
                MyVector2 getCenterOfMass();

                // Returns the number of sides of the triangle.
                int getSidesNum();

                // Returns the side of the triangle that corresponds to the given index.
                MySegment getSide(int index);

                // Returns the number of vertices of the triangle.
                int getVerticesNum();

                // Returns the vertex of the triangle that corresponds to the given index.
                MyVector2 getVertex(int index);

                // Draws the triangle in a raylib window.
                void draw(Color color);
        };

        // Rectangle structure that holds values for the origin point, width and height.
        class MyRectangle
        {        
            public:            
                // Attributes.
                MyVector2  origin;
                double     width, height;
                ShapeTypes type = ShapeTypes::RECTANGLE;

                // Constructor.
                MyRectangle();                                                          // Null rectangle.
                MyRectangle(MyVector2 new_origin, double new_width, double new_height); // Rectangle from origin, width and height.

                // Destroyer.
                ~MyRectangle();

                // ---------- RECTANGLE METHODS ---------- //

                // Returns the center of mass of the rectangle.
                MyVector2 getCenterOfMass();

                // Returns the number of sides of the rectangle.
                int getSidesNum();

                // Returns the side of the rectangle that corresponds to the given index.
                MySegment getSide(int index);

                // Returns the number of vertices of the rectangle.
                int getVerticesNum();

                // Returns the vertex of the rectangle that corresponds to the given index.
                MyVector2 getVertex(int index);

                // Draws the rectangle in a raylib window.
                void draw(Color color);

                // Converts the rectangle to a raylib rectangle.
                Rectangle toRayRec();
        };

        // Polygon structure that holds values for the origin point, the radius and the number of sides.
        class MyPolygon
        {
        public:
                // Attributes.
                MyVector2  origin;
                double     radius, rotation;
                int        sides;
                ShapeTypes type = ShapeTypes::POLYGON;

                // Constructor.
                MyPolygon();                                                                            // Null polygon.
                MyPolygon(MyVector2 new_origin, double new_radius, double new_rotation, int new_sides); // Polygon from origin, radius, rotation and side amount.

                // Destroyer.
                ~MyPolygon();

                // ---------- POLYGON METHODS ---------- //

                // Returns the center of mass of the polygon.
                MyVector2 getCenterOfMass();

                // Returns the number of sides of the polygon.
                int getSidesNum();

                // Returns the side of the polygon that corresponds to the given index.
                MySegment getSide(int index);

                // Returns the number of vertices of the polygon.
                int getVerticesNum();

                // Returns the vertex of the polygon that corresponds to the given index.
                MyVector2 getVertex(int index);

                // Draws the polygon in a raylib window.
                void draw(Color color);
        };

        // Circle structure that holds values for the origin point and radius.
        class MyCircle 
        {
            public :
                // Attributes.
                MyVector2  origin;
                double     radius;
                ShapeTypes type = ShapeTypes::CIRCLE;


                // Constructor.
                MyCircle();                                        // Null circle.
                MyCircle(MyVector2 new_origin, double new_radius); // Circle from origin and radius.

                // Destroyer.
                ~MyCircle();

                // ---------- CIRCLE METHODS ---------- //

                // Returns the center of mass of the circle.
                MyVector2 getCenterOfMass();

                // Returns the number of sides of the circle.
                int getSidesNum();

                // Does nothing and returns a null segment.
                MySegment getSide(int index);

                // Returns the number of vertices of the circle.
                int getVerticesNum();

                // Does nothing and returns a null vector.
                MyVector2 getVertex(int index);

                // Draws the circle in a raylib window.
                void draw(Color color);
        };
    }


    // ------------------- COLLISIONS -------------------- //

    namespace collisions
    {
        using namespace geometry;
        
        // Returns the smallest rectangle that contanins the given shape.
        template <typename T> MyRectangle getBoundingBox(T shape);

        // Returns an axis that passes through the center of the given circle and the center of the given shape.
        template <typename T> MySegment CircleGetAxis(MyCircle circle, T shape);

        // Returns the axis of the given shapes that corresponds to the given index.
        template <typename T1, typename T2>
        MySegment ShapesGetAxis(T1 shape1, T2 shape2, int index);

        // Returns true if the given point is colliding with the given circle.
        bool collisionCirclePoint(MyCircle c, MyVector2 p);

        // Returns true if the given circles are in collision.
        bool collisionCircles(MyCircle c1, MyCircle c2);

        // Checks for collision between two rectangles.
        bool collisionAABB(MyRectangle rec1, MyRectangle rec2);

        // Project a shape onto a given axis.
        template <typename T> MySegment projectShapeOnAxis(MySegment axis, T shape);

        // Returns true if the given point is colliding with the given segment.
        bool collisionSegmentPoint(MySegment segment, MyVector2 point);

        // Returns true if the given projections are colliding each others
        bool collisionProjections(MySegment projection1, MySegment projection2);

        // Checks for collision between two given shapes.
        template <typename T1, typename T2> bool collisionSAT(T1 shape1, T2 shape2);
    }
}