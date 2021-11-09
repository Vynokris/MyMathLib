#include "my_math.hpp"

using namespace MyMathLib;

// ------------------------------ ARITHMECTIC ------------------------------ //

// Rounds the given value to the nearest int.
int arithmetic::roundInt(double val)    { return (int)round(val); }

// Rounds down the given value.
int arithmetic::floorInt(double val)    { return (int)floor(val); }

// Rounds up the given value.
int arithmetic::ceilInt(double val)     { return (int)ceil(val); }

// Returns the sqare power of the given value.
double arithmetic::sqpow(double val)    { return val * val; }

// Returns 1 if the given value is positive or null, and -1 if it is negative.
int arithmetic::signOf(double val)      { if (val == 0) return 1; return val / abs((int)val); }

// Converts the given angle from degrees to radians.
double arithmetic::degToRad(double deg) { return deg * (PI / 180.0f); }

// Converts the given angle from radians to degrees.
double arithmetic::radToDeg(double rad) { return rad * (180.0f / PI); }

// Clamps the given value to be superior or equal to the minimum value and inferior or equal to the maximum value.
double arithmetic::clamp(double val, double min, double max)
{
    assert (min <= max); 
    if (val < min) val = min;
    if (val > max) val = max;
    return val;
}

// Clamps the given value to be inferior or equal to the maximum value.
double arithmetic::clampUnder(double val, double max) { if (val > max) val = max; return val; }

// Clamps the given value to be superior or equal to the minimum value.
double arithmetic::clampAbove(double val, double min) { if (val < min) val = min; return val; }

// Remaps the given value from one range to another.
double arithmetic::remap(double val, double inputStart, double inputEnd, double outputStart, double outputEnd)
{
    return outputStart + (val - inputStart) * (outputEnd - outputStart) / (inputEnd - inputStart);
}

// Return a random integer according to a seed and a state
int32_t arithmetic::getNextRandomNumber(int32_t state)
{
    uint64_t result = (((int64_t)1103515245) * state + 12345) % 2147483648;
    return result;
}

// Return a random value according to a range
int arithmetic::getRandomValue(int start, int end)
{
    assert(start > end);
    int r = arithmetic::getNextRandomNumber(0);
    for (int i = 0; r >= start && r <= end; i++) r = arithmetic::getNextRandomNumber(i);
    return r;
}



// ------------------------------ GEOMETRY ------------------------------ //



// -------------------- VECTOR 2 -------------------- //

using namespace geometry;

// Constructors.
MyVector2::MyVector2()                                                          : x(0),                 y(0)                 {}; // Null vector.
MyVector2::MyVector2(const double new_x, const double new_y)                    : x(new_x),             y(new_y)             {}; // Vector with 2 coordinates.
MyVector2::MyVector2(const MyVector2& p1, const MyVector2& p2)                  : x(p2.x - p1.x),       y(p2.y - p1.y)       {}; // Vector from 2 points.
MyVector2::MyVector2(const MySegment& seg)                                      : x(seg.b.x - seg.a.x), y(seg.b.y - seg.a.y) {}; // Vector from semgent.
MyVector2::MyVector2(const double rad, const double length, const bool isAngle) : x(cos(rad) * length), y(sin(rad) * length) {}; // Vector from angle (useless bool).

// Destroyer.
MyVector2::~MyVector2() {};

// ---------- VECTOR 2 OPERATORS ---------- //

// Copy constructor.
void MyVector2::operator=(const MyVector2& other) { x = other.x; y = other.y; }

// Vector2 addition.
template <typename T>
MyVector2 MyVector2::operator+(const T& val)
{
    assert (typeid(T) == typeid(const MyVector2) || isNumeral(val));

    if (typeid(T) == typeid(const MyVector2)) return MyVector2(x + val.x, y + val.y);
    // else if (isNumeral(val))                  return MyVector2(x + val,   y + val  );

    return MyVector2(); // Will never happen thanks to the assert.
}

// Vector2 substraction.
template <typename T>
MyVector2 MyVector2::operator-(const T& val)
{
    assert (typeid(T) == typeid(const MyVector2) || isNumeral(val));

    if (typeid(T) == typeid(const MyVector2)) return MyVector2(x - val.x, y - val.y);
    else if (isNumeral(val))                  return MyVector2(x - val,   y - val  );

    return MyVector2(); // Will never happen thanks to the assert.
}

// Vector2 multiplication.
template <typename T>
MyVector2 MyVector2::operator*(const T& val)
{
    assert (typeid(T) == typeid(const MyVector2) || isNumeral(val));

    if (typeid(T) == typeid(const MyVector2)) return MyVector2(x * val.x, y * val.y);
    else if (isNumeral(val))                  return MyVector2(x * val,   y * val  );

    return MyVector2(); // Will never happen thanks to the assert.
}

// Vector2 division.
template <typename T>
MyVector2 MyVector2::operator/(const T& val)
{
    assert (typeid(T) == typeid(const MyVector2) || isNumeral(val));

    if (typeid(T) == typeid(const MyVector2)) return MyVector2(x / val.x, y / val.y);
    else if (isNumeral(val))                  return MyVector2(x / val,   y / val  );

    return MyVector2(); // Will never happen thanks to the assert.
}

// Vector2 addition assignement.
template <typename T>
void MyVector2::operator+=(const T& val)
{
    if (typeid(T) == typeid(const MyVector2)) { x += val.x; y += val.y; }
    else if (isNumeral(val))                  { x += val;   y += val;   }
}

// Vector2 substraction assignement.
template<typename T>
void MyVector2::operator-=(const T& val)
{
    if(typeid(T) == typeid(const MyVector2)) { x -= val.x; y -= val.y; }
    else if (isNumeral(val))                 { x -= val;   y -= val;   }
}

// Vector2 multiplication assignement.
template <typename T>
void MyVector2::operator*=(const T& val)
{
    if (typeid(T) == typeid(const MyVector2)) { x *= val.x; y *= val.y; }
    else if (isNumeral(val))                  { x *= val;   y *= val;   }
}

// Vector2 division assignement.
template <typename T>
void MyVector2::operator/=(const T& val)
{
    if (typeid(T) == typeid(const MyVector2)) { x /= val.x; y /= val.y; }
    else if (isNumeral(val))                  { x /= val;   y /= val;   }
}

// Vector2 dot product.
double MyVector2::operator&(MyVector2 val)         { return (x * val.x) + (y * val.y); }

// Vector2 cross product.
double MyVector2::operator^(MyVector2 val)         { return (x * val.y) - (y * val.x); }

// ---------- VECTOR2 METHODS ---------- //

// Returns the length of the given vector.
double MyVector2::getLength()                              { return sqrt(arithmetic::sqpow(x) + arithmetic::sqpow(y)); }

// Returns the angle (in radians) of the given vector.
double MyVector2::getAngle()                               { return std::copysign(std::acos(this->getNormalized().x), std::asin(this->getNormalized().y)); }

// Returns the middle of the given vector.
MyVector2 MyVector2::getMiddle()                           { return MyVector2(x / 2, y / 2); }

// Normalizes the given vector so that its length is 1.
void MyVector2::normalize()                                { x /= this->getLength(); y /= this->getLength(); }

// Normalizes the given vector so that its length is 1.
MyVector2 MyVector2::getNormalized()                       { return MyVector2(x / this->getLength(), y / this->getLength()); }

// Modifies the length of the given vector to correspond to the given value.
void MyVector2::setLength(double length)                   { *(this) = MyVector2(this->getAngle(), length, true); }

// Negates both of the coordinates of the given vector.
void MyVector2::negate()                                   { *(this) = MyVector2(-x, -y); }

// Copies the signs from the source vector to the destination vector.
void MyVector2::copysign(MyVector2 source)                 { *(this) = MyVector2(std::copysign(x, source.x), std::copysign(y, source.y)); }

// Copies the signs from the source vector to the destination vector.
MyVector2 MyVector2::getCopiedSign(MyVector2 source)       { return MyVector2(std::copysign(x, source.x), std::copysign(y, source.y)); }

// Returns the normal of a given vector.
MyVector2 MyVector2::getNormal()                           { return MyVector2(-y, x); }

// Interprets the vector as a point and returns the distance to another point.
double MyVector2::getDistanceFromPoint(const MyVector2& p) { return MyVector2(*this, p).getLength(); }

// Rotates the given vector by the given angle.
void MyVector2::rotate(double angle)
{
    double this_length = this->getLength();
    double this_angle  = this->getAngle();
    *(this) = MyVector2(this_angle + angle, this_length, true);
}

// Returns the angle (in radians) between two vectors.
double MyVector2::getAngleWithVector2(MyVector2 v)
{
    double this_angle = this->getAngle();
    double v_angle    = v.getAngle();
    return (this_angle >= v_angle ? (this_angle - v_angle) : (v_angle - this_angle));
}

// Draws a vector at a certain origin point of a raylib window.
void MyVector2::draw(MyVector2 origin, Color color)
{
    DrawLine(origin.x, origin.y, origin.x + x, origin.y + y, color);
    MyVector2 tmp = origin + *this;
    DrawPoly(tmp.toRayVec(), 3, 4, arithmetic::radToDeg(this->getAngle() - PI / 2), color);
}

// Draws a point in a raylib window.
void MyVector2::drawAsPoint(Color color) { DrawCircle(x, y, 2, color); }

// Converts a my_math 2D vector to a raylib 2D vector.
Vector2 MyVector2::toRayVec() { return (Vector2){ (float)x, (float)y }; }



// -------------------- SEGMENT -------------------- //

// Constructors.
MySegment::MySegment()                                                     : a(MyVector2()), b(MyVector2())   {}; // Null segment.
MySegment::MySegment(const MyVector2& new_a, const MyVector2& new_b)       : a(new_a),       b(new_b)         {}; // Segment from points.
MySegment::MySegment(MyVector2& origin, MyVector2& vec, const bool vector) : a(origin),      b(origin + vec)  {}; // Segment from point and vector.

// Destroyer.
MySegment::~MySegment() {};

// ---------- SEGMENT METHODS ---------- //

// Returns the center of mass of the segment.
MyVector2 MySegment::getCenterOfMass()
{
    return MyVector2((a.x + b.x) / 2, (a.y +b.y) / 2);
}

// Returns the number of sides of the segment.
int MySegment::getSidesNum() { return 1; }

// Returns the side of the segment that corresponds to the given index.
MySegment MySegment::getSide(int index)
{
    assert (0 <= index && index < 1);
    return *this;
}

// Returns the number of vertices of the segment.
int MySegment::getVerticesNum() { return 2; }

// Returns the vertex of the segment that corresponds to the given index.
MyVector2 MySegment::getVertex(int index)
{
    assert (0 <= index && index < 2);

    switch (index)
    {
        case 0:  return a;
        case 1:  return b;
        default: break;
    }
    return MyVector2();
}

// Draws the segment in a raylib window.
void MySegment::draw(Color color)
{
    DrawLine(a.x, a.y, b.x, b.y, color);
    DrawCircle(a.x, a.y, 2, color);
    DrawCircle(b.x, b.y, 2, color);
}



// -------------------- TRIANGLE -------------------- //

// Constructor.
MyTriangle::MyTriangle()                                                  : a(MyVector2()), b(MyVector2()), c(MyVector2()) {}; // Null triangle.
MyTriangle::MyTriangle(MyVector2 new_a, MyVector2 new_b, MyVector2 new_c) : a(new_a),       b(new_b),       c(new_c)       {}; // Triangle from points.

// Destroyer.
MyTriangle::~MyTriangle() {};

// ---------- TRIANGLE METHODS ---------- //

// Returns the center of mass of the triangle.
MyVector2 MyTriangle::getCenterOfMass() { return MyVector2((a.x + b.x + c.x) / 3, (a.y + b.y + c.y) / 3); }

// Returns the number of sides of the triangle.
int MyTriangle::getSidesNum() { return 3; }

// Returns the side of the triangle that corresponds to the given index.
MySegment MyTriangle::getSide(int index) 
{
    assert (0<= index && index < 3);

    switch (index)
    {
        case 0:   return MySegment(a, b);
        case 1:   return MySegment(b, c);
        case 2:   return MySegment(c, a);
        default:  break;
    }

    return MySegment(MyVector2(), MyVector2());
}

// Returns the number of vertices of the triangle.
int MyTriangle::getVerticesNum() { return 3; }

// Returns the vertex of the triangle that corresponds to the given index.
MyVector2 MyTriangle::getVertex(int index)
{
    assert (0 <= index && index < 3);

    switch (index)
    {
    case 0:  return a;
    case 1:  return b;
    case 2:  return c;
    default: break;
    }

    return MyVector2();
}

// Draws the triangle in a raylib window.
void MyTriangle::draw(Color color) { DrawTriangleLines(a.toRayVec(), b.toRayVec(), c.toRayVec(), color); }



// -------------------- RECTANGLE -------------------- //

// Constructor.
MyRectangle::MyRectangle()                                                          : origin(MyVector2()), width(0),         height(0)          {}; // Null rectangle.
MyRectangle::MyRectangle(MyVector2 new_origin, double new_width, double new_height) : origin(new_origin),  width(new_width), height(new_height) {}; // Rectangle from origin, width and height.

// Destroyer.
MyRectangle::~MyRectangle() {};

// ---------- RECTANGLE METHODS ---------- //

// Returns the center of mass of the rectangle.
MyVector2 MyRectangle::getCenterOfMass() { return MyVector2(origin.x + width / 2, origin.y + height / 2); };

// Returns the number of sides of the rectangle.
int MyRectangle::getSidesNum() { return 4; }

// Returns the side of the rectangle that corresponds to the given index.
MySegment MyRectangle::getSide(int index)
{
    assert (0 <= index && index < 4);

    switch (index)
    {
        case 0:  return MySegment(MyVector2(origin.x + width, origin.y)         , origin);
        case 1:  return MySegment(origin                                        , MyVector2(origin.x, origin.y + height));
        case 2:  return MySegment(MyVector2(origin.x, origin.y + height)        , MyVector2(origin.x + width, origin.y + height));
        case 3:  return MySegment(MyVector2(origin.x + width, origin.y + height), MyVector2(origin.x + width, origin.y));
        default: break;
    }
    return MySegment(MyVector2(), MyVector2());
}

// Returns the number of vertices of the rectangle.
int MyRectangle::getVerticesNum() { return 4; }

// Returns the vertex of the rectangle that corresponds to the given index.
MyVector2 MyRectangle::getVertex(int index)
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

    return MyVector2();
}

// Draws the rectangle in a raylib window.
void MyRectangle::draw(Color color) { DrawRectangleLines(origin.x, origin.y, width, height, color); }

// Converts the rectangle to a raylib rectangle.
Rectangle MyRectangle::toRayRec() { return (Rectangle){ (float)origin.x, (float)origin.y, (float)width, (float)height }; }



// -------------------- POLYGON -------------------- //

// Constructor.
MyPolygon::MyPolygon()                                                                            : origin(MyVector2()), radius(0),          rotation(0),            sides(3)        {}; // Null polygon.
MyPolygon::MyPolygon(MyVector2 new_origin, double new_radius, double new_rotation, int new_sides) : origin(new_origin),  radius(new_radius), rotation(new_rotation), sides(new_sides){}; // Polygon from origin, radius, rotation and side amount.

// Destroyer.
MyPolygon::~MyPolygon() {};

// ---------- POLYGON METHODS ---------- //

// Returns the center of mass of the polygon.
MyVector2 MyPolygon::getCenterOfMass() { return origin; }

// Returns the number of sides of the polygon.
int MyPolygon::getSidesNum() { return sides; }

// Returns the side of the polygon that corresponds to the given index.
MySegment MyPolygon::getSide(int index)
{
    assert(0 <= index && index < sides);

    double corner_angle    = arithmetic::degToRad(360 / sides);
    double angle_offset    = PI / 2 + (index * corner_angle);
    MyVector2 poly_point_a = origin + MyVector2(angle_offset + rotation, radius, true);
    MyVector2 poly_point_b = origin + MyVector2(angle_offset + corner_angle + rotation, radius, true);

    return MySegment(poly_point_a, poly_point_b);
}

// Returns the number of vertices of the polygon.
int MyPolygon::getVerticesNum() { return sides; }

// Returns the vertex of the polygon that corresponds to the given index.
MyVector2 MyPolygon::getVertex(int index)
{
    assert (0 <= index && index < sides);
    return this->getSide(index).a;
}

// Draws the polygon in a raylib window.
void MyPolygon::draw(Color color) { DrawPolyLines(origin.toRayVec(), sides, radius, arithmetic::radToDeg(rotation), color); }



// -------------------- CIRCLE -------------------- //

// Constructor.
MyCircle::MyCircle()                                        : origin(MyVector2()), radius(0)          {}; // Null circle.
MyCircle::MyCircle(MyVector2 new_origin, double new_radius) : origin(new_origin),  radius(new_radius) {}; // Circle from origin and radius.

// Destroyer.
MyCircle::~MyCircle() {};

// ---------- CIRCLE METHODS ---------- //

// Returns the center of mass of the circle.
MyVector2 MyCircle::getCenterOfMass() { return origin; }

// Returns the number of sides of the circle.
int MyCircle::getSidesNum() { return 1; }

// Returns the number of vertices of the circle.
int MyCircle::getVerticesNum() { return 0; }

// Draws the circle in a raylib window.
void MyCircle::draw(Color color) { DrawCircleLines(origin.x, origin.y, radius, color); }





// ------------------------------ COLLISIONS ------------------------------ //

// Returns the smallest rectangle that contanins the given shape.
template <typename T> MyRectangle collisions::getBoundingBox(T shape)
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
MySegment collisions::CircleGetAxis(MyCircle circle, T shape)
{
    if (isShape(shape, true))
    {
        // Make a segment that starts at the center of the circle, goes in the direction of the center of the shape and is of length 1.
        return MySegment(circle.origin, MyVector2(circle.origin, shape.getCenterOfMass()).normalize(), true);
    }
}

// Returns the axis of the given shapes that corresponds to the given index.
template <typename T1, typename T2>
MySegment collisions::ShapesGetAxis(T1 shape1, T2 shape2, int index)
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
            if (typeid(T1) != typeid(MyCircle)) {
                side = shape1.ShapeGetSide(index);
                axis = MySegment(side.getCenterOfMass(),
                                    MyVector2(side).getNormalized().getNormal(),
                                    true);
            }
            // If the first shape is a circle, get its axis.
            else
                axis = CircleGetAxis(shape1, shape2);
        }
        // If the given index refers to an axis of the second shape...
        else
        {
            // If the second shape is not a circle, get the side pointed to by the index and calculate its normal.
            if (typeid(T2) != typeid(MyCircle)) {
                side = shape2.getSide(index - getSidesNum(shape1));
                axis = MySegment(side.getCenterOfMass(),
                                    MyVector2(side).getNormalized().getNormal(),
                                    true);
            }
            // If the second shape is a circle, get its axis.
            else
                axis = CircleGetAxis(shape2.data.circle, shape1);
        }

        //! Debug render.
        if (__debug_axes) (MyVector2(axis) * 100).draw(axis.a, BLUE);

        return axis;
    }
}

// Returns true if the given point is colliding with the given circle.
bool collisions::collisionCirclePoint(MyCircle c, MyVector2 p) { return (c.origin.getDistanceFromPoint(p) <= c.radius ? true : false); }

// Returns true if the given circles are in collision.
bool collisions::collisionCircles(MyCircle c1, MyCircle c2)    { return (c1.origin.getDistanceFromPoint(c2.origin) <= c1.radius + c2.radius ? true : false); }

// Checks for collision between two rectangles.
bool collisions::collisionAABB(MyRectangle rec1, MyRectangle rec2)
{
    return (rec1.origin.x + rec1.width  >= rec2.origin.x              &&
            rec1.origin.x               <= rec2.origin.x + rec2.width &&
            rec1.origin.y + rec1.height >= rec2.origin.y              &&
            rec1.origin.y               <= rec2.origin.y + rec2.height);
}

// Project a shape onto a given axis.
template <typename T>
MySegment collisions::projectShapeOnAxis(MySegment axis, T shape)
{
    if (isShape(shape, true))
    {
        // Get the axis' vector.
        MyVector2 axis_vec = MyVector2(axis);

        // Handle circles.
        if (typeid(T) == typeid(MyCircle))
        {
            // Project the circle's origin onto the axis.
            MyVector2 origin_projection = axis.a + axis_vec * (MyVector2(axis.a, shape.origin) & axis_vec);

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
            projected_points[i] = axis.a + axis_vec * (MyVector2(axis.a, vertex) & axis_vec);

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

        MyVector2 axis_orig_to_min_point = MyVector2(axis.a, min_point);
        MySegment projection = MySegment(axis.a + axis_orig_to_min_point, MyVector2(min_point, max_point), true);

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
bool collisions::collisionSegmentPoint(MySegment segment, MyVector2 point)
{
    if (arithmetic::roundInt(MyVector2(segment) ^ MyVector2(segment.a, point)) == 0)
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
bool collisions::collisionProjections(MySegment projection1, MySegment projection2)
{
    return (collisionSegmentPoint(projection1, projection2.a) ||
            collisionSegmentPoint(projection1, projection2.b) ||
            collisionSegmentPoint(projection2, projection1.a) ||
            collisionSegmentPoint(projection2, projection1.b));
}

// Checks for collision between two given shapes.
template<typename T1, typename T2>
bool collisions::collisionSAT(T1 shape1, T2 shape2)
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