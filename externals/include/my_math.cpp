#include "my_math.hpp"

using namespace MyMathLib;

// ------------------------------ ARITHMECTIC ------------------------------ //

// --------- ARITHMETIC FUNCTIONS -------- //

// Fast inverse square root from Quake III.
float arithmetic::Q_rsqrt(float number)
{
    long i;
    float x2, y;
    const float threehalfs = 1.5F;

    x2 = number * 0.5F;
    y = number;
    i = * ( long * ) &y; // evil floating point bit level hacking
    i = 0x5f3759df - ( i >> 1 ); // what the fuck?
    y = * ( float * ) &i;
    y = y * ( threehalfs - ( x2 * y * y ) ); // 1st iteration
//	y = y * ( threehalfs - ( x2 * y * y ) ); // 2nd iteration, this can be removed

    return y;
}

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

// Returns true if the given number is a power of 2.
bool arithmetic::is_power_of_two(int val)
{
    return val == (int)pow(2, (int)(log(val)/log(2)));
}

// Returns the closest power of 2 that is inferior or equal to val.
int arithmetic::get_power_of_two_under(int val)
{
    return (int)pow(2, (int)(log(val)/log(2)));
}

// Returns the closest power of 2 that is superior or equal to val.
int arithmetic::get_power_of_two_above(int val)
{
    if (is_power_of_two(val)) return (int)pow(2, (int)(log(val)/log(2)));
    else                      return (int)pow(2, (int)(log(val)/log(2)) + 1);
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

// --------------- VECTOR 2 -------------- //

using namespace geometry;

// Constructors.
MyVector2::MyVector2()                                                          : x(0),                 y(0)                 {}; // Null vector.
MyVector2::MyVector2(const double new_x, const double new_y)                    : x(new_x),             y(new_y)             {}; // Vector with 2 coordinates.
MyVector2::MyVector2(const MyVector2& p1, const MyVector2& p2)                  : x(p2.x - p1.x),       y(p2.y - p1.y)       {}; // Vector from 2 points.
MyVector2::MyVector2(const MySegment& seg)                                      : x(seg.b.x - seg.a.x), y(seg.b.y - seg.a.y) {}; // Vector from segment.
MyVector2::MyVector2(const double rad, const double length, const bool isAngle) : x(cos(rad) * length), y(sin(rad) * length) {}; // Vector from angle (useless bool).

// Destroyer.
MyVector2::~MyVector2() {};

// ---------- VECTOR 2 OPERATORS ---------- //

// Copy constructor.
void MyVector2::operator=(const MyVector2& other) { x = other.x; y = other.y; }

// Vector2 dot product.
double MyVector2::operator&(const MyVector2& val) const { return (x * val.x) + (y * val.y); }

// Vector2 cross product.
double MyVector2::operator^(const MyVector2& val) const { return (x * val.y) - (y * val.x); }

// ------------ VECTOR2 METHODS ----------- //

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
MySegment::MySegment()                                                                 : a(MyVector2()), b(MyVector2())   {}; // Null segment.
MySegment::MySegment(const MyVector2& new_a, const MyVector2& new_b)                   : a(new_a),       b(new_b)         {}; // Segment from points.
MySegment::MySegment(const MyVector2& origin, const MyVector2& vec, const bool vector) : a(origin),      b(origin + vec)  {}; // Segment from point and vector.

// Destroyer.
MySegment::~MySegment() {};

// ---------------- SEGMENT METHODS --------------- //

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

// Does nothing and returns a null segment.
MySegment getSide(int index) { return MySegment(); }

// Returns the number of vertices of the circle.
int MyCircle::getVerticesNum() { return 0; }

// Does nothing and returns a null vector.
MyVector2 getVertex(int index) { return MyVector2(); }

// Draws the circle in a raylib window.
void MyCircle::draw(Color color) { DrawCircleLines(origin.x, origin.y, radius, color); }

// ------------------------------ COLLISIONS ------------------------------ //

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