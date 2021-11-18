#pragma once
#include "my_math.hpp"

using namespace MyMathLib;
using namespace MyMathLib::geometry;

// ----------- MATRIX ---------- //

// Matrix addition.
template <typename T>
arithmetic::MyMatrix arithmetic::MyMatrix::operator+(const T& val) const
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < columns; j++)
            matrix[i][j] += val;
    return matrix;
}

template<>
MyMatrix MyMatrix::operator+<MyMatrix>(const MyMatrix& val) const
{
    assert(rows == val.rows && columns == val.columns); // Matrix must have the same dimension
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < columns; j++)
            matrix[i][j] += val[i][j];
    return matrix;
}

// Matrix substraction.
template <typename T>
MyMatrix MyMatrix::operator-(const T& val) const
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < columns; j++)
            matrix[i][j] -= val;
    return matrix;
}
template <>
MyMatrix MyMatrix::operator-<MyMatrix>(const MyMatrix& val) const
{
    assert(rows == val.rows && columns == val.columns); // Matrix must have the same dimension
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < columns; j++)
            matrix[i][j] -= val[i][j];
    return matrix;
}

// Matrix multiplication.
template <typename T>
MyMatrix MyMatrix::operator*(const T& val) const
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < columns; j++)
            matrix[i][j] *= val;
    return matrix;
}

template <typename T>
MyMatrix MyMatrix::operator*<MyMatrix>(const MyMatrix& m) const
{
    assert(columns == m.rows); // Size condition to calculate
    
    MyMatrix<T> result(columns > m.columns ? columns : m.rows, rows > m.rows ? rows : m.rows);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < m.columns; j++)
        {
            result[i][j] = 0;
            for (int k = 0; k < columns; k++) result[i][j] = matrix[i][k] + m[k][j];
        }
    return result;
}

// Vector2 division.
template <typename T>
MyMatrix MyMatrix::operator/(const T& val) const
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < columns; j++)
            matrix[i][j] /= val;
    return matrix;
}
template <>
MyMatrix MyMatrix::operator/<MyMatrix>(const MyMatrix& m) const
{
    //todo
}

// Vector2 addition assignement.
template <typename T>
void MyVector2::operator+=(const T& val)
{
    x += val;
    y += val;
}
template <>
void MyVector2::operator+=<MyVector2>(const MyVector2& val)
{
    x += val.x; 
    y += val.y;
}

// Vector2 substraction assignement.
template <typename T>
void MyVector2::operator-=(const T& val)
{
    x -= val;
    y -= val;
}
template <>
void MyVector2::operator-=<MyVector2>(const MyVector2& val)
{
    x -= val.x; 
    y -= val.y;
}

// Vector2 multiplication assignement.
template <typename T>
void MyVector2::operator*=(const T& val)
{
    x *= val;
    y *= val;
}
template <>
void MyVector2::operator*=<MyVector2>(const MyVector2& val)
{
    x *= val.x; 
    y *= val.y;
}

// Vector2 division assignement.
template <typename T>
void MyVector2::operator/=(const T& val)
{
    x /= val;
    y /= val;
}
template <>
void MyVector2::operator/=<MyVector2>(const MyVector2& val)
{
    x /= val.x; 
    y /= val.y;
}

// ---------- VECTOR2 ---------- //

// Vector2 addition.
template <typename T>
MyVector2 MyVector2::operator+(const T& val) const
{
    return MyVector2(x + val, y + val);
}
template<>
MyVector2 MyVector2::operator+<MyVector2>(const MyVector2& val) const
{
    return MyVector2(x + val.x, y + val.y);
}

// Vector2 substraction.
template <typename T>
MyVector2 MyVector2::operator-(const T& val) const
{
    return MyVector2(x - val,   y - val  );
}
template <>
MyVector2 MyVector2::operator-<MyVector2>(const MyVector2& val) const
{
    return MyVector2(x - val.x, y - val.y);
}

// Vector2 multiplication.
template <typename T>
MyVector2 MyVector2::operator*(const T& val) const
{
    return MyVector2(x * val,   y * val  );
}
template <>
MyVector2 MyVector2::operator*<MyVector2>(const MyVector2& val) const
{
    return MyVector2(x * val.x, y * val.y);
}

// Vector2 division.
template <typename T>
MyVector2 MyVector2::operator/(const T& val) const
{
    return MyVector2(x / val,   y / val  );
}
template <>
MyVector2 MyVector2::operator/<MyVector2>(const MyVector2& val) const
{
    return MyVector2(x / val.x, y / val.y);
}

// Vector2 addition assignement.
template <typename T>
void MyVector2::operator+=(const T& val)
{
    x += val;
    y += val;
}
template <>
void MyVector2::operator+=<MyVector2>(const MyVector2& val)
{
    x += val.x; 
    y += val.y;
}

// Vector2 substraction assignement.
template <typename T>
void MyVector2::operator-=(const T& val)
{
    x -= val;
    y -= val;
}
template <>
void MyVector2::operator-=<MyVector2>(const MyVector2& val)
{
    x -= val.x; 
    y -= val.y;
}

// Vector2 multiplication assignement.
template <typename T>
void MyVector2::operator*=(const T& val)
{
    x *= val;
    y *= val;
}
template <>
void MyVector2::operator*=<MyVector2>(const MyVector2& val)
{
    x *= val.x; 
    y *= val.y;
}

// Vector2 division assignement.
template <typename T>
void MyVector2::operator/=(const T& val)
{
    x /= val;
    y /= val;
}
template <>
void MyVector2::operator/=<MyVector2>(const MyVector2& val)
{
    x /= val.x; 
    y /= val.y;
}

// ---------- COLLISIONS ---------- //

// Returns the smallest rectangle that contanins the given shape.
template <typename T>
MyRectangle collisions::getBoundingBox(T shape)
{
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

// Bounding box specialization for a circle.
template <>
MyRectangle collisions::getBoundingBox<MyCircle>(MyCircle shape)
{
    MyRectangle bounding_box = MyRectangle(shape.origin - shape.radius, 
                                           shape.radius * 2, 
                                           shape.radius * 2);

    //! Debug render.
    if (__debug_bounding_boxes) bounding_box.draw(GRAY);

    return bounding_box;
}

// Returns an axis that passes through the center of the given circle and the center of the given shape.
template <typename T>
MySegment collisions::CircleGetAxis(MyCircle circle, T shape)
{
    // Make a segment that starts at the center of the circle, goes in the direction of the center of the shape and is of length 1.
    return MySegment(circle.origin, MyVector2(circle.origin, shape.getCenterOfMass()).getNormalized(), true);
}

// Returns the axis of the given shapes that corresponds to the given index.
template <typename T1, typename T2>
MySegment collisions::ShapesGetAxis(T1 shape1, T2 shape2, int index)
{
    assert (0 <= index && index < shape1.getSidesNum() + shape2.getSidesNum());

    MySegment side;
    MySegment axis;

    // If the given index refers to an axis of the first shape...
    if (index < shape1.getSidesNum())
    {
        // If the first shape is not a circle, get the side pointed to by the index and calculate its normal.
        side = shape1.getSide(index);
        axis = MySegment(side.getCenterOfMass(),
                            MyVector2(side).getNormalized().getNormal(),
                            true);
    }
    // If the given index refers to an axis of the second shape...
    else
    {
        // If the second shape is not a circle, get the side pointed to by the index and calculate its normal.
        side = shape2.getSide(index - shape1.getSidesNum());
        axis = MySegment(side.getCenterOfMass(),
                            MyVector2(side).getNormalized().getNormal(),
                            true);
    }

    //! Debug render.
    if (__debug_axes) (MyVector2(axis) * 100).draw(axis.a, BLUE);

    return axis;
}

// Shapes axis specialization for a circle and a shape.
template <>
MySegment collisions::ShapesGetAxis<MyCircle, MySegment>(MyCircle shape1, MySegment shape2, int index)
{
    assert (0 <= index && index < shape1.getSidesNum() + shape2.getSidesNum());

    MySegment side;
    MySegment axis;

    // If the given index refers to an axis of the first shape...
    if (index < shape1.getSidesNum())
    {
        // If the first shape is a circle, get its axis.
        axis = CircleGetAxis(shape1, shape2);
    }
    // If the given index refers to an axis of the second shape...
    else
    {
        // If the second shape is not a circle, get the side pointed to by the index and calculate its normal.
        side = shape2.getSide(index - shape1.getSidesNum());
        axis = MySegment(side.getCenterOfMass(),
                            MyVector2(side).getNormalized().getNormal(),
                            true);
    }

    //! Debug render.
    if (__debug_axes) (MyVector2(axis) * 100).draw(axis.a, BLUE);

    return axis;
}

// Shapes axis specialization for a circle and a shape.
template <>
MySegment collisions::ShapesGetAxis<MyCircle, MyTriangle>(MyCircle shape1, MyTriangle shape2, int index)
{
    assert (0 <= index && index < shape1.getSidesNum() + shape2.getSidesNum());

    MySegment side;
    MySegment axis;

    // If the given index refers to an axis of the first shape...
    if (index < shape1.getSidesNum())
    {
        // If the first shape is a circle, get its axis.
        axis = CircleGetAxis(shape1, shape2);
    }
    // If the given index refers to an axis of the second shape...
    else
    {
        // If the second shape is not a circle, get the side pointed to by the index and calculate its normal.
        side = shape2.getSide(index - shape1.getSidesNum());
        axis = MySegment(side.getCenterOfMass(),
                            MyVector2(side).getNormalized().getNormal(),
                            true);
    }

    //! Debug render.
    if (__debug_axes) (MyVector2(axis) * 100).draw(axis.a, BLUE);

    return axis;
}

// Shapes axis specialization for a circle and a shape.
template <>
MySegment collisions::ShapesGetAxis<MyCircle, MyRectangle>(MyCircle shape1, MyRectangle shape2, int index)
{
    assert (0 <= index && index < shape1.getSidesNum() + shape2.getSidesNum());

    MySegment side;
    MySegment axis;

    // If the given index refers to an axis of the first shape...
    if (index < shape1.getSidesNum())
    {
        // If the first shape is a circle, get its axis.
        axis = CircleGetAxis(shape1, shape2);
    }
    // If the given index refers to an axis of the second shape...
    else
    {
        // If the second shape is not a circle, get the side pointed to by the index and calculate its normal.
        side = shape2.getSide(index - shape1.getSidesNum());
        axis = MySegment(side.getCenterOfMass(),
                            MyVector2(side).getNormalized().getNormal(),
                            true);
    }

    //! Debug render.
    if (__debug_axes) (MyVector2(axis) * 100).draw(axis.a, BLUE);

    return axis;
}

// Shapes axis specialization for a circle and a shape.
template <>
MySegment collisions::ShapesGetAxis<MyCircle, MyPolygon>(MyCircle shape1, MyPolygon shape2, int index)
{
    assert (0 <= index && index < shape1.getSidesNum() + shape2.getSidesNum());

    MySegment side;
    MySegment axis;

    // If the given index refers to an axis of the first shape...
    if (index < shape1.getSidesNum())
    {
        // If the first shape is a circle, get its axis.
        axis = CircleGetAxis(shape1, shape2);
    }
    // If the given index refers to an axis of the second shape...
    else
    {
        // If the second shape is not a circle, get the side pointed to by the index and calculate its normal.
        side = shape2.getSide(index - shape1.getSidesNum());
        axis = MySegment(side.getCenterOfMass(),
                            MyVector2(side).getNormalized().getNormal(),
                            true);
    }

    //! Debug render.
    if (__debug_axes) (MyVector2(axis) * 100).draw(axis.a, BLUE);

    return axis;
}

// Shapes axis specialization for a shape and a circle.
template <>
MySegment collisions::ShapesGetAxis<MySegment, MyCircle>(MySegment shape1, MyCircle shape2, int index)
{
    assert (0 <= index && index < shape1.getSidesNum() + shape2.getSidesNum());

    MySegment side;
    MySegment axis;

    // If the given index refers to an axis of the first shape...
    if (index < shape1.getSidesNum())
    {
        // If the first shape is not a circle, get the side pointed to by the index and calculate its normal.
        side = shape1.getSide(index);
        axis = MySegment(side.getCenterOfMass(),
                            MyVector2(side).getNormalized().getNormal(),
                            true);
    }
    // If the given index refers to an axis of the second shape...
    else
    {
        // If the second shape is a circle, get its axis.
        axis = CircleGetAxis(shape2, shape1);
    }

    //! Debug render.
    if (__debug_axes) (MyVector2(axis) * 100).draw(axis.a, BLUE);

    return axis;
}

template <>
MySegment collisions::ShapesGetAxis<MyTriangle, MyCircle>(MyTriangle shape1, MyCircle shape2, int index)
{
    assert (0 <= index && index < shape1.getSidesNum() + shape2.getSidesNum());

    MySegment side;
    MySegment axis;

    // If the given index refers to an axis of the first shape...
    if (index < shape1.getSidesNum())
    {
        // If the first shape is not a circle, get the side pointed to by the index and calculate its normal.
        side = shape1.getSide(index);
        axis = MySegment(side.getCenterOfMass(),
                            MyVector2(side).getNormalized().getNormal(),
                            true);
    }
    // If the given index refers to an axis of the second shape...
    else
    {
        // If the second shape is a circle, get its axis.
        axis = CircleGetAxis(shape2, shape1);
    }

    //! Debug render.
    if (__debug_axes) (MyVector2(axis) * 100).draw(axis.a, BLUE);

    return axis;
}

template <>
MySegment collisions::ShapesGetAxis<MyRectangle, MyCircle>(MyRectangle shape1, MyCircle shape2, int index)
{
    assert (0 <= index && index < shape1.getSidesNum() + shape2.getSidesNum());

    MySegment side;
    MySegment axis;

    // If the given index refers to an axis of the first shape...
    if (index < shape1.getSidesNum())
    {
        // If the first shape is not a circle, get the side pointed to by the index and calculate its normal.
        side = shape1.getSide(index);
        axis = MySegment(side.getCenterOfMass(),
                            MyVector2(side).getNormalized().getNormal(),
                            true);
    }
    // If the given index refers to an axis of the second shape...
    else
    {
        // If the second shape is a circle, get its axis.
        axis = CircleGetAxis(shape2, shape1);
    }

    //! Debug render.
    if (__debug_axes) (MyVector2(axis) * 100).draw(axis.a, BLUE);

    return axis;
}


template <>
MySegment collisions::ShapesGetAxis<MyPolygon, MyCircle>(MyPolygon shape1, MyCircle shape2, int index)
{
    assert (0 <= index && index < shape1.getSidesNum() + shape2.getSidesNum());

    MySegment side;
    MySegment axis;

    // If the given index refers to an axis of the first shape...
    if (index < shape1.getSidesNum())
    {
        // If the first shape is not a circle, get the side pointed to by the index and calculate its normal.
        side = shape1.getSide(index);
        axis = MySegment(side.getCenterOfMass(),
                            MyVector2(side).getNormalized().getNormal(),
                            true);
    }
    // If the given index refers to an axis of the second shape...
    else
    {
        // If the second shape is a circle, get its axis.
        axis = CircleGetAxis(shape2, shape1);
    }

    //! Debug render.
    if (__debug_axes) (MyVector2(axis) * 100).draw(axis.a, BLUE);

    return axis;
}

// Shapes axis specialization for 2 circles.
template <>
MySegment collisions::ShapesGetAxis<MyCircle, MyCircle>(MyCircle shape1, MyCircle shape2, int index)
{
    assert (0 <= index && index < shape1.getSidesNum() + shape2.getSidesNum());

    MySegment side;
    MySegment axis;

    // If the given index refers to an axis of the first shape...
    if (index < shape1.getSidesNum())
    {
        // If the first shape is a circle, get its axis.
        axis = CircleGetAxis(shape1, shape2);
    }
    // If the given index refers to an axis of the second shape...
    else
    {
        // If the second shape is a circle, get its axis.
        axis = CircleGetAxis(shape2, shape1);
    }

    //! Debug render.
    if (__debug_axes) (MyVector2(axis) * 100).draw(axis.a, BLUE);

    return axis;
}

// Project a shape onto a given axis.
template <typename T>
MySegment collisions::projectShapeOnAxis(MySegment axis, T shape)
{
    // Get the axis' vector.
    MyVector2 axis_vec = MyVector2(axis);

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

// Projection specialization for a circle.
template <>
MySegment collisions::projectShapeOnAxis<MyCircle>(MySegment axis, MyCircle shape)
{
    // Get the axis' vector.
    MyVector2 axis_vec = MyVector2(axis);

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

// Checks for collision between two given shapes.
template<typename T1, typename T2>
bool collisions::collisionSAT(T1 shape1, T2 shape2)
{
    // Check for collisions on the shapes' bounding boxes to not have to check if they are not in collision.
    if (collisionAABB(getBoundingBox(shape1), getBoundingBox(shape2)))
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
    return false;
}

// SAT specialization for circles.
template<>
bool collisions::collisionSAT<MyCircle, MyCircle>(MyCircle shape1, MyCircle shape2)
{
    return collisionCircles(shape1, shape2);
}

// SAT specialization for rectangles.
template<>
bool collisions::collisionSAT<MyRectangle, MyRectangle>(MyRectangle shape1, MyRectangle shape2)
{
    return collisionAABB(shape1, shape2);
}