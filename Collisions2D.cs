using System;
using System.Numerics;
using static System.Math;
using System.Diagnostics;
using static MyMathLib.Geometry2D;
using static MyMathLib.Arithmetic;

// TODO.

namespace MyMathLib
{
    public static class Collisions2D
    {
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
            MyRectangle bounding_box = RectangleCreate(Vector2Create(xmin, ymin), xmax - xmin, ymax - ymin);

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
        static inline MySegment projectShapeOnAxis(MySegment axis, ShapeInfo shape)
        {
            // Get the axis' vector.
            MyVector2 axis_vec = Vector2FromSegment(axis);

            // Handle circles.
            if (shape.type == CIRCLE)
            {
                // Project the circle's origin onto the axis.
                MyVector2 origin_projection = Vector2Add(axis.a, Vector2MultiplyVal(axis_vec, Vector2DotProduct(Vector2FromPoints(axis.a, shape.data.circle.origin), axis_vec)));

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
        static inline bool collisionSAT(ShapeInfo shape1, ShapeInfo shape2)
        {
            // If both shapes are circles, don't use SAT.
            if (shape1.type == CIRCLE && shape2.type == CIRCLE)
                return collisionCircles(shape1.data.circle, shape2.data.circle);

            // If both shapes are rectangles, don't use SAT.
            else if (shape1.type == RECTANGLE && shape2.type == RECTANGLE)
                return collisionAABB(shape1.data.rectangle, shape2.data.rectangle);

            // Check for collisions on the shapes' bounding boxes to not have to check if they are not in collision.
            else if (collisionAABB(getBoundingBox(shape1), getBoundingBox(shape2)))
            {
                //! Debug render.
                if (__debug_shapes) {
                    DrawShape(shape1, Vector2Zero(), GREEN); DrawShape(shape2, Vector2Zero(), GREEN);
                }

                // Get the number of sides of both shapes.
                int sides = getSidesNum(shape1) + getSidesNum(shape2);

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
                            DrawMySegment(projection1, PINK); DrawMySegment(projection2, PINK);
                        }
                        return false;
                    }
                }
                return true;
            }
            //! Debug render.
            if (__debug_shapes) {
                DrawShape(shape1, Vector2Zero(), GREEN); DrawShape(shape2, Vector2Zero(), GREEN);
            }

            return false;
        }
    }
}
