#include <my_math.h>
#include <stdio.h>

int main(void)
{
    // ---------- INITIALIZATION ---------- //

    __debug_shapes = true;

    // Define the screen size.
    const int screenWidth = 830;
    const int screenHeight = 1060;

    // Initialize the game window.
    InitWindow(screenWidth, screenHeight, "My Math Lib Test");
    SetTargetFPS(60);

    // Variables for dynamic shapes.
    int current_shape = 0;
    int x[2] = { GetScreenWidth() / 2, GetScreenWidth() / 2 };
    int y[2] = { GetScreenHeight() / 2, GetScreenHeight() / 2 };
    int sizes[2] = { 200, 50 };
    double rotations[2] = { 0, 0 };
    int sides[2] = { 5, 3 };

    ShapeInfo shapes[2];
    shapes[0].type = POLYGON;
    shapes[1].type = POLYGON;

    // Static shapes.
    ShapeInfo circle;
    circle.type = CIRCLE;
    circle.data.circle = CircleCreate(Vector2Create(200, 200), 90);

    ShapeInfo triangle;
    triangle.type = TRIANGLE;
    triangle.data.triangle = TriangleCreate(Vector2Create(GetScreenWidth() - 100, GetScreenHeight() - 100), Vector2Create(GetScreenWidth() - 300, GetScreenHeight() - 100), Vector2Create(GetScreenWidth() - 200, GetScreenHeight() - 250));

    ShapeInfo rectangle;
    rectangle.type = RECTANGLE;
    rectangle.data.rectangle = RectangleCreate(Vector2Create(GetScreenWidth() - 300, 150), 200, 100);

    ShapeInfo segment;
    segment.type = SEGMENT;
    segment.data.segment = SegmentCreate(Vector2Create(100, GetScreenHeight() - 200), Vector2Create(300, GetScreenHeight() - 100));

    // ---------- GAME LOOP ---------- //
    while (!WindowShouldClose())
    {
        // ------- UPDATE ------- //

        // Switch shape.
        if (IsKeyPressed(KEY_SPACE))
            current_shape = !current_shape;
        // Move shape.
        if (IsKeyDown(KEY_W))
            y[current_shape] -= 2;
        if (IsKeyDown(KEY_S))
            y[current_shape] += 2;
        if (IsKeyDown(KEY_D))
            x[current_shape] += 2;
        if (IsKeyDown(KEY_A))
            x[current_shape] -= 2;
        // Change side number.
        if (IsKeyPressed(KEY_KP_ADD) && sides[current_shape] < 360)
            sides[current_shape]++;
        if (IsKeyPressed(KEY_KP_SUBTRACT) && sides[current_shape] > 3)
            sides[current_shape]--;
        // Change size.
        if (IsKeyDown(KEY_KP_5) && sizes[current_shape] <= 1000)
            sizes[current_shape] += 5;
        if (IsKeyDown(KEY_KP_4) && sizes[current_shape] > 10)
            sizes[current_shape] -= 5;
        // Rotate.
        if (IsKeyDown(KEY_KP_2))
            rotations[current_shape] += PI / sizes[current_shape];
        if (IsKeyDown(KEY_KP_1))
            rotations[current_shape] -= PI / sizes[current_shape];
        // Circle button.
        if (IsKeyPressed(KEY_KP_8))
            sides[current_shape] = 360;
        if (IsKeyPressed(KEY_KP_7))
            sides[current_shape] = 3;
        // Show/hide shapes.
        if (IsKeyPressed(KEY_KP_DECIMAL))
            __debug_shapes = !__debug_shapes;
        // Show/hide bounding boxes.
        if (IsKeyPressed(KEY_KP_3))
            __debug_bounding_boxes = !__debug_bounding_boxes;
        // Show/hide axes.
        if (IsKeyPressed(KEY_KP_6))
            __debug_axes = !__debug_axes;
        // Show/hide projections.
        if (IsKeyPressed(KEY_KP_9))
            __debug_projections = !__debug_projections;
        // Show/hide failed projections.
        if (IsKeyPressed(KEY_KP_MULTIPLY))
            __debug_failed_projections = !__debug_failed_projections;


        // ------- DISPLAY ------- //
        BeginDrawing();
        {
            ClearBackground(BLACK);

            // ---------- SANDBOX COLLISIONS ---------- //

            // Update the two dynamic polygons.
            for (int i = 0; i < 2; i++)
            {
                shapes[i].data.polygon = PolygonCreate(Vector2Create(x[i], y[i]), sizes[i], rotations[i], sides[i]);
            }

            // Test collisions using SAT.
            bool dynamic_shapes = collisionSAT(shapes[0], shapes[1]);
            bool circle_shape[2] = { collisionSAT(circle, shapes[0]), collisionSAT(circle, shapes[1]) };
            bool triangle_shape[2] = { collisionSAT(triangle, shapes[0]), collisionSAT(triangle, shapes[1]) };
            bool rectangle_shape[2] = { collisionSAT(rectangle, shapes[0]), collisionSAT(rectangle, shapes[1]) };
            bool segment_shape[2] = { collisionSAT(segment, shapes[0]), collisionSAT(segment, shapes[1]) };

            // Show the colliding shapes.
            if (__debug_shapes)
            {
                for (int i = 0; i < 2; i++) 
                    if (dynamic_shapes || circle_shape[i] || triangle_shape[i] || rectangle_shape[i] || segment_shape[i])
                        DrawShape(shapes[i], Vector2Zero(), RED);

                if (circle_shape[0] || circle_shape[1])
                    DrawShape(circle, Vector2Zero(), RED);
                if (triangle_shape[0] || triangle_shape[1])
                    DrawShape(triangle, Vector2Zero(), RED);
                if (rectangle_shape[0] || rectangle_shape[1])
                    DrawShape(rectangle, Vector2Zero(), RED);
                if (segment_shape[0] || segment_shape[1])
                    DrawShape(segment, Vector2Zero(), RED);
            }

            // FPS.
            DrawText(TextFormat("%d", GetFPS()), GetScreenWidth() - 50, 25, 15, WHITE);
        }
        EndDrawing();
    }

    // ---------- DE-INITIALIZATION ---------- //
    CloseWindow();

    return 0;
}