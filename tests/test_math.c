#include <my_math.h>
#include <stdio.h>

int main(void)
{
    // ---------- INITIALIZATION ---------- //

    // Define the screen size.
    const int screenWidth = 830;
    const int screenHeight = 1060;

    // Initialize the game window.
    InitWindow(screenWidth, screenHeight, "My Math Lib Test");
    SetTargetFPS(60);

    int current_shape = 0;
    int x[2] = { GetScreenWidth() / 2, GetScreenWidth() / 2 };
    int y[2] = { GetScreenHeight() / 2, GetScreenHeight() / 2 };
    int sizes[2] = { 200, 50 };
    double rotations[2] = { 0, 0 };
    int sides[2] = { 5, 3 };

    ShapeInfo shapes[2];
    shapes[0].type = POLYGON;
    shapes[1].type = POLYGON;
    MyRectangle bounding_boxes[2];

    // ---------- GAME LOOP ---------- //
    while (!WindowShouldClose())
    {
        // ------- UPDATE ------- //
        //! KEYBINDS AXIS AND PROJECTIONS DISPLAY
        if (IsKeyPressed(KEY_SPACE))
            current_shape = !current_shape;
        if (IsKeyDown(KEY_W))
            y[current_shape] -= 2;
        if (IsKeyDown(KEY_S))
            y[current_shape] += 2;
        if (IsKeyDown(KEY_D))
            x[current_shape] += 2;
        if (IsKeyDown(KEY_A))
            x[current_shape] -= 2;
        if (IsKeyPressed(KEY_KP_ADD) && sides[current_shape] <= 360)
            sides[current_shape]++;
        if (IsKeyPressed(KEY_KP_SUBTRACT) && sides[current_shape] > 3)
            sides[current_shape]--;
        if (IsKeyDown(KEY_KP_MULTIPLY) && sizes[current_shape] <= 1000)
            sizes[current_shape] += 5;
        if (IsKeyDown(KEY_KP_DIVIDE) && sizes[current_shape] > 10)
            sizes[current_shape] -= 5;
        if (IsKeyDown(KEY_KP_0))
            rotations[current_shape] += PI / sizes[current_shape];
        if (IsKeyDown(KEY_KP_DECIMAL))
            rotations[current_shape] -= PI / sizes[current_shape];
        if (IsKeyPressed(KEY_KP_9))
            sides[current_shape] = 360;
        if (IsKeyPressed(KEY_KP_8) || IsKeyPressed(KEY_KP_7))
            sides[current_shape] = 3;

        // ------- DISPLAY ------- //
        BeginDrawing();
        {
            ClearBackground(BLACK);

            //! ------ SANDBOX COLLISIONS ------ //

            for (int i = 0; i < 2; i++)
            {
                shapes[i].data.polygon = PolygonCreate(Vector2Create(x[i], y[i]), sizes[i], rotations[i], sides[i]);
                bounding_boxes[i] = getBoundingBox(shapes[i]);
            }

            //! Test collisions using SAT.
            bool AABB_result = collisionAABB(bounding_boxes[0], bounding_boxes[1]);
            bool SAT_result = collisionSAT(shapes[0], shapes[1]);

            // Show the bounding boxes and the shapes.
            for (int i = 0; i < 2; i++) 
            {
                DrawMyRectangle(bounding_boxes[i], (AABB_result ? RED : GREEN));
                DrawShape(shapes[i], Vector2Zero(), (SAT_result ? RED : GREEN));
            }

            //! Text indications.
            DrawText(TextFormat("%d", GetFPS()), GetScreenWidth() - 50, 25, 15, WHITE);
        }
        EndDrawing();
    }

    // ---------- DE-INITIALIZATION ---------- //
    CloseWindow();

    return 0;
}