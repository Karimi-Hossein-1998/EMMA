#include "AppState.hpp"
#include <cstring>
int main()
{
    SetConfigFlags(FLAG_WINDOW_RESIZABLE);
    AppState appState;
    appState.initW = 800;
    appState.initH = 600;
    appState.appTitle = "EMMA - Easy Mathematical Modelling App!";
    appState.BgColor = Color(15.0f,15.0f,15.0f);

    InitWindow(appState.initW,appState.initH,appState.appTitle.c_str());
    SetTargetFPS(60);

    // rlImGuiSetup(true);
    initUI(true);


    while(!WindowShouldClose())
    {
        BeginDrawing();

        	ClearBackground(appState.BgColor);

            rlImGuiBegin();
            	appState.RenderUI();
            rlImGuiEnd();

        EndDrawing();
    }

    rlImGuiShutdown();
    CloseWindow();

    return 0;
}
