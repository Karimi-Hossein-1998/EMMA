#include "AppState.hpp"
#include "implot.h"
#include <cstring>
int main()
{
    SetConfigFlags(FLAG_WINDOW_RESIZABLE);
    AppState appState;
    appState.initW = 800;
    appState.initH = 600;
    appState.appTitle = "EMMA - Easy Mathematical Modelling App!";
    appState.BgColor = Color(15.0f,15.0f,15.0f);
    appState.DarkTheme = true;

    InitWindow(appState.initW,appState.initH,appState.appTitle.c_str());
    SetTargetFPS(60);

    // rlImGuiSetup(true);
    initUI(appState.DarkTheme);
    ImPlot::CreateContext();


    while(!WindowShouldClose())
    {
        BeginDrawing();

        	ClearBackground(appState.BgColor);

            rlImGuiBegin();
            	appState.RenderUI();
            rlImGuiEnd();

        EndDrawing();
    }

    ImPlot::DestroyContext();
    rlImGuiShutdown();
    CloseWindow();

    return 0;
}
