#include "AppState.hpp"
#include "implot.h"
#include <cstring>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#endif

static AppState* g_appState = nullptr;

static void UpdateDrawFrame()
{
    if (!g_appState) return;

    BeginDrawing();
        ClearBackground(g_appState->BgColor);
        rlImGuiBegin();
            g_appState->RenderUI();
        rlImGuiEnd();
    EndDrawing();
}

int main()
{
    SetConfigFlags(FLAG_WINDOW_RESIZABLE);
    static AppState appState;
    g_appState = &appState;

    appState.initW = 800;
    appState.initH = 600;
    appState.appTitle = "EMMA - Easy Mathematical Modelling App!";
    appState.BgColor = Color{15,15,15,255};
    appState.DarkTheme = true;

    InitWindow(appState.initW,appState.initH,appState.appTitle.c_str());
    SetTargetFPS(60);

    // rlImGuiSetup(true);
    initUI(appState.DarkTheme);
    ImPlot::CreateContext();

#ifdef __EMSCRIPTEN__
    emscripten_set_main_loop(UpdateDrawFrame, 0, 1);
#else
    while(!WindowShouldClose())
    {
        UpdateDrawFrame();
    }
#endif

    ImPlot::DestroyContext();
    rlImGuiShutdown();
    CloseWindow();

    return 0;
}
