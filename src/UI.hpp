#pragma once
#include "raylib.h"
#include "rlImGui.h"
#include "imgui.h"
// Fonts
#include "InterFont.hpp"
#include "BrawlerRegular.hpp"
#include "BNSCRegular.hpp"
#include "CSRegular.hpp"

struct AppFonts
{
	ImFont* inter = nullptr;
	ImFont* bonanova = nullptr;
	ImFont* brawler = nullptr;
	ImFont* clearsans = nullptr;
}g_FONTs;

inline void ApplyCustomStyle()
{
    ImGuiStyle& style = ImGui::GetStyle();

    // --- Rounding (Floats) ---
    style.WindowRounding    = 6.0f; // Main window corners
    style.ChildRounding     = 6.0f; // Sub-panels
    style.FrameRounding     = 12.0f; // Buttons, input fields, checkboxes
    style.PopupRounding     = 6.0f; // Dropdown menus and popups
    style.ScrollbarRounding = 12.0f; // Scrollbars
    style.GrabRounding      = 6.0f; // Slider thumb handles
    style.TabRounding       = 9.0f; // Tab bars
    style.MenuItemRounding  = 9.0f;

    // --- Padding & Spacing (ImVec2) ---
    style.WindowPadding     = ImVec2(10.0f, 10.0f); // Space inside window edges
    style.FramePadding      = ImVec2(8.0f, 4.0f);   // Inner padding for buttons/inputs
    style.ItemSpacing       = ImVec2(9.0f, 6.0f);   // Gap between distinct UI elements
    style.ItemInnerSpacing  = ImVec2(6.0f, 4.0f);   // Gap inside composite widgets (e.g., checkbox + label)

    // --- Borders (Floats) ---
    style.WindowBorderSize  = 1.0f;
    style.FrameBorderSize   = 0.0f;
    style.PopupBorderSize   = 1.0f;

    // Separator
    style.SeparatorTextAlign   = ImVec2(0.5f,0.5f);
    style.SeparatorTextPadding = ImVec2(0.0f,20.0f);

	// Window
    style.WindowTitleAlign = ImVec2(0.5f,0.5f);

    ImVec4* colors = ImGui::GetStyle().Colors;
    colors[ImGuiCol_Text]                   = ImVec4(1.00f, 1.00f, 1.00f, 1.00f);
    colors[ImGuiCol_TextDisabled]           = ImVec4(0.50f, 0.50f, 0.50f, 1.00f);
    colors[ImGuiCol_WindowBg]               = ImVec4(0.06f, 0.06f, 0.06f, 1.00f);
    colors[ImGuiCol_ChildBg]                = ImVec4(0.00f, 0.00f, 0.00f, 1.00f);
    colors[ImGuiCol_PopupBg]                = ImVec4(0.04f, 0.04f, 0.04f, 0.94f);
    colors[ImGuiCol_Border]                 = ImVec4(0.43f, 0.43f, 0.64f, 0.60f);
    colors[ImGuiCol_BorderShadow]           = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
    colors[ImGuiCol_FrameBg]                = ImVec4(0.24f, 0.16f, 0.20f, 1.00f);
    colors[ImGuiCol_FrameBgHovered]         = ImVec4(0.63f, 0.31f, 0.39f, 1.00f);
    colors[ImGuiCol_FrameBgActive]          = ImVec4(0.51f, 0.20f, 0.31f, 1.00f);
    colors[ImGuiCol_TitleBg]                = ImVec4(0.04f, 0.04f, 0.04f, 1.00f);
    colors[ImGuiCol_TitleBgActive]          = ImVec4(0.10f, 0.10f, 0.10f, 1.00f);
    colors[ImGuiCol_TitleBgCollapsed]       = ImVec4(0.00f, 0.00f, 0.00f, 0.51f);
    colors[ImGuiCol_MenuBarBg]              = ImVec4(0.14f, 0.14f, 0.14f, 1.00f);
    colors[ImGuiCol_ScrollbarBg]            = ImVec4(0.02f, 0.02f, 0.02f, 0.53f);
    colors[ImGuiCol_ScrollbarGrab]          = ImVec4(0.31f, 0.31f, 0.31f, 1.00f);
    colors[ImGuiCol_ScrollbarGrabHovered]   = ImVec4(0.41f, 0.41f, 0.41f, 1.00f);
    colors[ImGuiCol_ScrollbarGrabActive]    = ImVec4(0.51f, 0.51f, 0.51f, 1.00f);
    colors[ImGuiCol_CheckMark]              = ImVec4(1.00f, 0.39f, 0.80f, 1.00f);
    colors[ImGuiCol_CheckboxSelectedBg]     = ImVec4(0.39f, 0.31f, 0.39f, 0.45f);
    colors[ImGuiCol_SliderGrab]             = ImVec4(0.82f, 0.67f, 0.08f, 1.00f);
    colors[ImGuiCol_SliderGrabActive]       = ImVec4(0.67f, 0.51f, 0.04f, 1.00f);
    colors[ImGuiCol_Button]                 = ImVec4(0.29f, 0.20f, 0.20f, 1.00f);
    colors[ImGuiCol_ButtonHovered]          = ImVec4(0.65f, 0.31f, 0.31f, 1.00f);
    colors[ImGuiCol_ButtonActive]           = ImVec4(0.22f, 0.14f, 0.14f, 1.00f);
    colors[ImGuiCol_Header]                 = ImVec4(0.35f, 0.24f, 0.24f, 0.31f);
    colors[ImGuiCol_HeaderHovered]          = ImVec4(0.29f, 0.20f, 0.20f, 1.00f);
    colors[ImGuiCol_HeaderActive]           = ImVec4(0.24f, 0.16f, 0.16f, 1.00f);
    colors[ImGuiCol_Separator]              = ImVec4(0.43f, 0.43f, 0.50f, 0.50f);
    colors[ImGuiCol_SeparatorHovered]       = ImVec4(0.10f, 0.40f, 0.75f, 0.78f);
    colors[ImGuiCol_SeparatorActive]        = ImVec4(0.10f, 0.40f, 0.75f, 0.00f);
    colors[ImGuiCol_ResizeGrip]             = ImVec4(0.26f, 0.59f, 0.98f, 0.20f);
    colors[ImGuiCol_ResizeGripHovered]      = ImVec4(0.26f, 0.59f, 0.98f, 0.67f);
    colors[ImGuiCol_ResizeGripActive]       = ImVec4(0.26f, 0.59f, 0.98f, 0.95f);
    colors[ImGuiCol_InputTextCursor]        = ImVec4(1.00f, 1.00f, 1.00f, 1.00f);
    colors[ImGuiCol_TabHovered]             = ImVec4(0.39f, 0.27f, 0.39f, 1.00f);
    colors[ImGuiCol_Tab]                    = ImVec4(0.24f, 0.00f, 0.24f, 1.00f);
    colors[ImGuiCol_TabSelected]            = ImVec4(0.31f, 0.16f, 0.31f, 1.00f);
    colors[ImGuiCol_TabSelectedOverline]    = ImVec4(0.59f, 0.27f, 0.59f, 1.00f);
    colors[ImGuiCol_TabDimmed]              = ImVec4(0.07f, 0.10f, 0.15f, 0.97f);
    colors[ImGuiCol_TabDimmedSelected]      = ImVec4(0.14f, 0.26f, 0.42f, 1.00f);
    colors[ImGuiCol_TabDimmedSelectedOverline]  = ImVec4(0.50f, 0.50f, 0.50f, 0.00f);
    colors[ImGuiCol_DockingPreview]         = ImVec4(0.26f, 0.59f, 0.98f, 0.70f);
    colors[ImGuiCol_DockingEmptyBg]         = ImVec4(0.20f, 0.20f, 0.20f, 1.00f);
    colors[ImGuiCol_PlotLines]              = ImVec4(0.61f, 0.61f, 0.61f, 1.00f);
    colors[ImGuiCol_PlotLinesHovered]       = ImVec4(1.00f, 0.43f, 0.35f, 1.00f);
    colors[ImGuiCol_PlotHistogram]          = ImVec4(0.90f, 0.70f, 0.00f, 1.00f);
    colors[ImGuiCol_PlotHistogramHovered]   = ImVec4(1.00f, 0.60f, 0.00f, 1.00f);
    colors[ImGuiCol_TableHeaderBg]          = ImVec4(0.19f, 0.19f, 0.20f, 1.00f);
    colors[ImGuiCol_TableBorderStrong]      = ImVec4(0.31f, 0.31f, 0.35f, 1.00f);
    colors[ImGuiCol_TableBorderLight]       = ImVec4(0.23f, 0.23f, 0.25f, 1.00f);
    colors[ImGuiCol_TableRowBg]             = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
    colors[ImGuiCol_TableRowBgAlt]          = ImVec4(1.00f, 1.00f, 1.00f, 0.06f);
    colors[ImGuiCol_TextLink]               = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
    colors[ImGuiCol_TextSelectedBg]         = ImVec4(0.26f, 0.59f, 0.98f, 0.35f);
    colors[ImGuiCol_TreeLines]              = ImVec4(0.43f, 0.43f, 0.50f, 0.50f);
    colors[ImGuiCol_DragDropTarget]         = ImVec4(1.00f, 1.00f, 0.00f, 0.90f);
    colors[ImGuiCol_DragDropTargetBg]       = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
    colors[ImGuiCol_UnsavedMarker]          = ImVec4(1.00f, 1.00f, 1.00f, 1.00f);
    colors[ImGuiCol_NavCursor]              = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
    colors[ImGuiCol_NavWindowingHighlight]  = ImVec4(1.00f, 1.00f, 1.00f, 0.70f);
    colors[ImGuiCol_NavWindowingDimBg]      = ImVec4(0.80f, 0.80f, 0.80f, 0.20f);
    colors[ImGuiCol_ModalWindowDimBg]       = ImVec4(0.80f, 0.80f, 0.80f, 0.35f);
}

inline void initUI(bool darkTheme)
{
    rlImGuiBeginInitImGui();
    ImGuiIO& io = ImGui::GetIO();

    if (darkTheme) ImGui::StyleColorsDark();
    else ImGui::StyleColorsLight();

    ImFontConfig config;
    config.FontDataOwnedByAtlas = false;
    config.OversampleH = 3;
    config.OversampleV = 1;
    config.PixelSnapH = true;

    g_FONTs.inter = io.Fonts->AddFontFromMemoryTTF(
        (void*)assets_fonts_Inter_opsz_wght__ttf,
        assets_fonts_Inter_opsz_wght__ttf_len,
        18.0f,
        &config
    );
    g_FONTs.bonanova = io.Fonts->AddFontFromMemoryTTF(
        (void*)assets_fonts_BonaNovaSC_Regular_ttf,
        assets_fonts_BonaNovaSC_Regular_ttf_len,
        18.0f,
        &config
    );
    g_FONTs.brawler = io.Fonts->AddFontFromMemoryTTF(
        (void*)assets_fonts_Brawler_Regular_ttf,
        assets_fonts_Brawler_Regular_ttf_len,
        18.0f,
        &config
    );
    g_FONTs.clearsans = io.Fonts->AddFontFromMemoryTTF(
        (void*)assets_fonts_ClearSans_Regular_ttf,
        assets_fonts_ClearSans_Regular_ttf_len,
        18.0f,
        &config
    );
    ImFont* customFont = g_FONTs.inter;
    if (customFont!=nullptr)
    {
        io.FontDefault = customFont;
        // rlImGuiReloadFonts();
        TraceLog(LOG_INFO,"UI: Modern font loaded and reloaded into GPU successfully!");
    }
    else
	{
        TraceLog(LOG_ERROR,"UI: Failed to parse font array from memory!");
	}
	ApplyCustomStyle();
    rlImGuiEndInitImGui();
}

inline void DrawFontMenu()
{
    ImGuiIO& io = ImGui::GetIO();
    if (ImGui::BeginMenu("Font"))
    {
        bool isInterActive = (io.FontDefault==g_FONTs.inter);
        if (ImGui::MenuItem("Inter Regular (Default)",nullptr,isInterActive)) io.FontDefault=g_FONTs.inter;
        bool isBonaNovaActive = (io.FontDefault==g_FONTs.bonanova);
        if (ImGui::MenuItem("BonaNovaSC Regular",nullptr,isBonaNovaActive)) io.FontDefault=g_FONTs.bonanova;
        bool isBrawlerActive = (io.FontDefault==g_FONTs.brawler);
        if (ImGui::MenuItem("Brawler Regular",nullptr,isBrawlerActive)) io.FontDefault=g_FONTs.brawler;
        bool isClearSansActive = (io.FontDefault==g_FONTs.clearsans);
        if (ImGui::MenuItem("ClearSans Regular",nullptr,isClearSansActive)) io.FontDefault=g_FONTs.clearsans;
        ImGui::EndMenu();
    }
}
