#pragma once
#include "MM/initializers/initials.hpp"
#include "MM/typedefs/header.hpp"
#include "MM/models/kuramoto.hpp"
#include "raylib.h"
#include "imgui.h"
#include "UI.hpp"
#include <climits>
#include <cstddef>
#include <cstring>
#include <string>

enum class ModelType
{
    // KuramotoGeneralParallel,
    // KuramotoSparseParallel,
    KuramotoGeneral=0,
    KuramotoSparce,
    KuramotoSpecial
};
enum class InitState
{
	RandomUniform=0,
    RandomNormal,
    RandomCauchy,
    RandomExponential,
    RandomCircle,
    Splay,
    SplayPerturbed,
    IdenticalModule
};
enum class NetworkTopology
{
	RandomUniform,
    RandomUniformSymmetric,
    ErdosRenyi,
    ErdosRenyiUniform,
    ErdosRenyiSymmetric,
    ErdosRenyiSymmetricUniform,
    SmallWorld,
    SmallWorldDirected,
    Modular,
    Hierarchical
};

struct GeneralModelParams
{
    AppFonts g_FONTs;
	size_t N = 50;
    size_t nModules = 2;
    size_t sModules = 25;
    double K = 1.0;
    double Q = 0.5;
    double alpha = 0.0;
    MathEngine::dVec iFrqnc;
    MathEngine::dVec iPhase;
    GeneralModelParams(size_t n=50) : N(n) {};
};
struct AppState
{
    double a=0.0,a1=0.0,a2=0.0,a3=1.0,a4=0.0,a5=0.0,b=1.0,b1=1.0,b2=1.0,b3=1.0;
    size_t d=41;
    double oa=0.0,oa1=0.0,oa2=0.0,oa3=1.0,oa4=0.0,oa5=0.0,ob=1.0,ob1=1.0,ob2=1.0,ob3=1.0;
    size_t od=41;
    ModelType modelType;
    // Panel Pos and Size
    ImVec2 StyleEditorPos=ImVec2(0.0f,0.0f);
    ImVec2 StyleEditorSize=ImVec2(0.0f,0.0f);
    ImVec2 SolverParametersPos=ImVec2(0.0f,0.0f);
    ImVec2 SolverParametersSize=ImVec2(0.0f,0.0f);
    ImVec2 ModelPickerPos=ImVec2(0.0f,0.0f);
    ImVec2 ModelPickerSize=ImVec2(0.0f,0.0f);
    ImVec2 InitialsPos=ImVec2(0.0f,0.0f);
    ImVec2 InitialsSize=ImVec2(0.0f,0.0f);
    InitState dsInitState;
    InitState otherInitState; // Used for things such as intrinsic frequencies and alike
    Color BgColor = Color(15.0f,15.0f,15.0f);
    size_t initW = 800;
    size_t initH = 600;
    std::string appTitle;
    bool showStyleEditor = false;
    bool showPhaseArray = false;
    bool showFrqncArray = false;
    bool DarkTheme = true;
    MathEngine::SolverParameters params;
    GeneralModelParams modelParams = GeneralModelParams(50);
    // MathEngine::dVec dsState=MathEngine::dVec(0.0,10); // Dynamical System State
    // MathEngine::dVec otherdsState;
    MathEngine::dMatrix adj; // Adjacency (for any system that might need it)
    MathEngine::SparsedMatrix sparseAdj = MathEngine::SparsedMatrix(modelParams.N); // Sparse adjacency
    inline void RenderUI()
    {
        drawTopMenuBar();

        if (showStyleEditor)
        {
            const ImGuiViewport* viewport = ImGui::GetMainViewport();

            float padding = 10.0f;
            ImVec2 tPos = ImVec2(viewport->WorkPos.x+padding,viewport->WorkPos.y+padding);
            ImGui::SetNextWindowPos(tPos,ImGuiCond_Always,ImVec2(0.0f,0.0f));
            ImGui::SetNextWindowSize(ImVec2(static_cast<size_t>(viewport->WorkSize.x-2*padding),static_cast<size_t>(viewport->WorkSize.y*0.45f)), ImGuiCond_Always);
            ImGui::Begin("ImGui Style Editor",&showStyleEditor);
	            ImGui::ShowStyleEditor();
                StyleEditorPos = ImGui::GetWindowPos();
                StyleEditorSize = ImGui::GetWindowSize();
            ImGui::End();
        }
        DrawSolverParametersPanel(params);
        DrawModelPanel();
        DrawInitials();
    }
    inline void drawTopMenuBar()
    {
        if (ImGui::BeginMainMenuBar())
        {
            if (ImGui::BeginMenu("Options"))
            {
                DrawFontMenu();

                ImGui::Separator();
                if (ImGui::MenuItem("Show Style Editor",nullptr,&showStyleEditor)) {}
                ImGui::EndMenu();
            }
            ImGui::EndMainMenuBar();
        }
    }
    inline void DrawSolverParametersPanel(MathEngine::SolverParameters& params);
	inline void DrawModelPanel();
    inline void DrawInitials();
};

inline void AppState::DrawSolverParametersPanel(MathEngine::SolverParameters& params)
{
    const ImGuiViewport* viewport = ImGui::GetMainViewport();
    float padding = 10.0f;
    ImVec2 tPos = ImVec2(viewport->WorkPos.x+padding,viewport->WorkPos.y+padding);
    if (showStyleEditor) tPos = ImVec2(StyleEditorPos.x,StyleEditorPos.y+StyleEditorSize.y+padding);
    ImGui::SetNextWindowPos(tPos,ImGuiCond_Always,ImVec2(0.0f,0.0f));
    ImGui::SetNextWindowSize(ImVec2(static_cast<size_t>(viewport->WorkSize.x*0.3f),static_cast<size_t>(viewport->WorkSize.y*0.45f)), ImGuiCond_Always);
    ImGui::Begin("Solver Parameters");
    if (ImGui::CollapsingHeader("Time & Basic Stepping", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputDouble("Start Time (t0)", &params.t0, 0.0001, 0.1, "%.4f");
        ImGui::InputDouble("End Time (t1)", &params.t1, 0.0001, 0.1, "%.4f");
        // Initial step size dt
        ImGui::InputDouble("Step Size (dt)", &params.dt, 0.000001, 0.01, "%.6f");
        ImGui::Separator();
        // Multi-step method controls (Adams-Bashforth / Adams-Moulton)
        ImGui::SliderInt("Method Order", &params.order, 1, 10);
        ImGui::SliderInt("ABM Iterations", &params.iterations, 1, 10);
    }
    if (ImGui::CollapsingHeader("Adaptive Step Control"))
    {
        ImGui::Checkbox("Estimate Error", &params.error_estimate);
        ImGui::SameLine();
        ImGui::Checkbox("Enable Variable Step Size", &params.variable_steps);
        if (params.error_estimate && params.variable_steps)
        {
            ImGui::Indent();
            ImGui::TextDisabled("Tolerances & Bounds");
            ImGui::InputDouble("Local Tolerance", &params.local_tol, 0.0, 0.0, "%.1e");
            ImGui::InputDouble("Absolute Tolerance", &params.absolute_tol, 0.0, 0.0, "%.1e");
            ImGui::InputDouble("Min dt", &params.min_dt, 0.0, 0.0, "%.6f");
            ImGui::InputDouble("Max dt", &params.max_dt, 0.0, 0.0, "%.4f");
            ImGui::Separator();
            ImGui::TextDisabled("Step Adaptation Factors");
            ImGui::InputDouble("Decrease Factor", &params.decrease_factor, 0.05, 0.1, "%.2f");
            ImGui::InputDouble("Increase Factor", &params.increase_factor, 0.1, 0.5, "%.2f");
            ImGui::InputDouble("Tol Error Ratio", &params.local_tol_error_ratio, 0.01, 0.05, "%.3f");
            // size_t cast to int for ImGui input
            int maxTrial = static_cast<int>(params.max_trial);
            if (ImGui::InputInt("Max Trials", &maxTrial))
                params.max_trial = static_cast<size_t>(maxTrial > 1 ? maxTrial : 1);
            ImGui::Unindent();
        }
    }
    if (ImGui::CollapsingHeader("Error Metrics & Flags"))
    {
        ImGui::Checkbox("Weighted Error Formula", &params.weighted_error);
        ImGui::Checkbox("Norm Error Formula", &params.norm_error);
        ImGui::Checkbox("Record Attempt History", &params.attempts_history);
    }
    if (ImGui::CollapsingHeader("Delay Differential Equations (DDE)"))
    {
        ImGui::Checkbox("Is DDE System", &params.is_dde);
        if (params.is_dde)
        {
            ImGui::Indent();
            int maxDelayOrder = static_cast<int>(params.max_delay_order);
            if (ImGui::InputInt("Max Delay Order", &maxDelayOrder))
            {
                params.max_delay_order = static_cast<size_t>(maxDelayOrder > 1 ? maxDelayOrder : 1);
            }
            ImGui::InputDouble("Interpolation Tol", &params.interpolation_tol, 0.0, 0.0, "%.1e");
            ImGui::InputDouble("dt Scale (Fine Step)", &params.dt_scale, 0.01, 0.05, "%.4f");
            ImGui::Text("Configured Delays: %zu", params.delay_times.size());
            ImGui::Unindent();
        }
    }
    SolverParametersPos = ImGui::GetWindowPos();
    SolverParametersSize = ImGui::GetWindowSize();
    ImGui::End();
}

inline void AppState::DrawModelPanel()
{
    const ImGuiViewport* viewport = ImGui::GetMainViewport();
    float padding = 10.0f;
    ImVec2 tPos = ImVec2(viewport->WorkPos.x+padding,viewport->WorkPos.y+padding);
    tPos = ImVec2(SolverParametersPos.x,SolverParametersPos.y+SolverParametersSize.y+padding);
    ImGui::SetNextWindowPos(tPos,ImGuiCond_Always,ImVec2(0.0f,0.0f));
    ImGui::SetNextWindowSize(ImVec2(static_cast<size_t>(viewport->WorkSize.x*0.3f),static_cast<size_t>(viewport->WorkSize.y*0.45f)), ImGuiCond_Always);
    if (ImGui::Begin("Model Selection"))
    {
        const char* modelNames[] = {"Kuramoto (General)","Kuramoto (Sparse Edition)","Kuramoto (Special Modular)"};
        int selectedIndex = static_cast<int>(modelType);
        if (ImGui::Combo("Model", &selectedIndex, modelNames, 3))
        {
            modelType = static_cast<ModelType>(selectedIndex);
        }
        ImGui::Separator();
        ImGui::Spacing();
        int n = 50;
        int nM = 2;
        int sM = 25;
        switch (modelType)
        {
            case ModelType::KuramotoGeneral:
                ImGui::TextDisabled("General Kuramoto Model Parameters");
                n = static_cast<int>(modelParams.N);
                if (ImGui::InputInt("Oscillators (N)", &n,1,50)) modelParams.N = static_cast<size_t>(n>1?n:1);
                ImGui::InputDouble("Coupling (K)", &modelParams.K, 0.0001, 0.01, "%.4f");
                ImGui::InputDouble("Phase Lag (alpha)", &modelParams.alpha, PI*0.001, PI*0.01, "%.6f rad");
                break;
            case ModelType::KuramotoSparce:
                ImGui::TextDisabled("Sparse Edition Parameters");
                n = static_cast<int>(modelParams.N);
                if (ImGui::InputInt("Oscillators (N)", &n,1,50)) modelParams.N = static_cast<size_t>(n>1?n:1);
                ImGui::InputDouble("Coupling (K)", &modelParams.K, 0.0001, 0.01, "%.4f");
                ImGui::InputDouble("Phase Lag (alpha)", &modelParams.alpha, PI*0.001, PI*0.01, "%.6f rad");
                break;
            case ModelType::KuramotoSpecial:
                ImGui::TextDisabled("Modular Model Parameters");
                nM = static_cast<size_t>(modelParams.nModules);
                if (ImGui::InputInt("Modules", &nM,1,5)) modelParams.nModules = static_cast<size_t>(nM>1?nM:1);
                sM = static_cast<size_t>(modelParams.sModules);
                if (ImGui::InputInt("Module Size", &sM,1,50)) modelParams.sModules = static_cast<size_t>(sM>1?sM:1);
                modelParams.N = static_cast<size_t>(sM*nM);
                ImGui::InputDouble("K Intra", &modelParams.K, 0.0001, 0.01, "%.4f");
                ImGui::InputDouble("K Inter", &modelParams.Q, 0.0001, 0.01, "%.4f");
                break;
        }
    }
    ModelPickerPos = ImGui::GetWindowPos();
    ModelPickerSize = ImGui::GetWindowSize();
    ImGui::End();
}

inline void AppState::DrawInitials()
{
    const ImGuiViewport* viewport = ImGui::GetMainViewport();
    float padding = 10.0f;
    ImVec2 tPos = ImVec2(viewport->WorkPos.x+padding,viewport->WorkPos.y+padding);
    tPos = ImVec2(ModelPickerPos.x,ModelPickerPos.y+ModelPickerSize.y+padding);
    ImGui::SetNextWindowPos(tPos,ImGuiCond_Always,ImVec2(0.0f,0.0f));
    ImGui::SetNextWindowSize(ImVec2(static_cast<size_t>(viewport->WorkSize.x*0.3f),static_cast<size_t>(viewport->WorkSize.y*0.45f)), ImGuiCond_Always);
    if (ImGui::Begin("Initial State"))
    {
        const char* dsStateNames[] = {"Random Uniform","Random Normal","Random Cauchy","Random Exponential","Random Circle",
            "Splay","Splay (Perturbed Randomly)","Identical Modules"
        };
        const char* moduleTypeNames[] = {"Random Uniform","Random Normal","Random Cauchy",
            "Random Exponential","Random Circle","Splay","Splay (Perturbed Randomly)"};
        const char* moduleTypeIds[] = {"uniform","normal","cauchy","exponential","circle","splay","splay_perturbed"};
        int selectedIndex = static_cast<int>(dsInitState);
        ImGui::Text("Initial Condition");
        if (modelType==static_cast<ModelType>(0) || modelType==static_cast<ModelType>(1) || modelType==static_cast<ModelType>(2))
        {
            ImGui::SeparatorText("Initial Phase");
            if (ImGui::Combo("Phase Distribution", &selectedIndex, dsStateNames,8))
            {
                dsInitState = static_cast<InitState>(selectedIndex);
            }
            int sEED = 41;
            std::string dist_type = "uniform";
            InitState subInitState=InitState::RandomUniform;
            int moduleTypeIndex = static_cast<int>(subInitState);
            bool shouldGenerate = false;
            switch (dsInitState)
            {
                case InitState::RandomUniform:
                    ImGui::InputDouble("Min Value",&a,0.0001f,0.01f,"%4f");
                    ImGui::InputDouble("Max Value",&b,0.0001f,0.01f,"%4f");
                    break;
                case InitState::RandomNormal:
                    ImGui::InputDouble("Mean Value",&a1,0.0001f,0.01f,"%4f");
                    ImGui::InputDouble("Stddev Value",&b1,0.0001f,0.01f,"%4f");
                    break;
                case InitState::RandomCauchy:
                    ImGui::InputDouble("Location",&a2,0.0001f,0.01f,"%4f");
                    ImGui::InputDouble("Scale",&b2,0.0001f,0.01f,"%4f");
                    break;
                case InitState::RandomExponential:
                    ImGui::InputDouble("Rate",&a3,0.0001f,0.01f,"%4f");
                    break;
                case InitState::RandomCircle:
                    break;
                case InitState::Splay:
                    break;
                case InitState::SplayPerturbed:
                    ImGui::InputDouble("Perturbation Amplitude",&a4,0.000001f,0.0001f,"%6f");
                    break;
                case InitState::IdenticalModule:
                    if (ImGui::Combo("Module State",&moduleTypeIndex,moduleTypeNames,7))
                    {
                        dist_type = static_cast<std::string>(moduleTypeIds[moduleTypeIndex]);
                    }
                    ImGui::InputDouble("Parameter 1",&a5,0.0001f,0.01f,"%4f");
                    ImGui::InputDouble("Parameter 2",&b3,0.0001f,0.01f,"%4f");
                    break;
            }
            sEED = static_cast<int>(d);
            if (ImGui::InputInt("Seed",&sEED,1,10)) d=static_cast<size_t>(sEED>1 ? sEED : 1);

            ImGui::Spacing();
            if(ImGui::Button("Generate Initial Phases"))
            {
                shouldGenerate=true;
            }
            if (shouldGenerate)
            {
                switch (dsInitState)
                {
                    case InitState::RandomUniform:
                        modelParams.iPhase=MathEngine::random_uniform(modelParams.N,a,b,d);
                        break;
                    case InitState::RandomNormal:
                        modelParams.iPhase=MathEngine::random_normal(modelParams.N,a1,b1,d);
                        break;
                    case InitState::RandomCauchy:
                        modelParams.iPhase=MathEngine::random_cauchy(modelParams.N,a2,b2,d);
                        break;
                    case InitState::RandomExponential:
                        modelParams.iPhase=MathEngine::random_exponential(modelParams.N,a3,d);
                        break;
                    case InitState::RandomCircle:
                        modelParams.iPhase=MathEngine::random_circle(modelParams.N,d);
                        break;
                    case InitState::Splay:
                        modelParams.iPhase=MathEngine::splay(modelParams.N);
                        break;
                    case InitState::SplayPerturbed:
                        modelParams.iPhase=MathEngine::splay_perturbed(modelParams.N,a4,d);
                        break;
                    case InitState::IdenticalModule:
                        modelParams.iPhase=MathEngine::identical_modules(modelParams.sModules,modelParams.nModules,dist_type,a5,b3,d);
                        break;
                }
            }
            // VIEW VALUES
            ImGui::SameLine();
            if (ImGui::Button("View Phase Array Values"))
            {
                showPhaseArray = true;
            }
            if (showPhaseArray)
            {
                ImVec2 center = viewport->GetCenter();
                ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
                ImGui::SetNextWindowSize(ImVec2(500, 480), ImGuiCond_Appearing);
                if (ImGui::Begin("Generated Values",&showPhaseArray))
                {
                    ImGui::Text("Total Size: %zu elements", modelParams.iPhase.size());
                    ImGui::Separator();
                    ImGui::Spacing();

                    // Scrollable child box for array elements
                    if (ImGui::BeginChild("ArrayList", ImVec2(0, 240), ImGuiChildFlags_Borders))
                    {
                        if (modelParams.iPhase.empty())
                        {
                            ImGui::TextDisabled("Array is empty.");
                        }
                        else
                        {
                            for (size_t i = 0; i < modelParams.iPhase.size(); ++i)
                            {
                                ImGui::Text("[%03zu]  %.6f", i, modelParams.iPhase[i]);
                            }
                        }
                    }
                    ImGui::EndChild();

                    ImGui::Spacing();
                    ImGui::Separator();
                    ImGui::Spacing();

                    // Close Button
                    if (ImGui::Button("Close", ImVec2(-1, 0)))
                    { // -1 fills full width
						showPhaseArray = false;
                    }
                    ImGui::End();
                }
            }
            // FREQUENCIES
            ImGui::SeparatorText("Intrinsic Frequency");
            ImGui::Spacing();
            int otherSelectedIndex = static_cast<int>(otherInitState);
            if (ImGui::Combo("Frequency Distribution", &otherSelectedIndex, dsStateNames,8))
            {
                otherInitState = static_cast<InitState>(otherSelectedIndex);
            }
            int osEED=41;
            std::string other_dist_type = "uniform";
            InitState otherSubInitState=InitState::RandomUniform;
            int otherModuleTypeIndex = static_cast<int>(otherSubInitState);
            bool shouldGenerateOther=false;
            switch (otherInitState)
            {
                case InitState::RandomUniform:
                    ImGui::InputDouble("Min Value *",&oa,0.0001f,0.01f,"%4f");
                    ImGui::InputDouble("Max Value *",&ob,0.0001f,0.01f,"%4f");
                    break;
                case InitState::RandomNormal:
                    ImGui::InputDouble("Mean Value *",&oa1,0.0001f,0.01f,"%4f");
                    ImGui::InputDouble("Stddev Value *",&ob1,0.0001f,0.01f,"%4f");
                    break;
                case InitState::RandomCauchy:
                    ImGui::InputDouble("Location *",&oa2,0.0001f,0.01f,"%4f");
                    ImGui::InputDouble("Scale *",&ob2,0.0001f,0.01f,"%4f");
                    break;
                case InitState::RandomExponential:
                    ImGui::InputDouble("Rate *",&oa3,0.0001f,0.01f,"%4f");
                    break;
                case InitState::RandomCircle:
                    break;
                case InitState::Splay:
                    break;
                case InitState::SplayPerturbed:
                    ImGui::InputDouble("Perturbation Amplitude *",&oa4,0.000001f,0.0001f,"%6f");
                    break;
                case InitState::IdenticalModule:
                    if (ImGui::Combo("Module State",&otherModuleTypeIndex,moduleTypeNames,7))
                    {
                        other_dist_type = static_cast<std::string>(moduleTypeIds[otherModuleTypeIndex]);
                    }
                    ImGui::InputDouble("Parameter 1 *",&oa5,0.0001f,0.01f,"%4f");
                    ImGui::InputDouble("Parameter 2 *",&ob3,0.0001f,0.01f,"%4f");
                    break;
            }
            osEED = static_cast<int>(od);
            if (ImGui::InputInt("Seed *",&osEED,1,10)) od=static_cast<size_t>(osEED>1 ? osEED : 1);

            ImGui::Spacing();
            if(ImGui::Button("Generate Intrinsic Frequencies"))
            {
                shouldGenerateOther=true;
            }
            if (shouldGenerateOther)
            {
                switch (otherInitState)
                {
                    case InitState::RandomUniform:
                        modelParams.iFrqnc=MathEngine::random_uniform(modelParams.N,oa,ob,od);
                        break;
                    case InitState::RandomNormal:
                        modelParams.iFrqnc=MathEngine::random_normal(modelParams.N,oa1,ob1,od);
                        break;
                    case InitState::RandomCauchy:
                        modelParams.iFrqnc=MathEngine::random_cauchy(modelParams.N,oa2,ob2,od);
                        break;
                    case InitState::RandomExponential:
                        modelParams.iFrqnc=MathEngine::random_exponential(modelParams.N,oa3,od);
                        break;
                    case InitState::RandomCircle:
                        modelParams.iFrqnc=MathEngine::random_circle(modelParams.N,od);
                        break;
                    case InitState::Splay:
                        modelParams.iFrqnc=MathEngine::splay(modelParams.N);
                        break;
                    case InitState::SplayPerturbed:
                        modelParams.iFrqnc=MathEngine::splay_perturbed(modelParams.N,oa4,od);
                        break;
                    case InitState::IdenticalModule:
                        modelParams.iFrqnc=MathEngine::identical_modules(modelParams.sModules,modelParams.nModules,dist_type,oa5,ob3,od);
                        break;
                }
            }
            // VIEW VALUES
            ImGui::SameLine();
            if (ImGui::Button("View Intrinsic Frequency Array Values"))
            {
                showFrqncArray = true;
            }
            if (showFrqncArray)
            {
                ImVec2 center = viewport->GetCenter();
                ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
                ImGui::SetNextWindowSize(ImVec2(500, 480), ImGuiCond_Appearing);
                if (ImGui::Begin("Generated Values *",&showFrqncArray))
                {
                    ImGui::Text("Total Size: %zu elements", modelParams.iFrqnc.size());
                    ImGui::Separator();
                    ImGui::Spacing();

                    // Scrollable child box for array elements
                    if (ImGui::BeginChild("ArrayList *", ImVec2(0, 240), ImGuiChildFlags_Borders))
                    {
                        if (modelParams.iFrqnc.empty())
                        {
                            ImGui::TextDisabled("Array is empty.");
                        }
                        else
                        {
                            for (size_t i = 0; i < modelParams.iFrqnc.size(); ++i)
                            {
                                ImGui::Text("[%03zu]  %.6f", i, modelParams.iFrqnc[i]);
                            }
                        }
                    }
                    ImGui::EndChild();

                    ImGui::Spacing();
                    ImGui::Separator();
                    ImGui::Spacing();

                    // Close Button
                    if (ImGui::Button("Close", ImVec2(-1, 0)))
                    { // -1 fills full width
                        showFrqncArray = false;
                    }
                    ImGui::End();
                }
            }
    	}
    }
    InitialsPos = ImGui::GetWindowPos();
    InitialsSize = ImGui::GetWindowSize();
    ImGui::End();
}
