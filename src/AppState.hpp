#pragma once
#include "MM/models/kuramoto/general.hpp"
#include "MM/models/kuramoto/sparse.hpp"
#include "MM/models/kuramoto/special.hpp"
#include "MM/solvers/ODE/rk/explicit/rk1-solver.hpp"
#include "MM/solvers/ODE/rk/explicit/rk2-solver.hpp"
#include "MM/solvers/ODE/rk/explicit/rk3-solver.hpp"
#include "MM/solvers/ODE/rk/explicit/rk4-solver.hpp"
#include "MM/solvers/ODE/rk/explicit/rk4-variants.hpp"
#include "MM/typedefs/header.hpp"
#include "MM/initializers/initials.hpp"
#include "MM/network/topology.hpp"
#include "MM/models/kuramoto/general.hpp"
#include "MM/models/kuramoto/sparse.hpp"
#include "MM/models/kuramoto/special.hpp"
#include "imgui_internal.h"
#include "raylib.h"
#include <atomic>
#ifdef PI
#undef PI
#endif
#include "imgui.h"
#include "implot.h"
#include "UI.hpp"
// #include <climits>
#include <cstddef>
#include <cstring>
#include <string>
#include <mutex>
#include <atomic>
using SolverFunc = std::function<MathEngine::SolverResults(const MathEngine::SolverParameters&)>;
inline SolverFunc rk1_wrapper()
{
    return [](const MathEngine::SolverParameters& Params) {return MathEngine::rk1_solver(Params);};
}
inline SolverFunc rk2_wrapper()
{
    return [](const MathEngine::SolverParameters& Params) {return MathEngine::rk2_solver(Params);};
}
inline SolverFunc rk3_wrapper()
{
    return [](const MathEngine::SolverParameters& Params) {return MathEngine::rk3_solver(Params);};
}
inline SolverFunc rk4_wrapper()
{
    return [](const MathEngine::SolverParameters& Params) {return MathEngine::rk4_solver(Params);};
}
inline SolverFunc rk4_38_wrapper()
{
    return [](const MathEngine::SolverParameters& Params) {return MathEngine::rk4_38_solver(Params);};
}
inline SolverFunc rk4_ralston_wrapper()
{
    return [](const MathEngine::SolverParameters& Params) {return MathEngine::rk4_ralston_solver(Params);};
}
inline SolverFunc rk4_gill_wrapper()
{
    return [](const MathEngine::SolverParameters& Params) {return MathEngine::rk4_gill_solver(Params);};
}

enum class ModelType
{
	Kuramoto=0
};
enum class KuramotoType
{
    KuramotoGeneral=0,
    KuramotoSparse,
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
    Modules
};
enum class NetworkTopology
{
	RandomUniform=0,
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
enum class SolverMethod
{
    RK1=0,
    RK2,
    RK3,
    RK4,
    RK4_38,
    RK4_Gill,
    RK4_Ralston
};
struct GeneralModelParams
{
    ModelType modelType = ModelType::Kuramoto;
    KuramotoType kuramotoType = KuramotoType::KuramotoGeneral;
	size_t N = 50;
    size_t nModules = 1;
    size_t sModules = 50;
    double K = 1.0;
    double Q = 0.5;
    double alpha = 0.0;
    MathEngine::dVec iFrqnc;
    MathEngine::dVec iPhase;
    GeneralModelParams(size_t n=50) : N(n) {};
    int modelSelectedIndex = 0;
    int kuramotoModelSelectedIndex = 0;
};
struct DistParams
{
    double minVal = 0.0, maxVal = 1.0;
    double mean = 0.0, stddev = 1.0;
    double location = 0.0, scale = 1.0;
    double rate = 1.0;
    double perturbation = 0.0001;
    double param1 = 0.0, param2 = 1.0;
    int seed = 41;
    int moduleTypeIndex = 0;
    int typeIndex = 0;
    InitState initState = InitState::RandomUniform;
    bool showArray = false;
    bool identical = false;
};
struct NetParams
{
    double weightMin = 0.0, weightMax = 1.0, weight = 0.5, weightIn = 0.0, weightOut = 1.0;
    double prob = 0.5, probIn = 0.5, probOut = 0.5;
    double decayRatio = 0.1;
    int meanDegree = 2;
    int seed = 41;
    size_t sModulesBase = 10, nModulesBase = 2, hLevels = 2;
    size_t sModulesM = 10, nModulesM = 2;
    NetworkTopology adjState = NetworkTopology::ErdosRenyi;
    bool showAdjMatrix = false;
    int adjSelectedIndex = 3;
};
struct SolverParams
{
    SolverMethod solverMethod = SolverMethod::RK4;
    int nDs = 1;
    int solverMethodSelectedIndex = 4;
    SolverFunc solverFunc=nullptr;
    MathEngine::SolverParameters solverParams;
    MathEngine::SolverResults solverResults;
};
struct PlotParams
{
    MathEngine::OneStepSolverResult liveSolverRes;
	std::mutex plotMutex;
    int Stride = 50;
    int trailCount = 5000;
    bool showPlot = false;
    bool showPlotSecond = false;
    bool showPlotThird = false;
    size_t offset = 0;
    MathEngine::Vec<ImVec4> plotColors;
    MathEngine::Vec<ImVec4> plotSecondColors;
    MathEngine::Vec<ImVec4> plotThirdColors;
    MathEngine::dVec liveTimePoints = {};
    MathEngine::dVec liveState = {};
    MathEngine::dVec plotX = {};
    MathEngine::dVec plotY = {};
    MathEngine::dVec plotXTrail = {};
    MathEngine::dVec plotYTrail = {};
    MathEngine::dMatrix plotYModules;
};
////////////////////////////////////
/////                          /////
/////     STATE OF THE APP     /////
/////                          /////
////////////////////////////////////
class AppState
{
	public:
        std::atomic<float> simProgress{0.0f};
        std::atomic<double> timeInv{0.0f};
        std::atomic<bool> isSimRunning{false};
        float padding = 10.0f;
        Color BgColor = Color(15.0f,15.0f,15.0f);
        size_t initW = 800;
        size_t initH = 600;
        std::string appTitle;
        bool showStyleEditor=false;
        bool showDelays=false;
        bool showPlot=false;
        bool runSimulation=false;
        bool DarkTheme=true;
        GeneralModelParams modelParams = GeneralModelParams(50);
        DistParams phaseParams;
        DistParams frqncParams;
        NetParams adjParams;
        SolverParams solverParams;
        PlotParams plotParams;

        MathEngine::dMatrix adj; // Adjacency (for any system that might need it)
        MathEngine::SparsedMatrix sparseAdj = MathEngine::SparsedMatrix(modelParams.N); // Sparse adjacency
        MathEngine::dVec delayTimes = {0.0};
        inline void RenderUI()
        {
            drawTopMenuBar();
            if (showStyleEditor)
            {
                const ImGuiViewport* viewport = ImGui::GetMainViewport();
                ImVec2 tPos = ImVec2(viewport->WorkPos.x+padding,viewport->WorkPos.y+padding);
                ImGui::SetNextWindowPos(tPos,ImGuiCond_Appearing,ImVec2(0.0f,0.0f));
                ImGui::SetNextWindowSize(ImVec2(static_cast<size_t>(viewport->WorkSize.x-2*padding),static_cast<size_t>(viewport->WorkSize.y*0.45f)), ImGuiCond_Appearing);
                ImGui::Begin("ImGui Style Editor",&showStyleEditor);
                    ImGui::ShowStyleEditor();
                ImGui::End();
            }
            const ImGuiViewport* viewport = ImGui::GetMainViewport();
            ImGui::SetNextWindowPos(ImVec2(viewport->WorkPos.x + padding, viewport->WorkPos.y + padding), ImGuiCond_Always);
            ImGui::SetNextWindowSize(ImVec2(viewport->WorkSize.x * 0.48f, viewport->WorkSize.y - (2 * padding)), ImGuiCond_Always);
            if (ImGui::Begin("Control Centre",nullptr,ImGuiWindowFlags_NoCollapse))
            {
                if (ImGui::BeginTabBar("Contral Tabs",ImGuiTabBarFlags_None))
                {
                    if (ImGui::BeginTabItem("Model Details"))
                    {
                        DrawModelPanelContent();
                        DrawTopologyPanelContent();
                        ImGui::EndTabItem();
                    }
                    if (ImGui::BeginTabItem("Initial Condition"))
                    {
                        DrawInitialsPanelContent();
                        ImGui::EndTabItem();
                    }
                    if (ImGui::BeginTabItem("Solver Parameters"))
                    {
                        DrawSolverParametersPanelContent();
                        ImGui::EndTabItem();
                    }
                    if (ImGui::BeginTabItem("Plot Parameters"))
                    {
                        DrawPlotPanelContent();
                        ImGui::EndTabItem();
                    }
                    if (ImGui::BeginTabItem("Run Simulation"))
                    {
                        DrawProgressBar();
                        ImGui::Spacing();
                        bool running = isSimRunning.load();
                        if (running) ImGui::BeginDisabled();
                        if (ImGui::Button(running?"Simulating...":"Begin Simulation",ImVec2(-1,0))) StartSimulation();
                        if (running) ImGui::EndDisabled();
                        ImGui::EndTabItem();
                    }
                    ImGui::EndTabBar();
                }
            }
            ImGui::End();
            RenderModals();
            DrawPlotWindow();
        }
        inline void drawTopMenuBar()
        {
            if (ImGui::BeginMainMenuBar())
            {
                if (ImGui::BeginMenu("Options"))
                {
                    DrawFontMenu();
                    ImGui::Separator();
                    if (ImGui::MenuItem("Show ImGui Style Editor",nullptr,&showStyleEditor)) {}
                    ImGui::EndMenu();
                }
                ImGui::EndMainMenuBar();
            }
        }
	private:
        static constexpr const char* modelNames[] = {"Kuramoto"};
        static constexpr const char* kuramotoModelNames[] = {"Kuramoto (General)", "Kuramoto (Sparse)", "Kuramoto (Modular)"};
        static constexpr const char* adjNames[] = {"Random (Uniform)", "Random (Uniform Symmetric)", "Erdos-Renyi",
            "Erdos-Renyi (True Count)","Erdos-Renyi (Symmetric)", "Erdos-Renyi (Symmetric True Count)",
            "Small World", "Small World (Directed)", "Modular", "Hierarchical"};
        static constexpr const char* dsStateNames[] = {"Random Uniform", "Random Normal", "Random Cauchy",
            "Random Exponential", "Random Circle", "Splay", "Splay Perturbed", "Modules (by type)"};
        static constexpr const char* moduleTypeNames[] = {"Random Uniform", "Random Normal", "Random Cauchy",
            "Random Exponential", "Random Circle", "Splay", "Splay Perturbed"};
        static constexpr const char* moduleTypeIds[] = {"uniform", "normal", "cauchy", "exponential",
            "circle", "splay", "splay_perturbed"};
        static constexpr const char* solverMethodNames[] = {"Euler (RK1)", "Midpoint (RK2)", "Runge-Kutta 3rd Order (RK3)",
            "Runge-Kutta 4th Order (Standard RK4)", "Runge-Kutta 4th Order (3/8 RK4 variant)", "Runge-Kutta 4th Order (Gill's RK4 variant)",
            "Runge-Kutta 4th Order (Ralston's RK4 variant)"
        };
        inline void DrawModelPanelContent();
        inline void DrawTopologyPanelContent();
		inline void DrawInitialsPanelContent();
		inline void DrawSolverParametersPanelContent();
		inline void RenderModals();
        inline void DrawPlotWindow();
        inline void StartSimulation();
        inline void DrawProgressBar();
        inline void DrawPlotPanelContent();
};

inline void AppState::DrawModelPanelContent()
{
    ImGui::SeparatorText("Model Configuration");
    if (ImGui::Combo("Model Type",&modelParams.modelSelectedIndex, modelNames,1))
    {
        modelParams.modelType = static_cast<ModelType>(modelParams.modelSelectedIndex);
    }
    if (modelParams.modelType==static_cast<ModelType>(0))
    {
        if (ImGui::Combo("Kuramoto Type",&modelParams.kuramotoModelSelectedIndex,kuramotoModelNames,3))
        {
            modelParams.kuramotoType=static_cast<KuramotoType>(modelParams.kuramotoModelSelectedIndex);
        }
        ImGui::Spacing();

        int n = static_cast<int>(modelParams.N);
        int nM = static_cast<int>(modelParams.nModules);
        int sM = static_cast<int>(modelParams.sModules);

        switch (modelParams.kuramotoType)
        {
            case KuramotoType::KuramotoGeneral:
            case KuramotoType::KuramotoSparse:
                if (ImGui::InputInt("Oscillators (N)", &n, 1, 50)) modelParams.N = static_cast<size_t>(std::max(1, n));
                ImGui::InputDouble("Coupling (K)", &modelParams.K, 0.0001, 0.01, "%.15g");
                ImGui::InputDouble("Phase Lag (alpha)", &modelParams.alpha, MathEngine::PI * 0.001, MathEngine::PI * 0.01, "%.15g rad");
                break;
            case KuramotoType::KuramotoSpecial:
                if (ImGui::InputInt("Modules", &nM, 1, 5)) modelParams.nModules = static_cast<size_t>(std::max(1, nM));
                if (ImGui::InputInt("Module Size", &sM, 1, 50)) modelParams.sModules = static_cast<size_t>(std::max(1, sM));
                modelParams.N = modelParams.sModules * modelParams.nModules;
                ImGui::TextDisabled("Total Oscillators (N): %zu", modelParams.N);
                ImGui::InputDouble("K Intra", &modelParams.K, 0.0001, 0.01, "%.15g");
                ImGui::InputDouble("K Inter", &modelParams.Q, 0.0001, 0.01, "%.15g");
                break;
        }
    }
}

inline void AppState::DrawTopologyPanelContent()
{
    if (modelParams.kuramotoType == KuramotoType::KuramotoSpecial) return;
    ImGui::Spacing();
    ImGui::SeparatorText("Network Topology");
    if (ImGui::Combo("Topology Type", &adjParams.adjSelectedIndex, adjNames, 10))
    {
        adjParams.adjState = static_cast<NetworkTopology>(adjParams.adjSelectedIndex);
    }
    switch (adjParams.adjState) {
        case NetworkTopology::RandomUniform:
        case NetworkTopology::RandomUniformSymmetric:
            ImGui::InputDouble("Min Weight", &adjParams.weightMin, 0.0001, 0.01, "%.15g");
            ImGui::InputDouble("Max Weight", &adjParams.weightMax, 0.0001, 0.01, "%.15g");
            break;
        case NetworkTopology::ErdosRenyi:
        case NetworkTopology::ErdosRenyiUniform:
        case NetworkTopology::ErdosRenyiSymmetric:
        case NetworkTopology::ErdosRenyiSymmetricUniform:
            ImGui::InputDouble("Min Weight", &adjParams.weightMin, 0.0001, 0.01, "%.15g");
            ImGui::InputDouble("Max Weight", &adjParams.weightMax, 0.0001, 0.01, "%.15g");
            if (ImGui::InputDouble("Connection Prob", &adjParams.prob, 0.0001, 0.01, "%.15g"))
                adjParams.prob = std::clamp(adjParams.prob, 0.0, 1.0);
	        break;
        case NetworkTopology::SmallWorld:
        case NetworkTopology::SmallWorldDirected:
            ImGui::InputDouble("Weight", &adjParams.weight, 0.0001, 0.01, "%.15g");
            if (ImGui::InputDouble("Rewiring Prob", &adjParams.prob, 0.0001, 0.01, "%.15g"))
                adjParams.prob = std::clamp(adjParams.prob, 0.0, 1.0);
            if (ImGui::InputInt("Mean Degree", &adjParams.meanDegree, 1, 10))
                adjParams.meanDegree=std::max(0,adjParams.meanDegree);
            break;
        case NetworkTopology::Modular: {
            int sM = static_cast<int>(adjParams.sModulesM);
            int nM = static_cast<int>(adjParams.nModulesM);
            if (ImGui::InputInt("Module Size", &sM, 1, 10)) adjParams.sModulesM = std::max(1, sM);
            if (ImGui::InputInt("Number of Modules", &nM, 1, 10)) adjParams.nModulesM = std::max(1, nM);
            ImGui::InputDouble("In Weight", &adjParams.weightIn, 0.0001, 0.01, "%.15g");
            ImGui::InputDouble("Out Weight", &adjParams.weightOut, 0.0001, 0.01, "%.15g");
            if (ImGui::InputDouble("Inner Prob", &adjParams.probIn, 0.0001, 0.01, "%.15g"))
                adjParams.probIn = std::clamp(adjParams.probIn, 0.0, 1.0);
            if (ImGui::InputDouble("Outer Prob", &adjParams.probOut, 0.0001, 0.01, "%.15g"))
                adjParams.probOut = std::clamp(adjParams.probOut, 0.0, 1.0);
            break;
        }
        case NetworkTopology::Hierarchical: {
            int sB = static_cast<int>(adjParams.sModulesBase);
            int nB = static_cast<int>(adjParams.nModulesBase);
            int hL = static_cast<int>(adjParams.hLevels);
            if (ImGui::InputInt("Module Size##h", &sB, 1, 10)) adjParams.sModulesBase = std::max(1, sB);
            if (ImGui::InputInt("Number of Modules##h", &nB, 1, 10)) adjParams.nModulesBase = std::max(1, nB);
            if (ImGui::InputInt("Hierarchy Levels##h", &hL, 1, 10)) adjParams.hLevels = std::max(1, hL);
            ImGui::InputDouble("In Weight##h", &adjParams.weightIn, 0.0001, 0.01, "%.15g");
            ImGui::InputDouble("Out Weight##h", &adjParams.weightOut, 0.0001, 0.01, "%.15g");
            if (ImGui::InputDouble("Inner Prob##h", &adjParams.probIn, 0.0001, 0.01, "%.15g"))
                adjParams.probIn = std::clamp(adjParams.probIn, 0.0, 1.0);
            if (ImGui::InputDouble("Outer Prob##h", &adjParams.probOut, 0.0001, 0.01, "%.15g"))
                adjParams.probOut = std::clamp(adjParams.probOut, 0.0, 1.0);
            if (ImGui::InputDouble("Decay Ratio", &adjParams.decayRatio, 0.001, 0.1, "%.15g"))
                adjParams.decayRatio = std::clamp(adjParams.decayRatio, 0.0, 1.0);
            break;
        }
    }
    ImGui::InputInt("Seed##Adj", &adjParams.seed, 1, 10);
    if (ImGui::Button("Generate Adjacency", ImVec2(-1, 0)))
    {
        size_t seedVal = static_cast<size_t>(std::max(1, adjParams.seed));
        switch (adjParams.adjState)
        {
            case NetworkTopology::RandomUniform:
                adj = MathEngine::random(modelParams.N,adjParams.weightMin,adjParams.weightMax,seedVal);
                break;
            case NetworkTopology::RandomUniformSymmetric:
                adj = MathEngine::random_symmetric(modelParams.N,adjParams.weightMin,adjParams.weightMax,seedVal);
                break;
            case NetworkTopology::ErdosRenyi:
                adj = MathEngine::erdos_renyi(modelParams.N,adjParams.prob,adjParams.weightMin,adjParams.weightMax,seedVal);
                break;
            case NetworkTopology::ErdosRenyiUniform:
                adj = MathEngine::erdos_renyi_uniform(modelParams.N,adjParams.prob,adjParams.weightMin,adjParams.weightMax,seedVal);
                break;
            case NetworkTopology::ErdosRenyiSymmetric:
                adj = MathEngine::erdos_renyi_symmetric(modelParams.N,adjParams.prob,adjParams.weightMin,adjParams.weightMax,seedVal);
                break;
            case NetworkTopology::ErdosRenyiSymmetricUniform:
                adj = MathEngine::erdos_renyi_symmetric_uniform(modelParams.N,adjParams.prob,adjParams.weightMin,adjParams.weightMax,seedVal);
                break;
            case NetworkTopology::SmallWorld:
                adj = MathEngine::small_world(modelParams.N,adjParams.meanDegree,adjParams.prob,adjParams.weight,seedVal);
                break;
            case NetworkTopology::SmallWorldDirected:
                adj = MathEngine::small_world_directed(modelParams.N,adjParams.meanDegree,adjParams.prob,adjParams.weight,seedVal);
                break;
            case NetworkTopology::Modular:
                modelParams.sModules = adjParams.sModulesM;
                modelParams.nModules = adjParams.nModulesM;
                modelParams.N = modelParams.sModules * modelParams.nModules;
                adj = MathEngine::modular(modelParams.sModules,modelParams.nModules,adjParams.probIn,adjParams.probOut,
                                          adjParams.weightIn,adjParams.weightOut,
                                          seedVal);
                break;
            case NetworkTopology::Hierarchical:
                modelParams.sModules = adjParams.sModulesBase;
                modelParams.nModules = adjParams.nModulesBase * static_cast<size_t>(std::pow(2, adjParams.hLevels - 1));
                modelParams.N = modelParams.sModules * modelParams.nModules;
                adj = MathEngine::hierarchical(adjParams.sModulesBase,adjParams.hLevels,adjParams.probIn, adjParams.probOut,
                                               adjParams.weightIn,adjParams.weightOut, adjParams.decayRatio, seedVal,
                                               adjParams.nModulesBase);
                break;
        }

        if (modelParams.kuramotoType == KuramotoType::KuramotoSparse)
        {
            sparseAdj = MathEngine::dense_to_sparse(adj);
        }
    }

    if (ImGui::Button("View Matrix Values", ImVec2(-1, 0)))
    {
        adjParams.showAdjMatrix = true;
    }
}

inline void AppState::DrawInitialsPanelContent()
{
    ImGui::SeparatorText("Initial Phases");
    if (ImGui::Combo("Phase Dist", &phaseParams.typeIndex, dsStateNames, 8))
        phaseParams.initState = static_cast<InitState>(phaseParams.typeIndex);

    switch (phaseParams.initState)
    {
        case InitState::RandomUniform:
            ImGui::InputDouble("Min##P", &phaseParams.minVal, 0.0001, 0.01, "%.15g");
            if (ImGui::InputDouble("Max##P", &phaseParams.maxVal, 0.0001, 0.01, "%.15g"))
                phaseParams.maxVal=std::max(phaseParams.minVal,phaseParams.maxVal);
            break;
        case InitState::RandomNormal:
            ImGui::InputDouble("Mean##P", &phaseParams.mean, 0.0001, 0.01, "%.15g");
            if (ImGui::InputDouble("Stddev##P", &phaseParams.stddev, 0.0001, 0.01, "%.15g"))
                phaseParams.stddev=std::max(1e-5,phaseParams.stddev);
            break;
        case InitState::RandomCauchy:
            ImGui::InputDouble("Location##P", &phaseParams.location, 0.0001, 0.01, "%.15g");
            if (ImGui::InputDouble("Scale##P", &phaseParams.scale, 0.0001, 0.01, "%.15g"))
                phaseParams.scale=std::max(1e-5,phaseParams.scale);
            break;
        case InitState::RandomExponential:
            if (ImGui::InputDouble("Rate##P", &phaseParams.rate, 0.0001, 0.01, "%.15g"))
                phaseParams.rate=std::max(1e-5,phaseParams.rate);
            break;
        case InitState::SplayPerturbed:
            ImGui::InputDouble("Perturbation##P", &phaseParams.perturbation, 0.000001, 0.0001, "%.15g");
            break;
        case InitState::Modules:
            ImGui::Combo("Module State##P", &phaseParams.moduleTypeIndex, moduleTypeNames, 7);
            ImGui::InputDouble("Param 1##P", &phaseParams.param1, 0.0001, 0.01, "%.15g");
            ImGui::InputDouble("Param 2##P", &phaseParams.param2, 0.0001, 0.01, "%.15g");
            ImGui::Checkbox("Identical Modules##P", &phaseParams.identical);
            break;
        default: break;
    }
    ImGui::InputInt("Seed##P", &phaseParams.seed, 1, 10);

    if (ImGui::Button("Generate Initial Phases", ImVec2(-1, 0)))
    {
        size_t pSeed = static_cast<size_t>(std::max(1, phaseParams.seed));
        switch (phaseParams.initState)
        {
            case InitState::RandomUniform:
                modelParams.iPhase = MathEngine::random_uniform(modelParams.N, phaseParams.minVal, phaseParams.maxVal, pSeed);
                break;
            case InitState::RandomNormal:
                modelParams.iPhase = MathEngine::random_normal(modelParams.N, phaseParams.mean, phaseParams.stddev, pSeed);
                break;
            case InitState::RandomCauchy:
                modelParams.iPhase = MathEngine::random_cauchy(modelParams.N, phaseParams.location, phaseParams.scale, pSeed);
                break;
            case InitState::RandomExponential:
                modelParams.iPhase = MathEngine::random_exponential(modelParams.N, phaseParams.rate, pSeed);
                break;
            case InitState::RandomCircle:
                modelParams.iPhase = MathEngine::random_circle(modelParams.N, pSeed);
                break;
            case InitState::Splay:
                modelParams.iPhase = MathEngine::splay(modelParams.N);
                break;
            case InitState::SplayPerturbed:
                modelParams.iPhase = MathEngine::splay_perturbed(modelParams.N, phaseParams.perturbation, pSeed);
                break;
            case InitState::Modules: {
                std::string dist_type = moduleTypeIds[phaseParams.moduleTypeIndex];
                modelParams.iPhase = MathEngine::modules(modelParams.sModules, modelParams.nModules, dist_type, phaseParams.param1, phaseParams.param2, pSeed, phaseParams.identical);
                break;
            }
        }
        solverParams.solverParams.initialConditions = modelParams.iPhase;
    }
    if (ImGui::Button("View Phase Array", ImVec2(-1, 0))) phaseParams.showArray = true;
	// ***************************************************************************** //
    ImGui::Spacing();
    ImGui::SeparatorText("Intrinsic Frequencies");
    if (ImGui::Combo("Freq Dist", &frqncParams.typeIndex, dsStateNames, 8))
        frqncParams.initState = static_cast<InitState>(frqncParams.typeIndex);
    switch (frqncParams.initState)
    {
        case InitState::RandomUniform:
            ImGui::InputDouble("Min##F", &frqncParams.minVal, 0.0001, 0.01, "%.15g");
            if (ImGui::InputDouble("Max##F", &frqncParams.maxVal, 0.0001, 0.01, "%.15g"))
                frqncParams.maxVal=std::max(frqncParams.minVal,frqncParams.maxVal);
            break;
        case InitState::RandomNormal:
            ImGui::InputDouble("Mean##F", &frqncParams.mean, 0.0001, 0.01, "%.15g");
            if (ImGui::InputDouble("Stddev##F", &frqncParams.stddev, 0.0001, 0.01, "%.15g"))
                frqncParams.stddev=std::max(1e-5,frqncParams.stddev);
        	break;
        case InitState::RandomCauchy:
            ImGui::InputDouble("Location##F", &frqncParams.location, 0.0001, 0.01, "%.15g");
            if (ImGui::InputDouble("Scale##F", &frqncParams.scale, 0.0001, 0.01, "%.15g"))
                frqncParams.scale=std::max(1e-5,frqncParams.scale);
        	break;
        case InitState::RandomExponential:
            if (ImGui::InputDouble("Rate##F", &frqncParams.rate, 0.0001, 0.01, "%.15g"))
                frqncParams.rate=std::max(1e-5,frqncParams.rate);
        	break;
        case InitState::SplayPerturbed:
            ImGui::InputDouble("Perturbation##F", &frqncParams.perturbation, 0.000001, 0.0001, "%.15g");
            break;
        case InitState::Modules:
            ImGui::Combo("Module State##F", &frqncParams.moduleTypeIndex, moduleTypeNames, 7);
            ImGui::InputDouble("Param 1##F", &frqncParams.param1, 0.0001, 0.01, "%.15g");
            ImGui::InputDouble("Param 2##F", &frqncParams.param2, 0.0001, 0.01, "%.15g");
            ImGui::Checkbox("Identical Modules##F", &frqncParams.identical);
            break;
        default: break;
    }
    ImGui::InputInt("Seed##F", &frqncParams.seed, 1, 10);
    if (ImGui::Button("Generate Frequencies", ImVec2(-1, 0)))
    {
        size_t fSeed = static_cast<size_t>(std::max(1, frqncParams.seed));
        switch (frqncParams.initState)
        {
            case InitState::RandomUniform:
                modelParams.iFrqnc = MathEngine::random_uniform(modelParams.N, frqncParams.minVal, frqncParams.maxVal, fSeed);
                break;
            case InitState::RandomNormal:
                modelParams.iFrqnc = MathEngine::random_normal(modelParams.N, frqncParams.mean, frqncParams.stddev, fSeed);
                break;
            case InitState::RandomCauchy:
                modelParams.iFrqnc = MathEngine::random_cauchy(modelParams.N, frqncParams.location, frqncParams.scale, fSeed);
                break;
            case InitState::RandomExponential:
                modelParams.iFrqnc = MathEngine::random_exponential(modelParams.N, frqncParams.rate, fSeed);
                break;
            case InitState::RandomCircle:
                modelParams.iFrqnc = MathEngine::random_circle(modelParams.N, fSeed);
                break;
            case InitState::Splay:
                modelParams.iFrqnc = MathEngine::splay(modelParams.N);
                break;
            case InitState::SplayPerturbed:
                modelParams.iFrqnc = MathEngine::splay_perturbed(modelParams.N, frqncParams.perturbation, fSeed);
                break;
            case InitState::Modules: {
                std::string dist_type = moduleTypeIds[frqncParams.moduleTypeIndex];
                modelParams.iFrqnc = MathEngine::modules(modelParams.sModules, modelParams.nModules, dist_type, frqncParams.param1, frqncParams.param2, fSeed, frqncParams.identical);
                break;
            }
        }
    }
    if (ImGui::Button("View Frequency Array", ImVec2(-1, 0))) frqncParams.showArray = true;

    ImGui::Spacing();
    ImGui::Separator();

    // MODEL DERIVATIVE WRAPPER COMPILATION
    if (ImGui::Button("Compile Model Function", ImVec2(-1, 35)))
    {
        if (modelParams.modelType==ModelType::Kuramoto)
        {
            if (modelParams.kuramotoType == KuramotoType::KuramotoGeneral)
            {
                MathEngine::KuramotoParams kParams;
                kParams.K = modelParams.K;
                kParams.N = modelParams.N;
                kParams.alpha = modelParams.alpha;
                kParams.omega = modelParams.iFrqnc;
                kParams.adj = adj;
                solverParams.solverParams.derivative = MathEngine::kuramoto_general_wrapper(kParams);
            }
            else if (modelParams.kuramotoType == KuramotoType::KuramotoSparse)
            {
                MathEngine::KuramotoSparseParams kParams;
                kParams.K = modelParams.K;
                kParams.N = modelParams.N;
                kParams.alpha = modelParams.alpha;
                kParams.omega = modelParams.iFrqnc;
                kParams.sparse_adj = sparseAdj;
                solverParams.solverParams.derivative = MathEngine::kuramoto_sparse_wrapper(kParams);
            }
            else if (modelParams.kuramotoType == KuramotoType::KuramotoSpecial)
            {
                MathEngine::KuramotoModularParams kParams;
                kParams.intra_K = modelParams.K;
                kParams.inter_K = modelParams.Q;
                kParams.N = modelParams.N;
                kParams.alpha = modelParams.alpha;
                kParams.omega = modelParams.iFrqnc;
                solverParams.solverParams.derivative = MathEngine::kuramoto_special_modular_wrapper(kParams);
            }
        }
    }
}

inline void AppState::DrawSolverParametersPanelContent()
{
	if (ImGui::Combo("Solver Method",&solverParams.solverMethodSelectedIndex,solverMethodNames,7))
        solverParams.solverMethod=static_cast<SolverMethod>(solverParams.solverMethodSelectedIndex);
    switch(solverParams.solverMethod)
    {
        case SolverMethod::RK1:
            solverParams.solverFunc = rk1_wrapper();
            break;
        case SolverMethod::RK2:
            solverParams.solverFunc = rk2_wrapper();
        case SolverMethod::RK3:
	        solverParams.solverFunc = rk3_wrapper();
        case SolverMethod::RK4:
            solverParams.solverFunc = rk4_wrapper();
        case SolverMethod::RK4_38:
            solverParams.solverFunc = rk4_38_wrapper();
        case SolverMethod::RK4_Gill:
            solverParams.solverFunc = rk4_gill_wrapper();
        case SolverMethod::RK4_Ralston:
            solverParams.solverFunc = rk4_ralston_wrapper();
        default:
            solverParams.solverFunc = rk1_wrapper();
            break;
    }
    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Time & Basic Stepping", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputDouble("Start Time (t0)", &solverParams.solverParams.t0, 0.0001, 0.1, "%.15g");
        ImGui::InputDouble("End Time (t1)", &solverParams.solverParams.t1, 0.0001, 0.1, "%.15g");
        // Initial step size dt
        ImGui::InputDouble("Step Size (dt)", &solverParams.solverParams.dt, 0.000001, 0.01, "%.15g");
        ImGui::Separator();
        // Multi-step method controls (Adams-Bashforth / Adams-Moulton)
        ImGui::SliderInt("Method Order", &solverParams.solverParams.order, 1, 10);
        ImGui::SliderInt("ABM Iterations", &solverParams.solverParams.iterations, 1, 10);
    }
    if (ImGui::CollapsingHeader("Adaptive Step Control"))
    {
        ImGui::Checkbox("Estimate Error", &solverParams.solverParams.errorEstimate);
        ImGui::SameLine();
        ImGui::Checkbox("Enable Variable Step Size", &solverParams.solverParams.variableSteps);
        if (solverParams.solverParams.errorEstimate && solverParams.solverParams.variableSteps)
        {
            ImGui::Indent();
            ImGui::TextDisabled("Tolerances & Bounds");
            ImGui::InputDouble("Local Tolerance", &solverParams.solverParams.localTol, 0.0, 0.0, "%.1e");
            ImGui::InputDouble("Absolute Tolerance", &solverParams.solverParams.absolute_tol, 0.0, 0.0, "%.1e");
            ImGui::InputDouble("Min dt", &solverParams.solverParams.minDt, 0.0, 0.0, "%.15g");
            ImGui::InputDouble("Max dt", &solverParams.solverParams.maxDt, 0.0, 0.0, "%.15g");
            ImGui::Separator();
            ImGui::TextDisabled("Step Adaptation Factors");
            ImGui::InputDouble("Decrease Factor", &solverParams.solverParams.decreaseFactor, 0.05, 0.1, "%.15g");
            ImGui::InputDouble("Increase Factor", &solverParams.solverParams.increaseFactor, 0.1, 0.5, "%.15g");
            ImGui::InputDouble("Tol Error Ratio", &solverParams.solverParams.localTolErrorRatio, 0.01, 0.05, "%.15g");
            // size_t cast to int for ImGui input
            int maxTrial = static_cast<int>(solverParams.solverParams.maxTrial);
            if (ImGui::InputInt("Max Trials", &maxTrial))
                solverParams.solverParams.maxTrial = static_cast<size_t>(maxTrial>1?maxTrial:1);
            ImGui::Unindent();
        }
    }
    if (ImGui::CollapsingHeader("Error Metrics & Flags"))
    {
        ImGui::Checkbox("Weighted Error Formula", &solverParams.solverParams.weightedError);
        ImGui::Checkbox("Norm Error Formula", &solverParams.solverParams.normError);
        ImGui::Checkbox("Record Attempt History", &solverParams.solverParams.attemptsHistory);
    }
    // double tau=0.0;
    if (ImGui::CollapsingHeader("Delay Differential Equations (DDE)"))
    {
        ImGui::Checkbox("Is DDE System", &solverParams.solverParams.isDDE);
        if (solverParams.solverParams.isDDE)
        {
            ImGui::Indent();
            int maxDelayOrder_ = static_cast<int>(solverParams.solverParams.maxDelayOrder);
            if (ImGui::InputInt("Max Delay Order", &maxDelayOrder_))
            {
                solverParams.solverParams.maxDelayOrder = static_cast<size_t>(maxDelayOrder_ > 1 ? maxDelayOrder_ : 1);
            }
            ImGui::InputDouble("Interpolation Tol", &solverParams.solverParams.interpolationTol, 1e-10, 1e-8, "%.1e");
            ImGui::InputDouble("dt Scale (Fine Step)", &solverParams.solverParams.dtScale, 0.01, 0.05, "%.15g");
            ImGui::SliderInt("Number of Delays",&solverParams.nDs,1,20);
            ImGui::SameLine();
            if (ImGui::Button("Set Delay Count"))
            {
                delayTimes.resize(solverParams.nDs);
            }
            if (ImGui::CollapsingHeader("Delays"))
            {
                // delayTimes = MathEngine::dVec(nDs,0.0);
                for (size_t i=0; i<delayTimes.size(); ++i)
                {
                    std::string label = "Delay number "+std::to_string(i+1);
                    ImGui::InputDouble(label.c_str(), &delayTimes[i],0.0001f,0.01f,"%.15g");
                    // delayTimes[i]=tau;
                }
            }
            if (ImGui::Button("Submit Delays"))
            {
                solverParams.solverParams.delayTimes=delayTimes;
            }
            ImGui::SameLine();
            if (ImGui::Button("View Delay Values"))
            {
                showDelays=true;
            }
            ImGui::Text("Configured Delays: %zu", solverParams.solverParams.delayTimes.size());
            ImGui::Unindent();
        }
    }
}

inline void AppState::RenderModals()
{
    const ImGuiViewport* viewport = ImGui::GetMainViewport();
    ImVec2 center = viewport->GetCenter();
    // Matrix View Popup
    if (adjParams.showAdjMatrix)
    {
        ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
        ImGui::SetNextWindowSize(ImVec2(550, 480), ImGuiCond_Appearing);
        if (ImGui::Begin("Adjacency Matrix Viewer", &adjParams.showAdjMatrix))
        {
            size_t nRows = adj.size();
            size_t nCols = nRows > 0 ? adj[0].size() : 0;
            ImGui::Text("Dimension: %zu x %zu", nRows, nCols);
            ImGui::Separator();

            if (nRows == 0 || nCols == 0)
            {
                ImGui::TextDisabled("Matrix is Empty. Generate an Adjacency Matrix first.");
            }
            else
            {
                ImGuiTableFlags tableFlags = ImGuiTableFlags_ScrollX | ImGuiTableFlags_ScrollY | ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingFixedFit;
                if (ImGui::BeginTable("MatrixGrid", static_cast<int>(nCols + 1), tableFlags, ImVec2(0, 360)))
                {
                    ImGui::TableSetupScrollFreeze(1, 1);
                    ImGui::TableSetupColumn("Row\\Col", ImGuiTableColumnFlags_NoHide);
                    for (size_t c = 0; c < nCols; ++c)
                    {
                        char colHeader[16];
                        snprintf(colHeader, sizeof(colHeader), "[%03zu]", c+1);
                        ImGui::TableSetupColumn(colHeader);
                    }
                    ImGui::TableHeadersRow();

                    for (size_t r = 0; r < nRows; ++r)
                    {
                        ImGui::TableNextRow();
                        ImGui::TableSetColumnIndex(0);
                        ImGui::TextDisabled("[%03zu]", r+1);
                        for (size_t c = 0; c < nCols; ++c)
                        {
                            ImGui::TableSetColumnIndex(static_cast<int>(c + 1));
                            double val = adj[r][c];
                            if (val >= 1e-5) ImGui::TextColored(ImVec4(0.4f, 0.9f, 1.0f, 1.0f), "%.15g", val);
                            else ImGui::TextDisabled("0.0000");
                        }
                    }
                    ImGui::EndTable();
                }
            }
        }
        ImGui::End();
    }

    if (phaseParams.showArray)
    {
        ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
        ImGui::SetNextWindowSize(ImVec2(400, 400), ImGuiCond_Appearing);
        if (ImGui::Begin("Phase Array Values", &phaseParams.showArray))
        {
            if (ImGui::BeginChild("PhaseList", ImVec2(0, 360), ImGuiChildFlags_Borders))
            {
                if (modelParams.iPhase.empty())
                {
                    ImGui::TextDisabled("Array is empty. Click Generate Initial Phases.");
                }
                else
                {
                    for (size_t i = 0; i < modelParams.iPhase.size(); ++i)
                    {
                        ImGui::Text("[%03zu]  %.15g", i + 1, modelParams.iPhase[i]);
                    }
                }
                ImGui::EndChild();
            }
        }
        ImGui::End();
    }

    if (frqncParams.showArray)
    {
        ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
        ImGui::SetNextWindowSize(ImVec2(400, 400), ImGuiCond_Appearing);
        if (ImGui::Begin("Frequency Array Values", &frqncParams.showArray))
        {
            if (ImGui::BeginChild("FreqList", ImVec2(0, 360), ImGuiChildFlags_Borders))
            {
                if (modelParams.iFrqnc.empty())
                {
                    ImGui::TextDisabled("Array is empty. Click Generate Frequencies.");
                }
                else
                {
                    for (size_t i = 0; i < modelParams.iFrqnc.size(); ++i)
                    {
                        ImGui::Text("[%03zu]  %.15g", i + 1, modelParams.iFrqnc[i]);
                    }
                }
                ImGui::EndChild();
            }
        }
        ImGui::End();
    }

    if (showDelays)
    {
        ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
        ImGui::SetNextWindowSize(ImVec2(500, 480), ImGuiCond_Appearing);
        if (ImGui::Begin("Submitted Delays",&showDelays))
        {
            ImGui::Text("Total Size: %zu elements", solverParams.solverParams.delayTimes.size());
            ImGui::Separator();
            ImGui::Spacing();

            // Scrollable child box for array elements
            if (ImGui::BeginChild("ArrayList", ImVec2(0, 360), ImGuiChildFlags_Borders))
            {
                if (solverParams.solverParams.delayTimes.empty())
                {
                    ImGui::TextDisabled("Array is empty.");
                }
                else
                {
                    for (size_t i = 0; i < solverParams.solverParams.delayTimes.size(); ++i)
                    {
                        ImGui::Text("[%03zu]  %.6f", i+1, solverParams.solverParams.delayTimes[i]);
                    }
                }
            }
            ImGui::EndChild();

            ImGui::Spacing();
            ImGui::Separator();
            ImGui::Spacing();
        }
        ImGui::End();
    }
}

inline void AppState::StartSimulation()
{
    simProgress.store(0.0f);
    isSimRunning.store(true);
    timeInv.store(static_cast<float>(1.0/std::abs(solverParams.solverParams.t1-solverParams.solverParams.t0)));
    // constexpr size_t Stride = 25;
    const size_t vectorSize = static_cast<size_t>((solverParams.solverParams.t1-solverParams.solverParams.t0)/solverParams.solverParams.dt);
    const size_t plotExpectedSize = static_cast<size_t>(vectorSize/plotParams.Stride)*2+100;
    {
		std::lock_guard<std::mutex> lock(plotParams.plotMutex);
        plotParams.liveTimePoints.clear();
        plotParams.liveState.clear();
        plotParams.plotX.clear();
        plotParams.plotY.clear();
        plotParams.plotX.reserve(plotExpectedSize);
        plotParams.plotY.reserve(plotExpectedSize);
        plotParams.plotXTrail.clear();
        plotParams.plotYTrail.clear();
        plotParams.plotYModules.clear();
        plotParams.plotYModules.clear();
        plotParams.plotYModules = MathEngine::dMatrix(modelParams.nModules,MathEngine::dVec());
        plotParams.offset = 0;
    }
    int stride = plotParams.Stride;
    bool condPlotThird = (modelParams.kuramotoType==KuramotoType::KuramotoSpecial || adjParams.adjState==NetworkTopology::Modular ||
    					adjParams.adjState==NetworkTopology::Hierarchical);
    solverParams.solverParams.onStep = [this, stepCount=0, stepCountCond=0, stride, condPlotThird](const MathEngine::OneStepSolverResult& res) mutable
    {
		float progress = (res.timePoint-solverParams.solverParams.t0)*timeInv;
        std::lock_guard<std::mutex> lock(plotParams.plotMutex);
        simProgress.store(progress);
        double rSine = 0.0, rCosine = 0.0, rho = 0.0;
        MathEngine::dVec rMSine(modelParams.nModules,0.0);
        MathEngine::dVec rMCosine(modelParams.nModules,0.0);
        MathEngine::dVec rhoM(modelParams.nModules,0.0);
        plotParams.liveTimePoints.push_back(res.timePoint);
        plotParams.liveState = res.sol;
        if (condPlotThird)
        {
            for (size_t i=0; i<modelParams.nModules; ++i)
            {
                for (size_t j=0; j<modelParams.sModules; ++j)
                {
                    rMSine[i] += sin(res.sol[i*modelParams.sModules+j]);
                    rMCosine[i] += cos(res.sol[i*modelParams.sModules+j]);
                }
                rSine += rMSine[i]; rCosine += rMCosine[i]; rMSine[i] /=modelParams.sModules; rMCosine[i] /= modelParams.sModules;
                rhoM[i] = sqrt(rMSine[i]*rMSine[i]+rMCosine[i]*rMCosine[i]);
            }
            rSine /= modelParams.N; rCosine /= modelParams.N;
            rho = sqrt(rSine*rSine+rCosine*rCosine);
            if (++stepCountCond%stride==0)
            {
                for (size_t i=0; i<modelParams.nModules; ++i)
                {
                    plotParams.plotYModules[i].push_back(rhoM[i]);
                }
            }
        }
        else
        {
            for (size_t i=0; i<res.sol.size(); ++i)
            {
                rSine += sin(res.sol[i]); rCosine += cos(res.sol[i]);
            }
            rSine /= modelParams.N; rCosine /= modelParams.N;
            rho = sqrt(rSine*rSine+rCosine*rCosine);
        }
        if (plotParams.plotXTrail.size()<plotParams.trailCount)
        {
            plotParams.plotXTrail.push_back(res.timePoint);
            plotParams.plotYTrail.push_back(rho);
        }
        else
        {
            plotParams.plotXTrail[plotParams.offset] = res.timePoint;
            plotParams.plotYTrail[plotParams.offset] = rho;
            plotParams.offset = static_cast<size_t>((plotParams.offset+1) % plotParams.trailCount);
        }
        if (++stepCount%stride==0)
        {
            plotParams.plotX.push_back(res.timePoint);
            plotParams.plotY.push_back(rho);
        }
    };
    std::thread([this]()
    {
        solverParams.solverResults = solverParams.solverFunc(solverParams.solverParams);
        simProgress.store(1.0);
        isSimRunning.store(false);
    }).detach();
}

inline void AppState::DrawProgressBar()
{
    float progress = simProgress.load();
    bool running = isSimRunning.load();
    char overlayBuf[64];
    if (running)
    {
        snprintf(overlayBuf,sizeof(overlayBuf),"Running... %.1f%%",progress*100.0f);
    }
    else if (progress>=1.0f)
    {
        snprintf(overlayBuf,sizeof(overlayBuf),"Completed 100%%");
    }
    else
    {
        snprintf(overlayBuf,sizeof(overlayBuf),"Idle 0.0%%");
    }
    ImGui::PushStyleVar(ImGuiStyleVar_FrameRounding,6.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize,1.0f);
    ImGui::PushStyleColor(ImGuiCol_PlotHistogram,ImVec4(0.15f,0.9f,0.6f,1.0f));
    ImGui::PushStyleColor(ImGuiCol_FrameBg,ImVec4(0.05f,0.3f,0.2f,1.0f));
    ImGui::PushStyleColor(ImGuiCol_Border,ImVec4(0.25f,0.25,0.3f,1.0f));
    ImGui::ProgressBar(progress,ImVec2(-1.0f,25.0f),overlayBuf);
    ImGui::PopStyleColor(3);
    ImGui::PopStyleVar(2);
}

inline void AppState::DrawPlotWindow()
{
    std::lock_guard<std::mutex> lock(plotParams.plotMutex);
    if (modelParams.modelType==ModelType::Kuramoto)
    {
        if (plotParams.showPlot)
        {
            const ImGuiViewport* viewport = ImGui::GetMainViewport();
            ImGui::SetNextWindowPos(ImVec2(viewport->WorkPos.x + 0.52*viewport->WorkSize.x - padding, viewport->WorkPos.y + padding),ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(viewport->WorkSize.x * 0.48f, 0.75*viewport->WorkSize.y - (2 * padding)), ImGuiCond_FirstUseEver);
            if (ImGui::Begin("Order Parameter",&plotParams.showPlot))
            {
                ImVec2 availableSpace = ImGui::GetContentRegionAvail();      // Set ImVec(-1,-1) to fill the whole window.
                if (ImPlot::BeginSubplots("##Plot-now", 2, 1, availableSpace))
                {
                    if (ImPlot::BeginPlot("(\U0001D70C-t) Plot"))
                    {
                        ImPlot::SetupAxes("Time (t)","Order (\U0001D70C)");
                        ImPlot::PlotLine("\U0001D70C",plotParams.plotX.data(),plotParams.plotY.data(),static_cast<int>(plotParams.plotX.size()));
                        // ImPlot::PlotScatter("Points 1",runtimeX.data(),runtimeY.data(),static_cast<int>(runtimeX.size()));
                        ImPlot::EndPlot();
                    }
                    if (ImPlot::BeginPlot("(\U0001D70C-t) Plot##trailing"))
                    {
                        ImPlotSpec spec;
                        spec.Offset = static_cast<int>(plotParams.offset);
                        ImPlot::SetupAxes("Time (t)","Order (\U0001D70C)");
                        ImPlot::PlotLine("\U0001D70C",plotParams.plotXTrail.data(),plotParams.plotYTrail.data(),static_cast<int>(plotParams.plotXTrail.size()),spec);
                        ImPlot::EndPlot();
                    }
                    ImPlot::EndSubplots();
                }
            }
            ImGui::End();
        }
        if (plotParams.showPlotSecond)
        {
            const ImGuiViewport* viewport = ImGui::GetMainViewport();
            ImGui::SetNextWindowPos(ImVec2(viewport->WorkPos.x+padding,viewport->WorkPos.y+0.25*viewport->WorkSize.y+padding),ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(viewport->WorkSize.x*0.48f, 0.75*viewport->WorkSize.y-(2*padding)), ImGuiCond_FirstUseEver);
            if (ImGui::Begin("Oscillators",&plotParams.showPlotSecond))
            {
                size_t plotN = plotParams.liveState.size();
                ImVec2 availableSpace = ImGui::GetContentRegionAvail();      // Set ImVec(-1,-1) to fill the whole window.
                MathEngine::dVec xTheta(plotN);
                MathEngine::dVec yTheta(plotN);
                for (size_t i=0; i<plotN; ++i)
                {
                    xTheta[i] = cos(plotParams.liveState[i]);
                    yTheta[i] = sin(plotParams.liveState[i]);
                }
                if (ImPlot::BeginPlot("\U0001D73D Plot",availableSpace))
                {
                    ImPlot::SetupAxes("x projection","y projection");
                    // ImPlot::PlotLine("\U0001D70C",plotParams.plotX.data(),plotParams.plotY.data(),static_cast<int>(plotParams.plotX.size()));
                    ImPlot::PlotScatter("\U0001D73D",xTheta.data(),yTheta.data(),static_cast<int>(plotN));
                    ImPlot::EndPlot();
                }
            }
            ImGui::End();
        }
        if (plotParams.showPlotThird)
        {
            const ImGuiViewport* viewport = ImGui::GetMainViewport();
            ImVec2 center = viewport->GetCenter();
            ImGui::SetNextWindowPos(center, ImGuiCond_FirstUseEver, ImVec2(0.5f, 0.5f));
            ImGui::SetNextWindowSize(ImVec2(400, 200), ImGuiCond_FirstUseEver);
            if (ImGui::Begin("Order Parameter (Modules)##mocules",&plotParams.showPlotThird))
            {
                ImVec2 availableSpace = ImGui::GetContentRegionAvail();      // Set ImVec(-1,-1) to fill the whole window.
                if (plotParams.plotYModules.empty() || plotParams.plotYModules.size()!=modelParams.nModules)
                {
                    if (ImPlot::BeginPlot("(\U0001D70C-t) Plot##modules",availableSpace))
                    {
                        ImPlot::EndPlot();
                    }
                }
                else
                {
                    if (ImPlot::BeginPlot("(\U0001D70C-t) Plot##modules",availableSpace))
                    {
                        ImPlot::SetupAxes("Time (t)","Order (\U0001D70C)##modules");
                        std::string label = "\U0001D70C";
                        for (size_t i=0; i<modelParams.nModules; ++i)
                        {
                            label = "\U0001D70C "+std::to_string(i);
                            ImPlot::PlotLine(label.c_str(),plotParams.plotX.data(),plotParams.plotYModules[i].data(),
                                            static_cast<int>(plotParams.plotX.size()));
                        }
                        ImPlot::EndPlot();
                    }
                }
            }
            ImGui::End();
        }
    }
}

inline void AppState::DrawPlotPanelContent()
{
    if (modelParams.modelType==ModelType::Kuramoto)
    {
        ImGui::SeparatorText("Plot Data Style");
        if (ImGui::CollapsingHeader("\U0001D70C-t Plot##main plot"))
        {
            if (ImGui::InputInt("Stride##main plot",&plotParams.Stride,1,10)) plotParams.Stride = std::max(plotParams.Stride,10);
            if (ImGui::InputInt("Trailing Data Count##main plot",&plotParams.trailCount,1,10)) plotParams.trailCount = std::clamp(plotParams.trailCount,100,10000);
            ImGui::Checkbox("Show \U0001D70C-t Plot##main plot", &plotParams.showPlot);
            ImGui::Spacing();
            ImGui::SeparatorText("Line Color##main plot");
            ImGui::Spacing();
            float colorR=0.2f, colorG=0.5f, colorB=0.5f, colorA=1.0f;
            ImGui::SliderFloat("R##main plot line", &colorR,0.0f,1.0f,"%.2f");
            ImGui::SliderFloat("G##main plot line", &colorG,0.0f,1.0f,"%.2f");
            ImGui::SliderFloat("B##main plot line", &colorB,0.0f,1.0f,"%.2f");
            ImGui::SliderFloat("A (opacity)##main plot line", &colorA,0.0f,1.0f,"%.2f");
            if (ImGui::Button("Set Color(s)##main plot",ImVec2(-1,0))) plotParams.plotColors.push_back(ImVec4(colorR,colorG,colorB,colorA));
        }
        if (ImGui::CollapsingHeader("\U0001D73D Plot##second plot"))
        {
            ImGui::Checkbox("Show \U0001D73D Plot##second plot", &plotParams.showPlotSecond);
            ImGui::Spacing();
            ImGui::SeparatorText("Line Color##second plot");
            ImGui::Spacing();
            float colorR=0.2f, colorG=0.5f, colorB=0.5f, colorA=1.0f;
            ImGui::SliderFloat("R##second plot line", &colorR,0.0f,1.0f,"%.2f");
            ImGui::SliderFloat("G##second plot line", &colorG,0.0f,1.0f,"%.2f");
            ImGui::SliderFloat("B##second plot line", &colorB,0.0f,1.0f,"%.2f");
            ImGui::SliderFloat("A (opacity)##second plot line", &colorA,0.0f,1.0f,"%.2f");
            if (ImGui::Button("Set Color(s)##second plot",ImVec2(-1,0)))
                plotParams.plotSecondColors.push_back(ImVec4(colorR,colorG,colorB,colorA));
        }
        if (adjParams.adjState==NetworkTopology::Modular || adjParams.adjState==NetworkTopology::Hierarchical || modelParams.kuramotoType==KuramotoType::KuramotoSpecial)
        {
            if (ImGui::CollapsingHeader("\U0001D73D Plot (Modules)##third plot"))
            {
                ImGui::Checkbox("Show \U0001D73D Plot##third plot", &plotParams.showPlotThird);
                ImGui::Spacing();
                ImGui::SeparatorText("Line Color##third plot");
                ImGui::Spacing();
                size_t nM = modelParams.nModules;
                MathEngine::Vec<float> colorAs(nM);MathEngine::Vec<float> colorRs(nM);
                MathEngine::Vec<float> colorGs(nM);MathEngine::Vec<float> colorBs(nM);
                for (size_t i=0; i<nM; ++i)
                {
                    std::string rLabel = "R ("+std::to_string(i+1)+")##third plot line";
                    std::string gLabel = "G ("+std::to_string(i+1)+")##third plot line";
                    std::string bLabel = "B ("+std::to_string(i+1)+")##third plot line";
                    std::string aLabel = "A ("+std::to_string(i+1)+")##third plot line";
                    ImGui::SliderFloat(rLabel.c_str(), &colorRs[i],0.0f,1.0f,"%.2f");
                    ImGui::SliderFloat(gLabel.c_str(), &colorGs[i],0.0f,1.0f,"%.2f");
                    ImGui::SliderFloat(bLabel.c_str(), &colorBs[i],0.0f,1.0f,"%.2f");
                    ImGui::SliderFloat(aLabel.c_str(), &colorAs[i],0.0f,1.0f,"%.2f");
					ImGui::Separator();
                }
                if (ImGui::Button("Set Color(s)##third plot",ImVec2(-1,0)))
                {
                    for (size_t i=0; i<nM; ++i)
                    {
                        plotParams.plotThirdColors.push_back(ImVec4(colorRs[i],colorGs[i],colorBs[i],colorAs[i]));
                    }
                }
            }
        }
    }
}
