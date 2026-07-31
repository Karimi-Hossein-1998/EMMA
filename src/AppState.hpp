#pragma once
#include "MM/models/kuramoto/general.hpp"
#include "MM/models/kuramoto/sparse.hpp"
#include "MM/models/kuramoto/special.hpp"
#include "MM/typedefs/header.hpp"
#include "MM/initializers/initials.hpp"
#include "MM/network/topology.hpp"
#include "MM/models/kuramoto.hpp"
#include "raylib.h"
#ifdef PI
#undef PI
#endif
#include "imgui.h"
#include "UI.hpp"
#include <climits>
#include <cstddef>
#include <cstring>
#include <string>

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
struct GeneralModelParams
{
    ModelType modelType = ModelType::Kuramoto;
    KuramotoType kuramotoType = KuramotoType::KuramotoGeneral;
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
};
struct WinPoses
{
    ImVec2 StyleEditorPos=ImVec2(0.0f,0.0f);
    ImVec2 StyleEditorSize=ImVec2(0.0f,0.0f);
};

class AppState
{
	public:
        float padding = 10.0f;
        Color BgColor = Color(15.0f,15.0f,15.0f);
        size_t initW = 800;
        size_t initH = 600;
        int nDs = 0;
        std::string appTitle;
        bool showStyleEditor=false;
        bool showDelays=false;
        bool DarkTheme=true;
        MathEngine::SolverParameters solverParams;
        GeneralModelParams modelParams = GeneralModelParams(50);
        DistParams phaseParams;
        DistParams frqncParams;
        NetParams adjParams;
        WinPoses wPoses;

        // MathEngine::dVec dsState=MathEngine::dVec(0.0,10); // Dynamical System State
        // MathEngine::dVec otherdsState;
        MathEngine::dMatrix adj; // Adjacency (for any system that might need it)
        MathEngine::SparsedMatrix sparseAdj = MathEngine::SparsedMatrix(modelParams.N); // Sparse adjacency
        MathEngine::dVec intermediaryState = {};
        MathEngine::dVec delayTimes = {0.0};
        inline void RenderUI()
        {
            drawTopMenuBar();
            if (showStyleEditor)
            {
                const ImGuiViewport* viewport = ImGui::GetMainViewport();
                ImVec2 tPos = ImVec2(viewport->WorkPos.x+padding,viewport->WorkPos.y+padding);
                ImGui::SetNextWindowPos(tPos,ImGuiCond_Always,ImVec2(0.0f,0.0f));
                ImGui::SetNextWindowSize(ImVec2(static_cast<size_t>(viewport->WorkSize.x-2*padding),static_cast<size_t>(viewport->WorkSize.y*0.45f)), ImGuiCond_Always);
                ImGui::Begin("ImGui Style Editor",&showStyleEditor);
                    ImGui::ShowStyleEditor();
                    wPoses.StyleEditorPos = ImGui::GetWindowPos();
                    wPoses.StyleEditorSize = ImGui::GetWindowSize();
                ImGui::End();
            }
            const ImGuiViewport* viewport = ImGui::GetMainViewport();
            ImGui::SetNextWindowPos(ImVec2(viewport->WorkPos.x + padding, viewport->WorkPos.y + padding), ImGuiCond_Always);
            ImGui::SetNextWindowSize(ImVec2(viewport->WorkSize.x * 0.35f, viewport->WorkSize.y - (2 * padding)), ImGuiCond_Always);
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
                    ImGui::EndTabBar();
                }
            }
            ImGui::End();
            RenderModals();
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
	private:
        static constexpr const char* modelNames[] = {"Kuramoto"};
        static constexpr const char* kuramotoModelNames[] = {"Kuramoto (General)", "Kuramoto (Sparse)", "Kuramoto (Modular)"};
        static constexpr const char* adjNames[] = {"Random (Uniform)", "Random (Uniform Symmetric)", "Erdos-Renyi",
            "Erdos-Renyi (True Count)","Erdos-Renyi (Symmetric)", "Erdos-Renyi (Symmetric True Count)",
            "Small World", "Small World (Directed)", "Modular", "Hierarchical"};
        static constexpr const char* dsStateNames[] = {"Random Uniform", "Random Normal", "Random Cauchy",
            "Random Exponential", "Random Circle", "Splay", "Splay Perturbed", "Identical Modules"};
        static constexpr const char* moduleTypeNames[] = {"Random Uniform", "Random Normal", "Random Cauchy",
            "Random Exponential", "Random Circle", "Splay", "Splay Perturbed"};
        static constexpr const char* moduleTypeIds[] = {"uniform", "normal", "cauchy", "exponential",
            "circle", "splay", "splay_perturbed"};
        int modelSelectedIndex = static_cast<int>(modelParams.modelType);
        int kuramotoModelSelectedIndex = static_cast<int>(modelParams.kuramotoType);
        int adjSelectedIndex = static_cast<int>(adjParams.adjState);
        inline void DrawModelPanelContent();
        inline void DrawTopologyPanelContent();
		inline void DrawInitialsPanelContent();
		inline void DrawSolverParametersPanelContent();
		inline void RenderModals();
};

inline void AppState::DrawModelPanelContent()
{
    ImGui::SeparatorText("Model Configuration");
    if (ImGui::Combo("Model Type",&modelSelectedIndex, modelNames,1))
    {
        modelParams.modelType = static_cast<ModelType>(modelSelectedIndex);
    }
    if (modelParams.modelType==static_cast<ModelType>(0))
    {
        if (ImGui::Combo("Kuramoto Type",&kuramotoModelSelectedIndex,kuramotoModelNames,3))
        {
            modelParams.kuramotoType=static_cast<KuramotoType>(kuramotoModelSelectedIndex);
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
    if (ImGui::Combo("Topology Type", &adjSelectedIndex, adjNames, 10))
    {
        adjParams.adjState = static_cast<NetworkTopology>(adjSelectedIndex);
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
                solverParams.derivative = MathEngine::kuramoto_general_wrapper(kParams);
            }
            else if (modelParams.kuramotoType == KuramotoType::KuramotoSparse)
            {
                MathEngine::KuramotoSparseParams kParams;
                kParams.K = modelParams.K;
                kParams.N = modelParams.N;
                kParams.alpha = modelParams.alpha;
                kParams.omega = modelParams.iFrqnc;
                kParams.sparse_adj = sparseAdj;
                solverParams.derivative = MathEngine::kuramoto_sparse_wrapper(kParams);
            }
            else if (modelParams.kuramotoType == KuramotoType::KuramotoSpecial)
            {
                MathEngine::KuramotoModularParams kParams;
                kParams.intra_K = modelParams.K;
                kParams.inter_K = modelParams.Q;
                kParams.N = modelParams.N;
                kParams.alpha = modelParams.alpha;
                kParams.omega = modelParams.iFrqnc;
                solverParams.derivative = MathEngine::kuramoto_special_modular_wrapper(kParams);
            }
        }
    }
}

inline void AppState::DrawSolverParametersPanelContent()
{
    const ImGuiViewport* viewport = ImGui::GetMainViewport();
    if (ImGui::CollapsingHeader("Time & Basic Stepping", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputDouble("Start Time (t0)", &solverParams.t0, 0.0001, 0.1, "%.15g");
        ImGui::InputDouble("End Time (t1)", &solverParams.t1, 0.0001, 0.1, "%.15g");
        // Initial step size dt
        ImGui::InputDouble("Step Size (dt)", &solverParams.dt, 0.000001, 0.01, "%.15g");
        ImGui::Separator();
        // Multi-step method controls (Adams-Bashforth / Adams-Moulton)
        ImGui::SliderInt("Method Order", &solverParams.order, 1, 10);
        ImGui::SliderInt("ABM Iterations", &solverParams.iterations, 1, 10);
    }
    if (ImGui::CollapsingHeader("Adaptive Step Control"))
    {
        ImGui::Checkbox("Estimate Error", &solverParams.error_estimate);
        ImGui::SameLine();
        ImGui::Checkbox("Enable Variable Step Size", &solverParams.variable_steps);
        if (solverParams.error_estimate && solverParams.variable_steps)
        {
            ImGui::Indent();
            ImGui::TextDisabled("Tolerances & Bounds");
            ImGui::InputDouble("Local Tolerance", &solverParams.local_tol, 0.0, 0.0, "%.1e");
            ImGui::InputDouble("Absolute Tolerance", &solverParams.absolute_tol, 0.0, 0.0, "%.1e");
            ImGui::InputDouble("Min dt", &solverParams.min_dt, 0.0, 0.0, "%.15g");
            ImGui::InputDouble("Max dt", &solverParams.max_dt, 0.0, 0.0, "%.15g");
            ImGui::Separator();
            ImGui::TextDisabled("Step Adaptation Factors");
            ImGui::InputDouble("Decrease Factor", &solverParams.decrease_factor, 0.05, 0.1, "%.15g");
            ImGui::InputDouble("Increase Factor", &solverParams.increase_factor, 0.1, 0.5, "%.15g");
            ImGui::InputDouble("Tol Error Ratio", &solverParams.local_tol_error_ratio, 0.01, 0.05, "%.15g");
            // size_t cast to int for ImGui input
            int maxTrial = static_cast<int>(solverParams.max_trial);
            if (ImGui::InputInt("Max Trials", &maxTrial))
                solverParams.max_trial = static_cast<size_t>(maxTrial>1?maxTrial:1);
            ImGui::Unindent();
        }
    }
    if (ImGui::CollapsingHeader("Error Metrics & Flags"))
    {
        ImGui::Checkbox("Weighted Error Formula", &solverParams.weighted_error);
        ImGui::Checkbox("Norm Error Formula", &solverParams.norm_error);
        ImGui::Checkbox("Record Attempt History", &solverParams.attempts_history);
    }
    // double tau=0.0;
    if (ImGui::CollapsingHeader("Delay Differential Equations (DDE)"))
    {
        ImGui::Checkbox("Is DDE System", &solverParams.is_dde);
        if (solverParams.is_dde)
        {
            ImGui::Indent();
            int maxDelayOrder = static_cast<int>(solverParams.max_delay_order);
            if (ImGui::InputInt("Max Delay Order", &maxDelayOrder))
            {
                solverParams.max_delay_order = static_cast<size_t>(maxDelayOrder > 1 ? maxDelayOrder : 1);
            }
            ImGui::InputDouble("Interpolation Tol", &solverParams.interpolation_tol, 1e-10, 1e-8, "%.1e");
            ImGui::InputDouble("dt Scale (Fine Step)", &solverParams.dt_scale, 0.01, 0.05, "%.15g");
            ImGui::SliderInt("Number of Delays",&nDs,1,20);
            ImGui::SameLine();
            if (ImGui::Button("Set Delay Count"))
            {
                delayTimes.resize(nDs);
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
                solverParams.delay_times=delayTimes;
            }
            ImGui::SameLine();
            if (ImGui::Button("View Delay Values"))
            {
                showDelays=true;
            }
            ImGui::Text("Configured Delays: %zu", solverParams.delay_times.size());
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
                if (ImGui::BeginTable("MatrixGrid", static_cast<int>(nCols + 1), tableFlags, ImVec2(0, 340)))
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
            if (ImGui::BeginChild("PhaseList", ImVec2(0, 300), ImGuiChildFlags_Borders))
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
            if (ImGui::BeginChild("FreqList", ImVec2(0, 300), ImGuiChildFlags_Borders))
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
            ImGui::Text("Total Size: %zu elements", solverParams.delay_times.size());
            ImGui::Separator();
            ImGui::Spacing();

            // Scrollable child box for array elements
            if (ImGui::BeginChild("ArrayList", ImVec2(0, 340), ImGuiChildFlags_Borders))
            {
                if (solverParams.delay_times.empty())
                {
                    ImGui::TextDisabled("Array is empty.");
                }
                else
                {
                    for (size_t i = 0; i < solverParams.delay_times.size(); ++i)
                    {
                        ImGui::Text("[%03zu]  %.6f", i+1, solverParams.delay_times[i]);
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
                showDelays = false;
            }
        }
        ImGui::End();
    }
}
