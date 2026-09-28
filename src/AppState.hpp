#pragma once
#include "MM/typedefs/header.hpp"
#include "MM/models/kuramoto/general.hpp"
#include "MM/models/kuramoto/sparse.hpp"
#include "MM/models/kuramoto/special.hpp"
#include "MM/models/OA-Ansatz.hpp"
#include "MM/models/molecular-dynamics.hpp"
#include "MM/solvers/ODE/rk/explicit/rk1-solver.hpp"
#include "MM/solvers/ODE/rk/explicit/rk2-solver.hpp"
#include "MM/solvers/ODE/rk/explicit/rk3-solver.hpp"
#include "MM/solvers/ODE/rk/explicit/rk4-solver.hpp"
#include "MM/solvers/ODE/rk/explicit/rk4-variants.hpp"
#include "MM/solvers/ODE/multistep/ab-solver.hpp"
#include "MM/solvers/ODE/multistep/abm-solver.hpp"
#include "MM/initializers/initials.hpp"
#include "MM/network/topology.hpp"
#include "MM/utility/write.hpp"
#include "MM/models/kuramoto/general.hpp"
#include "MM/models/kuramoto/sparse.hpp"
#include "MM/models/kuramoto/special.hpp"
// #include "imgui_internal.h"
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
#include <thread>
#include <atomic>
#include <chrono>

static constexpr const uint8_t MultiStepOrderMin{1};
static constexpr const uint8_t MultiStepOrderMax{10};
static constexpr const uint8_t MultiStepIterationsMin{1};
static constexpr const uint8_t MultiStepIterationsMax{10};



inline MathEngine::SolverFunc rk1_wrapper()
{
    return [](const MathEngine::ODESolverParameters& Params) {return MathEngine::rk1_solver_callback(Params);};
}
inline MathEngine::SolverFunc rk2_wrapper()
{
    return [](const MathEngine::ODESolverParameters& Params) {return MathEngine::rk2_solver_callback(Params);};
}
inline MathEngine::SolverFunc rk3_wrapper()
{
    return [](const MathEngine::ODESolverParameters& Params) {return MathEngine::rk3_solver_callback(Params);};
}
inline MathEngine::SolverFunc rk4_wrapper()
{
    return [](const MathEngine::ODESolverParameters& Params) {return MathEngine::rk4_solver_callback(Params);};
}
inline MathEngine::SolverFunc rk4_38_wrapper()
{
    return [](const MathEngine::ODESolverParameters& Params) {return MathEngine::rk4_38_solver_callback(Params);};
}
inline MathEngine::SolverFunc rk4_ralston_wrapper()
{
    return [](const MathEngine::ODESolverParameters& Params) {return MathEngine::rk4_ralston_solver_callback(Params);};
}
inline MathEngine::SolverFunc rk4_gill_wrapper()
{
    return [](const MathEngine::ODESolverParameters& Params) {return MathEngine::rk4_gill_solver_callback(Params);};
}
inline MathEngine::SolverFunc adams_bashforth_wrapper()
{
    return [](const MathEngine::ODESolverParameters& Params) {return MathEngine::adams_bashforth_solver_callback(Params);};
}
inline MathEngine::SolverFunc adams_bashforth_moulton_wrapper()
{
    return [](const MathEngine::ODESolverParameters& Params) {return MathEngine::adams_bashforth_moulton_solver_callback(Params);};
}

enum class ModelType
{
	Kuramoto=0,
	OttAntonsen,
	MolecularDynamics
};
enum class KuramotoType
{
    KuramotoGeneral=0,
    KuramotoSparse,
    KuramotoSpecial
};
enum class OAType
{
    OASingle=0,
    OAGeneral
};
enum class MolecularDynamicsType
{
    LennardJones=0,
    WCA,
    Morse
};
enum class MDIntegratorType
{
    VelocityVerlet=0,
    Leapfrog
};
enum class MDThermostatType
{
    None=0,
    Rescale,
    Berendsen,
    Andersen,
    Langevin,
    NoseHoover
};
enum class MDInitialConditionType
{
    SquareLattice=0,
    HexagonalLattice,
    Random,
    TwoPhaseSlab,
    BinaryMixture
};
enum class SolverMethod
{
    RK1=0,
    RK2,
    RK3,
    RK4,
    RK4_38,
    RK4_Gill,
    RK4_Ralston,
    AB,
    ABM
};
enum class SidebarTab : int
{
    Model=0,
    Topology,
    Initials,
    Solver,
    Plot,
    Run,
    Save,
    Count
};
struct SidebarTabInfo
{
    const char* icon;
    const char* title;
};
constexpr SidebarTabInfo SidebarTabs[] = {
    { ICON_FA_DIAGRAM_PROJECT, "Model"             },
    { ICON_FA_NETWORK_WIRED,   "Topology"          },
    { ICON_FA_WAVE_SQUARE,     "Initial Condition" },
    { ICON_FA_GEARS,           "Solver"            },
    { ICON_FA_CHART_LINE,      "Plot"              },
    { ICON_FA_PLAY,            "Run"               },
    { ICON_FA_FLOPPY_DISK,     "Save"              },
};
constexpr int SidebarTabCount = static_cast<int>(sizeof(SidebarTabs) / sizeof(SidebarTabs[0]));
constexpr float SidebarRailWidth = 52.0f;
static_assert(SidebarTabCount == static_cast<int>(SidebarTab::Count), "Sidebar icon table must match SidebarTab");

struct GeneralModelParams
{
    MathEngine::dVec iFrqnc;
    MathEngine::dVec iPhase;
    double K = 1.0;
    double Q = 0.5;
    double alpha = 0.0;
	size_t N = 50;
    size_t nModules = 1;
    size_t sModules = 50;
    ModelType modelType = ModelType::Kuramoto;
    KuramotoType kuramotoType = KuramotoType::KuramotoGeneral;
    OAType oaType = OAType::OASingle;
    int modelSelectedIndex = 0;
    int kuramotoModelSelectedIndex = 0;
    int oaModelSelectedIndex = 0;
    // Ott-Antonsen (single community): uses `K` (coupling) plus these Lorentzian params.
    double oaGamma = 1.0;
    double oaMu    = 0.0;
    // Ott-Antonsen (general / multi-community).
    size_t oaC = 2;              // number of communities
    double oaRho  = 0.5;         // fallback initial order-parameter magnitude
    double oaPhi  = 0.0;         // fallback initial order-parameter phase
    MathEngine::dVec     oaGammas;   // per-community gamma (size C)
    MathEngine::dVec     oaMus;      // per-community mu    (size C)
    MathEngine::dVec     oaEta;      // per-community population fractions (size C)
    MathEngine::dMatrix  oaK;        // C x C coupling strengths
    MathEngine::dVec     oaIC;       // initial state (interleaved Re/Im)
    GeneralModelParams(size_t n=50) : N(n) {};
};
struct MDParams
{
    size_t numParticles = 100;
    double width  = 800.0;
    double height = 600.0;

    double mass        = 1.0;
    double radius      = 0.1;
    double massRatio   = 2.0;
    double radiusRatio = 1.5;

    double sigma       = 1.0;
    double epsilon     = 1.0;
    double cutoffCoeff = 2.5;
    double morseAlpha  = 1.0;

    double temperature   = 1.0;
    double restitution   = 0.5;
    double minSeparation = 0.8;
    int    seed          = 41;

    bool periodicBoundaryCondition = false;
    bool bounce           = true;
    bool hardSphereCollisions = false;

    MolecularDynamicsType  potential       = MolecularDynamicsType::LennardJones;
    MDInitialConditionType initialCondition = MDInitialConditionType::SquareLattice;
    MDIntegratorType       integrator       = MDIntegratorType::VelocityVerlet;
    MDThermostatType       thermostat       = MDThermostatType::None;

    double thermostatT   = 1.0;   // thermostat target temperature
    double thermostatTau = 1.0;   // Berendsen / Nose-Hoover relaxation time
    double langevinGamma = 1.0;   // Langevin friction
    double andersenNu    = 5.0;   // Andersen collision frequency

    bool   barostat       = false;
    double targetPressure = 0.0;
    double barostatTau    = 10.0;

    double dt     = 0.001;
    double t1     = 10.0;
    int    stride = 50;

    int potentialIndex         = 0;
    int initialConditionIndex  = 0;
    int integratorIndex        = 0;
    int thermostatIndex        = 0;
};
struct MDRunState
{
    MathEngine::dVec posX, posY, velX, velY;      // latest particle state (for the live view)
    std::vector<double> time, temperature, kineticEnergy, potentialEnergy, totalEnergy,
                        pressure, psi4, psi6, msd; // observables time series
    void clear()
    {
        posX.clear(); posY.clear(); velX.clear(); velY.clear();
        time.clear(); temperature.clear(); kineticEnergy.clear(); potentialEnergy.clear();
        totalEnergy.clear(); pressure.clear(); psi4.clear(); psi6.clear(); msd.clear();
    }
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
    MathEngine::InitState initState = MathEngine::InitState::Uniform;
    MathEngine::InitType moduleType = MathEngine::InitType::Uniform;
    bool showArray = false;
    bool identical = false;
};
struct NetParams
{
    double weightMin = 0.0, weightMax = 1.0, weight = 0.5, weightIn = 0.0, weightOut = 1.0;
    double prob = 0.5, probIn = 0.5, probOut = 0.5;
    double decayRatio = 0.1;
    size_t sModulesBase = 10, nModulesBase = 2, hLevels = 2;
    size_t sModulesM = 10, nModulesM = 2;
    MathEngine::NetworkTopology adjState = MathEngine::NetworkTopology::ErdosRenyi;
    int adjSelectedIndex = 3;
    int meanDegree = 2;
    int seed = 41;
    bool showAdjMatrix = false;
};
struct SolverParams
{
    MathEngine::ODESolverParameters solverParams;
    MathEngine::SolverResults solverResults;
    MathEngine::SolverFunc solverFunc=nullptr;
    SolverMethod solverMethod = SolverMethod::RK4;
    int nDs = 1;
    int solverMethodSelectedIndex = 3;
};
struct PlotParams
{
    MathEngine::Matrix<double> plotYModules;
    MathEngine::OneStepSolverResult liveSolverRes;
    MathEngine::Vec<double> liveTimePoints = {};
    MathEngine::Vec<double> liveState      = {};
    MathEngine::Vec<double> plotX          = {};
    MathEngine::Vec<double> plotY          = {};
    MathEngine::Vec<double> plotXTrail     = {};
    MathEngine::Vec<double> plotYTrail     = {};
    MathEngine::Vec<ImVec4> plotColors;
    MathEngine::Vec<ImVec4> plotSecondColors;
    MathEngine::Vec<ImVec4> plotThirdColors;
	std::mutex              plotMutex;
    size_t offset                         = 0;
    int    Stride                         = 50;
    int    trailCount                     = 5000;
    bool   showPlot                       = false;
    bool   showPlotSecond                 = false;
    bool   showPlotThird                  = false;
};
struct SaveParams
{
    char outputDir[256]   = "EMMAOutput";
    bool saveAdjacency      = true;
    bool savePhases         = true;
    bool saveFrequencies    = true;
    bool saveSolution       = true;
    bool saveTimePoints     = true;
    bool saveOrderParameter = true;
    bool saveOAInitial      = true;
    bool saveOAGamma        = true;
    bool saveOAMu           = true;
    bool saveOAEta          = true;
    bool saveOACoupling     = true;
    bool saveMDFinalState   = true;
    bool saveMDObservables  = true;
    bool binary             = false;
    bool append             = false;
    int  precision          = 15;
    int  colWidth           = 20;
    int  fpFormatIndex      = 0; // 0=Scientific, 1=Fixed, 2=Default
    int  alignmentIndex     = 2; // 0=Left, 1=Right, 2=Center, 3=None
};

// Each savable artifact is tagged with the model it belongs to, so we only
// write (and offer in the UI) the data relevant to the active model.
// Artifacts marked universal (e.g. the solution trajectory) apply to every model.
enum class SaveArtifactKind : int
{
    Adjacency = 0,
    InitialPhases,
    IntrinsicFrequencies,
    Solution,
    TimePoints,
    OrderParameter,
    OAInitialOrder,
    OAGamma,
    OAMu,
    OAEta,
    OACoupling,
    MDFinalState,
    MDObservables,
    Count
};
struct SaveArtifact
{
    const char* label;
    const char* subDir;
    const char* fileName;
    ModelType   model;     // owning model (ignored when universal)
    bool        universal; // applies to every model type
    bool        SaveParams::* toggle;
};
constexpr SaveArtifact SaveArtifacts[] = {
    { "Adjacency Matrix",       "Topology",          "AdjacencyMatrix.csv",      ModelType::Kuramoto, false, &SaveParams::saveAdjacency      },
    { "Initial Phases",         "InitialConditions", "InitialPhases.csv",        ModelType::Kuramoto, false, &SaveParams::savePhases         },
    { "Intrinsic Frequencies",  "InitialConditions", "IntrinsicFrequencies.csv", ModelType::Kuramoto, false, &SaveParams::saveFrequencies    },
    { "Solution",               "Solution",          "Solution.csv",             ModelType::Kuramoto, true,  &SaveParams::saveSolution       },
    { "Time Points",            "Solution",          "TimePoints.csv",           ModelType::Kuramoto, true,  &SaveParams::saveTimePoints     },
    { "Order Parameter",        "Analysis",          "OrderParameter.csv",       ModelType::Kuramoto, true,  &SaveParams::saveOrderParameter },
    { "Initial Order Parameters","InitialConditions", "InitialOrderParameters.csv", ModelType::OttAntonsen, false, &SaveParams::saveOAInitial  },
    { "Lorentzian Width (gamma)","InitialConditions", "Gamma.csv",                 ModelType::OttAntonsen, false, &SaveParams::saveOAGamma    },
    { "Lorentzian Center (mu)",  "InitialConditions", "Mu.csv",                    ModelType::OttAntonsen, false, &SaveParams::saveOAMu       },
    { "Population Fractions (eta)","InitialConditions","Eta.csv",                  ModelType::OttAntonsen, false, &SaveParams::saveOAEta      },
    { "Coupling Matrix (K)",     "Topology",         "CouplingMatrix.csv",         ModelType::OttAntonsen, false, &SaveParams::saveOACoupling },
    { "Final Particle State",    "Solution",         "FinalState.csv",             ModelType::MolecularDynamics, false, &SaveParams::saveMDFinalState },
    { "MD Observables",          "Analysis",         "Observables.csv",            ModelType::MolecularDynamics, false, &SaveParams::saveMDObservables },
};
constexpr int SaveArtifactCount = static_cast<int>(sizeof(SaveArtifacts) / sizeof(SaveArtifacts[0]));
static_assert(SaveArtifactCount == static_cast<int>(SaveArtifactKind::Count), "SaveArtifact table must match SaveArtifactKind");

////////////////////////////////////
/////                          /////
/////     STATE OF THE APP     /////
/////                          /////
////////////////////////////////////
class AppState
{
	public:
        GeneralModelParams modelParams = GeneralModelParams(50);
        MDParams mdParams;
        MDRunState mdRunState;
        DistParams phaseParams;
        DistParams frqncParams;
        DistParams oaRhoParams{0.0, 0.5};               // initial order magnitude (rho)
        DistParams oaPhiParams{-MathEngine::PI, MathEngine::PI}; // initial order phase (phi)
        DistParams oaGammaParams{0.5, 1.5};             // Lorentzian half-width
        DistParams oaMuParams{-1.0, 1.0};               // Lorentzian center
        DistParams oaEtaParams;                         // population fractions (normalized)
        NetParams adjParams;
        SolverParams solverParams;
        PlotParams plotParams;
        SaveParams saveParams;
        std::string saveStatus;
        MathEngine::dMatrix adj; // Adjacency (for any system that might need it)
        MathEngine::SparsedMatrix sparseAdj = MathEngine::SparsedMatrix(modelParams.N); // Sparse adjacency
        // MathEngine::dVec delayTimes = {0.0};
        Color BgColor = Color{15,15,15,255};
        std::string appTitle;
        std::atomic<double> timeInv{0.0f};
        std::chrono::steady_clock::time_point processStartTime{};
        std::chrono::steady_clock::time_point processStopTime{};
        long milliSec{}, totalSec{}, sec{}, min{}, hour{}, processDuration{};
        std::atomic<float> simProgress{0.0f};
        float colorR=0.2,colorG=0.8,colorB=0.8,colorA=1.0;
        float padding = 10.0f;
        float cellWidthBase = 1.0;
        float cellHeightBase = 1.0;
        size_t initW = 800;
        size_t initH = 600;
        std::atomic<bool> isSimRunning{false};
        std::thread simThread;
        bool showStyleEditor=false;
        bool showDelays=false;
        bool showPlot=false;
        bool hasSimRan=false;
        bool DarkTheme=true;
        bool showAbout=false;
        int activeSidebarTab = -1; // -1 = drawer closed (no section selected)

        ~AppState()
        {
            // Join the simulation thread (if any) on shutdown to avoid std::terminate.
            if (simThread.joinable()) simThread.join();
        }

        inline void RenderUI()
        {
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

            // Activity bar: a seamless strip flush against the app border (always visible).
            ImGui::SetNextWindowPos(ImVec2(viewport->WorkPos.x, viewport->WorkPos.y), ImGuiCond_Always);
            ImGui::SetNextWindowSize(ImVec2(SidebarRailWidth, viewport->WorkSize.y), ImGuiCond_Always);
            ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
            ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
            ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.03f, 0.03f, 0.03f, 1.0f));
            if (ImGui::Begin("##ActivityBar", nullptr,
                ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
                ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoSavedSettings |
                ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus))
            {
                DrawActivityBar();
            }
            ImGui::End();
            ImGui::PopStyleColor();
            ImGui::PopStyleVar(2);

            // Drawer: a resizable content panel, shown only while a section is open.
            if (activeSidebarTab >= 0)
            {
                const char* title = SidebarTabs[activeSidebarTab].title;
                char drawerName[64];
                snprintf(drawerName, sizeof(drawerName), "%s##SidebarDrawer", title);
                ImGui::SetNextWindowPos(ImVec2(viewport->WorkPos.x + SidebarRailWidth, viewport->WorkPos.y), ImGuiCond_Always);
                ImGui::SetNextWindowSize(ImVec2(viewport->WorkSize.x * 0.46f, viewport->WorkSize.y), ImGuiCond_FirstUseEver);
                ImGui::SetNextWindowSizeConstraints(ImVec2(280.0f, 200.0f), ImVec2(FLT_MAX, FLT_MAX));
                ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(12.0f, 12.0f));
                if (ImGui::Begin(drawerName, nullptr, ImGuiWindowFlags_NoCollapse))
                {
                    DrawDrawerContent();
                }
                ImGui::End();
                ImGui::PopStyleVar();
            }

            RenderModals();
            DrawAboutPage(showAbout);
            DrawPlotWindow();
        }
	private:
        static constexpr const char* modelNames[] = {"Kuramoto", "Ott-Antonsen", "Molecular Dynamics"};
        static constexpr const char* kuramotoModelNames[] = {"Kuramoto (General)", "Kuramoto (Sparse)", "Kuramoto (Modular)"};
        static constexpr const char* oaModelNames[] = {"OA (Single Community)", "OA (Multi-community)"};
        static constexpr const char* mdPotentialNames[] = {"Lennard-Jones", "WCA", "Morse"};
        static constexpr const char* mdIntegratorNames[] = {"Velocity-Verlet", "Leapfrog"};
        static constexpr const char* mdThermostatNames[] = {"None", "Rescale", "Berendsen", "Andersen", "Langevin", "Nose-Hoover"};
        static constexpr const char* mdInitNames[] = {"Square Lattice", "Hexagonal Lattice", "Random", "Two-Phase Slab", "Binary Mixture"};
        static constexpr const char* adjNames[] = {"Random (Uniform)", "Random (Uniform Symmetric)", "Erdos-Renyi",
            "Erdos-Renyi (True Count)","Erdos-Renyi (Symmetric)", "Erdos-Renyi (Symmetric True Count)",
            "Small World", "Small World (Directed)", "Modular", "Hierarchical"};
        static constexpr const char* dsStateNames[] = {"Random Uniform", "Random Normal", "Random Cauchy",
            "Random Exponential", "Random Circle", "Splay", "Splay Perturbed", "Modules (by type)"};
        static constexpr const char* moduleTypeNames[] = {"Random Uniform", "Random Normal", "Random Cauchy",
            "Random Exponential", "Random Circle", "Splay", "Splay Perturbed"};
        static constexpr const char* solverMethodNames[] = {"Euler (RK1)", "Midpoint (RK2)", "Runge-Kutta 3rd Order (RK3)",
            "Runge-Kutta 4th Order (Standard RK4)", "Runge-Kutta 4th Order (3/8 RK4 variant)", "Runge-Kutta 4th Order (Gill's RK4 variant)",
            "Runge-Kutta 4th Order (Ralston's RK4 variant)","Adams-Bashforth (Multistep Predictor)", "Adams-Bashforth-Moulton (Predictor-Corrector)"
        };
        inline void DrawModelPanelContent();
        inline void DrawTopologyPanelContent();
        inline MathEngine::dMatrix GenerateOACouplingMatrix();
        inline void DrawInitDistPicker(DistParams& p, const char* label);
        inline MathEngine::dVec GenerateOAVector(const DistParams& p, size_t C, int seedOffset = 0);
		inline void DrawInitialsPanelContent();
		inline void DrawODESolverParametersPanelContent();
		inline void RenderModals();
        inline void DrawPlotWindow();
        inline void StartSimulation();
        inline void StartMolecularDynamics();
        inline void DrawProgressBar();
        inline void DrawPlotPanelContent();
        inline void RenderChrono();
        inline bool DrawSidebarIcon(int tabIndex, const char* icon, const char* title);
        inline bool DrawActivityButton(const char* icon, const char* title, bool active);
        inline void DrawActivityBar();
        inline void DrawDrawerContent();
        inline void DrawRunPanelContent();
        inline void DrawSavePanelContent();
        inline void SaveSimulationData(const std::filesystem::path& outputDir);
        inline std::filesystem::path BuildDefaultOutputPath();
        inline MathEngine::IO::WriteOptions MakeWriteOptions(const std::filesystem::path& filePath, std::string_view header = {});
        inline bool ArtifactApplies(const SaveArtifact& artifact) const;
        inline bool WriteArtifactData(SaveArtifactKind kind, const std::filesystem::path& filePath);
        inline void DrawVectorViewer(bool& open, const char* title, const MathEngine::dVec& data, const char* emptyText);
};

inline bool AppState::DrawActivityButton(const char* icon, const char* title, bool active)
{
    const float iconSize = 40.0f;
    bool clicked = false;

    ImGui::PushFont(g_FONTs.icons);
    // Flat, borderless activity-bar button (seamless look).
    ImGui::PushStyleColor(ImGuiCol_Button,        ImVec4(0.0f, 0.0f, 0.0f, 0.0f));
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(1.0f, 1.0f, 1.0f, 0.10f));
    ImGui::PushStyleColor(ImGuiCol_ButtonActive,  ImVec4(1.0f, 1.0f, 1.0f, 0.16f));
    ImGui::PushStyleColor(ImGuiCol_Text, active ? ImVec4(1.00f, 0.39f, 0.80f, 1.0f)
                                                : ImVec4(0.72f, 0.72f, 0.72f, 1.0f));

    if (ImGui::Button(icon, ImVec2(-1.0f, iconSize)))
    {
        clicked = true;
    }

    if (active)
    {
        const ImVec2 min = ImGui::GetItemRectMin();
        const ImVec2 max = ImGui::GetItemRectMax();
        ImGui::GetWindowDrawList()->AddRectFilled(
            ImVec2(min.x, min.y + 10.0f),
            ImVec2(min.x + 3.0f, max.y - 10.0f),
            IM_COL32(255, 100, 205, 255));
    }

    ImGui::PopStyleColor(4);
    ImGui::PopFont();

    if (ImGui::IsItemHovered())
    {
        ImGui::SetTooltip("%s", title);
    }
    return clicked;
}

inline bool AppState::DrawSidebarIcon(int tabIndex, const char* icon, const char* title)
{
    return DrawActivityButton(icon, title, activeSidebarTab == tabIndex);
}

inline void AppState::DrawActivityBar()
{
    ImGui::Dummy(ImVec2(0.0f, 6.0f));
    for (int i = 0; i < SidebarTabCount; ++i)
    {
        if (DrawSidebarIcon(i, SidebarTabs[i].icon, SidebarTabs[i].title))
        {
            activeSidebarTab = (activeSidebarTab == i) ? -1 : i;
        }
    }

    // Bottom section: Options + Help (moved here from the old top menu bar).
    const float bottomHeight = 2.0f * 40.0f + ImGui::GetStyle().ItemSpacing.y + 8.0f;
    ImGui::SetCursorPosY(ImGui::GetWindowHeight() - bottomHeight);

    if (DrawActivityButton(ICON_FA_GEAR, "Options", false))
    {
        ImGui::OpenPopup("##ActivityOptionsMenu");
    }
    const ImVec2 gearMin = ImGui::GetItemRectMin();
    ImGui::SetNextWindowPos(
        ImVec2(ImGui::GetWindowPos().x + SidebarRailWidth + 4.0f, ImGui::GetWindowPos().y + gearMin.y),
        ImGuiCond_Always);
    if (ImGui::BeginPopup("##ActivityOptionsMenu"))
    {
        DrawFontMenu();
        ImGui::Separator();
        if (ImGui::MenuItem("Show ImGui Style Editor", nullptr, &showStyleEditor)) {}
        ImGui::EndPopup();
    }

    if (DrawActivityButton(ICON_FA_CIRCLE_QUESTION, "Help", false))
    {
        showAbout = true;
    }
}

inline void AppState::DrawDrawerContent()
{
    switch (static_cast<SidebarTab>(activeSidebarTab))
    {
        case SidebarTab::Model:    DrawModelPanelContent();               break;
        case SidebarTab::Topology: DrawTopologyPanelContent();            break;
        case SidebarTab::Initials: DrawInitialsPanelContent();            break;
        case SidebarTab::Solver:   DrawODESolverParametersPanelContent(); break;
        case SidebarTab::Plot:     DrawPlotPanelContent();                break;
        case SidebarTab::Run:      DrawRunPanelContent();                 break;
        case SidebarTab::Save:     DrawSavePanelContent();                break;
        default: break;
    }
}

inline void AppState::DrawRunPanelContent()
{
    DrawProgressBar();
    ImGui::Spacing();
    RenderChrono();
    ImGui::Spacing();
    bool running = isSimRunning.load();
    if (running) ImGui::BeginDisabled();
    if (ImGui::Button(running ? "Running..." : "Begin Simulation", ImVec2(-1, 0))) StartSimulation();
    if (running) ImGui::EndDisabled();
}

inline MathEngine::IO::WriteOptions AppState::MakeWriteOptions(const std::filesystem::path& filePath, std::string_view header)
{
    MathEngine::IO::FPFormat fmt = MathEngine::IO::FPFormat::Scientific;
    switch (saveParams.fpFormatIndex)
    {
        case 1: fmt = MathEngine::IO::FPFormat::Fixed; break;
        case 2: fmt = MathEngine::IO::FPFormat::Default; break;
        default: fmt = MathEngine::IO::FPFormat::Scientific; break;
    }
    MathEngine::IO::Alignment align = MathEngine::IO::Alignment::Center;
    switch (saveParams.alignmentIndex)
    {
        case 0: align = MathEngine::IO::Alignment::Left; break;
        case 1: align = MathEngine::IO::Alignment::Right; break;
        case 3: align = MathEngine::IO::Alignment::None; break;
        default: align = MathEngine::IO::Alignment::Center; break;
    }
    return MathEngine::IO::WriteOptions{
        .path      = filePath,
        .separator = ",",
        .header    = header,
        .comment   = {},
        .footer    = {},
        .colWidth  = static_cast<size_t>(std::max(0, saveParams.colWidth)),
        .precision = saveParams.precision,
        .format    = fmt,
        .alignment = align,
        .append    = saveParams.append,
        .binary    = saveParams.binary
    };
}

inline std::filesystem::path AppState::BuildDefaultOutputPath()
{
    std::string kType = "General";
    switch (modelParams.kuramotoType)
    {
        case KuramotoType::KuramotoSparse:  kType = "Sparse";  break;
        case KuramotoType::KuramotoSpecial: kType = "Special"; break;
        default: kType = "General"; break;
    }

    std::string topo = "None";
    if (modelParams.kuramotoType == KuramotoType::KuramotoSpecial)
    {
        topo = "Modular";
    }
    else
    {
        switch (adjParams.adjState)
        {
            case MathEngine::NetworkTopology::Uniform:                    topo = "Uniform"; break;
            case MathEngine::NetworkTopology::UniformSymmetric:           topo = "UniformSymmetric"; break;
            case MathEngine::NetworkTopology::ErdosRenyi:                 topo = "ErdosRenyi"; break;
            case MathEngine::NetworkTopology::ErdosRenyiUniform:          topo = "ErdosRenyiUniform"; break;
            case MathEngine::NetworkTopology::ErdosRenyiSymmetric:        topo = "ErdosRenyiSymmetric"; break;
            case MathEngine::NetworkTopology::ErdosRenyiSymmetricUniform: topo = "ErdosRenyiSymmetricUniform"; break;
            case MathEngine::NetworkTopology::SmallWorld:                 topo = "SmallWorld"; break;
            case MathEngine::NetworkTopology::SmallWorldDirected:         topo = "SmallWorldDirected"; break;
            case MathEngine::NetworkTopology::Modular:                    topo = "Modular"; break;
            case MathEngine::NetworkTopology::Hierarchical:               topo = "Hierarchical"; break;
            default: topo = "None"; break;
        }
    }

    std::string solver = "RK4";
    switch (solverParams.solverMethod)
    {
        case SolverMethod::RK1:         solver = "RK1"; break;
        case SolverMethod::RK2:         solver = "RK2"; break;
        case SolverMethod::RK3:         solver = "RK3"; break;
        case SolverMethod::RK4:         solver = "RK4"; break;
        case SolverMethod::RK4_38:      solver = "RK438"; break;
        case SolverMethod::RK4_Gill:    solver = "RK4Gill"; break;
        case SolverMethod::RK4_Ralston: solver = "RK4Ralston"; break;
        case SolverMethod::AB:          solver = "AB"; break;
        case SolverMethod::ABM:         solver = "ABM"; break;
        default: solver = "RK4"; break;
    }

    std::string folderName;
    if (modelParams.modelType==ModelType::OttAntonsen)
    {
        const std::string oaType = (modelParams.oaType==OAType::OASingle) ? "Single" : "General";
        folderName = "OttAntonsen-" + oaType + "-C" + std::to_string(modelParams.oaC) + "-" + solver;
    }
    else if (modelParams.modelType==ModelType::MolecularDynamics)
    {
        folderName = "MolecularDynamics-" + std::string(mdPotentialNames[mdParams.potentialIndex])
                   + "-N" + std::to_string(mdParams.numParticles);
    }
    else
    {
        folderName = "Kuramoto-" + kType + "-" + topo + "-N" + std::to_string(modelParams.N) + "-" + solver;
    }
    return std::filesystem::path(saveParams.outputDir) / folderName;
}

inline bool AppState::ArtifactApplies(const SaveArtifact& artifact) const
{
    return artifact.universal || artifact.model == modelParams.modelType;
}

inline bool AppState::WriteArtifactData(SaveArtifactKind kind, const std::filesystem::path& filePath)
{
    switch (kind)
    {
        case SaveArtifactKind::Adjacency:
            if (adj.empty()) return false;
            MathEngine::IO::WriteMatrix(adj, MakeWriteOptions(filePath, "Adjacency Matrix"));
            return true;
        case SaveArtifactKind::InitialPhases:
            if (modelParams.iPhase.empty()) return false;
            MathEngine::IO::WriteVector(std::span<const double>(modelParams.iPhase), MakeWriteOptions(filePath, "Initial Phases"));
            return true;
        case SaveArtifactKind::IntrinsicFrequencies:
            if (modelParams.iFrqnc.empty()) return false;
            MathEngine::IO::WriteVector(std::span<const double>(modelParams.iFrqnc), MakeWriteOptions(filePath, "Intrinsic Frequencies"));
            return true;
        case SaveArtifactKind::Solution:
            if (solverParams.solverResults.solution.empty()) return false;
            MathEngine::IO::WriteMatrix(solverParams.solverResults.solution, MakeWriteOptions(filePath, "Solution (rows = time steps, cols = state)"));
            return true;
        case SaveArtifactKind::TimePoints:
            if (solverParams.solverResults.timePoints.empty()) return false;
            MathEngine::IO::WriteVector(std::span<const double>(solverParams.solverResults.timePoints), MakeWriteOptions(filePath, "Time Points"));
            return true;
        case SaveArtifactKind::OrderParameter:
        {
            MathEngine::dMatrix op;
            std::string_view header;
            {
                std::lock_guard<std::mutex> lock(plotParams.plotMutex);
                if (plotParams.plotYModules.Rows() > 0 && plotParams.plotYModules.Cols() > 0)
                {
                    // Per-module order parameter (modular / hierarchical systems).
                    const size_t nModules = plotParams.plotYModules.Rows();
                    const size_t nTime = std::min(plotParams.plotX.size(), static_cast<size_t>(plotParams.plotYModules.Cols()));
                    if (nTime > 0)
                    {
                        op = MathEngine::dMatrix(nTime, nModules + 1);
                        for (size_t t = 0; t < nTime; ++t)
                        {
                            op[t, 0] = plotParams.plotX[t];
                            for (size_t m = 0; m < nModules; ++m)
                                op[t, m + 1] = plotParams.plotYModules[m, t];
                        }
                        header = "Per-module order parameter (time, module 1..M)";
                    }
                }
                else
                {
                    // Global order parameter.
                    const size_t n = std::min(plotParams.plotX.size(), plotParams.plotY.size());
                    if (n > 0)
                    {
                        op = MathEngine::dMatrix(n, 2);
                        for (size_t i = 0; i < n; ++i)
                        {
                            op[i, 0] = plotParams.plotX[i];
                            op[i, 1] = plotParams.plotY[i];
                        }
                        header = "Order parameter (time, rho)";
                    }
                }
            }
            if (op.empty()) return false;
            MathEngine::IO::WriteMatrix(op, MakeWriteOptions(filePath, header));
            return true;
        }
        case SaveArtifactKind::OAInitialOrder:
        {
            if (modelParams.oaIC.empty()) return false;
            const size_t C = modelParams.oaIC.size() / 2;
            MathEngine::dMatrix op(C, 2, 0.0); // [rho, phi]
            for (size_t c = 0; c < C; ++c)
            {
                const double x = modelParams.oaIC[2*c + 0];
                const double y = modelParams.oaIC[2*c + 1];
                op[c, 0] = std::hypot(x, y);
                op[c, 1] = std::atan2(y, x);
            }
            MathEngine::IO::WriteMatrix(op, MakeWriteOptions(filePath, "Initial order parameters (rho, phi)"));
            return true;
        }
        case SaveArtifactKind::OAGamma:
            if (modelParams.oaGammas.empty()) return false;
            MathEngine::IO::WriteVector(std::span<const double>(modelParams.oaGammas), MakeWriteOptions(filePath, "Lorentzian half-width (gamma)"));
            return true;
        case SaveArtifactKind::OAMu:
            if (modelParams.oaMus.empty()) return false;
            MathEngine::IO::WriteVector(std::span<const double>(modelParams.oaMus), MakeWriteOptions(filePath, "Lorentzian center (mu)"));
            return true;
        case SaveArtifactKind::OAEta:
            if (modelParams.oaEta.empty()) return false;
            MathEngine::IO::WriteVector(std::span<const double>(modelParams.oaEta), MakeWriteOptions(filePath, "Population fractions (eta)"));
            return true;
        case SaveArtifactKind::OACoupling:
            if (modelParams.oaK.empty()) return false;
            MathEngine::IO::WriteMatrix(modelParams.oaK, MakeWriteOptions(filePath, "Community coupling matrix (K)"));
            return true;
        case SaveArtifactKind::MDFinalState:
        {
            if (mdRunState.posX.empty()) return false;
            const size_t N = mdRunState.posX.size();
            MathEngine::dMatrix state(N, 4, 0.0);
            for (size_t i = 0; i < N; ++i)
            {
                state[i, 0] = mdRunState.posX[i];
                state[i, 1] = mdRunState.posY[i];
                state[i, 2] = mdRunState.velX[i];
                state[i, 3] = mdRunState.velY[i];
            }
            MathEngine::IO::WriteMatrix(state, MakeWriteOptions(filePath, "Final state (x, y, vx, vy)"));
            return true;
        }
        case SaveArtifactKind::MDObservables:
        {
            if (mdRunState.time.empty()) return false;
            const size_t rows = mdRunState.time.size();
            const size_t cols = 9;
            MathEngine::dMatrix obs(rows, cols, 0.0);
            for (size_t r = 0; r < rows; ++r)
            {
                obs[r, 0] = mdRunState.time[r];
                obs[r, 1] = mdRunState.temperature[r];
                obs[r, 2] = mdRunState.kineticEnergy[r];
                obs[r, 3] = mdRunState.potentialEnergy[r];
                obs[r, 4] = mdRunState.totalEnergy[r];
                obs[r, 5] = mdRunState.pressure[r];
                obs[r, 6] = mdRunState.psi4[r];
                obs[r, 7] = mdRunState.psi6[r];
                obs[r, 8] = mdRunState.msd[r];
            }
            MathEngine::IO::WriteMatrix(obs, MakeWriteOptions(filePath, "time, temperature, kineticEnergy, potentialEnergy, totalEnergy, pressure, psi4, psi6, msd"));
            return true;
        }
        default: return false;
    }
}

inline void AppState::SaveSimulationData(const std::filesystem::path& outputDir)
{
    try
    {
        std::filesystem::create_directories(outputDir);
        size_t written = 0;

        for (int i = 0; i < SaveArtifactCount; ++i)
        {
            const SaveArtifact& artifact = SaveArtifacts[i];
            if (!ArtifactApplies(artifact)) continue;
            if (!(saveParams.*(artifact.toggle))) continue;

            if (WriteArtifactData(static_cast<SaveArtifactKind>(i), outputDir / artifact.subDir / artifact.fileName))
                ++written;
        }

        saveStatus = "Saved " + std::to_string(written) + " file(s) to: " + outputDir.string();
    }
    catch (const std::exception& e)
    {
        saveStatus = std::string("Save failed: ") + e.what();
    }
}

inline void AppState::DrawSavePanelContent()
{
    ImGui::SeparatorText("Save Simulation Data");

    if (ImGui::Button("Save (Default Layout)", ImVec2(-1, 0)))
    {
        SaveSimulationData(BuildDefaultOutputPath());
    }

    if (!saveStatus.empty())
    {
        ImGui::Spacing();
        ImGui::TextWrapped("%s", saveStatus.c_str());
        ImGui::Spacing();
    }

    if (ImGui::CollapsingHeader("Advanced Save Options"))
    {
        ImGui::InputText("Output Directory", saveParams.outputDir, sizeof(saveParams.outputDir));

        ImGui::Spacing();
        ImGui::SeparatorText("What to Save");
        for (int i = 0; i < SaveArtifactCount; ++i)
        {
            const SaveArtifact& artifact = SaveArtifacts[i];
            if (!ArtifactApplies(artifact)) continue;
            ImGui::Checkbox(artifact.label, &(saveParams.*(artifact.toggle)));
        }

        ImGui::Spacing();
        ImGui::SeparatorText("Format");
        static constexpr const char* fpFormatNames[] = { "Scientific", "Fixed", "Default" };
        ImGui::Combo("Number Format", &saveParams.fpFormatIndex, fpFormatNames, 3);
        static constexpr const char* alignNames[] = { "Left", "Right", "Center", "None" };
        ImGui::Combo("Alignment", &saveParams.alignmentIndex, alignNames, 4);
        if (ImGui::InputInt("Precision", &saveParams.precision, 1, 1))
            saveParams.precision = std::clamp(saveParams.precision, 0, 17);
        if (ImGui::InputInt("Column Width", &saveParams.colWidth, 1, 1))
            saveParams.colWidth = std::clamp(saveParams.colWidth, 0, 100);
        ImGui::Checkbox("Binary", &saveParams.binary);
        ImGui::Checkbox("Append", &saveParams.append);

        ImGui::Spacing();
        if (ImGui::Button("Save (Custom Layout)", ImVec2(-1, 0)))
        {
            SaveSimulationData(std::filesystem::path(saveParams.outputDir));
        }
    }
}

inline void AppState::DrawModelPanelContent()
{
    ImGui::SeparatorText("Model Configuration");
    if (ImGui::Combo("Model Type",&modelParams.modelSelectedIndex, modelNames,3))
    {
        modelParams.modelType = static_cast<ModelType>(modelParams.modelSelectedIndex);
    }
    if (modelParams.modelType==ModelType::Kuramoto)
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
    else if (modelParams.modelType==ModelType::OttAntonsen)
    {
        if (ImGui::Combo("OA Type",&modelParams.oaModelSelectedIndex,oaModelNames,2))
        {
            modelParams.oaType=static_cast<OAType>(modelParams.oaModelSelectedIndex);
        }
        ImGui::Spacing();
        switch (modelParams.oaType)
        {
            case OAType::OASingle:
                ImGui::InputDouble("Coupling (K)", &modelParams.K, 0.0001, 0.01, "%.15g");
                ImGui::InputDouble("Lorentzian Width (gamma)", &modelParams.oaGamma, 0.0001, 0.01, "%.15g");
                ImGui::InputDouble("Mean Frequency (mu)", &modelParams.oaMu, 0.0001, 0.01, "%.15g");
                break;
            case OAType::OAGeneral:
            {
                int C = static_cast<int>(modelParams.oaC);
                if (ImGui::InputInt("Communities (C)", &C, 1, 5)) modelParams.oaC = static_cast<size_t>(std::max(1, C));
                ImGui::TextDisabled("State dimension: %zu (interleaved Re/Im)", 2 * modelParams.oaC);
                ImGui::TextDisabled("Coupling matrix (K) is generated in the Topology tab.");
                break;
            }
        }
    }
    else if (modelParams.modelType==ModelType::MolecularDynamics)
    {
        if (ImGui::Combo("Potential",&mdParams.potentialIndex,mdPotentialNames,3))
            mdParams.potential = static_cast<MolecularDynamicsType>(mdParams.potentialIndex);
        ImGui::Spacing();

        int n = static_cast<int>(mdParams.numParticles);
        if (ImGui::InputInt("Particles (N)", &n, 1, 50)) mdParams.numParticles = static_cast<size_t>(std::max(1, n));
        ImGui::InputDouble("Mass", &mdParams.mass, 0.001, 0.1, "%.15g");
        ImGui::InputDouble("Radius", &mdParams.radius, 0.001, 0.01, "%.15g");
        ImGui::Spacing();
        ImGui::InputDouble("Sigma", &mdParams.sigma, 0.001, 0.1, "%.15g");
        ImGui::InputDouble("Epsilon", &mdParams.epsilon, 0.001, 0.1, "%.15g");
        ImGui::InputDouble("Cutoff Coefficient", &mdParams.cutoffCoeff, 0.01, 0.1, "%.15g");
        if (mdParams.potential==MolecularDynamicsType::Morse)
            ImGui::InputDouble("Morse Alpha", &mdParams.morseAlpha, 0.01, 0.1, "%.15g");
        ImGui::Spacing();
        ImGui::InputDouble("Temperature", &mdParams.temperature, 0.001, 0.01, "%.15g");
        ImGui::TextDisabled("Box & boundary conditions are set in the Topology tab.");
    }
}

inline void AppState::DrawTopologyPanelContent()
{
    if (modelParams.modelType==ModelType::MolecularDynamics)
    {
        ImGui::TextDisabled("Molecular dynamics has no adjacency matrix.");
        ImGui::Spacing();
        ImGui::SeparatorText("Box & Boundary Conditions");
        ImGui::InputDouble("Box Width", &mdParams.width, 1.0, 10.0, "%.15g");
        ImGui::InputDouble("Box Height", &mdParams.height, 1.0, 10.0, "%.15g");
        ImGui::Checkbox("Periodic Boundary", &mdParams.periodicBoundaryCondition);
        ImGui::Checkbox("Bounce (walls)", &mdParams.bounce);
        ImGui::Checkbox("Hard-Sphere Collisions", &mdParams.hardSphereCollisions);
        ImGui::InputDouble("Restitution", &mdParams.restitution, 0.01, 0.1, "%.15g");
        return;
    }

    const bool isOAGeneral = (modelParams.modelType==ModelType::OttAntonsen) && (modelParams.oaType==OAType::OAGeneral);

    if (modelParams.modelType==ModelType::OttAntonsen && modelParams.oaType==OAType::OASingle)
    {
        ImGui::TextDisabled("Single community: all-to-all coupling (no topology needed).");
        return;
    }
    if (!isOAGeneral && modelParams.kuramotoType == KuramotoType::KuramotoSpecial) return;

    ImGui::Spacing();
    ImGui::SeparatorText(isOAGeneral ? "Community Coupling (K)" : "Network Topology");
    if (isOAGeneral)
    {
        ImGui::TextDisabled("C = %zu communities", modelParams.oaC);
        ImGui::InputDouble("K Intra (diagonal)", &modelParams.K, 0.0001, 0.01, "%.15g");
    }
    if (ImGui::Combo("Topology Type", &adjParams.adjSelectedIndex, adjNames, 10))
    {
        adjParams.adjState = static_cast<MathEngine::NetworkTopology>(adjParams.adjSelectedIndex);
    }
    switch (adjParams.adjState) {
        case MathEngine::NetworkTopology::Uniform:
        case MathEngine::NetworkTopology::UniformSymmetric:
            ImGui::InputDouble("Min Weight", &adjParams.weightMin, 0.0001, 0.01, "%.15g");
            ImGui::InputDouble("Max Weight", &adjParams.weightMax, 0.0001, 0.01, "%.15g");
            break;
        case MathEngine::NetworkTopology::ErdosRenyi:
        case MathEngine::NetworkTopology::ErdosRenyiUniform:
        case MathEngine::NetworkTopology::ErdosRenyiSymmetric:
        case MathEngine::NetworkTopology::ErdosRenyiSymmetricUniform:
            ImGui::InputDouble("Min Weight", &adjParams.weightMin, 0.0001, 0.01, "%.15g");
            ImGui::InputDouble("Max Weight", &adjParams.weightMax, 0.0001, 0.01, "%.15g");
            if (ImGui::InputDouble("Connection Prob", &adjParams.prob, 0.0001, 0.01, "%.15g"))
                adjParams.prob = std::clamp(adjParams.prob, 0.0, 1.0);
	        break;
        case MathEngine::NetworkTopology::SmallWorld:
        case MathEngine::NetworkTopology::SmallWorldDirected:
            ImGui::InputDouble("Weight", &adjParams.weight, 0.0001, 0.01, "%.15g");
            if (ImGui::InputDouble("Rewiring Prob", &adjParams.prob, 0.0001, 0.01, "%.15g"))
                adjParams.prob = std::clamp(adjParams.prob, 0.0, 1.0);
            if (ImGui::InputInt("Mean Degree", &adjParams.meanDegree, 1, 10))
                adjParams.meanDegree=std::max(0,adjParams.meanDegree);
            break;
        case MathEngine::NetworkTopology::Modular: {
            if (isOAGeneral)
            {
                int C = static_cast<int>(modelParams.oaC);
                if (ImGui::InputInt("Communities (C)", &C, 1, 5)) modelParams.oaC = static_cast<size_t>(std::max(1, C));
                ImGui::InputDouble("Out Weight", &adjParams.weightOut, 0.0001, 0.01, "%.15g");
                if (ImGui::InputDouble("Outer Prob", &adjParams.probOut, 0.0001, 0.01, "%.15g"))
                    adjParams.probOut = std::clamp(adjParams.probOut, 0.0, 1.0);
            }
            else
            {
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
            }
            break;
        }
        case MathEngine::NetworkTopology::Hierarchical: {
            if (isOAGeneral)
            {
                int nB = static_cast<int>(adjParams.nModulesBase);
                int hL = static_cast<int>(adjParams.hLevels);
                if (ImGui::InputInt("Base Modules", &nB, 1, 10)) adjParams.nModulesBase = std::max(1, nB);
                if (ImGui::InputInt("Hierarchy Levels", &hL, 1, 10)) adjParams.hLevels = std::max(1, hL);
                ImGui::TextDisabled("C = %zu communities (derived)", adjParams.nModulesBase * (static_cast<size_t>(1) << (adjParams.hLevels - 1)));
                ImGui::InputDouble("Out Weight", &adjParams.weightOut, 0.0001, 0.01, "%.15g");
                if (ImGui::InputDouble("Outer Prob", &adjParams.probOut, 0.0001, 0.01, "%.15g"))
                    adjParams.probOut = std::clamp(adjParams.probOut, 0.0, 1.0);
                if (ImGui::InputDouble("Decay Ratio", &adjParams.decayRatio, 0.001, 0.1, "%.15g"))
                    adjParams.decayRatio = std::clamp(adjParams.decayRatio, 0.0, 1.0);
            }
            else
            {
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
            }
            break;
        }
    }
    ImGui::InputInt("Seed##Adj", &adjParams.seed, 1, 10);
    if (ImGui::Button(isOAGeneral ? "Generate Coupling Matrix" : "Generate Adjacency", ImVec2(-1, 0)))
    {
        if (isOAGeneral)
        {
            modelParams.oaK = GenerateOACouplingMatrix();
        }
        else
        {
            size_t seedVal = static_cast<size_t>(std::max(1, adjParams.seed));
            switch (adjParams.adjState)
            {
                case MathEngine::NetworkTopology::Uniform:
                    adj = MathEngine::random(modelParams.N,adjParams.weightMin,adjParams.weightMax,seedVal);
                    break;
                case MathEngine::NetworkTopology::UniformSymmetric:
                    adj = MathEngine::random_symmetric(modelParams.N,adjParams.weightMin,adjParams.weightMax,seedVal);
                    break;
                case MathEngine::NetworkTopology::ErdosRenyi:
                    adj = MathEngine::erdos_renyi(modelParams.N,adjParams.prob,adjParams.weightMin,adjParams.weightMax,seedVal);
                    break;
                case MathEngine::NetworkTopology::ErdosRenyiUniform:
                    adj = MathEngine::erdos_renyi_uniform(modelParams.N,adjParams.prob,adjParams.weightMin,adjParams.weightMax,seedVal);
                    break;
                case MathEngine::NetworkTopology::ErdosRenyiSymmetric:
                    adj = MathEngine::erdos_renyi_symmetric(modelParams.N,adjParams.prob,adjParams.weightMin,adjParams.weightMax,seedVal);
                    break;
                case MathEngine::NetworkTopology::ErdosRenyiSymmetricUniform:
                    adj = MathEngine::erdos_renyi_symmetric_uniform(modelParams.N,adjParams.prob,adjParams.weightMin,adjParams.weightMax,seedVal);
                    break;
                case MathEngine::NetworkTopology::SmallWorld:
                    adj = MathEngine::small_world(modelParams.N,adjParams.meanDegree,adjParams.prob,adjParams.weight,seedVal);
                    break;
                case MathEngine::NetworkTopology::SmallWorldDirected:
                    adj = MathEngine::small_world_directed(modelParams.N,adjParams.meanDegree,adjParams.prob,adjParams.weight,seedVal);
                    break;
                case MathEngine::NetworkTopology::Modular:
                    modelParams.sModules = adjParams.sModulesM;
                    modelParams.nModules = adjParams.nModulesM;
                    modelParams.N = modelParams.sModules * modelParams.nModules;
                    adj = MathEngine::modular(modelParams.sModules,modelParams.nModules,adjParams.probIn,adjParams.probOut,
                                              adjParams.weightIn,adjParams.weightOut,
                                              seedVal);
                    break;
                case MathEngine::NetworkTopology::Hierarchical:
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
            for (size_t i=0; i<modelParams.N; ++i)
            {
                for (size_t j=0; j<modelParams.N; ++j)
                {
                    char adjs[64];
                    snprintf(adjs,sizeof(adjs),"%.15g",adj[i][j]);
                    cellWidthBase = std::max(cellWidthBase,ImGui::CalcTextSize(adjs).x);
                    if (cellHeightBase ==1.0f) cellHeightBase = ImGui::CalcTextSize(adjs).y;
                }
            }
        }
    }

    if (ImGui::Button(isOAGeneral ? "View Coupling Matrix" : "View Matrix Values", ImVec2(-1, 0)))
    {
        adjParams.showAdjMatrix = true;
    }
}

inline MathEngine::dMatrix AppState::GenerateOACouplingMatrix()
{
    const size_t seedVal = static_cast<size_t>(std::max(1, adjParams.seed));
    MathEngine::dMatrix K;
    switch (adjParams.adjState)
    {
        case MathEngine::NetworkTopology::Uniform:
            K = MathEngine::random(modelParams.oaC, adjParams.weightMin, adjParams.weightMax, seedVal); break;
        case MathEngine::NetworkTopology::UniformSymmetric:
            K = MathEngine::random_symmetric(modelParams.oaC, adjParams.weightMin, adjParams.weightMax, seedVal); break;
        case MathEngine::NetworkTopology::ErdosRenyi:
            K = MathEngine::erdos_renyi(modelParams.oaC, adjParams.prob, adjParams.weightMin, adjParams.weightMax, seedVal); break;
        case MathEngine::NetworkTopology::ErdosRenyiUniform:
            K = MathEngine::erdos_renyi_uniform(modelParams.oaC, adjParams.prob, adjParams.weightMin, adjParams.weightMax, seedVal); break;
        case MathEngine::NetworkTopology::ErdosRenyiSymmetric:
            K = MathEngine::erdos_renyi_symmetric(modelParams.oaC, adjParams.prob, adjParams.weightMin, adjParams.weightMax, seedVal); break;
        case MathEngine::NetworkTopology::ErdosRenyiSymmetricUniform:
            K = MathEngine::erdos_renyi_symmetric_uniform(modelParams.oaC, adjParams.prob, adjParams.weightMin, adjParams.weightMax, seedVal); break;
        case MathEngine::NetworkTopology::SmallWorld:
            K = MathEngine::small_world(modelParams.oaC, adjParams.meanDegree, adjParams.prob, adjParams.weight, seedVal); break;
        case MathEngine::NetworkTopology::SmallWorldDirected:
            K = MathEngine::small_world_directed(modelParams.oaC, adjParams.meanDegree, adjParams.prob, adjParams.weight, seedVal); break;
        case MathEngine::NetworkTopology::Modular:
            K = MathEngine::modular(1, modelParams.oaC, adjParams.probIn, adjParams.probOut, adjParams.weightIn, adjParams.weightOut, seedVal); break;
        case MathEngine::NetworkTopology::Hierarchical:
            modelParams.oaC = adjParams.nModulesBase * (static_cast<size_t>(1) << (adjParams.hLevels - 1));
            K = MathEngine::hierarchical(1, adjParams.hLevels, adjParams.probIn, adjParams.probOut,
                                         adjParams.weightIn, adjParams.weightOut, adjParams.decayRatio,
                                         seedVal, adjParams.nModulesBase);
            break;
        default:
            K = MathEngine::random(modelParams.oaC, adjParams.weightMin, adjParams.weightMax, seedVal); break;
    }
    // Add the intra-community coupling on the diagonal.
    for (size_t c = 0; c < K.Rows(); ++c) K[c, c] = modelParams.K;
    return K;
}

inline void AppState::DrawInitDistPicker(DistParams& p, const char* label)
{
    char buf[128];
    if (ImGui::Combo(label, &p.typeIndex, moduleTypeNames, 7))
        p.initState = static_cast<MathEngine::InitState>(p.typeIndex);

    switch (p.initState)
    {
        case MathEngine::InitState::Uniform:
            snprintf(buf, sizeof(buf), "Min##%s", label);
            ImGui::InputDouble(buf, &p.minVal, 0.0001, 0.01, "%.15g");
            snprintf(buf, sizeof(buf), "Max##%s", label);
            if (ImGui::InputDouble(buf, &p.maxVal, 0.0001, 0.01, "%.15g"))
                p.maxVal = std::max(p.minVal, p.maxVal);
            break;
        case MathEngine::InitState::Normal:
            snprintf(buf, sizeof(buf), "Mean##%s", label);
            ImGui::InputDouble(buf, &p.mean, 0.0001, 0.01, "%.15g");
            snprintf(buf, sizeof(buf), "Stddev##%s", label);
            if (ImGui::InputDouble(buf, &p.stddev, 0.0001, 0.01, "%.15g"))
                p.stddev = std::max(1e-5, p.stddev);
            break;
        case MathEngine::InitState::Cauchy:
            snprintf(buf, sizeof(buf), "Location##%s", label);
            ImGui::InputDouble(buf, &p.location, 0.0001, 0.01, "%.15g");
            snprintf(buf, sizeof(buf), "Scale##%s", label);
            if (ImGui::InputDouble(buf, &p.scale, 0.0001, 0.01, "%.15g"))
                p.scale = std::max(1e-5, p.scale);
            break;
        case MathEngine::InitState::Exponential:
            snprintf(buf, sizeof(buf), "Rate##%s", label);
            if (ImGui::InputDouble(buf, &p.rate, 0.0001, 0.01, "%.15g"))
                p.rate = std::max(1e-5, p.rate);
            break;
        case MathEngine::InitState::SplayPerturbed:
            snprintf(buf, sizeof(buf), "Perturbation##%s", label);
            ImGui::InputDouble(buf, &p.perturbation, 0.000001, 0.0001, "%.15g");
            break;
        default: break;
    }
}

inline MathEngine::dVec AppState::GenerateOAVector(const DistParams& p, size_t C, int seedOffset)
{
    const size_t seedVal = static_cast<size_t>(std::max(1, p.seed + seedOffset));
    switch (p.initState)
    {
        case MathEngine::InitState::Uniform:        return MathEngine::random_uniform(C, p.minVal, p.maxVal, seedVal);
        case MathEngine::InitState::Normal:         return MathEngine::random_normal(C, p.mean, p.stddev, seedVal);
        case MathEngine::InitState::Cauchy:         return MathEngine::random_cauchy(C, p.location, p.scale, seedVal);
        case MathEngine::InitState::Exponential:    return MathEngine::random_exponential(C, p.rate, seedVal);
        case MathEngine::InitState::Circle:         return MathEngine::random_circle<double>(C, seedVal);
        case MathEngine::InitState::Splay:          return MathEngine::splay<double>(C);
        case MathEngine::InitState::SplayPerturbed: return MathEngine::splay_perturbed(C, p.perturbation, seedVal);
        default:                                    return MathEngine::random_uniform(C, p.minVal, p.maxVal, seedVal);
    }
}

inline void AppState::DrawInitialsPanelContent()
{
    if (modelParams.modelType==ModelType::MolecularDynamics)
    {
        ImGui::SeparatorText("Initial Configuration");
        if (ImGui::Combo("Configuration",&mdParams.initialConditionIndex,mdInitNames,5))
            mdParams.initialCondition = static_cast<MDInitialConditionType>(mdParams.initialConditionIndex);
        ImGui::InputInt("Seed##MD-IC", &mdParams.seed, 1, 10);
        if (mdParams.initialCondition==MDInitialConditionType::Random)
            ImGui::InputDouble("Min Separation (sigma)", &mdParams.minSeparation, 0.01, 0.1, "%.15g");
        if (mdParams.initialCondition==MDInitialConditionType::BinaryMixture)
        {
            ImGui::InputDouble("Mass Ratio", &mdParams.massRatio, 0.1, 0.5, "%.15g");
            ImGui::InputDouble("Radius Ratio", &mdParams.radiusRatio, 0.1, 0.5, "%.15g");
        }
        ImGui::Spacing();
        ImGui::TextDisabled("Velocities: Maxwell-Boltzmann at T = %.4g (zero centre-of-mass).", mdParams.temperature);
        return;
    }

    if (modelParams.modelType==ModelType::OttAntonsen)
    {
        ImGui::SeparatorText("Initial Order Parameters");
        DrawInitDistPicker(oaRhoParams, "rho (magnitude)");
        DrawInitDistPicker(oaPhiParams, "phi (phase)");
        ImGui::InputInt("Seed##OA-IC", &oaRhoParams.seed, 1, 10);
        oaPhiParams.seed = oaRhoParams.seed;
        if (ImGui::Button("Generate Initial Order Parameters", ImVec2(-1, 0)))
        {
            const size_t C = (modelParams.oaType==OAType::OASingle) ? 1 : modelParams.oaC;
            MathEngine::dVec rho = GenerateOAVector(oaRhoParams, C);
            MathEngine::dVec phi = GenerateOAVector(oaPhiParams, C, 7);
            for (double& r : rho) r = std::clamp(std::abs(r), 0.0, 1.0);
            modelParams.oaIC.resize(2 * C);
            for (size_t c = 0; c < C; ++c)
            {
                modelParams.oaIC[2*c + 0] = rho[c] * std::cos(phi[c]);
                modelParams.oaIC[2*c + 1] = rho[c] * std::sin(phi[c]);
            }
            solverParams.solverParams.initialConditions = modelParams.oaIC;
        }
        if (ImGui::Button("View Order Parameters", ImVec2(-1, 0))) oaRhoParams.showArray = true;

        if (modelParams.oaType==OAType::OAGeneral)
        {
            ImGui::Spacing();
            ImGui::SeparatorText("Lorentzian Parameters (gamma, mu)");
            DrawInitDistPicker(oaGammaParams, "gamma (width)");
            DrawInitDistPicker(oaMuParams, "mu (center)");
            ImGui::InputInt("Seed##OA-Lor", &oaGammaParams.seed, 1, 10);
            oaMuParams.seed = oaGammaParams.seed;
            if (ImGui::Button("Generate gamma/mu", ImVec2(-1, 0)))
            {
                const size_t C = modelParams.oaC;
                modelParams.oaGammas = GenerateOAVector(oaGammaParams, C);
                for (double& g : modelParams.oaGammas) g = std::abs(g) + 1e-6;
                modelParams.oaMus = GenerateOAVector(oaMuParams, C, 1);
            }
            if (ImGui::Button("View gamma", ImVec2(-1, 0))) oaGammaParams.showArray = true;
            if (ImGui::Button("View mu", ImVec2(-1, 0))) oaMuParams.showArray = true;

            ImGui::Spacing();
            ImGui::SeparatorText("Population Fractions (eta)");
            DrawInitDistPicker(oaEtaParams, "eta (fractions)");
            ImGui::InputInt("Seed##OA-Eta", &oaEtaParams.seed, 1, 10);
            if (ImGui::Button("Generate eta", ImVec2(-1, 0)))
            {
                const size_t C = modelParams.oaC;
                modelParams.oaEta = GenerateOAVector(oaEtaParams, C);
                double s = 0.0;
                for (double& e : modelParams.oaEta) { e = std::abs(e); s += e; }
                if (s > 0.0) { for (double& e : modelParams.oaEta) e /= s; }
                else modelParams.oaEta.assign(C, 1.0 / static_cast<double>(C));
            }
            if (ImGui::Button("View eta", ImVec2(-1, 0))) oaEtaParams.showArray = true;
        }

        ImGui::Spacing();
        ImGui::Separator();
        if (ImGui::Button("Compile Model Function", ImVec2(-1, 35)))
        {
            if (modelParams.oaType==OAType::OAGeneral && modelParams.oaK.Rows() != modelParams.oaC)
            {
                modelParams.oaK = GenerateOACouplingMatrix();
            }

            const size_t C = (modelParams.oaType==OAType::OASingle) ? 1 : modelParams.oaC;

            // Ensure initial conditions exist.
            if (modelParams.oaIC.size() != 2 * C)
            {
                modelParams.oaIC.resize(2 * C);
                for (size_t c = 0; c < C; ++c)
                {
                    modelParams.oaIC[2*c + 0] = modelParams.oaRho * std::cos(modelParams.oaPhi);
                    modelParams.oaIC[2*c + 1] = modelParams.oaRho * std::sin(modelParams.oaPhi);
                }
            }
            solverParams.solverParams.initialConditions = modelParams.oaIC;

            if (modelParams.oaType==OAType::OASingle)
            {
                MathEngine::OAParams p;
                p.gamma = modelParams.oaGamma;
                p.mu    = modelParams.oaMu;
                p.K     = modelParams.K;
                solverParams.solverParams.derivative = MathEngine::OA_wrapper(p);
            }
            else
            {
                if (modelParams.oaGammas.size() != C) modelParams.oaGammas = GenerateOAVector(oaGammaParams, C);
                if (modelParams.oaMus.size()    != C) modelParams.oaMus    = GenerateOAVector(oaMuParams, C, 1);
                if (modelParams.oaEta.size()    != C) modelParams.oaEta.assign(C, 1.0 / static_cast<double>(C));
                for (double& g : modelParams.oaGammas) g = std::abs(g) + 1e-6;

                MathEngine::OAGeneralParams p;
                p.gammas = modelParams.oaGammas;
                p.mus    = modelParams.oaMus;
                p.eta    = modelParams.oaEta;
                p.K      = modelParams.oaK;
                p.C      = static_cast<int>(C);
                solverParams.solverParams.derivative = MathEngine::OAGeneral_wrapper(p);
            }
        }
        return;
    }

    ImGui::SeparatorText("Initial Phases");
    if (ImGui::Combo("Phase Dist", &phaseParams.typeIndex, dsStateNames, 8))
        phaseParams.initState = static_cast<MathEngine::InitState>(phaseParams.typeIndex);

    switch (phaseParams.initState)
    {
        case MathEngine::InitState::Uniform:
            ImGui::InputDouble("Min##P", &phaseParams.minVal, 0.0001, 0.01, "%.15g");
            if (ImGui::InputDouble("Max##P", &phaseParams.maxVal, 0.0001, 0.01, "%.15g"))
                phaseParams.maxVal=std::max(phaseParams.minVal,phaseParams.maxVal);
            break;
        case MathEngine::InitState::Normal:
            ImGui::InputDouble("Mean##P", &phaseParams.mean, 0.0001, 0.01, "%.15g");
            if (ImGui::InputDouble("Stddev##P", &phaseParams.stddev, 0.0001, 0.01, "%.15g"))
                phaseParams.stddev=std::max(1e-5,phaseParams.stddev);
            break;
        case MathEngine::InitState::Cauchy:
            ImGui::InputDouble("Location##P", &phaseParams.location, 0.0001, 0.01, "%.15g");
            if (ImGui::InputDouble("Scale##P", &phaseParams.scale, 0.0001, 0.01, "%.15g"))
                phaseParams.scale=std::max(1e-5,phaseParams.scale);
            break;
        case MathEngine::InitState::Exponential:
            if (ImGui::InputDouble("Rate##P", &phaseParams.rate, 0.0001, 0.01, "%.15g"))
                phaseParams.rate=std::max(1e-5,phaseParams.rate);
            break;
        case MathEngine::InitState::SplayPerturbed:
            ImGui::InputDouble("Perturbation##P", &phaseParams.perturbation, 0.000001, 0.0001, "%.15g");
            break;
        case MathEngine::InitState::Modules:
            if (ImGui::Combo("Module State##P", &phaseParams.moduleTypeIndex, moduleTypeNames, 7))
                phaseParams.moduleType=static_cast<MathEngine::InitType>(phaseParams.moduleTypeIndex);
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
            case MathEngine::InitState::Uniform:
                modelParams.iPhase = MathEngine::random_uniform(modelParams.N, phaseParams.minVal, phaseParams.maxVal, pSeed);
                break;
            case MathEngine::InitState::Normal:
                modelParams.iPhase = MathEngine::random_normal(modelParams.N, phaseParams.mean, phaseParams.stddev, pSeed);
                break;
            case MathEngine::InitState::Cauchy:
                modelParams.iPhase = MathEngine::random_cauchy(modelParams.N, phaseParams.location, phaseParams.scale, pSeed);
                break;
            case MathEngine::InitState::Exponential:
                modelParams.iPhase = MathEngine::random_exponential(modelParams.N, phaseParams.rate, pSeed);
                break;
            case MathEngine::InitState::Circle:
                modelParams.iPhase = MathEngine::random_circle<double>(modelParams.N, pSeed);
                break;
            case MathEngine::InitState::Splay:
                modelParams.iPhase = MathEngine::splay<double>(modelParams.N);
                break;
            case MathEngine::InitState::SplayPerturbed:
                modelParams.iPhase = MathEngine::splay_perturbed(modelParams.N, phaseParams.perturbation, pSeed);
                break;
            case MathEngine::InitState::Modules: {
                modelParams.iPhase = MathEngine::modules(modelParams.sModules,modelParams.nModules,phaseParams.moduleType,
                                                         phaseParams.param1,phaseParams.param2,pSeed,phaseParams.identical);
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
        frqncParams.initState = static_cast<MathEngine::InitState>(frqncParams.typeIndex);
    switch (frqncParams.initState)
    {
        case MathEngine::InitState::Uniform:
            ImGui::InputDouble("Min##F", &frqncParams.minVal, 0.0001, 0.01, "%.15g");
            if (ImGui::InputDouble("Max##F", &frqncParams.maxVal, 0.0001, 0.01, "%.15g"))
                frqncParams.maxVal=std::max(frqncParams.minVal,frqncParams.maxVal);
            break;
        case MathEngine::InitState::Normal:
            ImGui::InputDouble("Mean##F", &frqncParams.mean, 0.0001, 0.01, "%.15g");
            if (ImGui::InputDouble("Stddev##F", &frqncParams.stddev, 0.0001, 0.01, "%.15g"))
                frqncParams.stddev=std::max(1e-5,frqncParams.stddev);
        	break;
        case MathEngine::InitState::Cauchy:
            ImGui::InputDouble("Location##F", &frqncParams.location, 0.0001, 0.01, "%.15g");
            if (ImGui::InputDouble("Scale##F", &frqncParams.scale, 0.0001, 0.01, "%.15g"))
                frqncParams.scale=std::max(1e-5,frqncParams.scale);
        	break;
        case MathEngine::InitState::Exponential:
            if (ImGui::InputDouble("Rate##F", &frqncParams.rate, 0.0001, 0.01, "%.15g"))
                frqncParams.rate=std::max(1e-5,frqncParams.rate);
        	break;
        case MathEngine::InitState::SplayPerturbed:
            ImGui::InputDouble("Perturbation##F", &frqncParams.perturbation, 0.000001, 0.0001, "%.15g");
            break;
        case MathEngine::InitState::Modules:
            if (ImGui::Combo("Module State##F", &frqncParams.moduleTypeIndex, moduleTypeNames, 7))
                frqncParams.moduleType=static_cast<MathEngine::InitType>(frqncParams.moduleTypeIndex);
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
            case MathEngine::InitState::Uniform:
                modelParams.iFrqnc = MathEngine::random_uniform(modelParams.N, frqncParams.minVal, frqncParams.maxVal, fSeed);
                break;
            case MathEngine::InitState::Normal:
                modelParams.iFrqnc = MathEngine::random_normal(modelParams.N, frqncParams.mean, frqncParams.stddev, fSeed);
                break;
            case MathEngine::InitState::Cauchy:
                modelParams.iFrqnc = MathEngine::random_cauchy(modelParams.N, frqncParams.location, frqncParams.scale, fSeed);
                break;
            case MathEngine::InitState::Exponential:
                modelParams.iFrqnc = MathEngine::random_exponential(modelParams.N, frqncParams.rate, fSeed);
                break;
            case MathEngine::InitState::Circle:
                modelParams.iFrqnc = MathEngine::random_circle<double>(modelParams.N, fSeed);
                break;
            case MathEngine::InitState::Splay:
                modelParams.iFrqnc = MathEngine::splay<double>(modelParams.N);
                break;
            case MathEngine::InitState::SplayPerturbed:
                modelParams.iFrqnc = MathEngine::splay_perturbed(modelParams.N, frqncParams.perturbation, fSeed);
                break;
            case MathEngine::InitState::Modules: {
                modelParams.iFrqnc = MathEngine::modules(modelParams.sModules,modelParams.nModules,frqncParams.moduleType,
                                                         frqncParams.param1,frqncParams.param2,fSeed,frqncParams.identical);
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
                kParams.intra_K     = modelParams.K;
                kParams.inter_K     = modelParams.Q;
                kParams.N           = modelParams.N;
                kParams.alpha       = modelParams.alpha;
                kParams.omega       = modelParams.iFrqnc;
                kParams.module_size = modelParams.sModules;
                kParams.num_modules = modelParams.nModules;
                solverParams.solverParams.derivative = MathEngine::kuramoto_special_modular_wrapper(kParams);
            }
        }
    }
}

inline void AppState::DrawODESolverParametersPanelContent()
{
    if (modelParams.modelType==ModelType::MolecularDynamics)
    {
        ImGui::SeparatorText("Integration");
        if (ImGui::Combo("Integrator",&mdParams.integratorIndex,mdIntegratorNames,2))
            mdParams.integrator = static_cast<MDIntegratorType>(mdParams.integratorIndex);
        ImGui::InputDouble("Step Size (dt)", &mdParams.dt, 0.000001, 0.01, "%.15g");
        ImGui::InputDouble("End Time (t1)", &mdParams.t1, 0.1, 1.0, "%.15g");
        ImGui::InputInt("Stride (samples)", &mdParams.stride, 1, 10);
        if (mdParams.stride < 1) mdParams.stride = 1;

        ImGui::Spacing();
        ImGui::SeparatorText("Thermostat");
        if (ImGui::Combo("Thermostat",&mdParams.thermostatIndex,mdThermostatNames,6))
            mdParams.thermostat = static_cast<MDThermostatType>(mdParams.thermostatIndex);
        if (mdParams.thermostat!=MDThermostatType::None)
        {
            ImGui::InputDouble("Target Temperature", &mdParams.thermostatT, 0.001, 0.01, "%.15g");
            switch (mdParams.thermostat)
            {
                case MDThermostatType::Berendsen:
                case MDThermostatType::NoseHoover:
                    ImGui::InputDouble("Relaxation Time (tau)", &mdParams.thermostatTau, 0.01, 0.1, "%.15g");
                    break;
                case MDThermostatType::Andersen:
                    ImGui::InputDouble("Collision Frequency (nu)", &mdParams.andersenNu, 0.1, 1.0, "%.15g");
                    break;
                case MDThermostatType::Langevin:
                    ImGui::InputDouble("Friction (gamma)", &mdParams.langevinGamma, 0.1, 1.0, "%.15g");
                    break;
                default: break;
            }
        }

        ImGui::Spacing();
        ImGui::SeparatorText("Barostat (NPT)");
        ImGui::Checkbox("Enable Barostat", &mdParams.barostat);
        if (mdParams.barostat)
        {
            ImGui::InputDouble("Target Pressure", &mdParams.targetPressure, 0.001, 0.01, "%.15g");
            ImGui::InputDouble("Relaxation Time (tauP)", &mdParams.barostatTau, 0.1, 1.0, "%.15g");
        }
        return;
    }

	if (ImGui::Combo("Solver Method",&solverParams.solverMethodSelectedIndex,solverMethodNames,9))
        solverParams.solverMethod=static_cast<SolverMethod>(solverParams.solverMethodSelectedIndex);
    switch(solverParams.solverMethod)
    {
        case SolverMethod::RK1:
            solverParams.solverFunc = rk1_wrapper();
            break;
        case SolverMethod::RK2:
            solverParams.solverFunc = rk2_wrapper();
            break;
        case SolverMethod::RK3:
	        solverParams.solverFunc = rk3_wrapper();
            break;
        case SolverMethod::RK4:
            solverParams.solverFunc = rk4_wrapper();
            break;
        case SolverMethod::RK4_38:
            solverParams.solverFunc = rk4_38_wrapper();
            break;
        case SolverMethod::RK4_Gill:
            solverParams.solverFunc = rk4_gill_wrapper();
            break;
        case SolverMethod::RK4_Ralston:
            solverParams.solverFunc = rk4_ralston_wrapper();
            break;
        case SolverMethod::AB:
            solverParams.solverFunc = adams_bashforth_wrapper();
            break;
        case SolverMethod::ABM:
            solverParams.solverFunc = adams_bashforth_moulton_wrapper();
            break;
        default:
            solverParams.solverFunc = rk4_wrapper();
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
        if (solverParams.solverMethod==SolverMethod::AB || solverParams.solverMethod==SolverMethod::ABM)
        {
            ImGui::SliderScalar("Method Order", ImGuiDataType_U8, &solverParams.solverParams.order, &MultiStepOrderMin, &MultiStepOrderMax, "%u");
            if (solverParams.solverMethod==SolverMethod::ABM) ImGui::SliderScalar("ABM PE(CE) Iterations", ImGuiDataType_U8, &solverParams.solverParams.iterations, &MultiStepIterationsMin, &MultiStepIterationsMax, "%u");
        }
    }
    // if (ImGui::CollapsingHeader("Adaptive Step Control"))
    // {
    //     ImGui::Checkbox("Estimate Error", &solverParams.solverParams.errorEstimate);
    //     ImGui::SameLine();
    //     ImGui::Checkbox("Enable Variable Step Size", &solverParams.solverParams.variableSteps);
    //     if (solverParams.solverParams.errorEstimate && solverParams.solverParams.variableSteps)
    //     {
    //         ImGui::Indent();
    //         ImGui::TextDisabled("Tolerances & Bounds");
    //         ImGui::InputDouble("Local Tolerance", &solverParams.solverParams.localTol, 0.0, 0.0, "%.1e");
    //         ImGui::InputDouble("Absolute Tolerance", &solverParams.solverParams.absolute_tol, 0.0, 0.0, "%.1e");
    //         ImGui::InputDouble("Min dt", &solverParams.solverParams.minDt, 0.0, 0.0, "%.15g");
    //         ImGui::InputDouble("Max dt", &solverParams.solverParams.maxDt, 0.0, 0.0, "%.15g");
    //         ImGui::Separator();
    //         ImGui::TextDisabled("Step Adaptation Factors");
    //         ImGui::InputDouble("Decrease Factor", &solverParams.solverParams.decreaseFactor, 0.05, 0.1, "%.15g");
    //         ImGui::InputDouble("Increase Factor", &solverParams.solverParams.increaseFactor, 0.1, 0.5, "%.15g");
    //         ImGui::InputDouble("Tol Error Ratio", &solverParams.solverParams.localTolErrorRatio, 0.01, 0.05, "%.15g");
    //         // size_t cast to int for ImGui input
    //         int maxTrial = static_cast<int>(solverParams.solverParams.maxTrial);
    //         if (ImGui::InputInt("Max Trials", &maxTrial))
    //             solverParams.solverParams.maxTrial = static_cast<size_t>(maxTrial>1?maxTrial:1);
    //         ImGui::Unindent();
    //     }
    // }
    // if (ImGui::CollapsingHeader("Error Metrics & Flags"))
    // {
    //     ImGui::Checkbox("Weighted Error Formula", &solverParams.solverParams.weightedError);
    //     ImGui::Checkbox("Norm Error Formula", &solverParams.solverParams.normError);
    //     ImGui::Checkbox("Record Attempt History", &solverParams.solverParams.attemptsHistory);
    // }
    // // double tau=0.0;
    // if (ImGui::CollapsingHeader("Delay Differential Equations (DDE)"))
    // {
    //     ImGui::Checkbox("Is DDE System", &solverParams.solverParams.isDDE);
    //     if (solverParams.solverParams.isDDE)
    //     {
    //         ImGui::Indent();
    //         int maxDelayOrder_ = static_cast<int>(solverParams.solverParams.maxDelayOrder);
    //         if (ImGui::InputInt("Max Delay Order", &maxDelayOrder_))
    //         {
    //             solverParams.solverParams.maxDelayOrder = static_cast<size_t>(maxDelayOrder_ > 1 ? maxDelayOrder_ : 1);
    //         }
    //         ImGui::InputDouble("Interpolation Tol", &solverParams.solverParams.interpolationTol, 1e-10, 1e-8, "%.1e");
    //         ImGui::InputDouble("dt Scale (Fine Step)", &solverParams.solverParams.dtScale, 0.01, 0.05, "%.15g");
    //         ImGui::SliderInt("Number of Delays",&solverParams.nDs,1,20);
    //         ImGui::SameLine();
    //         if (ImGui::Button("Set Delay Count"))
    //         {
    //             delayTimes.resize(solverParams.nDs);
    //         }
    //         if (ImGui::CollapsingHeader("Delays"))
    //         {
    //             // delayTimes = MathEngine::dVec(nDs,0.0);
    //             for (size_t i=0; i<delayTimes.size(); ++i)
    //             {
    //                 std::string label = "Delay number "+std::to_string(i+1);
    //                 ImGui::InputDouble(label.c_str(), &delayTimes[i],0.0001f,0.01f,"%.15g");
    //                 // delayTimes[i]=tau;
    //             }
    //         }
    //         if (ImGui::Button("Submit Delays"))
    //         {
    //             solverParams.solverParams.delayTimes=delayTimes;
    //         }
    //         ImGui::SameLine();
    //         if (ImGui::Button("View Delay Values"))
    //         {
    //             showDelays=true;
    //         }
    //         ImGui::Text("Configured Delays: %zu", solverParams.solverParams.delayTimes.size());
    //         ImGui::Unindent();
    //     }
    // }
}

inline void AppState::DrawVectorViewer(bool& open, const char* title, const MathEngine::dVec& data, const char* emptyText)
{
    const ImGuiViewport* viewport = ImGui::GetMainViewport();
    ImVec2 center = viewport->GetCenter();
    ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    ImGui::SetNextWindowSize(ImVec2(500, 480), ImGuiCond_Appearing);
    if (ImGui::Begin(title, &open))
    {
        if (ImGui::BeginChild("VecList", ImVec2(0, 330), ImGuiChildFlags_Borders))
        {
            if (data.empty())
                ImGui::TextDisabled("%s", emptyText);
            else
                for (size_t i = 0; i < data.size(); ++i)
                    ImGui::Text("[%03zu]  %.15g", i + 1, data[i]);
        }
        ImGui::EndChild();
    }
    ImGui::End();
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
        const MathEngine::dMatrix& m = (modelParams.modelType==ModelType::OttAntonsen) ? modelParams.oaK : adj;
        const char* viewerTitle = (modelParams.modelType==ModelType::OttAntonsen) ? "Coupling Matrix Viewer" : "Adjacency Matrix Viewer";
        if (ImGui::Begin(viewerTitle, &adjParams.showAdjMatrix))
        {
            size_t nRows = m.Rows();
            size_t nCols = m.Cols();
            ImGui::Text("Dimension: %zu x %zu", nRows, nCols);
            ImGui::Separator();

            if (nRows == 0 || nCols == 0)
            {
                ImGui::TextDisabled("Matrix is empty. Generate it first.");
            }
            else
            {
                if (nCols>500 || nRows>500)
                {
                    const float cellHeight = cellHeightBase+5.0f;
                    const float cellWidth = cellWidthBase+10.0f;
                    const ImVec2 totalCanvasSize = ImVec2((nCols+1)*cellWidth,(nRows+1)*cellHeight);
                    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding,ImVec2(0,0));
                    ImGuiWindowFlags windowFlags = ImGuiWindowFlags_HorizontalScrollbar | ImGuiWindowFlags_NoCollapse;
                    if (ImGui::BeginChild("Matrix Grid",ImVec2(0,330), true, windowFlags))
                    {
                        const float scrollX = ImGui::GetScrollX();
                        const float scrollY = ImGui::GetScrollY();
                        ImVec2 winSize = ImGui::GetWindowSize();
                        int minCol = std::max(0,static_cast<int>(scrollX/cellWidth)-1);
                        int maxCol = std::min(static_cast<int>(nCols),static_cast<int>((scrollX+winSize.x)/cellWidth)+1);
                        int minRow = std::max(0,static_cast<int>(scrollY/cellHeight)-1);
                        int maxRow = std::min(static_cast<int>(nRows),static_cast<int>((scrollY+winSize.y)/cellHeight)+1);
                        ImGui::SetCursorPos(totalCanvasSize);
                        ImGui::SetCursorPos(ImVec2(5.0f,0.0f));
                        ImGui::TextColored(ImVec4(1.0,1.0,1.0,1.0), "Row\\Col");
                        for (int c=minCol; c<maxCol; ++c)
                        {
                            ImGui::SetCursorPos(ImVec2((c+1)*cellWidth+5.0f,0.0f));
                            ImGui::TextColored(ImVec4(0.9,0.9,0.9,1.0), "[%03d]",c+1);
                        }
                        for (int r=minRow; r<maxRow; ++r)
                        {
                            ImGui::SetCursorPos(ImVec2(5.0f,(r+1)*cellHeight));
                            ImGui::TextColored(ImVec4(0.0,1.0,1.0,1.0), "[%03d]",r+1);
                            for (int c=minCol; c<maxCol; ++c)
                            {
                                ImGui::SetCursorPos(ImVec2((c+1)*cellWidth+5.0f,(r+1)*cellHeight));
                                double val = m[r][c];
                                if (val >= 1e-5) ImGui::TextColored(ImVec4(0.4f, 0.9f, 1.0f, 1.0f), "%.15g", val);
                                else ImGui::TextDisabled("0.0000");
                            }
                        }
                    }
                    ImGui::EndChild();
                    ImGui::PopStyleVar();
                }
                else
                {
                    ImGuiTableFlags tableFlags = ImGuiTableFlags_ScrollX | ImGuiTableFlags_ScrollY | ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingFixedFit;
                    if (ImGui::BeginTable("MatrixGrid", static_cast<int>(nCols + 1), tableFlags, ImVec2(0, 400)))
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
                            ImGui::TextColored(ImVec4(0.0,1.0,1.0,1.0),"[%03zu]", r+1);
                            for (size_t c = 0; c < nCols; ++c)
                            {
                                ImGui::TableSetColumnIndex(static_cast<int>(c + 1));
                                double val = m[r][c];
                                if (val >= 1e-5) ImGui::Text("%.15g", val);
                                else ImGui::TextDisabled("0.0000");
                            }
                        }
                        ImGui::EndTable();
                    }
                }
            }
        }
        ImGui::End();
    }

    if (phaseParams.showArray)
    {
        ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
        ImGui::SetNextWindowSize(ImVec2(500, 480), ImGuiCond_Appearing);
        if (ImGui::Begin("Phase Array Values", &phaseParams.showArray))
        {
            if (ImGui::BeginChild("PhaseList", ImVec2(0, 330), ImGuiChildFlags_Borders))
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
        ImGui::SetNextWindowSize(ImVec2(500, 480), ImGuiCond_Appearing);
        if (ImGui::Begin("Frequency Array Values", &frqncParams.showArray))
        {
            if (ImGui::BeginChild("FreqList", ImVec2(0, 330), ImGuiChildFlags_Borders))
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

    // ---- Ott-Antonsen view modals ----
    if (oaRhoParams.showArray)
    {
        ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
        ImGui::SetNextWindowSize(ImVec2(500, 480), ImGuiCond_Appearing);
        if (ImGui::Begin("Initial Order Parameters", &oaRhoParams.showArray))
        {
            if (ImGui::BeginChild("OAList", ImVec2(0, 330), ImGuiChildFlags_Borders))
            {
                if (modelParams.oaIC.empty())
                    ImGui::TextDisabled("Array is empty. Click Generate Initial Order Parameters.");
                else
                {
                    const size_t C = modelParams.oaIC.size() / 2;
                    for (size_t c = 0; c < C; ++c)
                    {
                        const double x = modelParams.oaIC[2*c + 0];
                        const double y = modelParams.oaIC[2*c + 1];
                        ImGui::Text("[%03zu]  rho=%.15g  phi=%.15g", c + 1, std::hypot(x, y), std::atan2(y, x));
                    }
                }
            }
            ImGui::EndChild();
        }
        ImGui::End();
    }
    if (oaGammaParams.showArray)
        DrawVectorViewer(oaGammaParams.showArray, "Lorentzian Width (gamma)", modelParams.oaGammas, "Array is empty. Click Generate gamma/mu.");
    if (oaMuParams.showArray)
        DrawVectorViewer(oaMuParams.showArray, "Lorentzian Center (mu)", modelParams.oaMus, "Array is empty. Click Generate gamma/mu.");
    if (oaEtaParams.showArray)
        DrawVectorViewer(oaEtaParams.showArray, "Population Fractions (eta)", modelParams.oaEta, "Array is empty. Click Generate eta.");

    // if (showDelays)
    // {
    //     ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    //     ImGui::SetNextWindowSize(ImVec2(500, 480), ImGuiCond_Appearing);
    //     if (ImGui::Begin("Submitted Delays",&showDelays))
    //     {
    //         ImGui::Text("Total Size: %zu elements", solverParams.solverParams.delayTimes.size());
    //         ImGui::Separator();
    //         ImGui::Spacing();

    //         // Scrollable child box for array elements
    //         if (ImGui::BeginChild("ArrayList", ImVec2(0, 330), ImGuiChildFlags_Borders))
    //         {
    //             if (solverParams.solverParams.delayTimes.empty())
    //             {
    //                 ImGui::TextDisabled("Array is empty.");
    //             }
    //             else
    //             {
    //                 for (size_t i = 0; i < solverParams.solverParams.delayTimes.size(); ++i)
    //                 {
    //                     ImGui::Text("[%03zu]  %.6f", i+1, solverParams.solverParams.delayTimes[i]);
    //                 }
    //             }
    //         }
    //         ImGui::EndChild();

    //         ImGui::Spacing();
    //         ImGui::Separator();
    //         ImGui::Spacing();
    //     }
    //     ImGui::End();
    // }
}

inline void AppState::StartMolecularDynamics()
{
    MathEngine::MDConfig cfg;
    cfg.numParticles = mdParams.numParticles;
    cfg.width  = mdParams.width;
    cfg.height = mdParams.height;
    cfg.mass         = mdParams.mass;
    cfg.radius       = mdParams.radius;
    cfg.massRatio    = mdParams.massRatio;
    cfg.radiusRatio  = mdParams.radiusRatio;
    cfg.sigma        = mdParams.sigma;
    cfg.epsilon      = mdParams.epsilon;
    cfg.cutoffCoeff  = mdParams.cutoffCoeff;
    cfg.morseAlpha   = mdParams.morseAlpha;
    cfg.temperature  = mdParams.temperature;
    cfg.restitution  = mdParams.restitution;
    cfg.minSeparation = mdParams.minSeparation;
    cfg.seed         = static_cast<size_t>(std::max(1, mdParams.seed));
    cfg.periodicBoundaryCondition = mdParams.periodicBoundaryCondition;
    cfg.bounce         = mdParams.bounce;
    cfg.hardSphereCollisions = mdParams.hardSphereCollisions;
    cfg.potential        = static_cast<MathEngine::PotentialType>(mdParams.potential);
    cfg.initialCondition = static_cast<MathEngine::InitialConditionType>(mdParams.initialCondition);

    const double dt = mdParams.dt;
    const double t1 = mdParams.t1;
    const int stride = std::max(1, mdParams.stride);
    const size_t numSteps = static_cast<size_t>(std::llround(t1 / dt));
    if (numSteps == 0) { isSimRunning.store(false); return; }

    timeInv.store(static_cast<float>(1.0 / std::abs(t1)));
    simProgress.store(0.0f);

    {
        std::lock_guard<std::mutex> lock(plotParams.plotMutex);
        plotParams.plotX.clear();
        plotParams.plotY.clear();
        plotParams.plotXTrail.clear();
        plotParams.plotYTrail.clear();
        plotParams.liveState.clear();
        mdRunState.clear();
    }

    const auto mdIntegrator = mdParams.integrator;
    const auto mdThermostat = mdParams.thermostat;
    const double thermostatT   = mdParams.thermostatT;
    const double thermostatTau = mdParams.thermostatTau;
    const double langevinGamma  = mdParams.langevinGamma;
    const double andersenNu     = mdParams.andersenNu;
    const bool   barostat       = mdParams.barostat;
    const double targetPressure = mdParams.targetPressure;
    const double barostatTau    = mdParams.barostatTau;

    auto runMD = [this, cfg, dt, stride, numSteps, mdIntegrator, mdThermostat,
                  thermostatT, thermostatTau, langevinGamma, andersenNu,
                  barostat, targetPressure, barostatTau]()
    {
        try
        {
            MathEngine::MolecularDynamics md(cfg);
            size_t rngSeed = cfg.seed;

            auto record = [&](double t)
            {
                const MathEngine::Observables o = MathEngine::CollectObservables(md, t, 1.4);
                std::lock_guard<std::mutex> lock(plotParams.plotMutex);
                plotParams.plotX.push_back(o.time);
                plotParams.plotY.push_back(o.psi6);
                mdRunState.posX = md.posX; mdRunState.posY = md.posY;
                mdRunState.velX = md.velX; mdRunState.velY = md.velY;
                mdRunState.time.push_back(o.time);
                mdRunState.temperature.push_back(o.temperature);
                mdRunState.kineticEnergy.push_back(o.kineticEnergy);
                mdRunState.potentialEnergy.push_back(o.potentialEnergy);
                mdRunState.totalEnergy.push_back(o.totalEnergy);
                mdRunState.pressure.push_back(o.pressure);
                mdRunState.psi4.push_back(o.psi4);
                mdRunState.psi6.push_back(o.psi6);
                mdRunState.msd.push_back(o.msd);
            };

            record(0.0);
            for (size_t s = 0; s < numSteps; ++s)
            {
                if (mdThermostat == MDThermostatType::NoseHoover)
                {
                    MathEngine::StepNoseHoover(md, dt, thermostatT, thermostatTau);
                }
                else
                {
                    if (mdIntegrator == MDIntegratorType::VelocityVerlet) md.Step(dt);
                    else md.StepLeapfrog(dt);

                    switch (mdThermostat)
                    {
                        case MDThermostatType::Rescale:   MathEngine::ApplyVelocityRescale(md, thermostatT); break;
                        case MDThermostatType::Berendsen: MathEngine::ApplyBerendsen(md, thermostatT, dt, thermostatTau); break;
                        case MDThermostatType::Andersen:  MathEngine::ApplyAndersen(md, thermostatT, dt, andersenNu, rngSeed); break;
                        case MDThermostatType::Langevin:  MathEngine::ApplyLangevin(md, thermostatT, dt, langevinGamma, rngSeed); break;
                        default: break;
                    }
                }
                if (barostat) MathEngine::ApplyBerendsenBarostat(md, targetPressure, dt, barostatTau);

                if ((s + 1) % static_cast<size_t>(stride) == 0 || s + 1 == numSteps)
                    record(static_cast<double>(s + 1) * dt);

                simProgress.store(static_cast<float>(static_cast<double>(s + 1) / static_cast<double>(numSteps)));
            }
        }
        catch (...) {}
        simProgress.store(1.0f);
        isSimRunning.store(false);
    };

    #ifndef __EMSCRIPTEN__
    simThread = std::thread(runMD);
    #else
    runMD();  // Web: synchronous run.
    #endif
    hasSimRan = true;
}

inline void AppState::StartSimulation()
{
    // Do not start a new simulation while one is still running.
    if (isSimRunning.load()) return;

    // Reap any previous (finished) simulation thread before launching a new one.
    if (simThread.joinable()) simThread.join();

    processStartTime = std::chrono::steady_clock::now();
    simProgress.store(0.0f);
    isSimRunning.store(true);
    timeInv.store(static_cast<float>(1.0/std::abs(solverParams.solverParams.t1-solverParams.solverParams.t0)));
    if (modelParams.modelType==ModelType::MolecularDynamics)
    {
        StartMolecularDynamics();
        return;
    }
    // constexpr size_t Stride = 25;
    const size_t vectorSize = static_cast<size_t>((solverParams.solverParams.t1-solverParams.solverParams.t0)/solverParams.solverParams.dt);
    const size_t plotExpectedSize = static_cast<size_t>(vectorSize/plotParams.Stride)*2+100;
    const bool isOA = (modelParams.modelType==ModelType::OttAntonsen);
    const bool isOAGeneral = isOA && (modelParams.oaType==OAType::OAGeneral);
    const size_t oaC = modelParams.oaC;
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
        plotParams.plotYModules = MathEngine::dMatrix(isOAGeneral ? oaC : modelParams.nModules, 0);
        plotParams.offset = 0;
    }
    int stride = plotParams.Stride;
    bool condPlotThird = isOAGeneral ||
        (modelParams.kuramotoType==KuramotoType::KuramotoSpecial || adjParams.adjState==MathEngine::NetworkTopology::Modular ||
         adjParams.adjState==MathEngine::NetworkTopology::Hierarchical);
    solverParams.solverParams.onStep = [this, isOA, isOAGeneral, oaC, stepCount=0, stepCountCond=0, stride, condPlotThird](const MathEngine::OneStepSolverResult& res) mutable
    {
		float progress = (res.timePoint-solverParams.solverParams.t0)*timeInv;
        std::lock_guard<std::mutex> lock(plotParams.plotMutex);
        simProgress.store(progress);
        double rSine = 0.0, rCosine = 0.0, rho = 0.0;
        MathEngine::dVec rhoM;
        plotParams.liveTimePoints.push_back(res.timePoint);
        plotParams.liveState = res.sol;
        if (isOA)
        {
            // OA: the state is already the order parameters (interleaved Re/Im).
            const size_t C = isOAGeneral ? oaC : 1;
            double sre = 0.0, sim = 0.0;
            rhoM.assign(C, 0.0);
            for (size_t c = 0; c < C; ++c)
            {
                const double x = res.sol[2*c + 0];
                const double y = res.sol[2*c + 1];
                rhoM[c] = std::hypot(x, y);
                const double w = (modelParams.oaEta.size() == C) ? modelParams.oaEta[c] : (1.0 / static_cast<double>(C));
                sre += w * x;
                sim += w * y;
            }
            rho = std::hypot(sre, sim);
            if (isOAGeneral && ++stepCountCond%stride==0)
            {
                plotParams.plotYModules.AppendCols(rhoM);
            }
        }
        else if (condPlotThird)
        {
            MathEngine::dVec rMSine(modelParams.nModules,0.0);
            MathEngine::dVec rMCosine(modelParams.nModules,0.0);
            rhoM.assign(modelParams.nModules, 0.0);
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
                plotParams.plotYModules.AppendCols(rhoM);
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
        else if (plotParams.trailCount>0)
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
#ifndef __EMSCRIPTEN__
    simThread = std::thread([this]()
    {
        try
        {
            solverParams.solverResults = solverParams.solverFunc(solverParams.solverParams);
        }
        catch (...)
        {
            // Keep the UI from getting stuck in "Running..." if the solver throws.
        }
        simProgress.store(1.0);
        isSimRunning.store(false);
    });
#else
    // Web: sync run (no pthreads).
    solverParams.solverResults = solverParams.solverFunc(solverParams.solverParams);
    simProgress.store(1.0);
    isSimRunning.store(false);
#endif
    hasSimRan = true;
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
    {
        if (plotParams.showPlot)
        {
            const ImGuiViewport* viewport = ImGui::GetMainViewport();
            ImVec2 center = viewport->GetCenter();
            ImGui::SetNextWindowPos(center, ImGuiCond_FirstUseEver, ImVec2(0.5f, 0.5f));
            ImGui::SetNextWindowSize(ImVec2(400, 480), ImGuiCond_FirstUseEver);
            if (ImGui::Begin("Order Parameter",&plotParams.showPlot))
            {
                ImVec2 availableSpace = ImGui::GetContentRegionAvail();      // Set ImVec(-1,-1) to fill the whole window.
                if (ImPlot::BeginSubplots("##Plot-now", 2, 1, availableSpace))
                {
                    if (ImPlot::BeginPlot("(\U0001D70C-t) Plot"))
                    {
                        ImPlotSpec spec;
                        std::pair<double,double> xmm;
                        std::pair<double,double> ymm;
                        if (plotParams.plotX.empty())
                        {
                            xmm = {0.0,1.0};
                            ymm = {0.0,1.0};
                        }
                        else
                        {
                            auto [xmin,xmax] = std::minmax_element(plotParams.plotX.begin(),plotParams.plotX.end());
                            xmm = {*xmin-0.01,*xmax+0.01};
                            auto [ymin,ymax] = std::minmax_element(plotParams.plotY.begin(),plotParams.plotY.end());
                            ymm = {*ymin-0.01,*ymax+0.01};
                        }
                        ImPlot::SetupAxesLimits(xmm.first,xmm.second,ymm.first,ymm.second,ImPlotCond_Always);
                        if (!plotParams.plotColors.empty())
	                        spec.LineColor = plotParams.plotColors[0];
                        ImPlot::SetupAxes("Time (t)","Order (\U0001D70C)");
                        ImPlot::PlotLine("\U0001D70C",plotParams.plotX.data(),plotParams.plotY.data(),static_cast<int>(plotParams.plotX.size()));
                        ImPlot::EndPlot();
                    }
                    if (ImPlot::BeginPlot("(\U0001D70C-t) Plot##trailing"))
                    {
                        std::pair<double,double> xmm;
                        std::pair<double,double> ymm;
                        if (plotParams.plotXTrail.empty())
                        {
                            xmm = {0.0,1.0};
                            ymm = {0.0,1.0};
                        }
                        else
                        {
                            auto [xmin,xmax] = std::minmax_element(plotParams.plotXTrail.begin(),plotParams.plotXTrail.end());
                            xmm = {*xmin-0.01,*xmax+0.01};
                            auto [ymin,ymax] = std::minmax_element(plotParams.plotYTrail.begin(),plotParams.plotYTrail.end());
                            ymm = {*ymin-0.01,*ymax+0.01};
                        }
                        ImPlot::SetupAxesLimits(xmm.first,xmm.second,ymm.first,ymm.second,ImPlotCond_Always);
                        ImPlotSpec spec;
                        spec.Offset = static_cast<int>(plotParams.offset);
                        if (!plotParams.plotColors.empty())
                            spec.LineColor = plotParams.plotColors[0];
                        ImPlot::SetupAxes("Time (t)","Order (\U0001D70C)");
                        ImPlot::PlotLine("\U0001D70C",plotParams.plotXTrail.data(),plotParams.plotYTrail.data(),static_cast<int>(plotParams.plotXTrail.size()),spec);
                        ImPlot::EndPlot();
                    }
                    ImPlot::EndSubplots();
                }
            }
            ImGui::End();
        }
        if (modelParams.modelType==ModelType::MolecularDynamics && plotParams.showPlotSecond)
        {
            const ImGuiViewport* viewport = ImGui::GetMainViewport();
            ImVec2 center = viewport->GetCenter();
            ImGui::SetNextWindowPos(center, ImGuiCond_FirstUseEver, ImVec2(0.5f, 0.5f));
            ImGui::SetNextWindowSize(ImVec2(500, 460), ImGuiCond_FirstUseEver);
            if (ImGui::Begin("Particles",&plotParams.showPlotSecond))
            {
                ImVec2 availableSpace = ImGui::GetContentRegionAvail();
                if (ImPlot::BeginPlot("Particle Positions",availableSpace))
                {
                    const double W = mdParams.width, H = mdParams.height;
                    ImPlot::SetupAxes("x","y");
                    ImPlot::SetupAxesLimits(0.0, std::max(1.0, W), 0.0, std::max(1.0, H), ImPlotCond_Always);
                    ImPlotSpec spec;
                    if (!plotParams.plotSecondColors.empty())
                        spec.LineColor = plotParams.plotSecondColors[0];
                    ImPlot::PlotScatter("particles", mdRunState.posX.data(), mdRunState.posY.data(),
                                        static_cast<int>(mdRunState.posX.size()), spec);
                    ImPlot::EndPlot();
                }
            }
            ImGui::End();
        }
        else if (modelParams.modelType==ModelType::Kuramoto && plotParams.showPlotSecond)
        {
            const ImGuiViewport* viewport = ImGui::GetMainViewport();
            ImVec2 center = viewport->GetCenter();
            ImGui::SetNextWindowPos(center, ImGuiCond_FirstUseEver, ImVec2(0.5f, 0.5f));
            ImGui::SetNextWindowSize(ImVec2(400, 250), ImGuiCond_FirstUseEver);
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
                    ImPlot::SetupAxesLimits(-1.01,1.01,-1.01,1.01,ImPlotCond_Always);
                    ImPlotSpec spec;
                    if (!plotParams.plotSecondColors.empty())
	                    spec.LineColor = plotParams.plotSecondColors[0];
                    ImPlot::PlotScatter("\U0001D73D",xTheta.data(),yTheta.data(),static_cast<int>(plotN), spec);
                    ImPlot::EndPlot();
                }
            }
            ImGui::End();
        }
        if (plotParams.showPlotThird)
        {
            const size_t nM = (modelParams.modelType==ModelType::OttAntonsen) ? modelParams.oaC : modelParams.nModules;
            const ImGuiViewport* viewport = ImGui::GetMainViewport();
            ImVec2 center = viewport->GetCenter();
            ImGui::SetNextWindowPos(center, ImGuiCond_FirstUseEver, ImVec2(0.5f, 0.5f));
            ImGui::SetNextWindowSize(ImVec2(400, 200), ImGuiCond_FirstUseEver);
            if (ImGui::Begin("Order Parameter (Modules)##mocules",&plotParams.showPlotThird))
            {
                ImVec2 availableSpace = ImGui::GetContentRegionAvail();      // Set ImVec(-1,-1) to fill the whole window.
                if (plotParams.plotYModules.empty() || plotParams.plotYModules.Rows()!=nM)
                {
                    if (ImPlot::BeginPlot("(\U0001D70C-t) Plot (empty modules)##modules",availableSpace))
                    {
                        ImPlot::EndPlot();
                    }
                }
                else
                {
                    if (ImPlot::BeginPlot("(\U0001D70C-t) Plot##modules",availableSpace))
                    {
                        std::pair<double,double> xmm;
                        std::pair<double,double> ymm;
                        if (plotParams.plotX.empty())
                        {
                            xmm = {0.0,1.0};
                            ymm = {0.0,1.0};
                        }
                        else
                        {
                            auto [xmin,xmax] = std::minmax_element(plotParams.plotX.begin(),plotParams.plotX.end());
                            xmm = {*xmin-0.01,*xmax+0.01};
                            MathEngine::dVec ymms;
                            for (size_t i=0; i<nM; ++i)
                            {
                                auto [ymin_,ymax_] = std::minmax_element(plotParams.plotYModules[i].begin(),plotParams.plotYModules[i].end());
                                ymms.push_back(*ymin_); ymms.push_back(*ymax_);
                            }
                            auto [ymin,ymax] = std::minmax_element(ymms.begin(),ymms.end());
                            ymm = {*ymin-0.01,*ymax+0.01};
                        }
                        ImPlot::SetupAxesLimits(xmm.first,xmm.second,ymm.first,ymm.second,ImPlotCond_Always);
                        ImPlot::SetupAxes("Time (t)","Order (\U0001D70C)");
                        std::string label = "\U0001D70C";
                        for (size_t i=0; i<nM; ++i)
                        {
                            label = "\U0001D70C "+std::to_string(i);
                            ImPlotSpec spec;
                            if (plotParams.plotThirdColors.size()>i)
                                spec.LineColor = plotParams.plotThirdColors[i];
                            ImPlot::PlotLine(label.c_str(),plotParams.plotX.data(),plotParams.plotYModules[i].data(),
                                            static_cast<int>(plotParams.plotX.size()),spec);
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
    ImGui::SeparatorText("Plot Data Style");
    if (ImGui::CollapsingHeader("\U0001D70C-t Plot##main plot"))
    {
        if (ImGui::InputInt("Stride##main plot",&plotParams.Stride,1,10)) plotParams.Stride = std::max(plotParams.Stride,10);
        if (ImGui::InputInt("Trailing Data Count##main plot",&plotParams.trailCount,1,10)) plotParams.trailCount = std::clamp(plotParams.trailCount,100,10000);
        ImGui::Checkbox("Show \U0001D70C-t Plot##main plot", &plotParams.showPlot);
        ImGui::Spacing();
        ImGui::SeparatorText("Line Color##main plot");
        ImGui::Spacing();
        if (plotParams.plotColors.empty())
            plotParams.plotColors.push_back(ImVec4(0.2f,0.8f,0.8f,1.0f));
        ImGui::ColorEdit4("\U0001D70C-t Colors (main)",&plotParams.plotColors[0].x);
    }
    if (modelParams.modelType==ModelType::Kuramoto || modelParams.modelType==ModelType::MolecularDynamics)
    {
        const bool isMD = (modelParams.modelType==ModelType::MolecularDynamics);
        if (ImGui::CollapsingHeader(isMD ? "Particles Plot##second plot" : "\U0001D73D Plot##second plot"))
        {
            ImGui::Checkbox(isMD ? "Show Particles Plot" : "Show \U0001D73D Plot", &plotParams.showPlotSecond);
            ImGui::Spacing();
            ImGui::SeparatorText("Line Color##second plot");
            ImGui::Spacing();
            if (plotParams.plotSecondColors.empty())
                plotParams.plotSecondColors.push_back(ImVec4(0.2f,0.8f,0.8f,1.0f));
            ImGui::ColorEdit4(isMD ? "Particle Color" : "\U0001D73D Colors",&plotParams.plotSecondColors[0].x);
        }
    }
    const bool showModules = (modelParams.modelType==ModelType::OttAntonsen && modelParams.oaType==OAType::OAGeneral)
        || adjParams.adjState==MathEngine::NetworkTopology::Modular
        || adjParams.adjState==MathEngine::NetworkTopology::Hierarchical
        || modelParams.kuramotoType==KuramotoType::KuramotoSpecial;
    if (showModules)
    {
        if (ImGui::CollapsingHeader("\U0001D73D Plot (Modules)##third plot"))
        {
            ImGui::Checkbox("Show \U0001D73D Plot##third plot", &plotParams.showPlotThird);
            ImGui::Spacing();
            ImGui::SeparatorText("Line Color##third plot");
            ImGui::Spacing();
            size_t nM = (modelParams.modelType==ModelType::OttAntonsen) ? modelParams.oaC : modelParams.nModules;
            if (plotParams.plotThirdColors.size()!=nM)
                plotParams.plotThirdColors.resize(nM,ImVec4(0.2f,0.8f,0.8f,1.0f));
            for (size_t i=0; i<nM; ++i)
            {
                std::string label = "\U0001D70C-t ("+std::to_string(i+1)+")##third plot line";
                ImGui::ColorEdit4(label.c_str(),&plotParams.plotThirdColors[i].x);
            }
        }
    }
}

inline void AppState::RenderChrono()
{
    bool running = isSimRunning.load();
    if (!hasSimRan) { ImGui::Text("00:00:00.000"); return; }
    if (running)
    {
        processStopTime = std::chrono::steady_clock::now();
        processDuration = std::chrono::duration_cast<std::chrono::milliseconds>(processStopTime-processStartTime).count();
        milliSec = processDuration % 1000;
        totalSec = processDuration / 1000;
        sec      = totalSec % 60;
        min      = (totalSec/60) % 60;
        hour     = totalSec / 3600;
    }
    ImGui::Text("%02ld:%02ld:%02ld.%03ld",hour,min,sec,milliSec);
}

