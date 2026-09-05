cmake_minimum_required(VERSION 3.16)
project(EMMA LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(POLICY CMP0169)
	cmake_policy(SET CMP0169 OLD)
endif()

include(FetchContent)

option(ENABLE_NATIVE_ARCH "Build with -march=native for the local CPU" ON)

set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
# set(USE_WAYLAND ON CACHE BOOL "Use Wayland" FORCE)
# set(USE_X11 OFF CACHE BOOL "Use X11" FORCE)
# set(GLFW_BUILD_WAYLAND ON CACHE BOOL "" FORCE)
# set(GLFW_BUILD_X11 OFF CACHE BOOL "" FORCE)

# Web platform
if(EMSCRIPTEN)
    set(PLATFORM Web CACHE STRING "" FORCE)
endif()

add_subdirectory(vendor/raylib)

add_library(imgui STATIC
	vendor/imgui/imgui.cpp
	vendor/imgui/imgui_draw.cpp
	vendor/imgui/imgui_tables.cpp
	vendor/imgui/imgui_widgets.cpp
	vendor/imgui/imgui_demo.cpp
)
target_include_directories(imgui PUBLIC vendor/imgui)

add_library(implot STATIC
	vendor/implot/implot.cpp
	vendor/implot/implot_items.cpp
)
target_include_directories(implot PUBLIC vendor/implot)
target_include_directories(implot PRIVATE vendor/imgui)

add_library(rlImGui STATIC
	vendor/rlImGui/rlImGui.cpp
)
target_include_directories(rlImGui PUBLIC
	vendor/imgui
	vendor/rlImGui
)
target_link_libraries(rlImGui PRIVATE raylib imgui)

add_executable(EMMA
	main.cpp
	src/AppState.hpp
	src/UI.hpp
	# src/MM/typedefs/header.hpp
)

target_compile_options(EMMA PRIVATE
	-O3
	-funroll-loops
)

if(ENABLE_NATIVE_ARCH)
	target_compile_options(EMMA PRIVATE -march=native)
endif()

target_include_directories(EMMA PRIVATE
	src
	src?MM
	src?C++Plot
	vendor/raylib
	vendor/raylib/src
	vendor/imgui
	vendor/implot
	vendor/rlImGui
)
target_link_libraries(EMMA PRIVATE raylib implot imgui rlImGui)

# Web: single-threaded (no COOP/COEP needed)
if(EMSCRIPTEN)
    set(CMAKE_EXECUTABLE_SUFFIX ".html")
    target_link_options(EMMA PRIVATE
        -sALLOW_MEMORY_GROWTH=1
        -sEXIT_RUNTIME=0
        -sSINGLE_FILE=0
    )
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/shell.html")
        target_link_options(EMMA PRIVATE "--shell-file=${CMAKE_CURRENT_SOURCE_DIR}/shell.html")
    endif()
endif()
