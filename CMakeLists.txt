cmake_minimum_required(VERSION 3.18)
project(sapling VERSION 0.1.5 DESCRIPTION "Sapling: a physics and chemistry sandbox built with Fizmo, Grove and Walnut" LANGUAGES C CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(SAPLING_WITH_WALNUT "Build Walnut alongside Grove so Sapling's Walnut scripts work" ON)
option(SAPLING_BUILD_TESTS "Build the Sapling tests" ON)
option(SAPLING_PACKAGE "Make release packages: sapling with its window and compiled script built in, linked so it runs without a compiler" OFF)
option(SAPLING_BUNDLE "Run sapling-host --bundle as part of the build to make sapling with its window built in (turn off when cross-compiling, then bundle on the target)" ON)
set(SAPLING_FIZMO_TAG "v0.3.5" CACHE STRING "Fizmo git tag to fetch (needs the numerics and chemistry components)")
set(SAPLING_GROVE_TAG "v0.2.6" CACHE STRING "Grove git tag to fetch")
set(SAPLING_WALNUT_TAG "v0.2.8" CACHE STRING "Walnut git tag Grove fetches")

if(ANDROID)
    set(SAPLING_WITH_WALNUT OFF CACHE BOOL "" FORCE)
    set(SAPLING_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(SAPLING_PACKAGE OFF CACHE BOOL "" FORCE)
    set(SAPLING_BUNDLE OFF CACHE BOOL "" FORCE)
endif()

if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
endif()

set(GROVE_BUILD_RUNTIME ON CACHE BOOL "" FORCE)
set(GROVE_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GROVE_WITH_WALNUT ${SAPLING_WITH_WALNUT} CACHE BOOL "" FORCE)
set(GROVE_WALNUT_TAG ${SAPLING_WALNUT_TAG} CACHE STRING "" FORCE)
set(GROVE_FIZMO_TAG ${SAPLING_FIZMO_TAG} CACHE STRING "" FORCE)

include(FetchContent)
set(SAPLING_FETCH_EXTRA "")
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.25)
    list(APPEND SAPLING_FETCH_EXTRA SYSTEM)
endif()

FetchContent_Declare(fizmo
    GIT_REPOSITORY https://github.com/ThatNewCoderJSCpp/Fizmo.git
    GIT_TAG        ${SAPLING_FIZMO_TAG}
    GIT_SHALLOW    TRUE
    ${SAPLING_FETCH_EXTRA})

FetchContent_Declare(grove
    GIT_REPOSITORY https://github.com/ThatNewCoderJSCpp/Grove.git
    GIT_TAG        ${SAPLING_GROVE_TAG}
    GIT_SHALLOW    TRUE
    ${SAPLING_FETCH_EXTRA})
FetchContent_MakeAvailable(fizmo grove)

if(NOT TARGET fizmo::numerics)
    message(FATAL_ERROR "Sapling needs Fizmo's numerics component (fizmo::numerics). Fizmo ${SAPLING_FIZMO_TAG} doesn't have it; use v0.2.38 or newer, or point FETCHCONTENT_SOURCE_DIR_FIZMO at a Fizmo checkout that does.")
endif()

function(sapling_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /utf-8)
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endfunction()

if(SAPLING_BUILD_TESTS)
    enable_testing()
endif()

add_library(sapling_core STATIC
    Source/canvas.cpp
    Source/canvas_shapes.cpp
    Source/canvas_text.cpp
    Source/simulation.cpp
    Source/chemistry/bonding_covalent.cpp
    Source/chemistry/bonding_lattice.cpp
    Source/chemistry/bonding_view.cpp
    Source/chemistry/molecule_draw.cpp
    Source/chemistry/orbital_views.cpp
    Source/chemistry/phase_view.cpp
    Source/chemistry/render.cpp
    Source/chemistry/storage.cpp
    Source/chemistry/structure_view.cpp
    Source/chemistry/tools.cpp
    Source/chemistry/view_util.cpp)
add_library(sapling::core ALIAS sapling_core)
target_compile_features(sapling_core PUBLIC cxx_std_17)
target_include_directories(sapling_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/Include" PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/Source")
target_link_libraries(sapling_core PUBLIC fizmo::numerics)
if(TARGET fizmo::chemistry)
    target_link_libraries(sapling_core PUBLIC fizmo::chemistry)
else()
    message(FATAL_ERROR "Sapling's chemistry section needs Fizmo's chemistry component (fizmo::chemistry). Use Fizmo ${SAPLING_FIZMO_TAG} or newer, or point FETCHCONTENT_SOURCE_DIR_FIZMO at a Fizmo checkout that has it.")
endif()
find_package(Threads REQUIRED)
target_link_libraries(sapling_core PUBLIC Threads::Threads)
sapling_warnings(sapling_core)

if(ANDROID)
    set(SAPLING_ANDROID_LIBS_DIR "${CMAKE_BINARY_DIR}/android-libs" CACHE PATH "Where the Android build puts <abi>/libsapling_android.so for grove build --runtime-lib")
    add_library(sapling_android SHARED Tools/sapling.cpp)
    target_link_libraries(sapling_android PRIVATE sapling_core grove::core grove::runtime)
    target_compile_definitions(sapling_android PRIVATE "SAPLING_VERSION_STRING=\"${PROJECT_VERSION}\"")
    target_link_options(sapling_android PRIVATE
        "LINKER:-z,max-page-size=16384"
        "LINKER:-u,ANativeActivity_onCreate"
        "LINKER:--gc-sections"
        "$<$<NOT:$<CONFIG:Debug>>:LINKER:--strip-debug>"
        "LINKER:--version-script=${grove_SOURCE_DIR}/cmake/android_exports.map")
    set_property(TARGET sapling_android APPEND PROPERTY LINK_DEPENDS "${grove_SOURCE_DIR}/cmake/android_exports.map")
    set_target_properties(sapling_android PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${SAPLING_ANDROID_LIBS_DIR}/${ANDROID_ABI}")
    sapling_warnings(sapling_android)
    return()
endif()

set(SAPLING_BUNDLE_FLAGS "")
if(SAPLING_PACKAGE)
    set(SAPLING_BUNDLE_FLAGS --portable)
endif()

add_executable(sapling-host Tools/sapling.cpp)
target_link_libraries(sapling-host PRIVATE sapling_core grove::core grove::runtime)
target_compile_definitions(sapling-host PRIVATE
    "SAPLING_VERSION_STRING=\"${PROJECT_VERSION}\""
    "SAPLING_SOURCE_UI_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/UI\"")
sapling_warnings(sapling-host)
if(WIN32)
    set_target_properties(sapling-host PROPERTIES WIN32_EXECUTABLE $<CONFIG:Release>)
endif()

add_custom_command(TARGET sapling-host POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/UI" "$<TARGET_FILE_DIR:sapling-host>/ui")

file(GLOB_RECURSE SAPLING_UI_FILES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/UI/*.grove" "${CMAKE_CURRENT_SOURCE_DIR}/UI/*.bark" "${CMAKE_CURRENT_SOURCE_DIR}/UI/*.wal")
set(SAPLING_BUNDLED "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/sapling${CMAKE_EXECUTABLE_SUFFIX}")
set(SAPLING_BUNDLE_DEPENDS sapling-host)
if(TARGET walnut)
    list(APPEND SAPLING_BUNDLE_DEPENDS walnut)
endif()
if(SAPLING_BUNDLE)
    add_custom_command(OUTPUT "${SAPLING_BUNDLED}"
        COMMAND ${CMAKE_COMMAND} -E env "SAPLING_UI_DIR=${CMAKE_CURRENT_SOURCE_DIR}/UI" $<TARGET_FILE:sapling-host> --bundle "${SAPLING_BUNDLED}" ${SAPLING_BUNDLE_FLAGS}
        DEPENDS ${SAPLING_BUNDLE_DEPENDS} ${SAPLING_UI_FILES}
        COMMENT "Building sapling with its window and compiled script built in"
        VERBATIM)
    add_custom_target(sapling ALL DEPENDS "${SAPLING_BUNDLED}")
endif()

if(SAPLING_PACKAGE)
    if(MSVC)
        message(FATAL_ERROR "SAPLING_PACKAGE needs a GCC or Clang toolchain, like Grove and Walnut")
    endif()

    string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" _sapling_cpu)
    if(_sapling_cpu MATCHES "^(x86_64|amd64|x64)$")
        if(CMAKE_SIZEOF_VOID_P EQUAL 4)
            set(_sapling_arch i686)
        else()
            set(_sapling_arch x86_64)
        endif()
    elseif(_sapling_cpu MATCHES "^(aarch64|arm64)$")
        set(_sapling_arch arm64)
    elseif(_sapling_cpu MATCHES "^(i[3-6]86|x86)$")
        set(_sapling_arch i686)
    else()
        set(_sapling_arch ${_sapling_cpu})
    endif()

    if(MINGW)
        target_link_options(sapling-host PRIVATE -static)
        include(CheckCXXSourceCompiles)
        check_cxx_source_compiles("#include <cstdlib>\n#ifndef _UCRT\n#error msvcrt\n#endif\nint main() { return 0; }" SAPLING_MINGW_UCRT)
        if(_sapling_arch STREQUAL "arm64")
            set(_sapling_platform windows-clangarm64)
        elseif(_sapling_arch STREQUAL "i686")
            set(_sapling_platform windows-i686)
        elseif(SAPLING_MINGW_UCRT)
            set(_sapling_platform windows-ucrt64)
        else()
            set(_sapling_platform windows-mingw64)
        endif()
        set(CPACK_GENERATOR ZIP)
    elseif(APPLE)
        set(_sapling_platform macos-${_sapling_arch})
        set(CPACK_GENERATOR TGZ)
        set(_sapling_mac_libs "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/lib")
        foreach(t sapling-host)
            add_custom_command(TARGET ${t} POST_BUILD
                COMMAND ${CMAKE_COMMAND} "-DFILES=$<TARGET_FILE:${t}>" "-DLIB_DIR=${_sapling_mac_libs}"
                    -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/macos_bundle_libs.cmake"
                COMMENT "Copying the libraries ${t} needs into ${_sapling_mac_libs}"
                VERBATIM)
        endforeach()
        install(DIRECTORY "${_sapling_mac_libs}/" DESTINATION lib COMPONENT sapling USE_SOURCE_PERMISSIONS)
    else()
        target_link_options(sapling-host PRIVATE -static-libstdc++ -static-libgcc)
        set(_sapling_platform linux-${_sapling_arch})
        set(CPACK_GENERATOR TGZ)
    endif()
    set(SAPLING_PACKAGE_PLATFORM ${_sapling_platform} CACHE STRING "Platform name used in the package file name")

    if(SAPLING_BUNDLE)
        install(PROGRAMS "${SAPLING_BUNDLED}" DESTINATION . COMPONENT sapling)
    endif()
    install(FILES README.md DESTINATION . COMPONENT sapling)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
        install(FILES LICENSE DESTINATION . COMPONENT sapling)
    endif()
    set(CPACK_INSTALL_CMAKE_PROJECTS "${CMAKE_BINARY_DIR};Sapling;sapling;/")

    set(CPACK_PACKAGE_NAME "Sapling")
    set(CPACK_PACKAGE_VENDOR "ThatNewCoderJSCpp")
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}")
    set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
    set(CPACK_PACKAGE_FILE_NAME "Sapling-${PROJECT_VERSION}-${SAPLING_PACKAGE_PLATFORM}")
    set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY ON)
    set(CPACK_PACKAGE_INSTALL_DIRECTORY "Sapling")
    include(CPack)
else()
    if(SAPLING_BUNDLE)
        install(PROGRAMS "${SAPLING_BUNDLED}" DESTINATION bin)
    endif()
endif()

if(SAPLING_BUILD_TESTS)
    add_executable(sapling_tests Tests/sapling_tests.cpp)
    target_link_libraries(sapling_tests PRIVATE sapling_core)
    sapling_warnings(sapling_tests)
    add_test(NAME sapling_tests COMMAND sapling_tests)
endif()
