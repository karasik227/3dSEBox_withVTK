cmake_minimum_required(VERSION 3.5)

project(3dSEBox_withVTK VERSION 0.1 LANGUAGES CXX)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Widgets)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets)
find_package(Qt6 REQUIRED COMPONENTS OpenGLWidgets)
find_package(Qt6 REQUIRED COMPONENTS Charts)

set(PROJECT_SOURCES
        main.cpp
        mainwindow.cpp
        mainwindow.h
        mainwindow.ui
)

set(VTK_DIR "E:/VTK-9.3.0-install-no-CASCADE/lib/cmake/vtk-9.3")

find_package(VTK COMPONENTS
  CommonColor
  CommonCore
  FiltersSources
  CommonDataModel
  InteractionStyle
  RenderingContextOpenGL2
  RenderingCore
  RenderingFreeType
  RenderingGL2PSOpenGL2
  RenderingOpenGL2
)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "Cube: Unable to find the VTK build folder.")
endif()

if(${QT_VERSION_MAJOR} GREATER_EQUAL 6)
    qt_add_executable(3dSEBox_withVTK
        MANUAL_FINALIZATION
        ${PROJECT_SOURCES}
        robcalculation.h
        robcalculation.cpp
    )
# Define target properties for Android with Qt 6 as:
#    set_property(TARGET 3dSEBox_withVTK APPEND PROPERTY QT_ANDROID_PACKAGE_SOURCE_DIR
#                 ${CMAKE_CURRENT_SOURCE_DIR}/android)
# For more information, see https://doc.qt.io/qt-6/qt-add-executable.html#target-creation
else()
    if(ANDROID)
        add_library(3dSEBox_withVTK SHARED
            ${PROJECT_SOURCES}
        )
# Define properties for Android with Qt 5 after find_package() calls as:
#    set(ANDROID_PACKAGE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/android")
    else()
        add_executable(3dSEBox_withVTK
            ${PROJECT_SOURCES}
        )
    endif()
endif()

target_link_libraries(3dSEBox_withVTK PRIVATE
    "E:/VTK-9.3.0-install-no-CASCADE/lib/vtkGUISupportQt-9.3.lib"
    "E:/VTK-9.3.0-install-no-CASCADE/lib/vtkInteractionWidgets-9.3.lib"
    "E:/VTK-9.3.0-install-no-CASCADE/lib/vtkRenderingAnnotation-9.3.lib"
    Qt6::Core
    Qt6::Gui
    Qt6::Charts
    Qt6::Widgets
    Qt6::OpenGL
    Qt6::OpenGLWidgets
        ${VTK_LIBRARIES})

# Qt for iOS sets MACOSX_BUNDLE_GUI_IDENTIFIER automatically since Qt 6.1.
# If you are developing for iOS or macOS you should consider setting an
# explicit, fixed bundle identifier manually though.
if(${QT_VERSION} VERSION_LESS 6.1.0)
  set(BUNDLE_ID_OPTION MACOSX_BUNDLE_GUI_IDENTIFIER com.example.3dSEBox_withVTK)
endif()
set_target_properties(3dSEBox_withVTK PROPERTIES
    ${BUNDLE_ID_OPTION}
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    MACOSX_BUNDLE TRUE
    WIN32_EXECUTABLE TRUE
)

include(GNUInstallDirs)
install(TARGETS 3dSEBox_withVTK
    BUNDLE DESTINATION .
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

if(QT_VERSION_MAJOR EQUAL 6)
    qt_finalize_executable(3dSEBox_withVTK)
endif()
