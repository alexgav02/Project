

add_library(PlotterBase PlotterBase.cpp)
target_include_directories(PlotterBase PUBLIC ${ROOT_INCLUDE_DIRS})

install(TARGETS PlotterBase DESTINATION lib)
