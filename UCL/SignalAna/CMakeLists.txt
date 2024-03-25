set(TARGETS SignalAna Overlay OverlayBackground WeightSum Merge Sum Normalise SumFrames MergeAll VertSum)


foreach(TARGET ${TARGETS})

  add_executable(UCL${TARGET} ${TARGET}.cpp)

  

  target_include_directories(UCL${TARGET} PUBLIC ../Plotter)
  target_include_directories(UCL${TARGET} PUBLIC ${ROOT_INCLUDE_DIRS})

  target_link_libraries(UCL${TARGET} PUBLIC PlotterBase)
  target_link_libraries(UCL${TARGET} PUBLIC ${ROOT_LIBRARIES})
  
  install(TARGETS UCL${TARGET} DESTINATION bin)

endforeach()
