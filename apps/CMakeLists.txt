add_executable(PrepareTagTree PrepareTagTree.cpp)

target_link_libraries(PrepareTagTree PUBLIC KKpipiStrongPhase)
target_link_libraries(PrepareTagTree PUBLIC ROOT::Physics ROOT::RIO)

install(TARGETS PrepareTagTree DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../bin)
