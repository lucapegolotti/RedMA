add_executable(RBNavierStokesSingleDomain RBNavierStokesSingleDomain.cpp ${SOURCES})
target_link_libraries(RBNavierStokesSingleDomain ${LINKED_LIBRARIES})
target_include_directories(RBNavierStokesSingleDomain PRIVATE ${HEADERS_DIRS})
