PROJECT(PowerOf2Tiling)

#INCLUDE_DIRECTORIES()

SET(PowerOf2Tiling_HDRS 
)

SET(PowerOf2Tiling_SRCS 
	main.cpp
)

ADD_EXECUTABLE(PowerOf2Tiling main.cpp ${PowerOf2Tiling_HDRS} ${PowerOf2Tiling_SRCS})
TARGET_LINK_LIBRARIES(PowerOf2Tiling vul ${ITK_LIBRARIES})
INSTALL( TARGETS PowerOf2Tiling RUNTIME DESTINATION ${INSTALL_BIN_DIR} )
