cmake_minimum_required(VERSION 2.4)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

FIND_PATH(NETCDF_C_INCLUDES netcdf.h
  /Users/cahrens/climate/netcdf/netcdfinstall/include
)

FIND_PATH(NETCDF_CPP_INCLUDES netcdfcpp.h
  /Users/cahrens/climate/netcdf/netcdfinstall/include
)

FIND_LIBRARY(NETCDF_C_LIBRARY netcdf
  /Users/cahrens/climate/netcdf/netcdfinstall/lib
)

FIND_LIBRARY(NETCDF_CPP_LIBRARY netcdf_c++
  /Users/cahrens/climate/netcdf/netcdfinstall/lib
)

INCLUDE_DIRECTORIES(${NETCDF_C_INCLUDES} ${NETCDF_CPP_INCLUDES})
LINK_LIBRARIES(${NETCDF_CPP_LIBRARY} ${NETCDF_C_LIBRARY})

#define the server-side portion of the reader plugin
ADD_PARAVIEW_PLUGIN(NewPopReaderSMPlugin "1.0"
  SERVER_MANAGER_SOURCES vtkNewPopReader.cxx
  SERVER_MANAGER_XML NewPopReaderSM.xml)

#define the client-side portion of the reader plugin
ADD_PARAVIEW_PLUGIN(NewPopReaderGUIPlugin "1.0"
  GUI_RESOURCE_FILES NewPopReaderGUI.xml)
