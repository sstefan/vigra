INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include)

INCLUDE(VigraAddTest)

IF(1)
ADD_SUBDIRECTORY(error)
ADD_SUBDIRECTORY(impex)
ADD_SUBDIRECTORY(utilities)
ADD_SUBDIRECTORY(pixeltypes)
ADD_SUBDIRECTORY(colorspaces)
ADD_SUBDIRECTORY(classifier)
ADD_SUBDIRECTORY(imgproc)
ADD_SUBDIRECTORY(functorexpression)
ADD_SUBDIRECTORY(convolution)
ADD_SUBDIRECTORY(morphology)
ADD_SUBDIRECTORY(fourier)
ADD_SUBDIRECTORY(math)
ADD_SUBDIRECTORY(optimization)
ADD_SUBDIRECTORY(tensorimaging)
#ADD_SUBDIRECTORY(simpleanalysis)
ADD_SUBDIRECTORY(image)
ADD_SUBDIRECTORY(multiarray)
ADD_SUBDIRECTORY(multiconvolution)
ADD_SUBDIRECTORY(voxelneighborhood)
#ADD_SUBDIRECTORY(volumelabeling)
#ADD_SUBDIRECTORY(watersheds3d)
ADD_SUBDIRECTORY(multidistance)
ADD_SUBDIRECTORY(multimorphology)
#ADD_SUBDIRECTORY(seededRegionGrowing3d)
ADD_SUBDIRECTORY(sampler)
ADD_SUBDIRECTORY(hdf5impex)
ADD_SUBDIRECTORY(sifImport)
ADD_SUBDIRECTORY(features)
ADD_SUBDIRECTORY(unsupervised)
ENDIF()

ADD_SUBDIRECTORY(coordinateiterator)
ADD_SUBDIRECTORY(objectfeatures)
# adapted for gridgraph (these have a gridgraph version):
ADD_SUBDIRECTORY(simpleanalysis)
ADD_SUBDIRECTORY(volumelabeling)
ADD_SUBDIRECTORY(watersheds3d)
ADD_SUBDIRECTORY(seededRegionGrowing3d)


OPTION(WITH_BOOST "Use Boost Graph library (if available)" ON)
IF(WITH_BOOST)
  FIND_PACKAGE(Boost COMPONENTS graph) # not strictly REQUIRED
ENDIF()
IF(Boost_FOUND AND WITH_BOOST)
  ADD_DEFINITIONS(-DWITH_BOOST_GRAPH)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
ENDIF()

# iterator experiments
#ADD_SUBDIRECTORY(multiiterator_benchmark1)
#ADD_SUBDIRECTORY(multiiterator_benchmark_localmin) # -> now in BGL_interface
ADD_SUBDIRECTORY(multiiterator_benchmark_watersheds)
#ADD_SUBDIRECTORY(multiiterator_benchmark_watersheds_simulatedPQ)
#ADD_SUBDIRECTORY(multiiterator_gridgraph1)

# BGL interface variants
ADD_SUBDIRECTORY(BGL_interface)

#ADD_SUBDIRECTORY(multiiterator_coupled)

