OPTION(USE_COMPILED_TEST_DATA "Simulated Watershed Test: Use compiled test data (large!)" OFF)
# SET(USE_COMPILED_TEST_DATA ON)

IF(USE_COMPILED_TEST_DATA)
  ADD_DEFINITIONS(-DUSE_COMPILED_TEST_DATA)
  VIGRA_ADD_TEST(test_multiiterator_benchmark_watersheds_simulatedPQ log.cxx test.cxx empty.cxx)
ELSE()
  VIGRA_COPY_TEST_DATA(log.bin)
  VIGRA_ADD_TEST(test_multiiterator_benchmark_watersheds_simulatedPQ test.cxx empty.cxx)
ENDIF()



IF(USE_COMPILED_TEST_DATA)
  ADD_DEFINITIONS(-DUSE_COMPILED_TEST_DATA)
  VIGRA_ADD_TEST(test_multiiterator_benchmark_watersheds_simulatedPQ2 log.cxx test2.cxx empty.cxx)
ELSE()
  VIGRA_COPY_TEST_DATA(log.bin)
  VIGRA_ADD_TEST(test_multiiterator_benchmark_watersheds_simulatedPQ2 test2.cxx empty.cxx)
ENDIF()


