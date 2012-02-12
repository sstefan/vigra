# FIND_PACKAGE(Boost COMPONENTS graph)
# INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

ADD_DEFINITIONS(-DVIGRA_HIRES_TIMING)
# ADD_DEFINITIONS(-DVIGRA_CPUTIME_TIMING)

SET(LINK_RT_LIB rt)

VIGRA_ADD_TEST(test_BGL1 test.cxx LIBRARIES vigraimpex ${LINK_RT_LIB})
VIGRA_COPY_TEST_DATA(test.png)

IF(0)
  # Code variants to test observed "compiler optimization noise":
  VIGRA_ADD_TEST(test_BGL_1 test1.cxx empty.cxx LIBRARIES vigraimpex ${LINK_RT_LIB})
  VIGRA_ADD_TEST(test_BGL_2 test2.cxx empty.cxx LIBRARIES vigraimpex ${LINK_RT_LIB})
  VIGRA_ADD_TEST(test_BGL_2a test2a.cxx empty.cxx LIBRARIES vigraimpex ${LINK_RT_LIB})
  VIGRA_ADD_TEST(test_BGL_2b test2b.cxx empty.cxx LIBRARIES vigraimpex ${LINK_RT_LIB})
  # with "doNothing" call
  VIGRA_ADD_TEST(test_BGL_2d test2d.cxx empty.cxx LIBRARIES vigraimpex ${LINK_RT_LIB})

  # put the "accelerating test" in a different compilation unit:
  VIGRA_ADD_TEST(test_BGL_2c test2c.cxx test2c_add.cxx empty.cxx LIBRARIES vigraimpex ${LINK_RT_LIB})
ENDIF()

VIGRA_ADD_TEST(testMST testMST.cxx LIBRARIES vigraimpex ${LINK_RT_LIB})
