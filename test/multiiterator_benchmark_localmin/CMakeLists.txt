# this one probably requires boost BGL includes.
VIGRA_ADD_TEST(test_multiiterator_benchmark_localmin test.cxx LIBRARIES vigraimpex)
VIGRA_COPY_TEST_DATA(test.png)



