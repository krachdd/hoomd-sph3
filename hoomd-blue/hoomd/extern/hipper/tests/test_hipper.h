#ifndef HIPPER_TESTS_TEST_HIPPER_H_
#define HIPPER_TESTS_TEST_HIPPER_H_

#include <catch2/catch.hpp>

#define REQUIRE_SUCCESS(f) REQUIRE(f == hipper::success);

#endif // HIPPER_TESTS_TEST_HIPPER_H_
