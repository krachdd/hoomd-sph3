// Copyright (c) 2020, Michael P. Howard
// Copyright (c) 2021, Auburn University
// This file is released under the Modified BSD License.

#ifndef HIPPER_TESTS_TEST_HIPPER_H_
#define HIPPER_TESTS_TEST_HIPPER_H_

#include <catch2/catch.hpp>

#define REQUIRE_SUCCESS(f) REQUIRE(f == hipper::success);

#endif // HIPPER_TESTS_TEST_HIPPER_H_
