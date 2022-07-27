#include <hipper/hipper_runtime.h>
#include <hipper/hipper_cub.h>
#include "test_hipper.h"

TEST_CASE("CUB operations", "[CUB]")
    {
    // input (array of integers)
    int* a;
    const int N = 2;
    REQUIRE_SUCCESS(hipper::mallocManaged(reinterpret_cast<void**>(&a), sizeof(int)*N));
    a[0] = 1;
    a[1] = 2;

    // output (sum of a)
    int* total;
    REQUIRE_SUCCESS(hipper::mallocManaged(reinterpret_cast<void**>(&total), sizeof(int)));
    *total = 0;

    // size temporary memory
    void *tmp = NULL;
    size_t tmp_bytes = 0;
    hipper::cub::DeviceReduce::Sum(tmp,tmp_bytes,a,total,N);
    REQUIRE_SUCCESS(hipper::mallocManaged(reinterpret_cast<void**>(&tmp), tmp_bytes));

    // take sum
    hipper::cub::DeviceReduce::Sum(tmp,tmp_bytes,a,total,N);
    REQUIRE_SUCCESS(hipper::deviceSynchronize());

    // check output
    REQUIRE(a[0] == 1);
    REQUIRE(a[1] == 2);
    REQUIRE(*total == 3);

    // free memory
    REQUIRE_SUCCESS(hipper::free(a));
    REQUIRE_SUCCESS(hipper::free(total));
    REQUIRE_SUCCESS(hipper::free(tmp));
    }
