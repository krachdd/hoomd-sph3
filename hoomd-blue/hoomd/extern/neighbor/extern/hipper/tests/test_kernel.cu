#include <hipper/hipper_runtime.h>
#include "test_hipper.h"

__global__ void set_thread(int* a, const int N)
    {
    int rank = hipper::threadRank<1,1>();
    if (rank >= N) return;

    a[rank] = rank+1;
    }

TEST_CASE("Basic kernel launch", "[kernel]")
    {
    int* a;
    const int N = 2;

    REQUIRE_SUCCESS(hipper::mallocManaged(reinterpret_cast<void**>(&a), sizeof(int)*N));

    REQUIRE_SUCCESS(hipper::memset(a, 0, sizeof(int)*N));
    REQUIRE(a[0] == 0);
    REQUIRE(a[1] == 0);

    hipper::KernelLauncher(1, 32)(set_thread, a, N);
    REQUIRE_SUCCESS(hipper::peekAtLastError());
    REQUIRE_SUCCESS(hipper::deviceSynchronize());
    REQUIRE(a[0] == 1);
    REQUIRE(a[1] == 2);

    REQUIRE_SUCCESS(hipper::free(a));
    }
