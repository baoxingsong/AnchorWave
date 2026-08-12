#include "ksw2.h"

#if defined(_MSC_VER)
#define KSW_THREAD_LOCAL __declspec(thread)
#elif defined(__GNUC__) || defined(__clang__)
#define KSW_THREAD_LOCAL __thread
#else
#define KSW_THREAD_LOCAL _Thread_local
#endif

static KSW_THREAD_LOCAL ksw_progress_callback_t ksw_progress_callback = 0;
static KSW_THREAD_LOCAL void *ksw_progress_arguments = 0;
static KSW_THREAD_LOCAL int ksw_progress_interval = 256;

void ksw_set_progress_callback(ksw_progress_callback_t callback,
                               void *arguments,
                               int antidiagonal_interval) {
    ksw_progress_callback = callback;
    ksw_progress_arguments = arguments;
    ksw_progress_interval = antidiagonal_interval > 0
            ? antidiagonal_interval : 1;
}

void ksw_clear_progress_callback(void) {
    ksw_progress_callback = 0;
    ksw_progress_arguments = 0;
    ksw_progress_interval = 256;
}

int ksw_progress_continue(int antidiagonal, int total_antidiagonals) {
    if (ksw_progress_callback == 0 ||
        (antidiagonal % ksw_progress_interval) != 0) {
        return 1;
    }
    return ksw_progress_callback(ksw_progress_arguments,
                                 antidiagonal,
                                 total_antidiagonals) != 0;
}
