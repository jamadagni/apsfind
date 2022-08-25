// based on https://stackoverflow.com/a/19898211/1503120

#include <time.h>

typedef struct timespec process_timer;

// call this function to start a nanosecond-resolution timer
process_timer process_timer_init()
{
    process_timer start_time;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time);
    return start_time;
}

// call this function to end a timer, returning nanoseconds elapsed as a long
double process_timer_read(process_timer start_time)
{
    process_timer cur_time;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cur_time);
    long elapsed_nanosecs = (cur_time.tv_sec - start_time.tv_sec) * (long)1e9 + (cur_time.tv_nsec - start_time.tv_nsec);
    return elapsed_nanosecs / 1e9;
}
