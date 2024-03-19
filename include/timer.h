#ifndef __TIMER__
#define __TIMER__

#include <math.h>

#ifndef __HELPER_BUFF_SIZE__
#define __HELPER_BUFF_SIZE__ 64
#endif

#ifdef _WIN32
#include <windows.h>

typedef struct timer_s
{
  LARGE_INTEGER frequency;
  LARGE_INTEGER start;
} timer;

void timer_start(timer *timer);
double timer_stop(timer *timer);
#endif

#endif // __TIMER__

#ifdef __TIMER_IMPLEMENTATION__

#ifdef _WIN32

void timer_start(timer *timer)
{
  QueryPerformanceFrequency(&(timer->frequency));
  QueryPerformanceCounter(&(timer->start));
  return;
}

double timer_stop(timer *timer)
{
  double interval;
  LARGE_INTEGER end;
  QueryPerformanceCounter(&end);
  interval = (double)(end.QuadPart - timer->start.QuadPart) / timer->frequency.QuadPart;
  return interval;
}

#else /* _WIN32 */
#error "timer is not (yet) defined for non windows systems (PR are welcome ;))"
#endif

#endif // __TIMER_IMPLEMENTATION__