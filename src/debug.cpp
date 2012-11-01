
#include <stdio.h>
#include <stdarg.h>

#include "debug.h"

using namespace PointPot;

/* "../stuff/name.c" -> "name.c" */
static const char *
fixFileName(const char *p)
{
  int i,l;

  l = strlen(p);
  for (i=l; i>0 && p[i-1] != '/'; i--);
  return p+i;
}

void PointPot::_Fatal(const char *file, const char *func, int line,
                      const char *frmt, ...)
{
  va_list va;

  fprintf(stderr," *** PointPot Fatal ***\n");
  fflush(stderr);
  fprintf(stderr,"In %.200s:%.200s:%i\n", fixFileName(file), func, line);
  fflush(stderr);
  va_start(va, frmt);
  vfprintf(stderr, frmt, va);
  va_end(va);
  fprintf(stderr,"\n");
  fflush(stderr);
  abort();
}

