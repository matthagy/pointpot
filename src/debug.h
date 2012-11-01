#ifndef DEBUG_H
#define DEBUG_H

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "opt.h"

#define Fatal(FRMT, ARGS...)                            \
        _Fatal(__FILE__, __FUNCTION__, __LINE__, FRMT, ## ARGS)

namespace PointPot {

void _Fatal(const char *file, const char *func, int line,
                 const char *frmt, ...)
        GCC_ATTRIBUTE((format (printf, 4, 5)))
        GCC_ATTRIBUTE((noreturn));

}

#endif /* DEBUG_H */
