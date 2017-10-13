#ifndef mfbvar__mfbvarConfig__h
#define mfbvar__mfbvarConfig__h
#if defined(WIN32) || defined(_WIN32)
   #define ARMA_USE_OPENMP
#else
   #include <mfbvarConfigGenerated.h>
#endif
#endif
