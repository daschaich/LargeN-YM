// -----------------------------------------------------------------
// Collect macros for preprocessor tweaks
// that accommodate differences in compilers, architecture and OS
#ifndef _CONFIG_H
#define _CONFIG_H
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compiler/Processor-dependent macros
// Specify the unsigned 32 bit integer base type for this compiler
// Run the script "getint.sh" to find out what to use
// One and only one of these should be defined
#define INT_IS_32BIT 1  // Most present systems
#undef SHORT_IS_32BIT   // Needed on T3E UNICOS, for example

// Define if you have a 64-byte cache line (if not, we assume 32 bytes)
// Processors that do: P4 (actually fetches 128), EV67, EV68
// Used only for prefetching, so it only affects performance
#define HAVE_64_BYTE_CACHELINE 1
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compiler/OS-dependent macros
// Define if you have the <ieeefp.h> header file.
// Systems that don't: T3E UNICOS, Exemplar, Linux gcc, SP AIX, HP True64
//#define HAVE_IEEEFP_H 1

// Define if you have the <unistd.h> header file.
// Systems that don't: NT
#define HAVE_UNISTD_H 1

// Define if you have the <sys/time.h> header file.
// Most systems do
#define HAVE_SYS_TIME_H 1

// Define if you have ANSI "fseeko"
// Systems that don't: T3E UNICOS
#define HAVE_FSEEKO 1

#endif
// -----------------------------------------------------------------
