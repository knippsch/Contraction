/* $Id: propagator_io.h,v 1.2 2007/11/24 14:37:24 urbach Exp $ */

#ifndef _PROPAGATOR_IO_HH
#define _PROPAGATOR_IO_HH

typedef unsigned int uint32_t;
typedef uint32_t DML_SiteRank;

typedef struct {
  uint32_t suma;
  uint32_t sumb;
} DML_Checksum;

typedef uint32_t uLong;            /* At least 32 bits */
typedef unsigned char Byte;
typedef Byte Bytef;
typedef uLong uLongf;
#define Z_NULL  0  /* for initializing zalloc, zfree, opaque */

int write_lime_spinor (double * const s, char * filename, const int append,
		const int prec, const unsigned int T, const unsigned int LX,
		const unsigned int LY, const unsigned int LZ);
int read_lime_spinor (double * const s, char * filename, const int position,
		const int number_data_points, const int size_data_point);

#endif
