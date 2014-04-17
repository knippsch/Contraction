/* $Id: propagator_io.c,v 1.3 2007/11/24 14:37:24 urbach Exp $ */
#define _FILE_OFFSET_BITS 64
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>
#include <sys/time.h>
#include <sys/types.h>
#include "lime.h"
#include "lime_fixed_types.h"
#include "io_utils.h"
#include "propagator_io.h"
#define local static
/* ========================================================================
 * Table of CRC-32's of all single-byte values (made by make_crc_table)
 */local uLongf crc_table[256] = { 0x00000000L, 0x77073096L, 0xee0e612cL,
		0x990951baL, 0x076dc419L, 0x706af48fL, 0xe963a535L, 0x9e6495a3L,
		0x0edb8832L, 0x79dcb8a4L, 0xe0d5e91eL, 0x97d2d988L, 0x09b64c2bL,
		0x7eb17cbdL, 0xe7b82d07L, 0x90bf1d91L, 0x1db71064L, 0x6ab020f2L,
		0xf3b97148L, 0x84be41deL, 0x1adad47dL, 0x6ddde4ebL, 0xf4d4b551L,
		0x83d385c7L, 0x136c9856L, 0x646ba8c0L, 0xfd62f97aL, 0x8a65c9ecL,
		0x14015c4fL, 0x63066cd9L, 0xfa0f3d63L, 0x8d080df5L, 0x3b6e20c8L,
		0x4c69105eL, 0xd56041e4L, 0xa2677172L, 0x3c03e4d1L, 0x4b04d447L,
		0xd20d85fdL, 0xa50ab56bL, 0x35b5a8faL, 0x42b2986cL, 0xdbbbc9d6L,
		0xacbcf940L, 0x32d86ce3L, 0x45df5c75L, 0xdcd60dcfL, 0xabd13d59L,
		0x26d930acL, 0x51de003aL, 0xc8d75180L, 0xbfd06116L, 0x21b4f4b5L,
		0x56b3c423L, 0xcfba9599L, 0xb8bda50fL, 0x2802b89eL, 0x5f058808L,
		0xc60cd9b2L, 0xb10be924L, 0x2f6f7c87L, 0x58684c11L, 0xc1611dabL,
		0xb6662d3dL, 0x76dc4190L, 0x01db7106L, 0x98d220bcL, 0xefd5102aL,
		0x71b18589L, 0x06b6b51fL, 0x9fbfe4a5L, 0xe8b8d433L, 0x7807c9a2L,
		0x0f00f934L, 0x9609a88eL, 0xe10e9818L, 0x7f6a0dbbL, 0x086d3d2dL,
		0x91646c97L, 0xe6635c01L, 0x6b6b51f4L, 0x1c6c6162L, 0x856530d8L,
		0xf262004eL, 0x6c0695edL, 0x1b01a57bL, 0x8208f4c1L, 0xf50fc457L,
		0x65b0d9c6L, 0x12b7e950L, 0x8bbeb8eaL, 0xfcb9887cL, 0x62dd1ddfL,
		0x15da2d49L, 0x8cd37cf3L, 0xfbd44c65L, 0x4db26158L, 0x3ab551ceL,
		0xa3bc0074L, 0xd4bb30e2L, 0x4adfa541L, 0x3dd895d7L, 0xa4d1c46dL,
		0xd3d6f4fbL, 0x4369e96aL, 0x346ed9fcL, 0xad678846L, 0xda60b8d0L,
		0x44042d73L, 0x33031de5L, 0xaa0a4c5fL, 0xdd0d7cc9L, 0x5005713cL,
		0x270241aaL, 0xbe0b1010L, 0xc90c2086L, 0x5768b525L, 0x206f85b3L,
		0xb966d409L, 0xce61e49fL, 0x5edef90eL, 0x29d9c998L, 0xb0d09822L,
		0xc7d7a8b4L, 0x59b33d17L, 0x2eb40d81L, 0xb7bd5c3bL, 0xc0ba6cadL,
		0xedb88320L, 0x9abfb3b6L, 0x03b6e20cL, 0x74b1d29aL, 0xead54739L,
		0x9dd277afL, 0x04db2615L, 0x73dc1683L, 0xe3630b12L, 0x94643b84L,
		0x0d6d6a3eL, 0x7a6a5aa8L, 0xe40ecf0bL, 0x9309ff9dL, 0x0a00ae27L,
		0x7d079eb1L, 0xf00f9344L, 0x8708a3d2L, 0x1e01f268L, 0x6906c2feL,
		0xf762575dL, 0x806567cbL, 0x196c3671L, 0x6e6b06e7L, 0xfed41b76L,
		0x89d32be0L, 0x10da7a5aL, 0x67dd4accL, 0xf9b9df6fL, 0x8ebeeff9L,
		0x17b7be43L, 0x60b08ed5L, 0xd6d6a3e8L, 0xa1d1937eL, 0x38d8c2c4L,
		0x4fdff252L, 0xd1bb67f1L, 0xa6bc5767L, 0x3fb506ddL, 0x48b2364bL,
		0xd80d2bdaL, 0xaf0a1b4cL, 0x36034af6L, 0x41047a60L, 0xdf60efc3L,
		0xa867df55L, 0x316e8eefL, 0x4669be79L, 0xcb61b38cL, 0xbc66831aL,
		0x256fd2a0L, 0x5268e236L, 0xcc0c7795L, 0xbb0b4703L, 0x220216b9L,
		0x5505262fL, 0xc5ba3bbeL, 0xb2bd0b28L, 0x2bb45a92L, 0x5cb36a04L,
		0xc2d7ffa7L, 0xb5d0cf31L, 0x2cd99e8bL, 0x5bdeae1dL, 0x9b64c2b0L,
		0xec63f226L, 0x756aa39cL, 0x026d930aL, 0x9c0906a9L, 0xeb0e363fL,
		0x72076785L, 0x05005713L, 0x95bf4a82L, 0xe2b87a14L, 0x7bb12baeL,
		0x0cb61b38L, 0x92d28e9bL, 0xe5d5be0dL, 0x7cdcefb7L, 0x0bdbdf21L,
		0x86d3d2d4L, 0xf1d4e242L, 0x68ddb3f8l, 0x1fda836eL, 0x81be16cdL,
		0xf6b9265bL, 0x6fb077e1L, 0x18b74777L, 0x88085ae6L, 0xff0f6a70L,
		0x66063bcaL, 0x11010b5cL, 0x8f659effL, 0xf862ae69L, 0x616bffd3L,
		0x166ccf45L, 0xa00ae278L, 0xd70dd2eeL, 0x4e048354L, 0x3903b3c2L,
		0xa7672661L, 0xd06016f7L, 0x4969474dL, 0x3e6e77dbL, 0xaed16a4aL,
		0xd9d65adcL, 0x40df0b66L, 0x37d83bf0L, 0xa9bcae53L, 0xdebb9ec5L,
		0x47b2cf7fL, 0x30b5ffe9L, 0xbdbdf21cL, 0xcabac28aL, 0x53b39330L,
		0x24b4a3a6L, 0xbad03605L, 0xcdd70693L, 0x54de5729L, 0x23d967bfL,
		0xb3667a2eL, 0xc4614ab8L, 0x5d681b02L, 0x2a6f2b94L, 0xb40bbe37L,
		0xc30c8ea1L, 0x5a05df1bL, 0x2d02ef8dL };
/* ========================================================================= */
#define DO1(buf) crc = crc_table[((int)crc ^ (*buf++)) & 0xff] ^ (crc >> 8);
#define DO2(buf)  DO1(buf); DO1(buf);
#define DO4(buf)  DO2(buf); DO2(buf);
#define DO8(buf)  DO4(buf); DO4(buf);
/* ========================================================================= */
uint32_t DML_crc32 (uint32_t crc, const unsigned char *buf, size_t len) {
	if(buf == Z_NULL) return 0L;
	crc = crc ^ 0xffffffffL;
	while(len >= 8){
		DO8(buf);
		len -= 8;
	}
	if(len) do{
		DO1(buf);
	}
	while(--len);
	return crc ^ 0xffffffffL;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void DML_checksum_init (DML_Checksum *checksum) {
	checksum->suma = 0;
	checksum->sumb = 0;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void DML_checksum_accum (DML_Checksum *checksum, DML_SiteRank rank, char *buf,
		size_t size) {

	DML_SiteRank rank29 = rank;
	DML_SiteRank rank31 = rank;
	uint32_t work = DML_crc32(0, (unsigned char*) buf, size);

	rank29 %= 29;
	rank31 %= 31;

	checksum->suma ^= work << rank29 | work >> (32 - rank29);
	checksum->sumb ^= work << rank31 | work >> (32 - rank31);
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
DML_Checksum write_binary_spinor_data (double * const s,
		LimeWriter * limewriter, const int prec, const unsigned int T,
		const unsigned int LX, const unsigned int LY, const unsigned int LZ) {

	int status = 0;
	double tmp[24];
	float tmp2[24];
	n_uint64_t bytes;
	//off_t bytes;
	DML_Checksum ans;
	DML_SiteRank rank;
	int words_bigendian;
	words_bigendian = big_endian();

	DML_checksum_init(&ans);
	rank = (DML_SiteRank) 0;

	if(prec == 32) bytes = 24 * sizeof(float);
	else bytes = 24 * sizeof(double);
	for(unsigned int t = 0; t < T; t++){
		for(unsigned int x = 0; x < LX; x++){
			for(unsigned int y = 0; y < LY; y++){
				for(unsigned int z = 0; z < LZ; z++){

					n_uint64_t ix = (t * LX * LY * LZ + x * LY * LZ + y * LZ + z)
							* (n_uint64_t) 12;
					rank = (DML_SiteRank) (((t * LZ + z) * LY + y) * LX + x);
					if(!words_bigendian){
						if(prec == 32){
							byte_swap_assign_double2single((float*) tmp2, &s[2 * ix], 24);
							DML_checksum_accum(&ans, rank, (char *) tmp2, bytes);
							status = limeWriteRecordData((void*) tmp2, &bytes, limewriter);
						}
						else{
							byte_swap_assign(tmp, &s[2 * ix], 24);
							DML_checksum_accum(&ans, rank, (char *) tmp, bytes);
							status = limeWriteRecordData((void*) tmp, &bytes, limewriter);
						}
					}
					else{
						if(prec == 32){
							double2single((float*) tmp2, &s[2 * ix], 24);
							DML_checksum_accum(&ans, rank, (char *) tmp2, bytes);
							status = limeWriteRecordData((void*) tmp2, &bytes, limewriter);
						}
						else{
							status = limeWriteRecordData((void*) &s[2 * ix], &bytes,
									limewriter);
							DML_checksum_accum(&ans, rank, (char *) &s[2 * ix], bytes);
						}
					}
				}
			}
		}
	}
//	printf("The final checksum is %#lx %#lx\n", (unsigned long) ans.suma,
//			(unsigned long) ans.sumb);
	return (ans);
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
int write_lime_spinor (double * const s, char * filename, const int append,
		const int prec, const unsigned int T, const unsigned int LX,
		const unsigned int LY, const unsigned int LZ) {

	FILE * ofs = NULL;
	LimeWriter * limewriter = NULL;
	LimeRecordHeader * limeheader = NULL;
	int status = 0;
	int ME_flag = 0, MB_flag = 0;
	n_uint64_t bytes;
	//off_t bytes;
	DML_Checksum checksum;

	if(append){
		ofs = fopen(filename, "a");
	}
	else{
		ofs = fopen(filename, "w");
	}
	if(ofs == (FILE*) NULL){
		fprintf(stderr, "Could not open file %s for writing!\n Aborting...\n",
				filename);
		exit(500);
	}
	limewriter = limeCreateWriter(ofs);
	if(limewriter == (LimeWriter*) NULL){
		fprintf(stderr, "LIME error in file %s for writing!\n Aborting...\n",
				filename);
		exit(500);
	}

	bytes = LX * LY * LZ * T * (n_uint64_t) 24 * sizeof(double) * prec / 64;
	MB_flag = 0;
	ME_flag = 1;
	limeheader = limeCreateHeader(MB_flag, ME_flag, "scidac-binary-data", bytes);
	status = limeWriteRecordHeader(limeheader, limewriter);
	if(status < 0){
		fprintf(stderr, "LIME write header (scidac-binary-data) error %d\n",
				status);
		exit(500);
	}
	limeDestroyHeader(limeheader);

	checksum = write_binary_spinor_data(s, limewriter, prec, T, LX, LY, LZ);
	printf("Final check sum is (%#lx  %#lx)\n", (unsigned long) checksum.suma,
			(unsigned long) checksum.sumb);
	if(ferror(ofs)){
		fprintf(stderr, "Warning! Error while writing to file %s \n", filename);
	}
	limeDestroyWriter(limewriter);
	fclose(ofs);
	fflush(ofs);
	return (0);
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
int read_binary_spinor_data (double * const s, LimeReader * limereader,
		const double prec, DML_Checksum &ans, const int number_data_points,
		const int size_data_point) {

	int status = 0;
	n_uint64_t bytes;
	//off_t bytes;
	double tmp[size_data_point];
	DML_SiteRank rank;
	float tmp2[24];
	int words_bigendian;
	words_bigendian = big_endian();

	DML_checksum_init(&ans);
	rank = (DML_SiteRank) 0;

	if(prec == 32) bytes = size_data_point * sizeof(float);
	else bytes = size_data_point * sizeof(double);

	for(int x = 0; x < number_data_points; x++){
		n_uint64_t ix = (n_uint64_t) x * (n_uint64_t) size_data_point;
		rank = (DML_SiteRank) x;
		if(prec == 32){
			status = limeReaderReadData(tmp2, &bytes, limereader);
			DML_checksum_accum(&ans, rank, (char *) tmp2, bytes);
		}
		else{
			status = limeReaderReadData(tmp, &bytes, limereader);
			DML_checksum_accum(&ans, rank, (char *) tmp, bytes);
		}
		if(!words_bigendian){
			if(prec == 32){
				byte_swap_assign_single2double(&s[ix], (float*) tmp2, size_data_point);
			}
			else{
				byte_swap_assign(&s[ix], tmp, size_data_point);
			}
		}
		else{
			if(prec == 32){
				single2double(&s[ix], (float*) tmp2, size_data_point);
			}
			else memcpy(&s[ix], tmp, bytes);
		}
		if(status < 0 && status != LIME_EOR ){
			return (-1);
		}
	}
//	printf("The final checksum is %#lx %#lx\n", (unsigned long) ans.suma,
//			(unsigned long) ans.sumb);
	return (0);
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
int read_lime_spinor (double * const s, char * filename, const int position,
		const int number_data_points, const int size_data_point) {

	FILE * ifs;
	int status = 0, getpos = -1;
	n_uint64_t bytes;
	char * header_type;
	LimeReader * limereader;
	n_uint64_t prec = 32;
	DML_Checksum checksum;

	if((ifs = fopen(filename, "r")) == (FILE*) NULL){
		fprintf(stderr, "Error opening file %s\n", filename);
		return (-1);
	}

	limereader = limeCreateReader(ifs);
	if(limereader == (LimeReader *) NULL){
		fprintf(stderr, "Unable to open LimeReader\n");
		return (-1);
	}
	while((status = limeReaderNextRecord(limereader)) != LIME_EOF ){
		if(status != LIME_SUCCESS ){
			fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n",
					status);
			status = LIME_EOF;
			break;
		}
		header_type = limeReaderType(limereader);

		if(strcmp("scidac-binary-data", header_type) == 0) getpos++;
//		printf("... found record of type %s pos = %d!\n", header_type, getpos);
		if(getpos == position) break;
	}
	if(status == LIME_EOF ){
		fprintf(stderr, "no scidac-binary-data record found in file %s\n",
				filename);
		limeDestroyReader(limereader);
		fclose(ifs);
		return (-1);
	}
	bytes = limeReaderBytes(limereader);
	if(bytes
			== number_data_points * (uint64_t) (size_data_point * sizeof(double))) prec =
			64;
	else if(bytes
			== number_data_points * (uint64_t) (size_data_point * sizeof(float))) prec =
			32;
	else{
		fprintf(stderr,
				"wrong length in eospinor: bytes = %lu, not %lu. Aborting read!\n",
				(unsigned long) bytes,
				(unsigned long) (number_data_points
						* (uint64_t) (size_data_point * sizeof(double))));
		return (-1);
	}
//	printf("# %lu Bit precision read\n", (unsigned long) prec);

	status = read_binary_spinor_data(s, limereader, prec, checksum,
			number_data_points, size_data_point);

	if(status < 0){
		fprintf(stderr,
				"LIME read error occured with status = %d while reading file %s!\n Aborting...\n",
				status, filename);
		exit(500);
	}

	limeDestroyReader(limereader);
	fclose(ifs);
	return (0);
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
//int write_checksum(DML_Checksum * c, char * filename) {
//
//  return(0);
//}
//int get_propagator_type(char * filename) {
//  FILE * ifs;
//  int status=0, ret=-1;
//  n_uint64_t bytes;
//  char * tmp;
//  LimeReader * limereader;
//
//  if((ifs = fopen(filename, "r")) == (FILE*)NULL) {
//    fprintf(stderr, "Error opening file %s\n", filename);
//    return(ret);
//  }
//
//  limereader = limeCreateReader( ifs );
//  if( limereader == (LimeReader *)NULL ) {
//    fprintf(stderr, "Unable to open LimeReader\n");
//    return(ret);
//  }
//  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
//    if(status != LIME_SUCCESS ) {
//      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
//      status = LIME_EOF;
//      break;
//    }
//    if(strcmp("etmc-propagator-type", limeReaderType(limereader)) == 0) break;
//  }
//  if(status == LIME_EOF) {
//    fprintf(stderr, "no etmc-propagator-type record found in file %s\n",filename);
//    limeDestroyReader(limereader);
//    fclose(ifs);
//    return(ret);
//  }
//  tmp = (char*) calloc(500, sizeof(char));
//  bytes = limeReaderBytes(limereader);
//  status = limeReaderReadData(tmp, &bytes, limereader);
//  limeDestroyReader(limereader);
//  fclose(ifs);
//  if(strcmp("DiracFermion_Sink", tmp) == 0) ret = 0;
//  else if(strcmp("DiracFermion_Source_Sink_Pairs", tmp) == 0) ret = 1;
//  else if(strcmp("DiracFermion_ScalarSource_TwelveSink", tmp) == 0) ret = 2;
//  else if(strcmp("DiracFermion_ScalarSource_FourSink", tmp) == 0) ret = 3;
//  free(tmp);
//  return(ret);
//}
//int read_lime_spinor(float * const s, char * filename, const int position,
//		     const int ts,
//		     const unsigned int T, const unsigned int LX, const unsigned int LY, const unsigned int LZ) {
//  FILE * ifs;
//  int status=0, getpos=-1;
//  n_uint64_t bytes;
//  char * header_type;
//  LimeReader * limereader;
//  n_uint64_t prec = 32;
//  DML_Checksum checksum;
//
//  if((ifs = fopen(filename, "r")) == (FILE*)NULL) {
//    fprintf(stderr, "Error opening file %s\n", filename);
//    return(-1);
//  }
//
//  limereader = limeCreateReader( ifs );
//  if( limereader == (LimeReader *)NULL ) {
//    fprintf(stderr, "Unable to open LimeReader\n");
//    return(-1);
//  }
//  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
//    if(status != LIME_SUCCESS ) {
//      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
//      status = LIME_EOF;
//      break;
//    }
//    header_type = limeReaderType(limereader);
//
//    if(strcmp("scidac-binary-data",header_type) == 0) getpos++;
//    printf("... found record of type %s pos = %d!\n", header_type, getpos);
//    if(getpos == position) break;
//  }
//  if(status == LIME_EOF) {
//    fprintf(stderr, "no scidac-binary-data record found in file %s\n",filename);
//    limeDestroyReader(limereader);
//    fclose(ifs);
//    return(-1);
//  }
//  bytes = limeReaderBytes(limereader);
//  if(bytes == LX*LY*LZ*T*(uint64_t)(24*sizeof(double))) prec = 64;
//  else if(bytes == LX*LY*LZ*T*(uint64_t)(24*sizeof(float))) prec = 32;
//  else {
//    fprintf(stderr, "wrong length in eospinor: bytes = %lu, not %lu. Aborting read!\n",
//	    (unsigned long)bytes, (unsigned long)(LX*LY*LZ*T*(uint64_t)(24*sizeof(double))));
//    return(-1);
//  }
//  printf("# %lu Bit precision read\n", (unsigned long)prec);
//
//  status = read_binary_spinor_data(s, limereader, prec, ts, checksum, &get_index, T, LX, LY, LZ);
//
//  if(status < 0) {
//    fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n",
//	    status, filename);
//    exit(500);
//  }
//
//  limeDestroyReader(limereader);
//  fclose(ifs);
//  return(0);
//}
//
//int write_propagator_type(const int type, char * filename) {
//
//  FILE * ofs = NULL;
//  LimeWriter * limewriter = NULL;
//  LimeRecordHeader * limeheader = NULL;
//  int status = 0;
//  int ME_flag=1, MB_flag=1;
//  char message[500];
//  n_uint64_t bytes;
//
//  ofs = fopen(filename, "w");
//
//  if(ofs == (FILE*)NULL) {
//    fprintf(stderr, "Could not open file %s for writing!\n Aboring...\n", filename);
//    exit(500);
//  }
//  limewriter = limeCreateWriter( ofs );
//  if(limewriter == (LimeWriter*)NULL) {
//    fprintf(stderr, "LIME error in file %s for writing!\n Aborting...\n", filename);
//    exit(500);
//  }
//
//  if(type == 0) {
//    sprintf(message,"DiracFermion_Sink");
//    bytes = strlen( message );
//  }
//  else if (type == 1) {
//    sprintf(message,"DiracFermion_Source_Sink_Pairs");
//    bytes = strlen( message );
//  }
//  else if (type == 2) {
//    sprintf(message,"DiracFermion_ScalarSource_TwelveSink");
//    bytes = strlen( message );
//  }
//  else if (type == 3) {
//    sprintf(message,"DiracFermion_ScalarSource_FourSink");
//    bytes = strlen( message );
//  }
//
//  limeheader = limeCreateHeader(MB_flag, ME_flag, "etmc-propagator-type", bytes);
//  status = limeWriteRecordHeader( limeheader, limewriter);
//  if(status < 0 ) {
//    fprintf(stderr, "LIME write header error %d\n", status);
//    exit(500);
//  }
//  limeDestroyHeader( limeheader );
//  limeWriteRecordData(message, &bytes, limewriter);
//
//  limeDestroyWriter( limewriter );
//  fclose(ofs);
//  fflush(ofs);
//  return(0);
//}
//
//int write_source_type(const int type, char * filename) {
//
//  FILE * ofs = NULL;
//  LimeWriter * limewriter = NULL;
//  LimeRecordHeader * limeheader = NULL;
//  int status = 0;
//  int ME_flag=1, MB_flag=1;
//  char message[500];
//  n_uint64_t bytes;
//
//  ofs = fopen(filename, "w");
//
//  if(ofs == (FILE*)NULL) {
//    fprintf(stderr, "Could not open file %s for writing!\n Aboring...\n", filename);
//    exit(500);
//  }
//  limewriter = limeCreateWriter( ofs );
//  if(limewriter == (LimeWriter*)NULL) {
//    fprintf(stderr, "LIME error in file %s for writing!\n Aborting...\n", filename);
//    exit(500);
//  }
//
//  sprintf(message,"DiracFermion_Source");
//  bytes = strlen( message );
//
//  limeheader = limeCreateHeader(MB_flag, ME_flag, "source-type", bytes);
//  status = limeWriteRecordHeader( limeheader, limewriter);
//  if(status < 0 ) {
//    fprintf(stderr, "LIME write header error %d\n", status);
//    exit(500);
//  }
//  limeDestroyHeader( limeheader );
//  limeWriteRecordData(message, &bytes, limewriter);
//
//  limeDestroyWriter( limewriter );
//  fclose(ofs);
//  fflush(ofs);
//  return(0);
//}
//
//int write_propagator_format(char * filename, const int prec, const int no_flavours,
//			    const int T, const int LX, const int LY, const int LZ) {
//  FILE * ofs = NULL;
//  LimeWriter * limewriter = NULL;
//  LimeRecordHeader * limeheader = NULL;
//  int status = 0;
//  int ME_flag=0, MB_flag=1;
//  char message[500];
//  n_uint64_t bytes;
//  /*   char * message; */
//
//
//  ofs = fopen(filename, "a");
//
//  if(ofs == (FILE*)NULL) {
//    fprintf(stderr, "Could not open file %s for writing!\n Aborting...\n", filename);
//    exit(500);
//  }
//  limewriter = limeCreateWriter( ofs );
//  if(limewriter == (LimeWriter*)NULL) {
//    fprintf(stderr, "LIME error in file %s for writing!\n Aborting...\n", filename);
//    exit(500);
//  }
//
//  sprintf(message, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<etmcFormat>\n<field>diracFermion</field>\n<precision>%d</precision>\n<flavours>%d</flavours>\n<lx>%d</lx>\n<ly>%d</ly>\n<lz>%d</lz>\n<lt>%d</lt>\n</etmcFormat>", prec, no_flavours, LX, LY, LZ, T);
//  bytes = strlen( message );
//  limeheader = limeCreateHeader(MB_flag, ME_flag, "etmc-propagator-format", bytes);
//  status = limeWriteRecordHeader( limeheader, limewriter);
//  if(status < 0 ) {
//    fprintf(stderr, "LIME write header error %d\n", status);
//    exit(500);
//  }
//  limeDestroyHeader( limeheader );
//  limeWriteRecordData(message, &bytes, limewriter);
//
//  limeDestroyWriter( limewriter );
//  fclose(ofs);
//  fflush(ofs);
//
//  return(0);
//}
//// **********
//// **********
//// **********
//// **********
//// **********
//
//// difference to read_lime_spinor: data is stored at the beginning of s.
//
//int read_lime_spinor_timeslice(double *const s, char *filename, const int position, const int ts, const unsigned int T, const unsigned int LX, const unsigned int LY, const unsigned int LZ) {
//  FILE * ifs;
//  int status=0, getpos=-1;
//  n_uint64_t bytes;
//  char * header_type;
//  LimeReader * limereader;
//  n_uint64_t prec = 32;
//  DML_Checksum checksum;
//
//  if((ifs = fopen(filename, "r")) == (FILE*)NULL) {
//    fprintf(stderr, "Error opening file %s\n", filename);
//    return(-1);
//  }
//
//  limereader = limeCreateReader( ifs );
//  if( limereader == (LimeReader *)NULL ) {
//    fprintf(stderr, "Unable to open LimeReader\n");
//    return(-1);
//  }
//  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
//    if(status != LIME_SUCCESS ) {
//      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
//      status = LIME_EOF;
//      break;
//    }
//    header_type = limeReaderType(limereader);
//    if(strcmp("scidac-binary-data",header_type) == 0) getpos++;
//    if(getpos == position) break;
//  }
//  if(status == LIME_EOF) {
//    fprintf(stderr, "no scidac-binary-data record found in file %s\n",filename);
//    limeDestroyReader(limereader);
//    fclose(ifs);
//    return(-1);
//  }
//  bytes = limeReaderBytes(limereader);
//  if(bytes == LX*LY*LZ*T*(uint64_t)(24*sizeof(double))) prec = 64;
//  else if(bytes == LX*LY*LZ*T*(uint64_t)(24*sizeof(float))) prec = 32;
//  else {
//    fprintf(stderr, "wrong length in eospinor: bytes = %llu, not %llu. Aborting read!\n",
//	    bytes, LX*LY*LZ*T*(uint64_t)(24*sizeof(double)));
//    return(-1);
//  }
//  fprintf(stderr, "# %llu Bit precision read\n", prec);
//
//  status = read_binary_spinor_data_timeslice(s, limereader, prec, ts, checksum, &get_index, T, LX, LY, LZ);
//
//  if(status < 0) {
//    fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n",
//	    status, filename);
//    exit(500);
//  }
//
//  limeDestroyReader(limereader);
//  fclose(ifs);
//  return(0);
//}
//
//// **********
//// **********
//// **********
//// **********
//// **********
//int write_source(double * const s, char * filename,
//		 const int append, const int prec,
//		 const int T, const int LX, const int LY, const int LZ) {
//
//  int err = 0;
//
//  write_source_type(0, filename);
//  write_propagator_format(filename, prec, 1, T, LX, LY, LZ);
//  err = write_lime_spinor(s, filename, 1, prec, T, LX, LY, LZ);
//  return(err);
//}
//
//int write_propagator(double * const s, char * filename,
//		     const int append, const int prec,
//		     const int T, const int LX, const int LY, const int LZ) {
//  int err = 0;
//
////   write_propagator_type(0, filename);
//  write_propagator_format(filename, prec, 1, T, LX, LY, LZ);
//  err = write_lime_spinor(s, filename, append, prec, T, LX, LY, LZ);
//  return(err);
//}
//
///* write two flavour operator to file */
//
//int write_double_propagator(double * const s, double * const r,
//			    char * filename, const int append, const int prec,
//			    const int T, const int LX, const int LY, const int LZ) {
//  int err = 0;
//
//  write_propagator_format(filename, prec, 2, T, LX, LY, LZ);
//  err = write_lime_spinor(s, filename, append, prec, T, LX, LY, LZ);
//  err += write_lime_spinor(r, filename, append, prec, T, LX, LY, LZ);
//  return(err);
//}
//
//DML_Checksum write_binary_spinor_data(double * const s, LimeWriter * limewriter,
//				      const int prec,
//				      const unsigned int T, const unsigned int LX,
//				      const unsigned int LY, const unsigned int LZ) {
//
//  int status=0;
//  double tmp[24];
//  float tmp2[24];
//  n_uint64_t bytes;
//  DML_Checksum ans;
//  DML_SiteRank rank;
//  int words_bigendian;
//  words_bigendian = big_endian();
//
//  DML_checksum_init(&ans);
//  rank = (DML_SiteRank) 0;
//
//  if(prec == 32) bytes = 24*sizeof(float);
//  else bytes = 24*sizeof(double);
//  for(unsigned int t = 0; t < T; t++) {
//    for(unsigned int z = 0; z < LZ; z++) {
//      for(unsigned int y = 0; y < LY; y++) {
//	for(unsigned int x = 0; x < LX; x++) {
//	  n_uint64_t ix = (t*LX*LY*LZ + x*LY*LZ + y*LZ + z)*(n_uint64_t)12;
//	  rank = (DML_SiteRank) (((t*LZ + z)*LY + y)*LX + x);
//	  if(!words_bigendian) {
//	    if(prec == 32) {
//	      byte_swap_assign_double2single((float*)tmp2, &s[2*ix], 24);
//	      DML_checksum_accum(&ans,rank,(char *) tmp2, bytes);
//	      status = limeWriteRecordData((void*)tmp2, &bytes, limewriter);
//	    }
//	    else {
//	      byte_swap_assign(tmp, &s[2*ix], 24);
//	      DML_checksum_accum(&ans,rank,(char *) tmp, bytes);
//	      status = limeWriteRecordData((void*)tmp, &bytes, limewriter);
//	    }
//	  }
//	  else {
//	    if(prec == 32) {
//	      double2single((float*)tmp2, &s[2*ix], 24);
//	      DML_checksum_accum(&ans,rank, (char *) tmp2, bytes);
//	      status = limeWriteRecordData((void*)tmp2, &bytes, limewriter);
//	    }
//	    else {
//	      status = limeWriteRecordData((void*) &s[2*ix], &bytes, limewriter);
//	      DML_checksum_accum(&ans,rank,(char *) &s[2*ix], bytes);
//	    }
//	  }
//	}
//      }
//    }
//  }
//  printf("The final checksum is %#lx %#lx\n", (unsigned long)ans.suma, (unsigned long)ans.sumb);
//  return(ans);
//}
//
///* if -1 < ts < T read only timeslice ts */
///* the other entries in s are untouched  */
//
//int read_binary_spinor_data(double * const s, LimeReader * limereader,
//			    const double prec, const int ts, DML_Checksum &ans,
//			    unsigned long int get_index(const int, const int, const int, const int, const int, const int),
//			    const unsigned int T, const unsigned int LX,
//			    const unsigned int LY, const unsigned int LZ) {
//  int status=0;
//  n_uint64_t bytes;
//  double tmp[24];
//  DML_SiteRank rank;
//  float tmp2[24];
//  int words_bigendian;
//  words_bigendian = big_endian();
//
//  DML_checksum_init(&ans);
//  rank = (DML_SiteRank) 0;
//
//  if(prec == 32) bytes = 24*sizeof(float);
//  else bytes = 24*sizeof(double);
//  for(unsigned int t = 0; t < T; t++){
//    if(ts > -1 && (unsigned int)abs(ts) < T) {
//      t = ts;
//      limeReaderSeek(limereader,(n_uint64_t)
//		     (t*LZ*LY*LX)*bytes,
//		     SEEK_SET);
//    }
//    for(unsigned int z = 0; z < LZ; z++){
//      for(unsigned int y = 0; y < LY; y++){
//	for(unsigned int x = 0; x < LX; x++){
//	  n_uint64_t ix = get_index(t, x, y, z, T, LX)*(n_uint64_t)12;
//	  rank = (DML_SiteRank) (((t*LZ + z)*LY + y)*LX + x);
//	  if(prec == 32) {
//	    status = limeReaderReadData(tmp2, &bytes, limereader);
//	    DML_checksum_accum(&ans,rank,(char *) tmp2, bytes);
//	  }
//	  else {
//	    status = limeReaderReadData(tmp, &bytes, limereader);
//	    DML_checksum_accum(&ans,rank,(char *) tmp, bytes);
//	  }
//	  if(!words_bigendian) {
//	    if(prec == 32) {
//	      byte_swap_assign_single2double(&s[2*ix], (float*)tmp2, 24);
//	    }
//	    else {
//	      byte_swap_assign(&s[2*ix], tmp, 24);
//	    }
//	  }
//	  else {
//	    if(prec == 32) {
//	      single2double(&s[2*ix], (float*)tmp2, 24);
//	    }
//	    else memcpy(&s[2*ix], tmp, bytes);
//	  }
//	  if(status < 0 && status != LIME_EOR) {
//	    return(-1);
//	  }
//	}
//      }
//    }
//    if(ts > -1 && ts < (int)T) {
//      t = T;
//    }
//  }
//  printf("The final checksum is %#lx %#lx\n", (unsigned long)ans.suma, (unsigned long)ans.sumb);
//  return(0);
//}
//
//
//int read_binary_spinor_data(float * const s, LimeReader * limereader,
//			    const double prec, const int ts, DML_Checksum &ans,
//			    unsigned long int get_index(const int, const int, const int, const int, const int, const int),
//			    const unsigned int T, const unsigned int LX,
//			    const unsigned int LY, const unsigned int LZ) {
//  int status=0;
//  n_uint64_t bytes;
//  double tmp[24];
//  DML_SiteRank rank;
//  float tmp2[24];
//  int words_bigendian;
//  words_bigendian = big_endian();
//
//  DML_checksum_init(&ans);
//  rank = (DML_SiteRank) 0;
//
//  if(prec == 32) bytes = 24*sizeof(float);
//  else return(-1);
//  for(unsigned int t = 0; t < T; t++){
//    if(ts > -1 && (unsigned int)abs(ts) < T) {
//      t = ts;
//      limeReaderSeek(limereader,(n_uint64_t)
//		     (t*LZ*LY*LX)*bytes,
//		     SEEK_SET);
//    }
//    for(unsigned int z = 0; z < LZ; z++){
//      for(unsigned int y = 0; y < LY; y++){
//	for(unsigned int x = 0; x < LX; x++){
//	  n_uint64_t ix = get_index(t, x, y, z, T, LX)*(n_uint64_t)12;
//	  rank = (DML_SiteRank) (((t*LZ + z)*LY + y)*LX + x);
//	  status = limeReaderReadData(tmp2, &bytes, limereader);
//	  DML_checksum_accum(&ans,rank,(char *) tmp2, bytes);
//	  if(!words_bigendian) {
//	    byte_swap_assign_singleprec(&s[2*ix], tmp2, 24);
//	  }
//	  else {
//	    memcpy(&s[2*ix], tmp2, bytes);
//	  }
//	  if(status < 0 && status != LIME_EOR) {
//	    return(-1);
//	  }
//	}
//      }
//    }
//    if(ts > -1 && ts < (int)T) {
//      t = T;
//    }
//  }
//  printf("The final checksum is %#lx %#lx\n", (unsigned long)ans.suma, (unsigned long)ans.sumb);
//  return(0);
//}
//
//
//
//// **********
//// **********
//// **********
//// **********
//// **********
//
//int read_binary_spinor_data_timeslice(double * const s, LimeReader * limereader,
//			    const double prec, const int ts, DML_Checksum &ans,
//			    unsigned long int get_index(const int, const int, const int, const int, const int, const int),
//			    const unsigned int T, const unsigned int LX,
//			    const unsigned int LY, const unsigned int LZ) {
//  int status=0;
//  n_uint64_t bytes;
//  double tmp[24];
//  DML_SiteRank rank;
//  float tmp2[24];
//  int words_bigendian;
//  words_bigendian = big_endian();
//
//  DML_checksum_init(&ans);
//  rank = (DML_SiteRank) 0;
//
//  if(prec == 32) bytes = 24*sizeof(float);
//  else bytes = 24*sizeof(double);
//  for(unsigned int t = 0; t < T; t++){
//    if(ts > -1 && abs(ts) < T) {
//      t = ts;
//      limeReaderSeek(limereader,(n_uint64_t)
//		     (t*LZ*LY*LX)*bytes,
//		     SEEK_SET);
//    }
//    else
//      {
//	fprintf(stderr, "Error: int read_binary_spinor_data_timeslice(...\n");
//	exit(EXIT_SUCCESS);
//      }
//    for(unsigned int z = 0; z < LZ; z++){
//      for(unsigned int y = 0; y < LY; y++){
//	for(unsigned int x = 0; x < LX; x++){
//
//
//
//	  // n_uint64_t ix = get_index(t, x, y, z, T, LX)*(n_uint64_t)12;
//	  n_uint64_t ix = get_index(0, x, y, z, T, LX)*(n_uint64_t)12;
//
//
//
//	  rank = (DML_SiteRank) (((t*LZ + z)*LY + y)*LX + x);
//	  if(prec == 32) {
//	    status = limeReaderReadData(tmp2, &bytes, limereader);
//	    DML_checksum_accum(&ans,rank,(char *) tmp2, bytes);
//	  }
//	  else {
//	    status = limeReaderReadData(tmp, &bytes, limereader);
//	    DML_checksum_accum(&ans,rank,(char *) tmp, bytes);
//	  }
//	  if(!words_bigendian) {
//	    if(prec == 32) {
//	      byte_swap_assign_single2double(&s[2*ix], (float*)tmp2, 24);
//	    }
//	    else {
//	      byte_swap_assign(&s[2*ix], tmp, 24);
//	    }
//	  }
//	  else {
//	    if(prec == 32) {
//	      single2double(&s[2*ix], (float*)tmp2, 24);
//	    }
//	    else memcpy(&s[2*ix], tmp, bytes);
//	  }
//	  if(status < 0 && status != LIME_EOR) {
//	    return(-1);
//	  }
//	}
//      }
//    }
//    if(ts > -1 && ts < T) {
//      t = T;
//    }
//  }
//  fprintf(stderr, "The final checksum is %#lx %#lx\n", ans.suma, ans.sumb);
//  return(0);
//}
// **********
// **********
// **********
// **********
// **********
