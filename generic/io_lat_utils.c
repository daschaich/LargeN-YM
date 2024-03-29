// -----------------------------------------------------------------
// Routines for gauge configuration I/O
// Works for most machines

// Wrappers for parallel I/O are in io_ansi.c and io_piofs.c

#include "generic_includes.h"
#include "../include/io_lat.h"
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#ifdef HAVE_QIO
#include <qio.h>
#endif

#define EPS 1e-6

#define PARALLEL 1   // Must evaluate to true
#define SERIAL 0     // Must evaluate to false

#define NODE_DUMP_ORDER 1
#define NATURAL_ORDER 0

#undef MAX_BUF_LENGTH
#define MAX_BUF_LENGTH 4096

/* Checksums
   The dataset from which each checksum is computed is the full gauge
   configuration for lattice files and for propagator files, the
   propagator for a single source spin-color combination.  Data in these
   files appear as a series of 32-bit floating point numbers.  We treat
   the 32-bit values as unsigned integers v(i) where i = 0, ..., N-1 ranges
   over the values in the order of appearance on the file and N is the
   total count of values in the data set.  The checksum is obtained by a
   combination of bit-rotations and exclusive or operations.  It is
   designed to be commutative and associative, unlike the BSD sum
   operation, so that the data set can be read in parallel with checksum
   contributions computed for each portion read, and then combined
   afterwards.  The sum29 checksum does a left bit rotation through i mod
   29 bits and forms an exclusive or with the accumulated checksum.  The
   sum31 checksum does the same thing, but with i mod 31 bits.

   In writing the file the bit rotation is done on the number as
   represented on the architecture and consequently as written on the
   file.  In reading and checking file integrity on an architecure with a
   relatively byte-reversed representation, byte reversal of the data
   must be done before doing the bit rotation and the resulting checksum
   must be compared with the checksum recorded on the file after
   byte-reversal.
*/
// For checksums we want a 32 bit unsigned int, for which
// we have u_int32type defined in include/int32type.h which is
// included in include/io_lat.h

#define SUCCESS  0
#define FAILURE -1
#define MAX_LINE_LENGTH 1024
#define MAX_TOKENS 512

#define TOL 0.0000001       // Tolerance for floating point checks
// -----------------------------------------------------------------


/*----------------------------------------------------------------------
   Routines for archive I/O
   -----------------------------------------------------------------------*/

int qcdhdr_get_str(char *s, QCDheader *hdr, char **q) {
  /* find a token and return the value */
  int i;
  for (i = 0; i<(char)(*hdr).ntoken; i++) {
    if (strcmp(s, (char *)(*hdr).token[i]) == 0) {
      *q = (*hdr).value[i];
      return(SUCCESS);
    }
  }
  *q = NULL;
  return FAILURE;
}

int qcdhdr_get_int(char *s, QCDheader *hdr, int *q) {
  char *p;
  qcdhdr_get_str(s, hdr,&p);
  if (p == NULL)
    return FAILURE;
  sscanf(p, "%d", q);
  return SUCCESS;
}
int qcdhdr_get_int32x(char *s, QCDheader *hdr, u_int32type *q) {
  char *p;
  int r;
  qcdhdr_get_str(s, hdr,&p);
  if (p==NULL) return (FAILURE);
  sscanf(p, "%x",&r);
  *q = r;
  return (SUCCESS);
}
int qcdhdr_get_float(char *s, QCDheader *hdr, Real *q) {
  char *p;
  qcdhdr_get_str(s, hdr,&p);
  if (p==NULL) return (FAILURE);
#if PRECISION == 1
  sscanf(p, "%f", q);
#else
  sscanf(p, "%lf", q);
#endif
  return (SUCCESS);
}

void error_exit(char *s) { printf("%s\n", s); terminate(1);}

void complete_U(float *u) {
  u[12] = u[ 2]*u[10] - u[ 4]*u[ 8] - u[ 3]*u[11] + u[ 5]*u[ 9];
  u[13] = u[ 4]*u[ 9] - u[ 2]*u[11] + u[ 5]*u[ 8] - u[ 3]*u[10];
  u[14] = u[ 4]*u[ 6] - u[ 0]*u[10] - u[ 5]*u[ 7] + u[ 1]*u[11];
  u[15] = u[ 0]*u[11] - u[ 4]*u[ 7] + u[ 1]*u[10] - u[ 5]*u[ 6];
  u[16] = u[ 0]*u[ 8] - u[ 2]*u[ 6] - u[ 1]*u[ 9] + u[ 3]*u[ 7];
  u[17] = u[ 2]*u[ 7] - u[ 0]*u[ 9] + u[ 3]*u[ 6] - u[ 1]*u[ 8];
}


int big_endian() {
  union  {
    long l;
    char c[sizeof (long)];
  } u;
  u.l = 1;
  return (u.c[sizeof (long) - 1] == 1);
}

QCDheader * qcdhdr_get_hdr(FILE *in)
{
#define MAX_LINE_LENGTH 1024
#define MAX_TOKENS 512
  char line[MAX_LINE_LENGTH];
  int n, len;
  QCDheader *hdr;
  char **tokens, **values;
  char *p, *q;

  /* Begin reading, and check for "BEGIN_HEADER" token */
  fgets(line, MAX_LINE_LENGTH, in);
  /*
  if (strcmp(line, "BEGIN_HEADER\n")!=0)
    error_exit("qcdhdr_get_hdr: Missing \"BEGIN_HEADER\"; punting \n");
  */
  /* Allocate space for QCDheader and its pointers */
  tokens = (char **) malloc(MAX_TOKENS*sizeof(char *));
  values = (char **) malloc(MAX_TOKENS*sizeof(char *));
  hdr = (QCDheader *) malloc(sizeof(QCDheader));
  (*hdr).token = tokens;
  (*hdr).value = values;

  /* Begin loop on tokens */
  n = 0;
  printf("reading Archive header:\n");
  while (1) {
    fgets(line, MAX_LINE_LENGTH, in);
    printf("%s", line);

    if (strcmp(line, "END_HEADER\n") == 0) break;

    /* Tokens are terminated by a space */
    q = strchr(line, (int)' ');

    /* Overwrite space with a terminating null */
    *q = '\0';
    len = (int)(q -  line);

    /* allocate space and copy the token in to it */
    p = (char *)malloc(len+1);
    (*hdr).token[n] = p;
    strcpy(p, line);

    q = strchr(++q, (int)'='); q++;
    len = strlen(q);
    q[len-1] = 0;
    p = (char *)malloc(len);
    (*hdr).value[n] = p;
    strcpy(p, q);
    n++;
  }
  (*hdr).ntoken = n;
  return (hdr);
}

/* Destroy header - for freeing up storage */
void qcdhdr_destroy_hdr(QCDheader *hdr) {
  int i;

  if (hdr == NULL)return;

  for (i = 0; i < hdr->ntoken; i++) {
    free(hdr->value[i]);
    free(hdr->token[i]);
  }

  free(hdr->token);
  free(hdr->value);
  free(hdr);
}

/*---------------------------------------------------------------------------*/
/* Convert (or copy) four single precision matrices to generic precision */
void f2d_4mat(fmatrix *a, matrix *b) {
  int dir, i, j;

  for (dir = 0; dir < 4; dir++) {
    for (i = 0; i < NCOL; i++)for (j = 0; j < NCOL; j++) {
      b[dir].e[i][j].real = a[dir].e[i][j].real;
      b[dir].e[i][j].imag = a[dir].e[i][j].imag;
    }
  }
}

/* Convert (or copy) four generic precision matrices to single precision */
void d2f_4mat(matrix *a, fmatrix *b) {
  int dir, i, j;

  for (dir = 0; dir < 4; dir++) {
    for (i = 0; i < NCOL; i++)for (j = 0; j < NCOL; j++) {
      b[dir].e[i][j].real = a[dir].e[i][j].real;
      b[dir].e[i][j].imag = a[dir].e[i][j].imag;
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void swrite_data(FILE* fp, void *src, size_t size,
                 char *myname, char *descrip) {

  if (fwrite(src, size, 1, fp) != 1) {
    printf("%s: Node %d %s write error %d\n", myname,
           this_node, descrip, errno);
    fflush(stdout);
    terminate(1);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void pwrite_data(FILE* fp, void *src, size_t size,
                 char *myname, char *descrip) {

  if (g_write(src, size, 1, fp) != 1) {
    printf("%s: Node %d %s descrip, write error %d\n", myname,
           this_node, descrip, errno);
    fflush(stdout);
    terminate(1);
  }
}
/*---------------------------------------------------------------------------*/
void pswrite_data(int parallel, FILE* fp, void *src, size_t size,
     char *myname, char *descrip) {

  if (parallel)pwrite_data(fp, src, size, myname, descrip);
  else        swrite_data(fp, src, size, myname, descrip);
}

/*---------------------------------------------------------------------------*/

int sread_data(FILE* fp, void *src, size_t size, char *myname, char *descrip) {
  if (g_read(src, size, 1, fp) != 1) {
    printf("%s: Node %d %s read error %d\n",
           myname, this_node, descrip, errno);
    fflush(stdout);
    return 1;
  }
  return 0;
}

/*---------------------------------------------------------------------------*/

int pread_data(FILE* fp, void *src, size_t size, char *myname, char *descrip) {
  if (g_read(src, size, 1, fp) != 1) {
    printf("%s: Node %d %s read error %d\n",
           myname, this_node, descrip, errno);
    fflush(stdout);
    return 1;
  }
  return 0;
}
/*---------------------------------------------------------------------------*/
int pread_byteorder(int byterevflag, FILE* fp, void *src, size_t size, char *myname, char *descrip)
{
  int status;

  status = pread_data(fp, src, size, myname, descrip);
  if (byterevflag==1)
    byterevn((int32type *)src, size/sizeof(int32type));
  return status;
}
/*---------------------------------------------------------------------------*/
int sread_byteorder(int byterevflag, FILE* fp, void *src, size_t size, char *myname, char *descrip)
{
  int status;

  status = sread_data(fp, src, size, myname, descrip);
  if (byterevflag==1)
    byterevn((int32type *)src, size/sizeof(int32type));
  return status;
}
/*---------------------------------------------------------------------------*/
int psread_data(int parallel, FILE* fp, void *src, size_t size,
     char *myname, char *descrip)
{
  if (parallel)return pread_data(fp, src, size, myname, descrip);
  else        return sread_data(fp, src, size, myname, descrip);
}
/*---------------------------------------------------------------------------*/
int psread_byteorder(int byterevflag, int parallel, FILE* fp,
          void *src, size_t size,
          char *myname, char *descrip)
{
  if (parallel)return pread_byteorder(byterevflag, fp, src, size, myname, descrip);
  else        return sread_byteorder(byterevflag, fp, src, size, myname, descrip);
}
/*---------------------------------------------------------------------------*/

/* This subroutine writes the gauge configuration header structure */
/* Parallel access version */
/* While the procedures for serial and parallel writing are
   identical, (the header is written only by node 0, no matter what),
   the file which is accessed can be opened either by all
   nodes in w_parallel or one node in w_serial.  We have to
   distinguish between these modes when writing */

void pwrite_gauge_hdr(FILE *fp, gauge_header *gh)
{

  char myname[] = "pwrite_gauge_hdr";

  pwrite_data(fp, (void *)&gh->magic_number, sizeof(gh->magic_number),
        myname, "magic_number");
  pwrite_data(fp, (void *)gh->dims, sizeof(gh->dims),
        myname, "dimensions");
  pwrite_data(fp, (void *)gh->time_stamp, sizeof(gh->time_stamp),
        myname, "time_stamp");
  pwrite_data(fp,&gh->order, sizeof(gh->order),
        myname, "order");

  /* Header byte length */

  gh->header_bytes = sizeof(gh->magic_number) + sizeof(gh->dims) +
    sizeof(gh->time_stamp) + sizeof(gh->order);

}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Write gauge configuration header structure -- serial version
void swrite_gauge_hdr(FILE *fp, gauge_header *gh) {
  char myname[] = "swrite_gauge_hdr";

  swrite_data(fp, (void *)&gh->magic_number, sizeof(gh->magic_number),
              myname, "magic_number");
  swrite_data(fp, (void *)gh->dims, sizeof(gh->dims),
              myname, "dimensions");
  swrite_data(fp, (void *)gh->time_stamp, sizeof(gh->time_stamp),
              myname, "time_stamp");
  swrite_data(fp, &gh->order, sizeof(gh->order),
              myname, "order");

  // Header byte length
  gh->header_bytes = sizeof(gh->magic_number) + sizeof(gh->dims)
                   + sizeof(gh->time_stamp) + sizeof(gh->order);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Write a data item to the gauge info file
int write_gauge_info_item(FILE *fpout,    /* ascii file pointer */
           char *keyword,   /* keyword */
           char *fmt,       /* output format -
                must use s, d, e, f, lu, or g */
           char *src,       /* address of starting data
             floating point data must be
             of type (Real) */
           int count,       /* number of data items if > 1 */
           int stride)      /* byte stride of data if
                                           count > 1 */
{

  int i, k, n;
  char *data;
  float tt;

  /* Check for valid keyword */
  for (i = 0; strlen(gauge_info_keyword[i])>0 &&
      strcmp(gauge_info_keyword[i], keyword) != 0; i++);
  if (strlen(gauge_info_keyword[i]) == 0)
    printf("write_gauge_info_item: WARNING: keyword %s not in table\n",
      keyword);

  /* Write keyword */

  fprintf(fpout, "%s =", keyword);

  /* Write count if more than one item */
  if (count > 1)
    fprintf(fpout, "[%d]", count);

  n = count; if (n == 0)n = 1;

  /* Write data */
  for (k = 0, data = (char *)src; k < n; k++, data += stride) {
    fprintf(fpout, " ");
    if (strstr(fmt, "s") != NULL)
      fprintf(fpout, fmt, data);
    else if (strstr(fmt, "d") != NULL)
      fprintf(fpout, fmt,*(int *)data);
    else if (strstr(fmt, "lu") != NULL)
      fprintf(fpout, fmt,*(unsigned long *)data);
    else if (strstr(fmt, "e") != NULL ||
        strstr(fmt, "f") != NULL ||
        strstr(fmt, "g") != NULL) {
      tt = *(Real *)data;
      fprintf(fpout, fmt, tt);
    }
    else {
      printf("write_gauge_info_item: Unrecognized data type %s\n", fmt);
      return 1;
    }
  }
  fprintf(fpout, "\n");
  return 0;
}

/*------------------------------------------------------------------------*/

/* Write a data item to a character string */
int sprint_gauge_info_item(
  char *string,    /* character string */
  size_t nstring,     /* string length */
  char *keyword,   /* keyword */
  char *fmt,       /* output format -
          must use s, d, e, f, or g */
  char *src,       /* address of starting data
          floating point data must be
          of type (Real) */
  int count,       /* number of data items if > 1 */
  int stride)      /* byte stride of data if
          count > 1 */
{

  int i, k, n;
  size_t bytes;
  char *data;
  float tt;

  /* Check for valid keyword */

  for (i = 0; strlen(gauge_info_keyword[i]) > 0 &&
      strcmp(gauge_info_keyword[i], keyword) != 0; i++);
  if (strlen(gauge_info_keyword[i]) == 0)
    printf("sprint_gauge_info_item: WARNING: keyword %s not in table\n",
      keyword);

  /* Write keyword */
  bytes = 0;

  snprintf(string, nstring-bytes, "%s =", keyword);
  bytes = strlen(string);
  if (bytes >= nstring)return 1;

  /* Write count if more than one item */
  if (count > 1) {
    snprintf(string+bytes, nstring-bytes, "[%d]", count);
    bytes = strlen(string);
    if (bytes >= nstring)return 1;
  }

  n = count; if (n == 0)n = 1;

  // Write data
  for (k = 0, data = (char *)src; k < n; k++, data += stride) {
    snprintf(string + bytes, nstring - bytes, " ");
    bytes = strlen(string);
    if (bytes >= nstring)
      return 1;

    if (strstr(fmt, "s") != NULL) {
      snprintf(string + bytes, nstring - bytes, fmt, data);
      bytes = strlen(string);
      if (bytes >= nstring)
        return 1;
    }
    else if (strstr(fmt, "d") != NULL) {
      snprintf(string + bytes, nstring - bytes, fmt,*(int *)data);
      bytes = strlen(string);
      if (bytes >= nstring)
        return 1;
    }
    else if (strstr(fmt, "lu") != NULL) {
      snprintf(string + bytes, nstring - bytes, fmt,*(unsigned long *)data);
      bytes = strlen(string);
      if (bytes >= nstring)
        return 1;
    }
    else if (strstr(fmt, "e") != NULL ||
        strstr(fmt, "f") != NULL ||
        strstr(fmt, "g") != NULL) {
      tt = *(Real *)data;
      snprintf(string + bytes, nstring - bytes, fmt, tt);
      bytes = strlen(string);
      if (bytes >= nstring)
        return 1;
    }
    else {
      printf("sprint_gauge_info_item: Unrecognized data type %s\n", fmt);
      return 1;
    }
  }
  snprintf(string+bytes, nstring-bytes, "\n");
  bytes = strlen(string);
  if (bytes >= nstring)return 1;

  return 0;
}

/*----------------------------------------------------------------------*/
/* Open, write, and close the ASCII info file */
void write_gauge_info_file(gauge_file *gf) {
  FILE *info_fp;
  gauge_header *gh;
  char info_filename[256];
  char sums[20];

  gh = gf->header;

  /* Construct header file name from lattice file name
   by adding filename extension to lattice file name */

  strcpy(info_filename, gf->filename);
  strcat(info_filename, ASCII_GAUGE_INFO_EXT);

  /* Open header file */
  if ((info_fp = g_open(info_filename, "w")) == NULL) {
      printf("write_gauge_info_file: Can't open ascii info file %s\n", info_filename);
      return;
    }

  // Write required information
  write_gauge_info_item(info_fp, "magic_number", "%d", (char *)&gh->magic_number, 0, 0);
  write_gauge_info_item(info_fp, "time_stamp", "\"%s\"", gh->time_stamp, 0, 0);
  sprintf(sums, "%x %x", gf->check.sum29, gf->check.sum31);
  write_gauge_info_item(info_fp, "checksums", "\"%s\"", sums, 0, 0);
  write_gauge_info_item(info_fp, "nx", "%d", (char *)&nx, 0, 0);
  write_gauge_info_item(info_fp, "ny", "%d", (char *)&ny, 0, 0);
  write_gauge_info_item(info_fp, "nz", "%d", (char *)&nz, 0, 0);
  write_gauge_info_item(info_fp, "nt", "%d", (char *)&nt, 0, 0);

  write_appl_gauge_info(info_fp);

  g_close(info_fp);

  printf("Wrote info file %s\n", info_filename);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set up the input gauge file and gauge header structures
gauge_file *setup_input_gauge_file(char *filename) {
  char myname[] = "setup_input_gauge_file";
  gauge_file *gf = malloc(sizeof *gf);
  gauge_header *gh = malloc(sizeof *gh);

  // Check that memory allocations succeeded
  if (gf == NULL) {
    printf("%s: Can't malloc gf\n", myname);
    terminate(1);
  }
  if (gh == NULL) {
    printf("%s: Can't malloc gh\n", myname);
    terminate(1);
  }

  gf->filename = filename;

  // Make sure compilation gave us a 32 bit integer type
  assert(sizeof(int32type) == 4);

  gf->header = gh;
  gf->rank2rcv = NULL;
  gf->check.sum29 = 0;
  gf->check.sum31 = 0;

  return gf;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set up the output gauge file and gauge header structure
gauge_file *setup_output_gauge_file() {
  char myname[] = "setup_output_gauge_file";
  gauge_file *gf = malloc(sizeof *gf);
  gauge_header *gh = malloc(sizeof *gh);
  time_t time_stamp;
  int i;

  // Check that memory allocations succeeded
  if (gf == NULL) {
    printf("%s: Can't malloc gf\n", myname);
    terminate(1);
  }
  if (gh == NULL) {
    printf("%s: Can't malloc gh\n", myname);
    terminate(1);
  }

  // Make sure compilation gave us a 32 bit integer type
  assert(sizeof(int32type) == 4);

  /* Load header pointer and file name */
  gf->header = gh;

  /* Initialize */
  gf->check.sum29 = 0;
  gf->check.sum31 = 0;

  /* Load header values */
  gh->magic_number = GAUGE_VERSION_NUMBER;

  gh->dims[0] = nx;
  gh->dims[1] = ny;
  gh->dims[2] = nz;
  gh->dims[3] = nt;

  // Get date and time stamp using local time on node 0
  if (this_node == 0) {
    time(&time_stamp);
    strcpy(gh->time_stamp, ctime(&time_stamp));
    /* For aesthetic reasons, don't leave trailing junk bytes here to be
       written to the file */
    for (i = strlen(gh->time_stamp) + 1; i < (int)sizeof(gh->time_stamp); i++)
      gh->time_stamp[i] = '\0';

    // Remove trailing end-of-line character
    if (gh->time_stamp[strlen(gh->time_stamp) - 1] == '\n')
      gh->time_stamp[strlen(gh->time_stamp) - 1] = '\0';
  }

  // Broadcast to all nodes
  broadcast_bytes(gh->time_stamp, sizeof(gh->time_stamp));

  return gf;
}

/*---------------------------------------------------------------------------*/
/* Read checksum and compare.  It is assumed that the file is already
   correctly positioned.

   Should be called only by one node */
void read_checksum(int parallel, gauge_file *gf, gauge_check *test_gc) {
  char myname[] = "read_checksum";

  /* Read checksums with byte reversal */
  if (psread_byteorder(gf->byterevflag, parallel, gf->fp,
      &gf->check.sum29, sizeof(gf->check.sum29), myname, "checksum")!=0)
    terminate(1);
  if (psread_byteorder(gf->byterevflag, parallel, gf->fp,
      &gf->check.sum31, sizeof(gf->check.sum31), myname, "checksum")!=0)
    terminate(1);

  if (gf->check.sum29 != test_gc->sum29 ||
      gf->check.sum31 != test_gc->sum31)
    printf("%s: Checksum violation. Computed %x %x.  Read %x %x.\n",
           myname, test_gc->sum29, test_gc->sum31,
           gf->check.sum29, gf->check.sum31);
  else {
    printf("Checksums %x %x OK\n", gf->check.sum29, gf->check.sum31);
    fflush(stdout);
  }
}

/*---------------------------------------------------------------------------*/
/* Write checksum to lattice file.  It is assumed that the file
   is already correctly positioned.

   Should be called only by one node */
void write_checksum(int parallel, gauge_file *gf) {
  char myname[] = "write_checksum";

  pswrite_data(parallel, gf->fp, &gf->check.sum29, sizeof(gf->check.sum29),
               myname, "checksum");
  pswrite_data(parallel, gf->fp, &gf->check.sum31, sizeof(gf->check.sum31),
               myname, "checksum");
  printf("Checksums %x %x\n", gf->check.sum29, gf->check.sum31);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
/* Subroutine for reading site list from gauge configuration file */
/* Only node 0 reads this list, so same for parallel and serial reading */
void read_site_list(int parallel, gauge_file *gf) {
  /* All nodes allocate space for site list table, if file is not in
     natural order */

  if (gf->header->order != NATURAL_ORDER) {
    gf->rank2rcv = (int32type *)malloc(volume*sizeof(int32type));
    if (gf->rank2rcv == NULL) {
      printf("read_site_list: Can't malloc rank2rcv table\n");
      terminate(1);
    }

    // Only node 0 reads the site list
    if (this_node == 0) {
      /* Reads receiving site coordinate if file is not in natural order */
      if (parallel) {
        if ((int)g_read(gf->rank2rcv, sizeof(int32type), volume, gf->fp) != volume) {
          printf("read_site_list: Node %d site list read error %d\n",
                 this_node, errno);
          terminate(1);
        }
      }
      else {
        if ((int)g_read(gf->rank2rcv, sizeof(int32type), volume, gf->fp) != volume) {
          printf("read_site_list: Node %d site list read error %d\n",
                 this_node, errno);
          terminate(1);
        }
      }

      if (gf->byterevflag == 1)
        byterevn(gf->rank2rcv, volume);
    }

    // Broadcast result to all nodes
    broadcast_bytes((char *)gf->rank2rcv, volume * sizeof(int32type));
  }
  else
    gf->rank2rcv = NULL;  /* If no site list */
}

int read_gauge_hdr(gauge_file *gf, int parallel) {
  /* parallel = 1 (TRUE) if all nodes are accessing the file */
  /*            0        for access from node 0 only */

  FILE *fp = gf->fp;
  gauge_header *gh = gf->header;
  int32type tmp, btmp;
  int j, stat, byterevflag = 0;
  char myname[] = "read_gauge_hdr";

  // Read header, do byte reversal if necessary, and check consistency
  // Read and verify magic number
  stat = psread_data(parallel, fp,&gh->magic_number, sizeof(gh->magic_number),
                     myname, "magic number");
  if (stat != 0)
    terminate(1);

  tmp = gh->magic_number;
  btmp = gh->magic_number;
  byterevn((int32type *)&btmp, 1);

  /* See if header chunk is BEGI = 1111836489 for big endian
     or the byte reverse 1229407554 for little endian **/
  if (tmp == GAUGE_VERSION_NUMBER)
    byterevflag = 0;
  else if (btmp == GAUGE_VERSION_NUMBER) {
    byterevflag = 1;
    gh->magic_number = btmp;
//    printf("Reading with byte reversal\n");
    if (sizeof(float) != sizeof(int32type)) {
      printf("%s: Can't byte reverse\n", myname);
      printf("requires size of int32type(%d) = size of float(%d)\n",
             (int)sizeof(int32type), (int)sizeof(float));
      terminate(1);
    }
  }
  else if (tmp == LIME_MAGIC_NO || btmp == LIME_MAGIC_NO) {
    // LIME format suggests a SciDAC file
    // Print error, set flag and return
    printf("%s: Reading as a SciDAC formatted file\n", myname);
    gh->magic_number = LIME_MAGIC_NO;
    return 0;
  }
  else {
    // End of the road
    printf("%s: Unrecognized magic number in gauge header\n", myname);
    printf("Expected %x but read %x\n", GAUGE_VERSION_NUMBER, tmp);
    printf("Expected %s but read %d\n", (char *)GAUGE_VERSION_NUMBER, tmp);
    terminate(1);
    return byterevflag;
  }

  // Read and process header information
  // Get lattice dimensions
  if (psread_byteorder(byterevflag, parallel, fp, gh->dims, sizeof(gh->dims),
                       myname, "dimensions") != 0)
    terminate(1);

  // Check lattice dimensions for consistency
  if (gh->dims[0] != nx ||
      gh->dims[1] != ny ||
      gh->dims[2] != nz ||
      gh->dims[3] != nt) {
    /* So we can use this routine to discover the dimensions,
       we provide that if nx = ny = nz = nt = -1 initially
       we don't die */
    if (nx != -1 || ny != -1 || nz != -1 || nt != -1) {
      printf("%s: Incorrect lattice dimensions ", myname);
      for (j = 0; j < 4; j++)
        printf("%d ", gh->dims[j]);
      printf("\n");
      fflush(stdout);
      terminate(1);
    }
    else {
      nx = gh->dims[0];
      ny = gh->dims[1];
      nz = gh->dims[2];
      nt = gh->dims[3];
      volume = nx * ny * nz * nt;
    }
  }

  // Read date and time stamp
  if (psread_data(parallel, fp, gh->time_stamp, sizeof(gh->time_stamp),
                  myname, "time stamp") != 0)
    terminate(1);

  // Read header byte length
  gh->header_bytes = sizeof(gh->magic_number) + sizeof(gh->dims)
                   + sizeof(gh->time_stamp) + sizeof(gh->order);

  // Read data order
  if (psread_byteorder(byterevflag, parallel, fp,&gh->order, sizeof(gh->order),
                       myname, "order parameter") != 0)
    terminate(1);

  return byterevflag;
}

/*---------------------------------------------------------------------------*/

/* Write site list - only for checkpoint files */
void write_site_list(FILE *fp, gauge_header *gh) {
  off_t offset;
  int i, buf_length;
  register site *s;
  int32type coords, *cbuf;

  /* All nodes write their site coordinate list in sequential
     blocks after the header.  The list is in the order of appearance
     in the lattice array.  Node 0 writes to the first block
     followed by node 1, etc.  The result is a contiguous table
     that can be used to remap the data to the corresponding space-time
     coordinate */

  /* Location of site list for this node */
  offset = gh->header_bytes + sizeof(int32type) * sites_on_node * this_node;

  cbuf = malloc(sites_on_node * sizeof(int32type));
  if (cbuf == NULL) {
    printf("write_site_list: node %d can't malloc cbuf\n", this_node);
    fflush(stdout);terminate(1);
  }

  if (g_seek(fp, offset, SEEK_SET) < 0) {
    printf("write_site_list: node %d g_seek %ld failed errno %d\n",
        this_node, (long)offset, errno);
    fflush(stdout);terminate(1);
  }

  buf_length = 0;

  FORALLSITES(i, s) {
    /* Encode the space-time coordinate vector as a 32-bit integer */
    coords = nx*(ny*(nz*s->t + s->z) + s->y) + s->x;
    cbuf[buf_length] = coords;
    buf_length++;
  }

  if ((int)g_write(cbuf, sizeof(int32type), sites_on_node, fp) != sites_on_node) {
    printf("write_site_list: Node %d coords write error %d\n",
        this_node, errno);fflush(stdout);terminate(1);
  }

  free(cbuf);
}

/*---------------------------------------------------------------------------*/

/* Open a file for parallel writing */
gauge_file *parallel_open(int order, char *filename) {
  /* All nodes open the same filename */
  /* Returns a file structure describing the opened file */

  /* order = NATURAL_ORDER for coordinate natural order
           = NODE_DUMP_ORDER for node-dump order */

  FILE *fp;
  gauge_file *gf;
  gauge_header *gh;

  /* Set up gauge file and gauge header structures and load header values */
  gf = setup_output_gauge_file();
  gh = gf->header;

  gh->order = order;

  /* All nodes open the requested file */
  fp = g_open(filename, "wb");
  g_sync();     // Make sure everyone has opened before attempting to write
  if (fp == NULL) {
    printf("parallel_open: Node %d can't open file %s, error %d\n",
           this_node, filename, errno);fflush(stdout);terminate(1);
  }

  /* Node 0 writes the header */
  if (this_node == 0)
    pwrite_gauge_hdr(fp, gh);

  broadcast_bytes((char *)&gh->header_bytes, sizeof(gh->header_bytes));

  /* All nodes write site list to file if order is not natural */
  if (order != NATURAL_ORDER)
    write_site_list(fp, gh);

  /* Assign values to file structure */
  gf->fp          = fp;
  gf->filename    = filename;
  gf->byterevflag = 0;            /* Not used for writing */
  gf->parallel    = 1;            /* File opened in parallel */

  return gf;
}

/*---------------------------------------------------------------------------*/

/* Position gauge configuration file for writing in parallel */
/* Returns pointer to malloc'ed write buffer */
/* gf = file descriptor as opened by w_checkpoint_i */
fmatrix *w_parallel_setup(gauge_file *gf, off_t *checksum_offset) {
  FILE *fp;
  fmatrix *lbuf = malloc(MAX_BUF_LENGTH * 4 * sizeof(*lbuf));

  off_t offset ;           /* File stream pointer */
  off_t gauge_node_size;   /* Size of a gauge configuration block for
                              all sites on one node */
  off_t coord_list_size;   /* Size of coordinate list in bytes */
  off_t head_size;         /* Size of header plus coordinate list */
  off_t gauge_check_size;  /* Size of checksum */
  char myname[] = "w_parallel_setup";

  if (!gf->parallel)
    printf("%s: Attempting parallel write to serial file.\n", myname);

  if (lbuf == NULL) {
    printf("%s: Node %d can't malloc lbuf\n", myname, this_node);
    fflush(stdout);
    terminate(1);
  }

  fp = gf->fp;

  gauge_node_size = sites_on_node*4*sizeof(fmatrix);

  if (gf->header->order == NATURAL_ORDER)coord_list_size = 0;
  else coord_list_size = sizeof(int32type)*volume;
  head_size = gf->header->header_bytes + coord_list_size;
  *checksum_offset = head_size;
  gauge_check_size = sizeof(gf->check.sum29) + sizeof(gf->check.sum31);

  offset = head_size + gauge_check_size;

  /* Each node writes its gauge configuration values */

  offset += gauge_node_size*this_node;

  if (g_seek(fp, offset, SEEK_SET) < 0) {
      printf("%s: Node %d g_seek %ld failed error %d file %s\n",
       myname, this_node, (long)offset, errno, gf->filename);
      fflush(stdout);terminate(1);
    }

  return lbuf;
}

/*-----------------------------------------------------------------------*/

/* Open a file for parallel writing in natural order */
gauge_file *w_parallel_i(char *filename) {
  /* All nodes open the same filename */
  /* Returns a file structure describing the opened file */
  return parallel_open(NATURAL_ORDER, filename);
}


/*---------------------------------------------------------------------------*/
/* Open a file for parallel writing in node-dump order */
gauge_file *w_checkpoint_i(char *filename) {
  /* All nodes open the same filename */
  /* Returns a file structure describing the opened file */

  return parallel_open(NODE_DUMP_ORDER, filename);

} /* w_checkpoint_i */

/*---------------------------------------------------------------------------*/

/* Close the file and free associated structures */
void w_serial_f(gauge_file *gf) {
  g_sync();
  if (this_node == 0) {
    if (gf->parallel)
      printf("w_serial_f: Attempting serial close on parallel file \n");

    g_close(gf->fp);
  }

  /* Node 0 writes ascii info file */
  if (this_node == 0)
    write_gauge_info_file(gf);

  /* Do not free gf and gf->header so calling program can use them */
}

/*---------------------------------------------------------------------------*/

gauge_file *r_serial_i(char *filename) {
  /* Returns file descriptor for opened file */

  gauge_header *gh;
  gauge_file *gf;
  FILE *fp;
  int byterevflag;
  char editfilename[513];

  /* All nodes set up a gauge file and gauge header structure for reading */

  gf = setup_input_gauge_file(filename);
  gh = gf->header;

  /* File opened for serial reading */
  gf->parallel = 0;

  /* Node 0 alone opens the file and reads the header */
  g_sync();

  if (this_node == 0) {
    fp = g_open(filename, "rb");
    if (fp == NULL) {
      /* If this is a partition format SciDAC file the node 0 name
         has an extension ".vol0000".  So try again. */
      printf("r_serial_i: Node %d can't open file %s, error %d\n",
          this_node, filename, errno);fflush(stdout);
      strncpy(editfilename, filename, 504);
      editfilename[504] = '\0';  /* Just in case of truncation */
      strcat(editfilename, ".vol0000");
      printf("r_serial_i: Trying SciDAC partition volume %s\n", editfilename);
      fp = g_open(editfilename, "rb");
      if (fp == NULL) {
        printf("r_serial_i: Node %d can't open file %s, error %d\n",
               this_node, editfilename, errno);
        fflush(stdout);
        terminate(1);
      }
      printf("r_serial_i: Open succeeded\n");
    }

    gf->fp = fp;

    byterevflag = read_gauge_hdr(gf, SERIAL);
  }

  else gf->fp = NULL;  /* The other nodes don't know about this file */

  /* Broadcast the byterevflag from node 0 to all nodes */

  broadcast_bytes((char *)&byterevflag, sizeof(byterevflag));
  gf->byterevflag = byterevflag;

  /* Node 0 broadcasts the header structure to all nodes */

  broadcast_bytes((char *)gh, sizeof(gauge_header));

  /* No further processing here if this is a SciDAC file */
  if (gh->magic_number == LIME_MAGIC_NO)
    return gf;

  /* Read site list and broadcast to all nodes */

  read_site_list(SERIAL, gf);

  return gf;
}

/*----------------------------------------------------------------------*/

/* Close the file and free associated structures */
void r_serial_f(gauge_file *gf) {
  g_sync();
  if (this_node == 0) {
    if (gf->parallel)
      printf("r_serial_f: Attempting serial close on parallel file \n");

    g_close(gf->fp);
  }

  if (gf->rank2rcv != NULL)
    free(gf->rank2rcv);

  /* Do not free gf and gf->header so calling program can use them */
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Close file (if still active) and release header and file structures
void w_parallel_f(gauge_file *gf) {

  g_sync();
  if (gf->fp != NULL) {
    if (!gf->parallel)
      printf("w_parallel_f: Attempting parallel close on serial file.\n");

    g_close(gf->fp);
    gf->fp = NULL;
  }

  /* Node 0 writes ascii info file */
  if (this_node == 0)
    write_gauge_info_file(gf);

  /* Do not free gf and gf->header so calling program can use them */
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void r_parallel_f(gauge_file *gf) {
  /* Close file (if active) and release header and file structures */

  g_sync();
  if (gf->fp != NULL) {
    if (!gf->parallel)
      printf("r_parallel_f: Attempting parallel close on serial file.\n");
    g_close(gf->fp);
    gf->fp = NULL;
  }

  // Do not free gf and gf->header so calling program can use them
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read lattice dimensions from a binary file and close the file
void read_lat_dim_gf(char *filename, int *ndim, int dims[]) {
  gauge_file *gf;
  int i;

  // Only four dimensions here
  *ndim = 4;

  // Open the file
  nx = -1;
  ny = -1;
  nz = -1;
  nt = -1;
  gf = r_serial_i(filename);

  for (i = 0; i < *ndim; i++)
    dims[i] = gf->header->dims[i];

  r_serial_f(gf);
}
// -----------------------------------------------------------------
