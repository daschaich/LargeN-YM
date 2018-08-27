// -----------------------------------------------------------------
// Functions for reading and writing through QIO
#include "generic_includes.h"
#include <qio.h>
#include "../include/io_lat.h"
#include "../include/io_scidac.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Conversion of types between prevailing and fixed precision
// Real type
static void f2p_real(Real *dest, float *src) {
  *dest = *src;
}

static void p2f_real(float *dest, Real *src) {
  *dest = *src;
}

// Complex type
static void f2p_complex(complex *dest, fcomplex *src) {
  dest->real = src->real;
  dest->imag = src->imag;
}

static void p2f_complex(fcomplex *dest, complex *src) {
  dest->real = src->real;
  dest->imag = src->imag;
}

static void d2p_complex(complex *dest, dcomplex *src) {
  dest->real = src->real;
  dest->imag = src->imag;
}

static void p2d_complex(dcomplex *dest, complex *src) {
  dest->real = src->real;
  dest->imag = src->imag;
}

// Color vector
static void f2p_vec(vector_f *dest, fvector_f *src) {
  int i;

  for (i = 0; i < NCOL; i++) {
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

static void p2f_vec(fvector_f *dest, vector_f *src) {
  int i;

  for (i = 0; i < NCOL; i++) {
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

static void d2p_vec(vector_f *dest, dvector_f *src) {
  int i;

  for (i = 0; i < NCOL; i++) {
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

static void p2d_vec(dvector_f *dest, vector_f *src) {
  int i;

  for (i = 0; i < NCOL; i++) {
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

// Color matrix
static void f2p_mat(matrix_f *dest, fmatrix_f *src) {
  int i, j;

  for (i = 0; i < NCOL; i++) {
    for(j = 0; j < NCOL; j++) {
      dest->e[i][j].real = src->e[i][j].real;
      dest->e[i][j].imag = src->e[i][j].imag;
    }
  }
}

static void p2f_mat(fmatrix_f *dest, matrix_f *src) {
  int i, j;

  for (i = 0; i < NCOL; i++) {
    for(j = 0; j < NCOL; j++) {
      dest->e[i][j].real = src->e[i][j].real;
      dest->e[i][j].imag = src->e[i][j].imag;
    }
  }
}

static void d2p_mat(matrix_f *dest, dmatrix_f *src) {
  int i, j;

  for (i = 0; i < NCOL; i++) {
    for(j = 0; j < NCOL; j++) {
      dest->e[i][j].real = src->e[i][j].real;
      dest->e[i][j].imag = src->e[i][j].imag;
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Factory function for moving data from site structure to output buffer
// arg contains pointer to field offset value
#define make_vget_from_site(P, C, T, FIXTYPE, VARTYPE, FUNC) \
void vget_##P##N_##T##_from_site(char *buf, size_t index, int count, \
                                 void *arg) { \
 \
  int i; \
  FIXTYPE *dest = (FIXTYPE *)buf; \
  field_offset src_off = *((field_offset *)arg); \
  site *s = &lattice[index]; \
  VARTYPE *src = (VARTYPE *)F_PT(s,src_off); \
  for (i = 0; i < count; i++) \
    FUNC(dest + i, src + i); \
}

// Factory function for moving data from a field to the output buffer
// arg contains pointer to the data array
#define make_vget_from_field(P, C, T, FIXTYPE, VARTYPE, FUNC) \
void vget_##P##N_##T##_from_field(char *buf, size_t index, int count, \
                                  void *arg) { \
 \
  int i; \
  FIXTYPE *dest = (FIXTYPE *)buf; \
  VARTYPE *src_pt = (VARTYPE *)arg; \
  VARTYPE *src = src_pt + index * count; \
  for (i = 0; i < count; i++) \
    FUNC(dest + i, src + i); \
}

#define make_vget(P, C,  T, FIXTYPE, VARTYPE, FUNC) \
 make_vget_from_site(P, C,  T, FIXTYPE, VARTYPE, FUNC); \
 make_vget_from_field(P, C,  T, FIXTYPE, VARTYPE, FUNC);

// Factory function for moving data from input buffer to site structure
// arg contains pointer to field offset value
#define make_vput_to_site(P, C, T, FIXTYPE, VARTYPE, FUNC) \
void vput_##P##N_##T##_to_site(char *buf, size_t index, int count, \
                               void *arg) { \
 \
  int i; \
  FIXTYPE *src = (FIXTYPE *)buf; \
  field_offset dest_off = *((field_offset *)arg); \
  site *s = &lattice[index]; \
  VARTYPE *dest = (VARTYPE *)F_PT(s, dest_off); \
  for (i = 0; i < count; i++) \
    FUNC(dest + i, src + i); \
}

// Factory function for moving data from input buffer to field
// arg contains pointer to the data array
#define make_vput_to_field(P, C, T, FIXTYPE, VARTYPE, FUNC) \
void vput_##P##N_##T##_to_field(char *buf, size_t index, int count, \
                                void *arg) { \
 \
  int i; \
  FIXTYPE *src = (FIXTYPE *)buf; \
  VARTYPE *dest_pt = (VARTYPE *)arg; \
  VARTYPE *dest = dest_pt + index * count; \
  for (i = 0; i < count; i++) \
    FUNC(dest + i, src + i); \
}

#define make_vput(P, C, T, FIXTYPE, VARTYPE, FUNC) \
 make_vput_to_site(P, C, T, FIXTYPE, VARTYPE, FUNC); \
 make_vput_to_field(P, C, T, FIXTYPE, VARTYPE, FUNC);

// Single precision
make_vget(F,     , R, float,     Real,     p2f_real);
make_vget(F,     , C, fcomplex,  complex,  p2f_complex);
make_vget(F, NCOL, V, fvector_f, vector_f, p2f_vec);
make_vget(F, NCOL, M, fmatrix_f, matrix_f, p2f_mat);

make_vput(F,     , R, float,     Real,     f2p_real);
make_vput(F,     , C, fcomplex,  complex,  f2p_complex);
make_vput(F, NCOL, V, fvector_f, vector_f, f2p_vec);
make_vput(F, NCOL, M, fmatrix_f, matrix_f, f2p_mat);

// Double precision
make_vget(D,     , C, dcomplex,  complex,  p2d_complex);
make_vget(D, NCOL, V, dvector_f, vector_f, p2d_vec);

make_vput(D,     , C, dcomplex,  complex,  d2p_complex);
make_vput(D, NCOL, V, dvector_f, vector_f, d2p_vec);
make_vput(D, NCOL, M, dmatrix_f, matrix_f, d2p_mat);

// Write site structure data, assuming output precision is single
#define make_write_all_from_site(P, PSTRING, C, CVAL, SVAL, T, TYPESTRING, \
                                 FIXTYPE, VARTYPE, MYREAL) \
int write_##P##N_##T##_from_site(QIO_Writer *outfile, QIO_String *xml_out, \
                                 field_offset src, int count) { \
 \
  int status, datum_size = sizeof(FIXTYPE), word_size = sizeof(MYREAL); \
  char qdptype[] = TYPESTRING, prec[] = PSTRING; \
  QIO_RecordInfo *rec_info; \
 \
  /* Create the record info for the field */ \
  rec_info = QIO_create_record_info(QIO_FIELD, 0, 0, 0, qdptype, prec, \
                                    CVAL, SVAL, datum_size, count); \
 \
  /* Write the record for the field */ \
  status = QIO_write(outfile, rec_info, xml_out, \
                     vget_##P##N_##T##_from_site, \
                     count * datum_size, word_size, (void *)&src); \
  if (status != QIO_SUCCESS) \
    return 1; \
 \
  QIO_destroy_record_info(rec_info); \
  return 0; \
}

// Write a time slice of site structure data
// Assume output precision is single
#define make_write_tslice_from_site(P, PSTRING, C, CVAL, SVAL, T, TYPESTRING, \
                                    FIXTYPE, VARTYPE, MYREAL) \
int write_##P##N_##T##_timeslice_from_site(QIO_Writer *outfile, \
                                           QIO_String *xml_out, \
                                           field_offset src, \
                                           int count, int t0) { \
 \
  char qdptype[] = TYPESTRING, prec[] = PSTRING; \
  int status, datum_size = sizeof(FIXTYPE), word_size = sizeof(MYREAL); \
  int lower[4] = {0, 0, 0, t0}, upper[4] = {nx - 1, ny - 1, nz - 1, t0}; \
  QIO_RecordInfo *rec_info; \
 \
  /* Create the record info for the field */ \
  rec_info = QIO_create_record_info(QIO_HYPER, lower, upper, 4, qdptype, \
                                    prec, CVAL, SVAL, datum_size, count); \
 \
  /* Write the record for the field */ \
  status = QIO_write(outfile, rec_info, xml_out, \
                     vget_##P##N_##T##_from_site, \
                     count * datum_size, word_size, (void *)&src); \
  if (status != QIO_SUCCESS) \
    return 1; \
 \
  QIO_destroy_record_info(rec_info); \
  return 0; \
}

// Write field data
#define make_write_all_from_field(P, PSTRING, C, CVAL, SVAL, T, TYPESTRING, \
                                  FIXTYPE, VARTYPE, MYREAL) \
int write_##P##N_##T##_from_field(QIO_Writer *outfile, QIO_String *xml_out, \
                                  VARTYPE *src, int count) { \
 \
  char qdptype[] = TYPESTRING, prec[] = PSTRING; \
  int status, datum_size = sizeof(FIXTYPE), word_size = sizeof(MYREAL); \
  QIO_RecordInfo *rec_info; \
 \
  /* Create the record info for the field */ \
  rec_info = QIO_create_record_info(QIO_FIELD, 0, 0, 0, qdptype, prec, \
                                    CVAL, SVAL, datum_size, count); \
 \
  /* Write the record for the field */ \
  status = QIO_write(outfile, rec_info, xml_out, \
                     vget_##P##N_##T##_from_field, \
                     count * datum_size, word_size, (void *)src); \
  if (status != QIO_SUCCESS) \
    return 1; \
 \
  QIO_destroy_record_info(rec_info); \
  return 0; \
}

// Write a time slice of field data
#define make_write_tslice_from_field(P, PSTRING, C, CVAL, SVAL, T, \
           TYPESTRING, FIXTYPE, VARTYPE, MYREAL) \
int write_##P##N_##T##_timeslice_from_field(QIO_Writer *outfile, \
                                            QIO_String *xml_out, \
                                            VARTYPE *src, int count, \
                                            int t0) { \
 \
  char qdptype[] = TYPESTRING, prec[] = PSTRING; \
  int status, datum_size = sizeof(FIXTYPE), word_size = sizeof(MYREAL); \
  int lower[4] = {0, 0, 0, t0}, upper[4] = {nx - 1, ny - 1, nz - 1, t0}; \
  QIO_RecordInfo *rec_info; \
 \
  /* Create the record info for the field */ \
  rec_info = QIO_create_record_info(QIO_HYPER, lower, upper, 4, qdptype, \
                                    prec, CVAL, SVAL, datum_size, count); \
 \
  /* Write the record for the field */ \
  status = QIO_write(outfile, rec_info, xml_out, \
                     vget_##P##N_##T##_from_field, \
                     count * datum_size, word_size, (void *)src); \
  if (status != QIO_SUCCESS) \
    return 1; \
 \
  QIO_destroy_record_info(rec_info); \
 \
  return 0; \
}

#define make_write_all(P, PSTRING, C, CVAL, SVAL, T, TYPESTRING, \
                       FIXTYPE, VARTYPE, MYREAL) \
  make_write_all_from_site(P, PSTRING, C, CVAL, SVAL, T, TYPESTRING, \
                           FIXTYPE, VARTYPE, MYREAL); \
  make_write_all_from_field(P, PSTRING, C, CVAL, SVAL, T, TYPESTRING, \
                            FIXTYPE, VARTYPE, MYREAL);

#define make_write_tslice(P, PSTRING, C, CVAL, SVAL, T, TYPESTRING, \
                          FIXTYPE, VARTYPE, MYREAL) \
  make_write_tslice_from_site(P, PSTRING, C, CVAL, SVAL, T, TYPESTRING, \
                              FIXTYPE, VARTYPE, MYREAL); \
  make_write_tslice_from_field(P, PSTRING, C, CVAL, SVAL, T, TYPESTRING, \
                               FIXTYPE, VARTYPE, MYREAL);

// Single precision
make_write_all(F, "F",     , 0,    0, R, "QLA_F_Real",           float,     Real,     float);
make_write_all(F, "F",     , 0,    0, C, "QLA_F_Complex",        fcomplex,  complex,  float);
make_write_all(F, "F", NCOL, NCOL, 0, V, "USQCD_FN_ColorVector", fvector_f, vector_f, float);
make_write_all(F, "F", NCOL, NCOL, 0, M, "USQCD_FN_ColorMatrix", fmatrix_f, matrix_f, float);

make_write_tslice(F, "F",     , 0,    0, R, "QLA_F_Real",           float,     Real,     float);
make_write_tslice(F, "F",     , 0,    0, C, "QLA_F_Complex",        fcomplex,  complex,  float);
make_write_tslice(F, "F", NCOL, NCOL, 0, V, "USQCD_FN_ColorVector", fvector_f, vector_f, float);

// Double precision
make_write_all(D, "D",     , 0,    0, C, "QLA_D_Complex",        dcomplex,  complex,  double);
make_write_all(D, "D", NCOL, NCOL, 0, V, "USQCD_DN_ColorVector", dvector_f, vector_f, double);
make_write_tslice(D, "D",     , 0,    0, C, "QLA_D_Complex",        dcomplex,  complex,  double);
make_write_tslice(D, "D", NCOL, NCOL, 0, V, "USQCD_DN_ColorVector", dvector_f, vector_f, double);

// Read site structure data
#define make_read_to_site(P, C, T, FIXTYPE, VARTYPE, MYREAL) \
int read_##P##N_##T##_to_site(QIO_Reader *infile, QIO_String *xml_in, \
                              field_offset dest, int count) { \
 \
  int status, datum_size = sizeof(FIXTYPE), word_size = sizeof(MYREAL); \
  QIO_RecordInfo rec_info; \
 \
  /* Read the field record */ \
  status = QIO_read(infile, &rec_info, xml_in, \
                    vput_##P##N_##T##_to_site, datum_size * count, \
                    word_size, (void *)&dest); \
  node0_printf("Record info \n\"%s\"\n", QIO_string_ptr(xml_in)); \
  if (status != QIO_SUCCESS) \
    return 1; \
 \
  node0_printf("Checksums %x %x\n", QIO_get_reader_last_checksuma(infile), \
                                    QIO_get_reader_last_checksumb(infile)); \
 \
  return 0; \
}

// Read field data
#define make_read_to_field(P, C, T, FIXTYPE, VARTYPE, MYREAL) \
int read_##P##N_##T##_to_field(QIO_Reader *infile, QIO_String *xml_in, \
                               VARTYPE *dest, int count) { \
 \
  int status, datum_size = sizeof(FIXTYPE), word_size = sizeof(MYREAL); \
  QIO_RecordInfo rec_info; \
 \
  /* Read the field record */ \
  status = QIO_read(infile, &rec_info, xml_in, \
                    vput_##P##N_##T##_to_field, datum_size * count, \
                    word_size, (void *)dest); \
  node0_printf("Record info \n\"%s\"\n", QIO_string_ptr(xml_in)); \
  if (status != QIO_SUCCESS) \
    return 1; \
 \
  node0_printf("Checksums %x %x\n", \
               QIO_get_reader_last_checksuma(infile), \
               QIO_get_reader_last_checksumb(infile)); \
 \
  return 0; \
}

#define make_read(P, C, T, FIXTYPE, VARTYPE, MYREAL) \
  make_read_to_site(P, C, T, FIXTYPE, VARTYPE, MYREAL); \
  make_read_to_field(P, C, T, FIXTYPE, VARTYPE, MYREAL);

// Single precision
make_read(F,     , R, float,     Real,     float);
make_read(F,     , C, fcomplex,  complex,  float);
make_read(F, NCOL, V, fvector_f, vector_f, float);
make_read(F, NCOL, M, fmatrix_f, matrix_f, float);

// Double precision
make_read(D,     , C, dcomplex,  complex,  double);
make_read(D, NCOL, V, dvector_f, vector_f, double);
make_read(D, NCOL, M, dmatrix_f, matrix_f, double);

// Factory function for moving RNG state from site structure to output
// arg contains pointer to field offset value
void vget_S_from_site(char *buf, size_t index, int count, void *arg) {
  char *dest = buf;
  field_offset src = *((field_offset *)arg);
  site *s = &lattice[index];
  char *src_prn = (char *)F_PT(s,src);

  memcpy(dest, src_prn, sizeof(double_prn));
}

int write_S_from_site(QIO_Writer *outfile, QIO_String *xml_out,
                      field_offset src) {

  char qdptype[] = "RandomState", prec[] = "";
  int datum_size = sizeof(double_prn), word_size = sizeof(float);
  int count = 1, status;
  QIO_RecordInfo *rec_info;

  // Create the record info for the field
  rec_info = QIO_create_record_info(QIO_FIELD, 0, 0, 0, qdptype, prec, 0,
                                    0, datum_size, count);
  // Write the record for the field
  status = QIO_write(outfile, rec_info, xml_out, vget_S_from_site,
                     count * datum_size, word_size, (void *)&src);
  if (status != QIO_SUCCESS)
    return 1;

  QIO_destroy_record_info(rec_info);
  return 0;
}

// Factory function for moving RNG state from input to site structure
// Assume input lattice is single precision, NCOL colors
void vput_S_to_site(char *buf, size_t index, int count, void *arg) {
  char *src = buf, *dest_prn;
  field_offset dest = *((field_offset *)arg);
  site *s = &lattice[index];

  dest_prn = (char *)F_PT(s, dest);
  memcpy(dest_prn, src, sizeof(double_prn));
}

// Read random number state
int read_S_to_site(QIO_Reader *infile, QIO_String *xml_in,
                   field_offset dest) {

  int status, count = 1;
  int datum_size = sizeof(double_prn), word_size = sizeof(float);
  QIO_RecordInfo rec_info;

  // Read the field record
  status = QIO_read(infile, &rec_info, xml_in,
                    vput_S_to_site, datum_size * count,
                    word_size, (void *)&dest);
  node0_printf("Record info \n\"%s\"\n", QIO_string_ptr(xml_in));
  if (status != QIO_SUCCESS)
    return 1;

  node0_printf("Checksums %x %x\n",
               QIO_get_reader_last_checksuma(infile),
               QIO_get_reader_last_checksumb(infile));

  return 0;
}
// -----------------------------------------------------------------
