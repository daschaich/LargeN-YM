/********************** io_helpers.c **********************************/
/* MIMD version 7 */
/* DT 8/97
     General purpose high level routines, to be used by any application
     that wants them.
*/

#include "generic_includes.h"
#include "../include/io_lat.h"
#include "../include/file_types.h"
#ifdef HAVE_QIO
#include <qio.h>
#endif

/**static file_type gauge_list[N_GAUGE_TYPES] =
  { {FILE_TYPE_GAUGE_V1,      GAUGE_VERSION_NUMBER_V1},
    {FILE_TYPE_GAUGE_V5,      GAUGE_VERSION_NUMBER},
    {FILE_TYPE_GAUGE_1996,    GAUGE_VERSION_NUMBER_1996},
    {FILE_TYPE_GAUGE_FNAL,    GAUGE_VERSION_NUMBER_FNAL},
    {FILE_TYPE_GAUGE_ARCHIVE, GAUGE_VERSION_NUMBER_ARCHIVE},
    {FILE_TYPE_GAUGE_SCIDAC,  LIME_MAGIC_NO}
    };**/

/* save a lattice in any of the formats:
    SAVE_ASCII, SAVE_SERIAL, SAVE_PARALLEL, SAVE_CHECKPOINT
*/
gauge_file *save_lattice( int flag, char *filename, char *stringLFN){
    double dtime;
    gauge_file *gf = NULL;

#ifndef NOLINKS
    d_plaquette(&g_ssplaq,&g_stplaq);
    d_linktrsum(&linktrsum);
    nersc_checksum = nersc_cksum();
#endif

    dtime = -dclock();
    switch( flag ){
	case SAVE_SERIAL:
	    gf = save_serial(filename);
	    break;
        case FORGET:
            gf = NULL;
            break;
/*	other cases not implemented yet -bqs 1/07
	case SAVE_ASCII:
	    gf = save_ascii(filename);
	    break;
	case SAVE_PARALLEL:
	    gf = save_parallel(filename);
	    break;
	case SAVE_CHECKPOINT:
	    gf = save_checkpoint(filename);
	    break;
        case SAVE_SERIAL_FM:
 	    printf("Save serial FNAL format not implemented\n");
	    break;
        case SAVE_SERIAL_ILDG:
#ifdef HAVE_QIO
 	    gf = save_serial_ildg(filename,stringLFN);
#else
	    node0_printf("save_serial_ildg requires QIO compilation\n");
	    terminate(1);
#endif
            break;
        case SAVE_PARALLEL_ILDG:
#ifdef HAVE_QIO
 	    gf = save_parallel_ildg(filename,stringLFN);
#else
	    node0_printf("save_parallel_ildg requires QIO compilation\n");
	    terminate(1);
#endif
            break;
        case SAVE_PARTITION_ILDG:
#ifdef HAVE_QIO
 	    gf = save_partition_ildg(filename,stringLFN);
#else
	    node0_printf("save_partition_ildg requires QIO compilation\n");
	    terminate(1);
#endif
            break;
        case SAVE_MULTIFILE_ILDG:
#ifdef HAVE_QIO
 	    gf = save_multifile_ildg(filename,stringLFN);
#else
	    node0_printf("save_multifile_ildg requires QIO compilation\n");
	    terminate(1);
#endif
            break;
        case SAVE_SERIAL_SCIDAC:
#ifdef HAVE_QIO
 	    gf = save_serial_scidac(filename);
#else
	    node0_printf("save_serial_scidac requires QIO compilation\n");
	    terminate(1);
#endif
            break;
        case SAVE_PARALLEL_SCIDAC:
#ifdef HAVE_QIO
 	    gf = save_parallel_scidac(filename);
#else
	    node0_printf("save_parallel_scidac requires QIO compilation\n");
	    terminate(1);
#endif
            break;
        case SAVE_MULTIFILE_SCIDAC:
#ifdef HAVE_QIO
	    gf = save_multifile_scidac(filename);
#else
	    node0_printf("save_multifile_scidac requires QIO compilation\n");
	    terminate(1);
#endif
            break;
        case SAVE_PARTITION_SCIDAC:
#ifdef HAVE_QIO
 	    gf = save_partition_scidac(filename);
#else
	    node0_printf("save_partition_scidac requires QIO compilation\n");
	    terminate(1);
#endif
            break;
	case SAVE_SERIAL_ARCHIVE:
	    gf = save_serial_archive(filename);
	    break;
*/
	default:
	    node0_printf("\nsave_lattice: ERROR: unknown type for saving lattice\n");
	    terminate(1);
    }
    dtime += dclock();
    if(flag != FORGET)
      node0_printf("Time to save = %e\n",dtime);
#if (PRECISION==1)
    node0_printf("CHECK PLAQ: %e %e\n",g_ssplaq,g_stplaq);
    node0_printf("CHECK NERSC LINKTR: %e CKSUM: %x\n",
		 linktrsum.real/(Real)NCOL,nersc_checksum);
#else
    /* Double precision */
    node0_printf("CHECK PLAQ: %.16e %.16e\n",g_ssplaq,g_stplaq);
    node0_printf("CHECK NERSC LINKTR: %.16e CKSUM: %x\n",
		 linktrsum.real/(Real)NCOL,nersc_checksum);
#endif
    return gf;
}

/* reload a lattice in any of the formats, or cold lattice, or keep
     current lattice:
    FRESH, CONTINUE,
    RELOAD_ASCII, RELOAD_SERIAL, RELOAD_PARALLEL
*/
void coldlat();

gauge_file *reload_lattice( int flag, char *filename){
    double dtime;
    gauge_file *gf = NULL;
    Real max_deviation;
#if PRECISION == 2
    Real max_deviation2;
#endif

    dtime = -dclock();
    switch(flag){
	case CONTINUE:	/* do nothing */
            gf = NULL;
	    break;
	case FRESH:	/* cold lattice */
	    coldlat();
            gf = NULL;
	    break;
	case RELOAD_SERIAL:	/* read binary lattice serially */
	    gf = restore_serial(filename);
	    break;
/*	other cases not implemented yet -bqs 1/07 */
	case RELOAD_ASCII:	/* read Ascii lattice */
/*	    gf = restore_ascii(filename);
	    break;*/
	case RELOAD_PARALLEL:	/* read binary lattice in parallel */
/*	    gf = restore_parallel(filename);
	    break;*/
	default:
	    if(this_node==0)printf("reload_lattice: Bad startflag %d\n",flag);
	    terminate(1);
    }
    dtime += dclock();
    if(flag != FRESH && flag != CONTINUE)
      node0_printf("Time to reload gauge configuration = %e\n",dtime);
/*
SCHROED_FUN call to  set_boundary_fields() moved to setup.c
*/
#ifndef NOLINKS
    d_plaquette(&g_ssplaq,&g_stplaq);
    d_linktrsum(&linktrsum);
    nersc_checksum = nersc_cksum();
#endif
    if(this_node==0){
#if (PRECISION==1)
    node0_printf("CHECK PLAQ: %e %e\n",g_ssplaq,g_stplaq);
    node0_printf("CHECK NERSC LINKTR: %e CKSUM: %x\n",
		 linktrsum.real/(Real)NCOL,nersc_checksum);
    fflush(stdout);
#else
    /* Double precision */
    node0_printf("CHECK PLAQ: %.16e %.16e\n",g_ssplaq,g_stplaq);
    node0_printf("CHECK NERSC LINKTR: %.16e CKSUM: %x\n",
		 linktrsum.real/(Real)NCOL,nersc_checksum);
    fflush(stdout);
#endif
    }
    dtime = -dclock();
    max_deviation = check_unitarity();
    g_floatmax(&max_deviation);
#if (PRECISION==1)
    if(this_node==0)printf("Unitarity checked.  Max deviation %.2e\n",
			   max_deviation); fflush(stdout);
#else
    reunitarize();
    max_deviation2 = check_unitarity();
    g_floatmax(&max_deviation2);
    if(this_node==0)
      printf("Reunitarized for double precision. Max deviation %.2e changed to %.2e\n",
                       max_deviation,max_deviation2); fflush(stdout);
#endif
    dtime += dclock();
    if(this_node==0)printf("Time to check unitarity = %e\n",dtime);
    return gf;
}

/* find out what kind of starting lattice to use, and lattice name if
   necessary.  This routine is only called by node 0.
*/
int ask_starting_lattice( FILE *fp, int prompt, int *flag, char *filename ){
    char savebuf[256];
    int status;

    if (prompt!=0) printf(
        "enter 'continue', 'fresh', 'reload_ascii', 'reload_serial', or 'reload_parallel'\n");
    status=fscanf(fp,"%s",savebuf);
    if (status == EOF){
      printf("ask_starting_lattice: EOF on STDIN.\n");
      return 1;
    }
    if(status !=1) {
        printf("\nask_starting_lattice: ERROR IN INPUT: can't read starting lattice command\n");
        return 1;
    }

    printf("%s ",savebuf);
    if(strcmp("fresh",savebuf) == 0 ){
       *flag = FRESH;
    printf("\n");
    }
    else if(strcmp("continue",savebuf) == 0 ) {
        *flag = CONTINUE;
	printf("\n");
    }
    else if(strcmp("reload_ascii",savebuf) == 0 ) {
       *flag = RELOAD_ASCII;
    }
    else if(strcmp("reload_serial",savebuf) == 0 ) {
       *flag = RELOAD_SERIAL;
    }
    else if(strcmp("reload_parallel",savebuf) == 0 ) {
       *flag = RELOAD_PARALLEL;
    }
    else{
    	printf(" is not a valid starting lattice command. INPUT ERROR.\n");
	return 1;
    }

    /*read name of file and load it */
    if( *flag != FRESH && *flag != CONTINUE ){
        if(prompt!=0)printf("enter name of file containing lattice\n");
        status=fscanf(fp,"%s",filename);
        if(status !=1) {
	    printf("\nask_starting_lattice: ERROR IN INPUT: error reading file name\n"); return 1;
        }
	printf("%s\n",filename);
    }
    return 0;
}

/* find out what do to with lattice at end, and lattice name if
   necessary.  This routine is only called by node 0.
*/
int ask_ending_lattice(FILE *fp, int prompt, int *flag, char *filename ){
    char savebuf[256];
    int status;

    if (prompt!=0) printf(
        "'forget' lattice at end,  'save_ascii', 'save_serial', 'save_parallel', 'save_checkpoint', 'save_serial_fm', 'save_serial_scidac', 'save_parallel_scidac', 'save_multifile_scidac', 'save_partition_scidac', 'save_serial_archive', 'save_serial_ildg', 'save_parallel_ildg', 'save_partition_ildg', or 'save_multifile_ildg'\n");
    status=fscanf(fp,"%s",savebuf);
    if(status !=1) {
        printf("\nask_ending_lattice: ERROR IN INPUT: error reading ending lattice command\n");
        return 1;
    }
    printf("%s ",savebuf);
    if(strcmp("save_ascii",savebuf) == 0 )  {
        *flag=SAVE_ASCII;
    }
    else if(strcmp("save_serial",savebuf) == 0 ) {
        *flag=SAVE_SERIAL;
    }
    else if(strcmp("save_parallel",savebuf) == 0 ) {
      *flag=SAVE_PARALLEL;
    }
    else if(strcmp("save_checkpoint",savebuf) == 0 ) {
        *flag=SAVE_CHECKPOINT;
    }
    else if(strcmp("save_serial_fm",savebuf) == 0 ) {
        *flag=SAVE_SERIAL_FM;
    }
    else if(strcmp("save_serial_scidac",savebuf) == 0 ) {
#ifdef HAVE_QIO
        *flag=SAVE_SERIAL_SCIDAC;
#else
	node0_printf("requires QIO compilation!\n");
	terminate(1);
#endif
    }
    else if(strcmp("save_parallel_scidac",savebuf) == 0 ) {
#ifdef HAVE_QIO
        *flag=SAVE_PARALLEL_SCIDAC;
#else
	node0_printf("requires QIO compilation!\n");
	terminate(1);
#endif
    }
    else if(strcmp("save_multifile_scidac",savebuf) == 0 ) {
#ifdef HAVE_QIO
        *flag=SAVE_MULTIFILE_SCIDAC;
#else
	node0_printf("requires QIO compilation!\n");
	terminate(1);
#endif
    }
    else if(strcmp("save_partition_scidac",savebuf) == 0 ) {
#ifdef HAVE_QIO
        *flag=SAVE_PARTITION_SCIDAC;
#else
	node0_printf("requires QIO compilation!\n");
	terminate(1);
#endif
    }
    else if(strcmp("save_serial_archive",savebuf) == 0 ) {
        *flag=SAVE_SERIAL_ARCHIVE;
    }
    else if(strcmp("save_serial_ildg",savebuf) == 0 ) {
#ifdef HAVE_QIO
        *flag=SAVE_SERIAL_ILDG;
#else
	node0_printf("requires QIO compilation!\n");
	terminate(1);
#endif
    }
    else if(strcmp("save_parallel_ildg",savebuf) == 0 ) {
#ifdef HAVE_QIO
        *flag=SAVE_PARALLEL_ILDG;
#else
	node0_printf("requires QIO compilation!\n");
	terminate(1);
#endif
    }
    else if(strcmp("save_partition_ildg",savebuf) == 0 ) {
#ifdef HAVE_QIO
        *flag=SAVE_PARTITION_ILDG;
#else
	node0_printf("requires QIO compilation!\n");
	terminate(1);
#endif
    }
    else if(strcmp("save_multifile_ildg",savebuf) == 0 ) {
#ifdef HAVE_QIO
        *flag=SAVE_MULTIFILE_ILDG;
#else
	node0_printf("requires QIO compilation!\n");
	terminate(1);
#endif
    }
    else if(strcmp("forget",savebuf) == 0 ) {
        *flag=FORGET;
	printf("\n");
    }
    else {
      printf("is not a save lattice command. INPUT ERROR\n");
      return 1;
    }

    if( *flag != FORGET ){
        if(prompt!=0)printf("enter filename\n");
        status=fscanf(fp,"%s",filename);
        if(status !=1){
    	    printf("\nask_ending_lattice: ERROR IN INPUT: error reading filename\n"); return 1;
        }
	printf("%s\n",filename);

    }
    return 0;
}

int ask_ildg_LFN(FILE *fp, int prompt, int flag, char *stringLFN){
  int status = 0;

  /* For ILDG output formats we require a logical file name next */
  if( flag == SAVE_SERIAL_ILDG ||
      flag == SAVE_PARALLEL_ILDG ||
      flag == SAVE_PARTITION_ILDG ||
      flag == SAVE_MULTIFILE_ILDG ){
    status = get_s(fp, prompt, "ILDG_LFN", stringLFN);
  }
  else
    stringLFN[0] = '\0';
  return status;
}

void coldlat(){
    /* sets link matrices to unit matrices       */
    /* linkf!  -bqs                              */
    /* SF: Set to initial bkgd field config -bqs */
    register int i,j,k,dir;
    register site *sit;
#ifdef SF
    double lr,li,angle,angle0[NCOL],dangle[NCOL];
#endif

    for(dir=XUP;dir<=ZUP;dir++){
#ifdef SF
        /* calculate electric field from SF bndry values */
	for(j=0;j<NCOL;j++){
	    lr=linkf_bndr_dn[dir].e[j][j].real;
	    li=linkf_bndr_dn[dir].e[j][j].imag;
	    angle0[j]=atan2(li,lr);
	    lr=linkf_bndr_up[dir].e[j][j].real;
	    li=linkf_bndr_up[dir].e[j][j].imag;
	    angle=atan2(li,lr);
	    dangle[j]=(angle-angle0[j])/nt;
	    /* printf("Angles %f %f\n",angle0[j],dangle[j]); */
	}
#endif
	FORALLDYNLINKS(i,sit,dir){
	    for(j=0; j<NCOL; j++)  {
		for(k=0; k<NCOL; k++)  {
		    if (j != k)  {
			sit->linkf[dir].e[j][k] = cmplx(0.0,0.0);
		    }
		    else  {
#ifdef SF
			angle=angle0[j]+sit->t*dangle[j];
			sit->linkf[dir].e[j][k].real=cos(angle);
			sit->linkf[dir].e[j][k].imag=sin(angle);
#else
			sit->linkf[dir].e[j][k] = cmplx(1.0,0.0);
#endif
		    }
		}
	    }
	} /* end loop for dir<=TUP */
    }     /* end loop over dir */

    FORALLSITES(i,sit){
	for(j=0; j<NCOL; j++)  {
	    for(k=0; k<NCOL; k++)  {
		if (j != k)  {
		    sit->linkf[TUP].e[j][k] = cmplx(0.0,0.0);
		}
		else  {
		    sit->linkf[TUP].e[j][k] = cmplx(1.0,0.0);
		}
	    }
	}
    } /* end loop for dir=TUP */

#ifndef SF
    node0_printf("unit gauge configuration loaded\n");
#else
    node0_printf("cold SF gauge configuration loaded\n");
#endif
}

void funnylat()  {
    /* sets link matrices to funny matrices for debugging */
    /* linkf!  -bqs	*/
    register int i,j,k,dir;
    register site *sit;

    FORALLSITES(i,sit){
	for(dir=XUP;dir<=TUP;dir++){
	    for(j=0; j<NCOL; ++j)  {
		for(k=0; k<NCOL; ++k)  {
		    sit->linkf[dir].e[j][k] = cmplx(0.0,0.0);
		}
	    }
	    sit->linkf[dir].e[0][0].real = dir;
	    sit->linkf[dir].e[1][1].real = 10*sit->x;
	    sit->linkf[dir].e[2][2].real = 100*sit->y;
	    sit->linkf[dir].e[0][0].imag = dir;
	    sit->linkf[dir].e[1][1].imag = 10*sit->z;
	    sit->linkf[dir].e[2][2].imag = 100*sit->t;
	}
    }
}

/* Read and echo the next tag.  Echo any intervening comments */
/* Comments begin with # and apply to the rest of the line */
/* Verify that the input tag agrees with the expected tag */

static int get_tag(FILE *fp, char *tag, char *myname){
  static char checktag[80];
  char line[512];
  int s;

  while(1){
    s = fscanf(fp,"%s",checktag);
    if (s == EOF){
      printf("%s(%d): EOF on input.\n",myname,this_node);
      return 1;
    }
    if(s == 0){
      printf("%s(%d) Error reading %s\n",myname,this_node,tag);
      return 1;
    }
    if(strchr(checktag,'#')!=NULL){
      printf("%s",checktag);
      if(fgets(line,512,fp)==NULL){
	printf("%s(%d) EOF on input.\n",myname,this_node);
	return 1;
      }
      printf("%s",line);
    }
    else{
      if(strcmp(checktag,tag) != 0){
	printf("\n%s: ERROR IN INPUT: expected %s but found %s\n",
	       myname,tag,checktag);
	return 1;
      }
      printf("%s ",tag);
      return 0;
    }
  }
}

/* Check return value of scanf */
static int check_read(int s, char *myname, char *tag){

  if (s == EOF){
    printf("\n%s: Expecting value for %s but found EOF.\n",
	   myname,tag);
    return 1;
  }
  else if(s == 0){
    printf("\n%s: Format error reading value for %s\n",
	   myname,tag);
    return 1;
  }
  else
    return 0;
}

/* get_f is used to get a floating point number.  If prompt is non-zero,
it will prompt for the input value with the variable_name_string.  If
prompt is zero, it will require that variable_name_string precede the
input value.  get_i gets an integer.
get_i and get_f return the values, and exit on error */

int get_f( FILE *fp, int prompt, char *tag, Real *value ){
    int s;
    char checkvalue[80];
    char myname[] = "get_f";

    if(prompt)  {
	s = 0;
	while(s != 1){
	  printf("enter %s ",tag);
	  fscanf(fp,"%s",checkvalue);
#if PRECISION == 1
	  s=sscanf(checkvalue,"%e",value);
#else
	  s=sscanf(checkvalue,"%le",value);
#endif
	  if(s==EOF)return 1;
	  if(s==0)printf("Data format error.\n");
	  else printf("%s %g\n",tag,*value);
	}
    }
    else  {
      if(get_tag(fp, tag, myname) == 1)return 1;

#if PRECISION == 1
      s = fscanf(fp,"%e",value);
#else
      s = fscanf(fp,"%le",value);
#endif
      if(check_read(s,myname,tag) == 1)return 1;
      printf("%g\n",*value);
    }

    return 0;
}

int get_i( FILE *fp, int prompt, char *tag, int *value ){
    int s;
    char checkvalue[80];
    char myname[] = "get_i";

    if(prompt)  {
      s = 0;
      while(s != 1){
    	printf("enter %s ",tag);
	fscanf(fp,"%s",checkvalue);
    	s=sscanf(checkvalue,"%d",value);
	if(s==EOF)return 1;
	if(s==0)printf("Data format error.\n");
	else printf("%s %d\n",tag,*value);
      }
    }
    else  {
      if(get_tag(fp, tag, myname) == 1)return 1;

      s = fscanf(fp,"%d",value);
      if(check_read(s,myname,tag) == 1)return 1;
      printf("%d\n",*value);
    }

    return 0;

}

#ifdef SF
/*  modifed from get_s(..), lists possible choices according to
requested input flag */
int ask_sf_flag( FILE *fp, int prompt, char *tag, int *value ){
    int s;
    char value_code[80];
    char myname[] = "ask_sf_flag";

    if(prompt)  {
      s = 0;
      while(s != 1){
        if(strcmp("bc_flag",tag) == 0 )  {
          printf("enter option for boundary values: 'trivial', or SU(2): 'alpha', or SU(3): 'alpha_one', or 'alpha_two', or SU(4): 'sym' \n");
        }
        else if(strcmp("ct_flag",tag) == 0 )  {
          printf("enter option for c_t: 'tree_level', or 'one_loop'\n");
        }
        else {
          printf("Unexpected flag!\n");
          terminate(1);
        }
    	s=fscanf(fp,"%s",value_code);
	if(s==EOF)return 1;
	else if(s==0)printf("Data format error.\n");
        else printf("%s ",tag);
      }
    }
    else  {
      if(get_tag(fp, tag, myname) == 1)return 1;

      s=fscanf(fp,"%s",value_code);
      if(check_read(s,myname,tag) == 1)return 1;

    }

    /* assign numerical codes for options */

    if(strcmp("bc_flag",tag) == 0 ) {
       if(strcmp("trivial",value_code) == 0 ) {
          *value=TRIVIAL;
/*        printf("= %d (%s)\n",*value,value_code); */
          printf(" %s\n", value_code);
       }
       else if(strcmp("alpha",value_code) == 0 ) {
          *value=ALPHA;
          printf(" %s\n", value_code);
       }
       else if(strcmp("alpha_one",value_code) == 0 ) {
          *value=ALPHA_ONE;
          printf(" %s\n", value_code);
       }
       else if(strcmp("alpha_two",value_code) == 0 ) {
          *value=ALPHA_TWO;
          printf(" %s\n", value_code);
       }
       else if(strcmp("sf_sym",value_code) == 0 ) {
          *value=SF_SYM;
          printf(" %s\n", value_code);
       }
       else {
          printf("UNKOWN OPTION: %s\n",value_code);
          terminate(1);
       }
    }
    else if(strcmp("ct_flag",tag) == 0 ) {
       if(strcmp("tree_level",value_code) == 0 ) {
          *value=TREE_LEVEL;
          printf(" %s\n", value_code);
       }
       else if(strcmp("one_loop",value_code) == 0 ) {
          *value=ONE_LOOP;
          printf(" %s\n", value_code);
       }
       else if(strcmp("two_loop",value_code) == 0 ) {
          *value=TWO_LOOP;
          printf(" %s\n", value_code);
       }
       else {
          printf("UNKOWN OPTION: %s\n",value_code);
          terminate(1);
       }
    }
    else {
       printf("UNKOWN FLAG: %s\n",tag);
       terminate(1);
    }

    return 0;

}
#endif /* SF */

/***  NOT IN USE!!!  ***  NOT IN USE!!! ***  NOT IN USE!!!  ***  NOT IN USE!!!
modifed from get_i(..), lists possible choices according to
requested input flag, with numerical codes for values

int ask_sf_flag( FILE *fp, int prompt, char *tag, int *value ){
    int s;
    char checkvalue[80];
    char myname[] = "ask_sf_flag";

    if(prompt)  {
      s = 0;
      while(s != 1){
        if(strcmp("bc_flag",tag) == 0 )  {
          printf("enter option for boundary values:\n '0' boundary value = identity\n '1' boundary value = alpha collab.\n");
        }
        else if(strcmp("ct_flag",tag) == 0 )  {
          printf("enter option for perturbative c_t:\n '0' c_t = tree-level\n '1' c_t = one-loop\n");
        }
        else {
          printf("Unexpected flag!\n");
          terminate(1);
        }
	fscanf(fp,"%s",checkvalue);
    	s=sscanf(checkvalue,"%d",value);
	if(s==EOF)return 1;
	if(s==0)printf("Data format error.\n");
	else printf("%s %d\n",tag,*value);
      }
    }
    else  {
      if(get_tag(fp, tag, myname) == 1)return 1;

      s = fscanf(fp,"%d",value);
      if(check_read(s,myname,tag) == 1)return 1;
      printf("%d\n",*value);
    }

    return 0;

}
END OF ALTERNATIVE-FORMAT SF INPUT */

/* Read a single word as a string */

int get_s( FILE *fp, int prompt, char *tag, char *value ){
    int s;
    char myname[] = "get_s";

    if(prompt)  {
      s = 0;
      while(s != 1){
    	printf("enter %s ",tag);
    	s=fscanf(fp,"%s",value);
	if(s==EOF)return 1;
	if(s==0)printf("Data format error.\n");
	else printf("%s %s\n",tag,value);
      }
    }
    else  {
      if(get_tag(fp, tag, myname) == 1)return 1;

      s = fscanf(fp,"%s",value);
      if(check_read(s,myname,tag) == 1)return 1;
      printf("%s\n",value);
    }
    return 0;
}

/* Read a vector of integers */
int get_vi( FILE* fp, int prompt, char *tag,
	    int *value, int nvalues ){
    int s,i;
    char myname[] = "get_vi";

    if(prompt)  {
      s = 0;
      printf("enter %s with %d values",tag, nvalues);
      for(i = 0; i < nvalues; i++){
	while(s != 1){
	  printf("\n[%d] ",i);
	  s=fscanf(fp,"%d",value+i);
	  if(s==EOF)return 1;
	  if(s==0)printf("Data format error.\n");
	  printf("%s %d\n",tag,value[i]);
	}
      }
    }
    else  {
      if(get_tag(fp, tag, myname) == 1)return 1;

      for(i = 0; i < nvalues; i++){
	s = fscanf(fp,"%d",value + i);
	if(check_read(s,myname,tag) == 1)return 1;
	printf("%d ",value[i]);
      }
      printf("\n");
    }

    return 0;

}

/* Read a vector of reals */
int get_vf( FILE* fp, int prompt, char *tag,
	    Real *value, int nvalues ){
    int s,i;
    char myname[] = "get_vf";

    if(prompt)  {
      s = 0;
      printf("enter %s with %d values",tag, nvalues);
      for(i = 0; i < nvalues; i++){
	while(s != 1){
	  printf("\n[%d] ",i);
#if PRECISION == 1
	  s=scanf("%e",value+i);
#else
	  s=scanf("%le",value+i);
#endif
	  if(s==EOF)return 1;
	  if(s==0)printf("Data format error.\n");
	  printf("%s %g\n",tag,*(value+i));
	}
      }
    }
    else  {
      if(get_tag(fp, tag, myname) == 1)return 1;

      for(i = 0; i < nvalues; i++){
#if PRECISION == 1
	s = fscanf(fp,"%e",value + i);
#else
	s = fscanf(fp,"%le",value + i);
#endif
	if(check_read(s,myname,tag) == 1)return 1;
	printf("%g ",value[i]);
      }
      printf("\n");
    }

    return 0;

}

/* get_prompt gets the initial value of prompt */
/* 0 for reading from file, 1 prompts for input from terminal */
/* should be called only by node 0 */
/* return 0 if sucessful, 1 if failure */
int get_prompt( FILE *fp, int *prompt ){
    char initial_prompt[512];
    int status;
    char myname[] = "get_prompt";

    *prompt = -1;
    printf( "type 0 for no prompts  or 1 for prompts\n");
    while(1){
      status = fscanf(fp, "%s",initial_prompt);
      if(status != 1){
	printf("\n%s: Can't read input\n",myname);
	terminate(1);
      }
      if(strchr(initial_prompt,'#')==NULL)break;
      /* Provide for comment lines with # before "prompt" */
      else{
	printf("%s",initial_prompt);
	if(fgets(initial_prompt,512,fp)==NULL){
	  printf("%s(%d) EOF on input.\n",myname,this_node);
	  return 1;
	}
	printf("%s",initial_prompt);
      }
    }
    if(strcmp(initial_prompt,"prompt") == 0)  {
      fscanf(fp, "%d",prompt);
    }
    else if(strcmp(initial_prompt,"0") == 0) *prompt=0;
    else if(strcmp(initial_prompt,"1") == 0) *prompt=1;

    if( *prompt==0 || *prompt==1 )return 0;
    else{
        printf("\n%s: ERROR IN INPUT: initial prompt\n",myname);
        return 1;
    }
}
