#include "generic_wilson_includes.h"
#include "../include/file_types.h"
#include "../include/io_lat.h"
#include "../include/io_wprop.h"
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#undef MAX_BUF_LENGTH
#define MAX_BUF_LENGTH 4096

#ifndef PAD_SEND_BUF
#define PAD_SEND_BUF 8
#endif
#define NATURAL_ORDER 0

int read_w_fm_prop_hdr(w_prop_file *wpf)
{
  w_prop_header *wph = wpf->header;
  int *dims = wph->dims;

  int32type   tmp, elements_per_site; 
  int32type   size_of_element, order; 
  int32type   t_stamp;
  int byterevflag = 0;
  char myname[] = "read_w_fm_prop_hdr";

  if( sizeof(float) != sizeof(int32type)) {
    printf("%s: Can't byte reverse\n", wpf->filename);
    printf("requires size of int32type(%d) = size of float(%d)\n",
	   (int)sizeof(int32type),(int)sizeof(float));
    terminate(1);
  }

  if(sread_data(wpf->fp,&wph->magic_number,sizeof(int32type),
		myname,"magic_no")) terminate(1);

  tmp = wph->magic_number;
  if(wph->magic_number == IO_UNI_MAGIC) byterevflag = 0;
  else 
    {
      byterevn((int32type *)&wph->magic_number,1);
      if(wph->magic_number == IO_UNI_MAGIC) 
	{
	  byterevflag = 1; 
	  /** printf("Reading with byte reversal\n"); **/
	}
      else
	{
	  /* Restore magic number as originally read */
	  wph->magic_number = tmp;
	  
	  /* End of the road. */
	  printf("%s: Unrecognized magic number in prop file header.\n",
		 wpf->filename);
	  printf("Expected %x but read %x\n", IO_UNI_MAGIC,tmp);
	  terminate(1);
	}
    };

  if(sread_byteorder(byterevflag,wpf->fp,&t_stamp,
	     sizeof(t_stamp),myname,"t_stamp"))terminate(1);
  if(sread_byteorder(byterevflag,wpf->fp,&size_of_element,
	     sizeof(int32type),myname,"size_of_element"))terminate(1);
  if(sread_byteorder(byterevflag,wpf->fp,&elements_per_site,
	     sizeof(int32type),myname,"elements_per_site"))terminate(1);
  if(sread_byteorder(byterevflag,wpf->fp,dims,
	     4*sizeof(int32type),myname,"dimensions"))terminate(1);

  if( dims[0]!=nx || dims[1]!=ny || dims[2]!=nz || dims[3]!=nt )
    {
      /* So we can use this routine to discover the dimensions,
	 we provide that if nx = ny = nz = nt = -1 initially
	 we don't die */
      
      if(nx != -1 || ny != -1 || nz != -1 || nt != -1)
	{
	  printf(" Incorrect lattice size %d,%d,%d,%d\n",
		 dims[0], dims[1], dims[2], dims[3]);
	  terminate(1);
	}
    }
  else
    {
      nx = dims[0]; ny = dims[1]; nz = dims[2]; nt = dims[3];
      volume = nx*ny*nz*nt;
    }
  
  /* Consistency checks */
  
  if( size_of_element != sizeof(float) ||
      elements_per_site != 288 /* wilson_propagator */)
    {	
      printf(" File %s is not a wilson propagator",
	     wpf->filename);
      printf(" got size_of_element %d and elements_per_site %d\n",
	     size_of_element, elements_per_site);
      terminate(1);
    }
  
  /* The site order parameter is ignored */
  
  if(sread_byteorder(byterevflag,wpf->fp,&order,sizeof(int32type),
		     myname,"order parameter"))terminate(1);
  
  return byterevflag;
}
 
w_prop_file *r_serial_w_fm_i(char *filename)
{
  /* Returns file descriptor for opened file */

  w_prop_file *wpf;
  w_prop_header *wph;
  int byterevflag;

  wpf = setup_input_w_prop_file(filename);
  wph = wpf->header;

  if(this_node==0){
    wpf->fp = fopen(filename,"rb");
    if(wpf->fp==NULL){
      printf("Can't open propagator file %s, error %d\n",filename,errno);
      terminate(1);
    }
    byterevflag = read_w_fm_prop_hdr(wpf);

  }
  else wpf->fp = NULL;

  /* Broadcast the byterevflag from node 0 to all nodes */
  broadcast_bytes((char *)&byterevflag,sizeof(byterevflag));
  wpf->byterevflag = byterevflag;
  
  /* Node 0 broadcasts the header structure to all nodes */
  broadcast_bytes((char *)wph,sizeof(w_prop_header));

  return wpf;

} /* r_serial_w_fm_i */


void r_serial_w_fm(w_prop_file *wpf, field_offset dest_site, 
		   wilson_propagator *dest_field)
{
  int rcv_rank, rcv_coords, status;
  int destnode;
  int x,y,z,t,i=0, byterevflag, c0,s0,c1,s1,a;
  struct {
    fwilson_matrix q;
    char pad[PAD_SEND_BUF];    /* Introduced because some switches
				  perform better if message lengths
				  are longer */
  } msg;

  int buf_length=0, where_in_buf=0;
  fwilson_matrix *pbuff=NULL;
  wilson_propagator *qp;
  w_prop_check test_wpc;
  u_int32type *val;
  int rank29=0,rank31=0;
  int k;

  byterevflag = wpf->byterevflag;

  if(this_node == 0)
    {
      pbuff = (fwilson_matrix *)malloc(MAX_BUF_LENGTH*sizeof(fwilson_matrix));
      if(pbuff == NULL)
	{
	  printf("Node %d can't malloc pbuff\n",this_node);
	  fflush(stdout);
	  terminate(1);
	}
      
      buf_length = 0;
      where_in_buf = 0;

    } /* end of if(this_node == 0)*/

  g_sync();

  /* Node 0 reads and deals out the values */
  status = 0;
  for(rcv_rank=0; rcv_rank<volume; rcv_rank++)
    {
      /* We do only natural order here */
      rcv_coords = rcv_rank;
      
      x = rcv_coords % nx;   rcv_coords /= nx;
      y = rcv_coords % ny;   rcv_coords /= ny;
      z = rcv_coords % nz;   rcv_coords /= nz;
      t = rcv_coords % nt;

      /* All nodes initialize timeslice checksums at the beginning of
	 a time slice */
      if(x == 0 && y == 0 && z == 0)
	{
	  test_wpc.sum31 = 0;
	  test_wpc.sum29 = 0;
	  /* counts 32-bit words mod 29 and mod 31 in order of appearance
	     on file */
	  /* Here all nodes see the same sequence because we read serially */
	  rank29 = 0;
	  rank31 = 0;
	}

      destnode=node_number(x,y,z,t);

      if(this_node==0){
	/* Node 0 fills its buffer, if necessary */
	if(where_in_buf == buf_length)
	  {  /* get new buffer */
	    /* new buffer length  = remaining sites, but never bigger 
	       than MAX_BUF_LENGTH */
	    buf_length = volume - rcv_rank;
	    if(buf_length > MAX_BUF_LENGTH) buf_length = MAX_BUF_LENGTH;
	    /* then do read */
	    a=(int)fread(pbuff,sizeof(fwilson_matrix),buf_length,wpf->fp);
	    
	    if( a  != buf_length)
	      {
		
		if(status == 0)
		  printf(" node %d propagator read error %d file %s\n",
			 this_node, errno, wpf->filename); 
		fflush(stdout); 
		status = 1;
	      }
	    where_in_buf = 0;  /* reset counter */
	  }  /*** end of the buffer read ****/
	
	/* Save Wilson matrix in msg.q for further processing */
	msg.q = pbuff[where_in_buf];

	if(destnode==0){	/* just copy su3_matrix */
	  i = node_index(x,y,z,t);
	}
	else {		        /* send to correct node */
	  send_field((char *)&msg, sizeof(msg), destnode);
	}
	where_in_buf++;
      }
       /* if(this_node==0) */
      else {	/* for all nodes other than node 0 */
	if(this_node==destnode){
	  i = node_index(x,y,z,t);
	  get_field((char *)&msg, sizeof(msg),0);
	}
      }

      /* The receiving node does the byte reversal.  At this point msg
	 contains the input vectors and i points to the destination
	 site structure */
      
      if(this_node==destnode)
	{
	  if(byterevflag==1)
	    byterevn((int32type *)&msg.q, 
		     sizeof(fwilson_matrix)/sizeof(int32type));
	  
	  /* Accumulate checksums */
	  for(k = 0, val = (u_int32type *)(&msg.q); 
	      k < (int)sizeof(fwilson_matrix)/(int)sizeof(int32type); 
	      k++, val++)
	    {
	      test_wpc.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
	      test_wpc.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
	      rank29++; if(rank29 >= 29)rank29 = 0;
	      rank31++; if(rank31 >= 31)rank31 = 0;
	    }

	  /* Now copy the site data into the site structure or field
	     switching to propagator storage and converting to generic
	     precision if needed */

	  if(dest_site == (field_offset)(-1))
	    qp = dest_field + i;
	  else
	    qp = (wilson_propagator *)F_PT(&lattice[i],dest_site);

	  for(s0=0;s0<4;s0++)for(c0=0;c0<DIMF;c0++)
	    for(s1=0;s1<4;s1++)for(c1=0;c1<DIMF;c1++)
	      {
		qp->c[c0].d[s0].d[s1].c[c1].real 
		  = msg.q.d[s0].c[c0].d[s1].c[c1].real;
		qp->c[c0].d[s0].d[s1].c[c1].imag 
		  = msg.q.d[s0].c[c0].d[s1].c[c1].imag;
	      }
	}
      else
	{
	  rank29 += sizeof(fwilson_matrix)/sizeof(int32type);
	  rank31 += sizeof(fwilson_matrix)/sizeof(int32type);
	  rank29 %= 29;
	  rank31 %= 31;
	}

      /* Accumulate and print checksum at the end of each time slice */
      if(x == nx - 1 && y == ny - 1 && z == nz - 1)
	{
	  /* Combine node checksum contributions with global exclusive or */
	  g_xor32(&test_wpc.sum29);
	  g_xor32(&test_wpc.sum31);
	  
	  node0_printf("quark.t[%d].checksum  \"%0x %0x\"\n",t,
		       test_wpc.sum29, test_wpc.sum31);
	}
    } /* rcv_rank */

  /**  if(this_node==0)
    {
      printf("Read Wilson prop serially from file %s\n", wpf->filename);
      }**/
}

/*----------------------------------------------------------------------*/

void r_serial_w_fm_f(w_prop_file *wpf)

/* Close the file and free associated structures */
{
  g_sync();
  if(this_node==0)
    {
      if(wpf->fp != NULL)fclose(wpf->fp);
      fflush(stdout);
    }
  
  free(wpf->header);
  free(wpf);
  
} /* r_serial_w_fm_f */

/*--------------------------------------------------------------------*/
void r_prop_w_fm(char *filename, field_offset dest_site, 
		 wilson_propagator *dest_field)
{
  w_prop_file *wpf;
  wpf = r_serial_w_fm_i(filename);
  r_serial_w_fm(wpf, dest_site, dest_field);
  r_serial_w_fm_f(wpf);
  
}

/*--------------------------------------------------------------------*/
void r_prop_w_fm_to_site(char *filename, field_offset dest_site)
{
  r_prop_w_fm(filename, dest_site, NULL);
}

/*--------------------------------------------------------------------*/
void r_prop_w_fm_to_field(char *filename, wilson_propagator *dest_field)
{
  r_prop_w_fm(filename, (field_offset)(-1), dest_field);
}

/*--------------------------------------------------------------------*/
