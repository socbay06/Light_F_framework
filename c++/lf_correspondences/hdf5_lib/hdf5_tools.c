#include "hdf5_tools.h"






float *hdf5_read_lf( const char* file_name, const char* dset_name,
			   size_t &W, size_t &H, size_t &S, size_t &T, size_t &C )
{
  // Open the light field data sets
  // HDF5 library differs across Ubuntu versions.
  // Try to uncomment if build fails for you.
  //hid_t dset = H5Dopen (file, dset_name.c_str(), H5P_DEFAULT );


  /* Open file for reading */
  hid_t file_id = H5Fopen(file_name, H5F_ACC_RDWR, H5P_DEFAULT);
  /* open dataset */
  hid_t dset_id = H5Dopen (file_id, dset_name, H5P_DEFAULT );
  if ( dset_id == 0 ) {
    fprintf(stderr, "could not open light field data set %s \n",dset_name);
    return NULL;
  }
  fprintf(stdout,"dataset %s opened.\n",dset_name);

  // Retrieve dimension attributes
  int ndims;
  H5LTget_dataset_ndims(file_id, dset_name, &ndims );
  if ( ndims != 4 && ndims!= 5 ) {
    fprintf(stderr, "light field data set %s is not four or five dimensional.",dset_name );
    return NULL;
  }

  hsize_t dims[5];
  H5LTget_dataset_info(file_id, dset_name, dims, NULL,NULL);
  if(ndims==5){
    C = dims[0];
    W = dims[1];
    H = dims[2];
    S = dims[3];
    T = dims[4];
  }else{
    C = 0;
    W = dims[0];
    H = dims[1];
    S = dims[2];
    T = dims[3];
  }
 fprintf(stdout, "Dataset %s with views %d x %d and resolution %d x %d, no color %d \n",dset_name,T,S,H,W,C);

  // create buffer for light field data
  float * data = new float[(C==0?1:C)*S*T*W*H ];
  H5Dread( dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,H5P_DEFAULT, data );
  H5Dclose( dset_id );
  H5Fclose(file_id);
  return data;
}


/*
* Check if group exists and create it recursively.
*/
void group_create_recursive(hid_t file_id, char* group_path, int pos_end){
	char group_name[80];
	char current_group_path[512];
	int pos, group_id;
	int status;
	// Seperate group name.
	pos=-1;
	for(int i=pos_end; i>=0; i--){
		if(group_path[i]=='/'){
			pos = i;
			break;
		}
	}
	if(pos<=-1)
		return;

	strcpy(group_name, &group_path[pos+1]);
	group_name[pos_end-pos]=NULL;
	if(pos>0) // there still be parent group
		group_create_recursive(file_id, group_path, pos-1);
	//reach root now. check and create group for current group.
	strcpy(current_group_path, group_path);
	current_group_path[pos_end+1]=NULL;
	status = H5Gget_objinfo(file_id,current_group_path,0,NULL);
	printf("Group %s : ", current_group_path);
	if(status == 0) printf(" exists\n");
	else{ //let create it
		//create group now
		group_id = H5Gcreate(file_id, current_group_path, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
		printf(" not exists but created\n");
	}

}

/*
 * very simple function to write data into Hdf5
 *
 *
*/
int hdf5_write_data_simple(const char* file_name, const char* dset_path,int* dimarray, int dimsize, float* data){
	hid_t file_id, dataspace_id, dataset_id,group_id; /* identifiers*/
	int status;
	char group_name[512];
	char dset_name[80];
	size_t pathlen;
	hsize_t  *dims;
	dims = new hsize_t[dimsize];



		// try to open h5 file
	file_id = H5Fopen(file_name, H5F_ACC_RDWR,H5P_DEFAULT);
	if(file_id<0){//file not exist try to create it.
		printf("File not found, try to create %s \n",file_name);
		file_id = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	}else
		printf("Existing %s HDF5 file.\n",file_name);

	//copy dimensions
	for(int i =0; i<dimsize;i++){
		dims[i] = dimarray[i];
	}

	// Seperate group dataset.
	pathlen = strlen(dset_path);
	int pos_dset=-1;
	for(int i=pathlen-1; i>=0; i--){
		if(dset_path[i]=='/'){
			pos_dset = i;
			break;
		}
	}
	if(pos_dset!=-1){
		strcpy(dset_name, &dset_path[pos_dset+1]);
		strcpy(group_name, dset_path);
		group_name[pos_dset] =NULL;
	}
	printf(" Got dset name : %s group name: %s\n", dset_name,group_name );
	printf(" *** Start checking and create group \n");
	group_create_recursive(file_id, group_name,pos_dset-1);
	//printf("finished create group %s \n",group_name);

	//create data space
	dataspace_id = H5Screate_simple(dimsize,dims,NULL);
	printf(" finish create dataspace_id %d dim %d \n",dataspace_id,dimsize);
	//create dataset.
	dataset_id = H5Dcreate2(file_id,dset_path,H5T_NATIVE_FLOAT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT );
	printf(" finish create dataset_id %d \n",dataset_id);
	// start to write data into set
	status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	printf(" finish write data \n");
	free(dims);

	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);
	H5Fclose(file_id);
	return 0;
}


/*
 * Process a group and all it's members
 *
 *   This can be used as a model to implement different actions and
 *   searches.
 */

void scan_group(hid_t gid) {
	int i;
	ssize_t len;
	hsize_t nobj;
	herr_t err;
	int otype;
	hid_t grpid, mytypeid, dsid;
	char group_name[MAX_NAME];
	char memb_name[MAX_NAME];

        /*
         * Information about the group:
         *  Name and attributes
         *
         *  Other info., not shown here: number of links, object id
         */
	len = H5Iget_name(gid, group_name, MAX_NAME  );

	printf("Group Name: %s\n",group_name);

	/*
	 *  process the attributes of the group, if any.
	 */
        scan_attrs(gid);

        /*
         *  Get all the members of the groups, one at a time.
         */
	err = H5Gget_num_objs(gid, &nobj);
	for (i = 0; i < nobj; i++) {
                /*
                 *  For each object in the group, get the name and
                 *   what type of object it is.
                 */
		printf("  Member: %d ",i);fflush(stdout);
		len = H5Gget_objname_by_idx(gid, (hsize_t)i,
			memb_name, (size_t)MAX_NAME );
		printf("   %ld ",len);fflush(stdout);
		printf("  Member: %s ",memb_name);fflush(stdout);
		otype =  H5Gget_objtype_by_idx(gid, (size_t)i );

                /*
                 * process each object according to its type
                 */
		switch(otype) {
			case H5G_LINK:
				printf(" SYM_LINK:\n");
				do_link(gid,memb_name);
				break;
			case H5G_GROUP:
				printf(" GROUP:\n");
				grpid = H5Gopen(gid,memb_name,H5P_DEFAULT);
				scan_group(grpid);
				H5Gclose(grpid);
				break;
			case H5G_DATASET:
				printf(" DATASET:\n");
				dsid = H5Dopen(gid,memb_name,H5P_DEFAULT);
				do_dset(dsid);
				H5Dclose(dsid);
				break;
			case H5G_TYPE:
				printf(" DATA TYPE:\n");
				mytypeid = H5Topen(gid,memb_name,H5P_DEFAULT);
				do_dtype(mytypeid);
				H5Tclose(mytypeid);
				break;
			default:
				printf(" unknown?\n");
				break;
			}

		}
}

/*
 *  Retrieve information about a dataset.
 *
 *  Many other possible actions.
 *
 *  This example does not read the data of the dataset.
 */
void
do_dset(hid_t did)
{
	hid_t tid;
	hid_t pid;
	hid_t sid;
	hsize_t size;
	char ds_name[MAX_NAME];

        /*
         * Information about the group:
         *  Name and attributes
         *
         *  Other info., not shown here: number of links, object id
         */
	H5Iget_name(did, ds_name, MAX_NAME  );
	printf("Dataset Name : ");
	puts(ds_name);
	printf("\n");

	/*
	 *  process the attributes of the dataset, if any.
	 */
	scan_attrs(did);

	/*
	 * Get dataset information: dataspace, data type
	 */
	sid = H5Dget_space(did); /* the dimensions of the dataset (not shown) */
	tid = H5Dget_type(did);
	printf(" DATA TYPE:\n");
	do_dtype(tid);

	/*
	 * Retrieve and analyse the dataset properties
	 */
	pid = H5Dget_create_plist(did); /* get creation property list */
	do_plist(pid);
	size = H5Dget_storage_size(did);
	printf("Total space currently written in file: %d\n",(int)size);

	/*
	 * The datatype and dataspace can be used to read all or
	 * part of the data.  (Not shown in this example.)
	 */

	  /* ... read data with H5Dread, write with H5Dwrite, etc. */

	H5Pclose(pid);
	H5Tclose(tid);
	H5Sclose(sid);
}

/*
 *  Analyze a data type description
 */
void
do_dtype(hid_t tid) {

	H5T_class_t t_class;
	t_class = H5Tget_class(tid);
	if(t_class < 0){
		puts(" Invalid datatype.\n");
	} else {
		/*
		 * Each class has specific properties that can be
		 * retrieved, e.g., size, byte order, exponent, etc.
		 */
		if(t_class == H5T_INTEGER) {
		      puts(" Datatype is 'H5T_INTEGER'.\n");
			/* display size, signed, endianess, etc. */
		} else if(t_class == H5T_FLOAT) {
		      puts(" Datatype is 'H5T_FLOAT'.\n");
			/* display size, endianess, exponennt, etc. */
		} else if(t_class == H5T_STRING) {
		      puts(" Datatype is 'H5T_STRING'.\n");
			/* display size, padding, termination, etc. */
		} else if(t_class == H5T_BITFIELD) {
		      puts(" Datatype is 'H5T_BITFIELD'.\n");
			/* display size, label, etc. */
		} else if(t_class == H5T_OPAQUE) {
		      puts(" Datatype is 'H5T_OPAQUE'.\n");
			/* display size, etc. */
		} else if(t_class == H5T_COMPOUND) {
		      puts(" Datatype is 'H5T_COMPOUND'.\n");
			/* recursively display each member: field name, type  */
		} else if(t_class == H5T_ARRAY) {
		      puts(" Datatype is 'H5T_COMPOUND'.\n");
			/* display  dimensions, base type  */
		} else if(t_class == H5T_ENUM) {
		      puts(" Datatype is 'H5T_ENUM'.\n");
			/* display elements: name, value   */
		} else  {
		      puts(" Datatype is 'Other'.\n");
		      /* and so on ... */
		}
	}
}


/*
 *  Analyze a symbolic link
 *
 * The main thing you can do with a link is find out
 * what it points to.
 */
void
do_link(hid_t gid, char *name) {
	herr_t status;
	char target[MAX_NAME];

	status = H5Gget_linkval(gid, name, MAX_NAME, target  ) ;
	printf("Symlink: %s points to: %s\n", name, target);
}

/*
 *  Run through all the attributes of a dataset or group.
 *  This is similar to iterating through a group.
 */
void
scan_attrs(hid_t oid) {
	int na;
	hid_t aid;
	int i;

	na = H5Aget_num_attrs(oid);

	for (i = 0; i < na; i++) {
		aid =	H5Aopen_idx(oid, (unsigned int)i );
		do_attr(aid);
		H5Aclose(aid);
	}
}

/*
 *  Process one attribute.
 *  This is similar to the information about a dataset.
 */
void do_attr(hid_t aid) {
	ssize_t len;
	hid_t atype;
	hid_t aspace;
	char buf[MAX_NAME];

	/*
	 * Get the name of the attribute.
	 */
	len = H5Aget_name(aid, MAX_NAME, buf );
	printf("    Attribute Name : %s\n",buf);

	/*
	 * Get attribute information: dataspace, data type
	 */
	aspace = H5Aget_space(aid); /* the dimensions of the attribute data */

	atype  = H5Aget_type(aid);
	do_dtype(atype);

	/*
	 * The datatype and dataspace can be used to read all or
	 * part of the data.  (Not shown in this example.)
	 */

	  /* ... read data with H5Aread, write with H5Awrite, etc. */

	H5Tclose(atype);
	H5Sclose(aspace);
}

/*
 *   Example of information that can be read from a Dataset Creation
 *   Property List.
 *
 *   There are many other possibilities, and there are other property
 *   lists.
 */
void
do_plist(hid_t pid) {
  unsigned int *filter_config;
	hsize_t chunk_dims_out[2];
	int  rank_chunk;
	int nfilters;
	H5Z_filter_t  filtn;
	int i;
	unsigned int   filt_flags;
	size_t cd_nelmts;
	unsigned int cd_values[32] ;
	char f_name[MAX_NAME];
	H5D_fill_time_t ft;
	H5D_alloc_time_t at;
	H5D_fill_value_t fvstatus;
	unsigned int szip_options_mask;
	unsigned int szip_pixels_per_block;

	/* zillions of things might be on the plist */
        /*  here are a few... */

	/*
	 * get chunking information: rank and dimensions.
	 *
	 *  For other layouts, would get the relevant information.
	 */
	if(H5D_CHUNKED == H5Pget_layout(pid)){
		rank_chunk = H5Pget_chunk(pid, 2, chunk_dims_out);
		printf("chunk rank %d, dimensions %lu x %lu\n", rank_chunk,
		   (unsigned long)(chunk_dims_out[0]),
		   (unsigned long)(chunk_dims_out[1]));
	} /* else if contiguous, etc. */

	/*
	 *  Get optional filters, if any.
	 *
	 *  This include optional checksum and compression methods.
	 */

	nfilters = H5Pget_nfilters(pid);
	for (i = 0; i < nfilters; i++)
	{
		/* For each filter, get
		 *   filter ID
		 *   filter specific parameters
		 */
		cd_nelmts = 32;
		filtn = H5Pget_filter(pid, (unsigned)i,
			&filt_flags, &cd_nelmts, cd_values,
			(size_t)MAX_NAME, f_name,filter_config);
  		/*
		 *  These are the predefined filters
		 */
		switch (filtn) {
			case H5Z_FILTER_DEFLATE:  /* AKA GZIP compression */
				printf("DEFLATE level = %d\n", cd_values[0]);
				break;
			case H5Z_FILTER_SHUFFLE:
				printf("SHUFFLE\n"); /* no parms */
				break;
		       case H5Z_FILTER_FLETCHER32:
				printf("FLETCHER32\n"); /* Error Detection Code */
				break;
		       case H5Z_FILTER_SZIP:
				szip_options_mask=cd_values[0];;
				szip_pixels_per_block=cd_values[1];

				printf("SZIP COMPRESSION: ");
				printf("PIXELS_PER_BLOCK %d\n",
					szip_pixels_per_block);
				 /* print SZIP options mask, etc. */
				break;
			default:
				printf("UNKNOWN_FILTER\n" );
				break;
	       }
      }

	/*
	 *  Get the fill value information:
	 *    - when to allocate space on disk
	 *    - when to fill on disk
	 *    - value to fill, if any
	 */

	printf("ALLOC_TIME ");
	H5Pget_alloc_time(pid, &at);

	switch (at)
	{
		case H5D_ALLOC_TIME_EARLY:
			printf("EARLY\n");
			break;
		case H5D_ALLOC_TIME_INCR:
			printf("INCR\n");
			break;
		case H5D_ALLOC_TIME_LATE:
			printf("LATE\n");
			break;
		default:
			printf("unknown allocation policy");
			break;
	}

	printf("FILL_TIME: ");
	H5Pget_fill_time(pid, &ft);
	switch ( ft )
	{
		case H5D_FILL_TIME_ALLOC:
			printf("ALLOC\n");
			break;
		case H5D_FILL_TIME_NEVER:
			printf("NEVER\n");
			break;
		case H5D_FILL_TIME_IFSET:
			printf("IFSET\n");
			break;
		default:
			printf("?\n");
		break;
	}


	H5Pfill_value_defined(pid, &fvstatus);

	if (fvstatus == H5D_FILL_VALUE_UNDEFINED)
	{
		printf("No fill value defined, will use default\n");
	} else {
		/* Read  the fill value with H5Pget_fill_value.
		 * Fill value is the same data type as the dataset.
		 * (detailse not shown)
		 **/
	}

	/* ... and so on for other dataset properties ... */
}
