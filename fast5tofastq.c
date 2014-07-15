/*
    fast5tofastq
    Copyright (C) 2014 German Tischler
    Copyright (C) 2014 Genome Research Limited

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include <assert.h>

int main(int argc, char * argv[])
{
	/* file name */
	char * fn = NULL;
	/* set of hdf paths we check for fastq data */
        char const * fastq_paths[] =
        {
	        "/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq",
	        "/Analyses/Basecall_2D_000/BaseCalled_template/Fastq",
	        "/Analyses/Basecall_2D_000/BaseCalled_complement/Fastq",
	        0
        };
	/* iterator for fastq_paths */
	char const ** fqp = NULL;
	/* still searching for fastq */
        int running = 1;
        /* file handle */
        hid_t file_id;
        /* op status */
        herr_t status;

	/* check for input file name */
	if ( !(1<argc) )
	{
		fprintf(stderr,"usage: %s <in.fast5>\n", argv[0]);
		return EXIT_FAILURE;
	}

	/* set input file name */
	fn = argv[1];
	
	/* open file with default properties */
	file_id = H5Fopen(fn, H5F_ACC_RDONLY, H5P_DEFAULT);
        
        /* check if opening file was succeful */
        if ( file_id < 0 )
        {
        	fprintf(stderr,"Failed to open/init input file %s\n", fn);
        	return EXIT_FAILURE;
        }

        /* try list of paths to fastq data */
        for ( fqp = &fastq_paths[0]; running && *fqp; ++fqp )
        {
        	/* fastq path */
        	char const * path = *fqp;

        	/* get root of hdf5 */
	        hid_t group_id = H5Gopen(file_id, "/",H5P_DEFAULT);
        
	        /* this should work */
	        if ( group_id < 0 )
	        {
	        	fprintf(stderr,"failed to open group /\n");
	        	break;
	        }
	 
	 	/* if path is valid/contained in file */       
	        if ( H5LTpath_valid ( group_id, path, 1) )
	        {
			hsize_t dims;
			H5T_class_t class_id;
			size_t type_size;
			
			/* get info on data set (length of string) */
			herr_t status = H5LTget_dataset_info( group_id, path, &dims, &class_id, &type_size );
						
			if ( status >= 0 && class_id == H5T_STRING )
			{
				/* allocate memory for string */
				char * p = (char*)calloc(type_size,1);
		        	assert ( p );

		        	/* get data */
				status = H5LTread_dataset_string (group_id, path, p);
			
				/* write data if we have received it */
				if ( status >= 0 )
					fwrite(p,type_size,1,stdout);
			
				/* free memory */
				free(p);
				
				/* break loop */
				running = 0;
			}
	        }
	        
	        /* close group handle */
	        H5Gclose(group_id);
	}

	/* close file */
	status = H5Fclose(file_id);
	
	if ( status < 0 )
	{
		fprintf(stderr,"Failed to close HDF file %s\n", fn);
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
