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
#include <ctype.h>
#include <string.h>

#define _XOPEN_SOURCE       /* See feature_test_macros(7) */
#define __USE_XOPEN
#include <time.h>
#include <locale.h>

/* 2014-Aug-22 09:30:27 */
static char const * timeformatstring = "%Y-%b-%d %H:%M:%S";

static int getSequenceLength(char const * fq, int len)
{
	int seqlen = 0;
	
	if ( ! len )
		return -1;
	if ( *fq != '@' )
		return -1;
	while ( len && *fq != '\n' )
		++fq, --len;
	while ( len && isspace(*fq) )
		++fq, --len;
	while ( len && *fq != '+' )
	{
		if ( ! isspace(*fq) )
		{
			++seqlen;
		}
		--len;
		++fq;
	}
	
	return seqlen;
}

static long getTimeStamp(hid_t root_id)
{
        char const * basecall_path = "/Analyses/Basecall_2D_000";

        if ( H5LTpath_valid(root_id, basecall_path, 1) )
        {
		hid_t g1 = -1;
		hid_t g2 = -1;
		hid_t a1 = -1;
		herr_t e1 = -1;
		H5A_info_t info;
		hid_t t1 = -1;
		/* H5T_class_t c1 = -1; */
		/* size_t s1 = -1; */
		htri_t v1 = -1;
		long int v = -1;

		g1 = H5Gopen(root_id,"Analyses",H5P_DEFAULT);
		
		if ( g1 < 0 )
		{
			fprintf(stderr,"Failed to open group Analyses: %d\n", g1);
			goto cleanup;
		}
		
		g2 = H5Gopen(g1,"Basecall_2D_000",H5P_DEFAULT);
		
		if ( g2 < 0 )
		{
			fprintf(stderr,"Failed to open group Basecall_2D_000: %d\n", g2);
			goto cleanup;	
		}

		a1 = H5Aopen_by_name(g2,".","time_stamp",H5P_DEFAULT,0);

		if ( a1 < 0 )
		{
			fprintf(stderr,"Failed to open attribute time_stamp: %d\n", a1);
			goto cleanup;	
		}

		memset(&info,0,sizeof(info));
		
		e1 = H5Aget_info(a1,&info);
		
		if ( e1 < 0 )
		{
			fprintf(stderr,"Failed to get info for attribute time_stamp: %d\n", e1);
			goto cleanup;			
		}

		t1 = H5Aget_type(a1);

		if ( t1 < 0 )
		{
			fprintf(stderr,"Failed to get type for attribute time_stamp: %d\n", t1);
			goto cleanup;			
		}

		/* c1 = H5Tget_class(t1); */
		
		/* s1 = H5Tget_size(t1); */
		
		v1 = H5Tis_variable_str(t1);

		if ( v1 >= 0 && v1 )
		{
			char const * tsmem = 0;
			herr_t const er = H5Aread(a1,t1,&tsmem);
			size_t tslen = 0;
			
			if ( er < 0 )
			{
				fprintf(stderr,"Failed to read attribute time_stamp: %d\n", er);
				goto cleanup;
			}
			
			tslen = strlen(tsmem);

			if ( tslen > 0 )
			{
				struct tm tm;
				/* char * nonproc = */ strptime(tsmem, timeformatstring, &tm);
				char secmem[256];
				char * secmemendptr = 0;
				size_t const secmemfilled = strftime(&secmem[0], sizeof(secmem), "%s",&tm);
				secmem[secmemfilled] = 0;
				
				v = strtoul(secmem,&secmemendptr,10);
			}
		}
		else
		{
			fprintf(stderr,"/Analyses/Basecall_2D_000 attribute time_stamp is not a variable length strength or invalid.\n");
		}

		cleanup:
		if ( t1 >= 0 )
			H5Tclose(t1);
		if ( a1 >= 0 )
			H5Aclose(a1);
		if ( g2 >= 0 )
			H5Gclose(g2);		
		if ( g1 >= 0 )
			H5Gclose(g1);
				        	
		return v;
        }
        else
        {
        	fprintf(stderr,"[E] no /Analyses/Basecall_2D_000 -> no time stamp\n");
        	return -1;
        }
}

static long getExpStartTime(hid_t root_id)
{
	hid_t g1 = -1;
	hid_t g2 = -1;
	hid_t a1 = -1;
	herr_t e1 = -1;
	hid_t t1 = -1;
	long r = -1;
	H5A_info_t info;
	/* H5T_class_t c1 = -1; */
	/* size_t s1 = -1; */
	htri_t v1 = -1;
	char const * s = 0;
	
	g1 = H5Gopen(root_id,"UniqueGlobalKey",H5P_DEFAULT);
	
	if ( g1 < 0 )
	{
		fprintf(stderr,"Failed to open group UniqueGlobalKey: %d\n", g1);
		goto cleanup;
	}
	
	g2 = H5Gopen(g1,"tracking_id",H5P_DEFAULT);
	
	if ( g2 < 0 )
	{
		fprintf(stderr,"Failed to open group tracking_id: %d\n", g2);
		goto cleanup;	
	}
	
	a1 = H5Aopen_by_name(g2,".","exp_start_time",H5P_DEFAULT,0);

	if ( a1 < 0 )
	{
		fprintf(stderr,"Failed to open attribute exp_start_time: %d\n", a1);
		goto cleanup;	
	}
	
	memset(&info,0,sizeof(info));
	
	e1 = H5Aget_info(a1,&info);
	
	if ( e1 < 0 )
	{
		fprintf(stderr,"Failed to get info for attribute exp_start_time: %d\n", e1);
		goto cleanup;			
	}
	
	t1 = H5Aget_type(a1);

	if ( t1 < 0 )
	{
		fprintf(stderr,"Failed to get type for attribute exp_start_time: %d\n", t1);
		goto cleanup;			
	}
	
	/* c1 = H5Tget_class(t1); */
	
	/* s1 = H5Tget_size(t1); */
	
	v1 = H5Tis_variable_str(t1);
	
	if ( v1 >= 0 && v1 )
	{
		herr_t const er = H5Aread(a1,t1,&s);
		size_t tslen = 0;
		char * endptr = 0;
		
		if ( er < 0 )
		{
			fprintf(stderr,"Failed to read attribute exp_start_time: %d\n", er);
			goto cleanup;
		}
		
		/* fprintf(stderr,"ok: %d %d %d %s\n", (int)(info.data_size), (int)c1, (int)s1, s); */
		
		tslen = strlen(s);

		r = strtoul(s,&endptr,10);
			
		if ( endptr != s + tslen )
		{
			fprintf(stderr,"Failed to parse string %s as number\n", s);
			goto cleanup;
		}

	}
	
	cleanup:
	if ( t1 >= 0 )
		H5Tclose(t1);
	if ( a1 >= 0 )
		H5Aclose(a1);
	if ( g2 >= 0 )
		H5Gclose(g2);		
	if ( g1 >= 0 )
		H5Gclose(g1);
	return r;
}

double getEventsStartTime(hid_t root_id)
{
	static char const * events_path = "/Analyses/Basecall_2D_000/BaseCalled_template/Events";

        if ( H5LTpath_valid(root_id, events_path, 1) )
        {
        	hsize_t attrdims = 0;
		H5T_class_t attrclass_id;
		size_t attr_type_size;
		herr_t attrerr = H5LTget_attribute_info(root_id, events_path, "start_time",  &attrdims, &attrclass_id, &attr_type_size );
				        	
		if ( attrerr >= 0 && attrclass_id == H5T_FLOAT )
		{
			herr_t err = -1;
			double data = -1;
			
			err = H5LTget_attribute_double(root_id,events_path,"start_time",&data);
			
			if ( err < 0 )
				return -1;
			
			return data;			
		}
		else
		{
			return -1;
		}
        }
        else
        {
        	return -1;
        }
}

static long getChannelAsChannelId(hid_t root_id)
{
	static char const * events_path = "/UniqueGlobalKey/channel_id";
	/* /Analyses/Basecall_2D_000/Configuration/general@channel */

        if ( H5LTpath_valid(root_id, events_path, 1) )
        {
        	hsize_t attrdims = 0;
		H5T_class_t attrclass_id;
		size_t attr_type_size;
		herr_t attrerr = H5LTget_attribute_info(root_id, events_path, "channel_number",  &attrdims, &attrclass_id, &attr_type_size );
				        	
		if ( attrerr >= 0 && attrclass_id == H5T_INTEGER )
		{
			herr_t err = -1;
			unsigned long data = 0;
			
			err = H5LTget_attribute_ulong(root_id,events_path,"channel_number",&data);
			
			if ( err < 0 )
				return -1;
			
			return data;			
		}
		else
		{
			return -1;
		}
        }
        else
        {
        	return -1;
        }
}

static long getChannelAsChannel(hid_t root_id)
{
	static char const * events_path = "/Analyses/Basecall_2D_000/Configuration/general";

        if ( H5LTpath_valid(root_id, events_path, 1) )
        {
		hid_t g1 = -1;
		hid_t g2 = -1;
		hid_t g3 = -1;
		hid_t g4 = -1;
		long int r = -1;
		
		hid_t a1 = -1;
		herr_t e1 = -1;
		hid_t t1 = -1;
		H5A_info_t info;
		/* H5T_class_t c1 = -1; */
		/* size_t s1 = -1; */
		htri_t v1 = -1;
		/* char const * s = 0; */
		
		g1 = H5Gopen(root_id,"Analyses",H5P_DEFAULT);
		
		if ( g1 < 0 )
		{
			fprintf(stderr,"Failed to open group Analyses: %d\n", g1);
			goto cleanup;
		}
		
		g2 = H5Gopen(g1,"Basecall_2D_000",H5P_DEFAULT);
		
		if ( g2 < 0 )
		{
			fprintf(stderr,"Failed to open group Basecall_2D_000: %d\n", g2);
			goto cleanup;	
		}

		g3 = H5Gopen(g2,"Configuration",H5P_DEFAULT);
		
		if ( g3 < 0 )
		{
			fprintf(stderr,"Failed to open group Configuration: %d\n", g3);
			goto cleanup;	
		}
		
		g4 = H5Gopen(g3,"general",H5P_DEFAULT);
		
		if ( g4 < 0 )
		{
			fprintf(stderr,"Failed to open group general: %d\n", g4);
			goto cleanup;	
		}

		a1 = H5Aopen_by_name(g4,".","channel",H5P_DEFAULT,0);

		if ( a1 < 0 )
		{
			fprintf(stderr,"Failed to open attribute channel: %d\n", a1);
			goto cleanup;	
		}
		
		memset(&info,0,sizeof(info));
		
		e1 = H5Aget_info(a1,&info);
		
		if ( e1 < 0 )
		{
			fprintf(stderr,"Failed to get info for attribute channel: %d\n", e1);
			goto cleanup;			
		}
		
		t1 = H5Aget_type(a1);

		if ( t1 < 0 )
		{
			fprintf(stderr,"Failed to get type for attribute channel: %d\n", t1);
			goto cleanup;			
		}
		
		/* c1 = H5Tget_class(t1); */
		
		/* s1 = H5Tget_size(t1); */
		
		v1 = H5Tis_variable_str(t1);

		if ( v1 >= 0 && v1 )
		{
			char const * s = 0;
			herr_t const er = H5Aread(a1,t1,&s);
			size_t tslen = 0;
			char * endptr = 0;
			
			if ( er < 0 )
			{
				fprintf(stderr,"Failed to read attribute channel: %d\n", er);
				goto cleanup;
			}
			
			/* fprintf(stderr,"ok: %d %d %d %s\n", (int)(info.data_size), (int)c1, (int)s1, s); */
			
			tslen = strlen(s);

			r = strtoul(s,&endptr,10);
				
			if ( endptr != s + tslen )
			{
				fprintf(stderr,"Failed to parse string %s as number\n", s);
				goto cleanup;
			}
		}
		else
		{
			fprintf(stderr,"Attribute channel in /Analyses/Basecall_2D_000/Configuration/general is not a variable length string.\n");		
			goto cleanup;
		}

		cleanup:
		if ( g4 >= 0 )
			H5Gclose(g4);		
		if ( g3 >= 0 )
			H5Gclose(g3);		
		if ( g2 >= 0 )
			H5Gclose(g2);		
		if ( g1 >= 0 )
			H5Gclose(g1);
		
		return r;		
        }
        else
        {
        	return -1;
        }
}

static long getChannel(hid_t root_id)
{
	long channel = -1;
	
	channel = getChannelAsChannelId(root_id);
	
	if ( channel >= 0 )
		return channel;

	channel = getChannelAsChannel(root_id);
	
	if ( channel >= 0 )
		return channel;
	
	return -1;
}

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
        /* char const * basecall_path = "/Analyses/Basecall_2D_000"; */
        /* char const * trackingid_path = "/UniqueGlobalKey/tracking_id"; */
	/* iterator for fastq_paths */
	char const ** fqp = NULL;
	/* still searching for fastq */
        int running = 1;
        /* file handle */
        hid_t file_id;
        /* op status */
        herr_t status;
        /* time stamp */
        char * tsmem = 0; 
        /* int tslen = -1; */
        /* sequence length */
        int seqlen = -1;

	/* check for input file name */
	if ( !(1<argc) )
	{
		fprintf(stderr,"usage: %s <in.fast5>\n", argv[0]);
		return EXIT_FAILURE;
	}

	setlocale(LC_ALL, "C");

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
				{
					fwrite(p,type_size,1,stdout);
					
					seqlen = getSequenceLength(p,type_size);

					fprintf(stderr,"[V] extracted read of type %s from file %s time stamp %ld sequence length %d exp start time %ld events start %f channel %d\n", path, fn, getTimeStamp(group_id), seqlen, getExpStartTime(group_id), getEventsStartTime(group_id), (int)getChannel(group_id));
				}
			
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
	
	if ( tsmem )
	{
		free(tsmem);
	}
	
	return EXIT_SUCCESS;
}
