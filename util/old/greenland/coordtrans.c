
#include <stdlib.h>
#include <stdio.h>
#include <proj_api.h>
#include <netcdf.h>


void nc_check( int status ) {
  if (status != NC_NOERR) {
    printf( "%s\n" , nc_strerror(status) );
    exit(1);
  }
}

// no actual input arguments... should just say void
main( int argc , char **argv ) {
  projPJ  pj_velxy , pj_latlon , pj_bxy;
  int     nx = 3010;
  int     ny = 5460;
  double  x0 = -645000;
  double  y0 = -3370000;
  double  *x , *y;
  int     i , j , p;
  int     ncid, xdim, ydim, xvar, yvar, xvar2, yvar2;
  int     ncerr;
  int     dims[2];

  x = (double *) malloc( nx*ny*sizeof(double) );
  y = (double *) malloc( nx*ny*sizeof(double) );

  for ( j=0 ; j < ny ; j++ ) {
    for ( i=0 ; i < nx ; i++ ) {
      x[j*nx+i] = x0+i*500.;
      y[j*nx+i] = y0+j*500.;
    }
  }

  if ( !( pj_velxy  = pj_init_plus("+proj=stere +ellps=WGS84 +datum=WGS84 +lat_ts=70 +lat_0=90 +lon_0=-45") ) ) { 
      printf("bad velxy \n"); 
      exit(1); 
  }
  if ( !( pj_bxy    = pj_init_plus("+proj=stere +ellps=WGS84 +datum=WGS84 +lat_ts=71 +lat_0=90 +lon_0=-39") ) ) { 
      printf("bad bxy   \n"); 
      exit(1); 
  }
  if ( !( pj_latlon = pj_init_plus("+proj=latlong +ellps=WGS84"                                           ) ) ) { 
      printf("bad latlon\n"); 
      exit(1); 
  }
 
  p = pj_transform( pj_velxy  , pj_latlon , nx*ny , 1 , x , y , NULL );  
  if ( p != 0) { 
      printf( "bad velxy2latlon \n"); 
      printf( "%s\n" , pj_strerrno(p) ); 
      exit(1); 
  }
  p = pj_transform( pj_latlon , pj_bxy    , nx*ny , 1 , x , y , NULL );  
  if ( p != 0) { 
      printf( "bad latlon2bxy   \n"); 
      printf( "%s\n" , pj_strerrno(p) ); 
      exit(1); 
  }

  nc_check( nc_create( "coords.nc" , NC_CLOBBER , &ncid ) );
  nc_check( nc_def_dim( ncid, "x" , nx , &xdim ) );
  nc_check( nc_def_dim( ncid, "y" , ny , &ydim ) );
  dims[0] = ydim;
  dims[1] = xdim;
  nc_check( nc_def_var( ncid , "x" , NC_DOUBLE , 2 , dims , &xvar ) );
  nc_check( nc_def_var( ncid , "y" , NC_DOUBLE , 2 , dims , &yvar ) );
  nc_check( nc_def_var( ncid , "xvar" , NC_DOUBLE , 2 , dims , &xvar2 ) );
  nc_check( nc_def_var( ncid , "yvar" , NC_DOUBLE , 2 , dims , &yvar2 ) );
  nc_check( nc_enddef( ncid ) );
  nc_check( nc_put_var_double( ncid , xvar , x ) );
  nc_check( nc_put_var_double( ncid , yvar , y ) );
  nc_check( nc_put_var_double( ncid , xvar2 , x ) );
  nc_check( nc_put_var_double( ncid , yvar2 , y ) );
  nc_check( nc_close( ncid ) );

  exit( EXIT_SUCCESS );

}


