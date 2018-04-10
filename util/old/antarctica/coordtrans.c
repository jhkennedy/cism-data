
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


void print( double *x , int nx , int ny ) {
  int i,j;
  for (i=0; i<nx; i++) {
    for (j=0; j<ny; j++) {
      printf("%lf\n",x[j*nx+i]);
    }
  }
}


main( int argc , char **argv ) {
  projPJ  pj_velxy , pj_latlon , pj_bxy;
  int     nx = 262;
  int     ny = 240;
  double  *x , *y;
  int     i , j , p;
  int     ncid, xdim, ydim, xvar, yvar, lonvar, latvar;
  int     ncerr;
  int     dims[2];

  x = (double *) malloc( nx*ny*sizeof(double) );
  y = (double *) malloc( nx*ny*sizeof(double) );

  nc_check( nc_open( "./originals/racmo/1979-2010-average-SMB_AIS.nc" , NC_NOWRITE , &ncid ) );
  nc_check( nc_inq_varid( ncid, "lon" , &xvar) );
  nc_check( nc_inq_varid( ncid, "lat" , &yvar) );
  nc_check( nc_get_var_double( ncid, xvar , x ) );
  nc_check( nc_get_var_double( ncid, yvar , y ) );
  nc_check( nc_close( ncid ) );

  nc_check( nc_create( "coords.nc" , NC_CLOBBER , &ncid ) );
  nc_check( nc_def_dim( ncid, "x" , nx , &xdim ) );
  nc_check( nc_def_dim( ncid, "y" , ny , &ydim ) );
  dims[0] = ydim;
  dims[1] = xdim;
  nc_check( nc_def_var( ncid , "xvar" , NC_DOUBLE , 2 , dims , &xvar ) );
  nc_check( nc_def_var( ncid , "yvar" , NC_DOUBLE , 2 , dims , &yvar ) );
  nc_check( nc_def_var( ncid , "lonvar" , NC_DOUBLE , 2 , dims , &lonvar ) );
  nc_check( nc_def_var( ncid , "latvar" , NC_DOUBLE , 2 , dims , &latvar ) );
  nc_check( nc_enddef( ncid ) );

  nc_check( nc_put_var_double( ncid , lonvar , x ) );
  nc_check( nc_put_var_double( ncid , latvar , y ) );

  for (i=0; i<nx; i++) {
    for (j=0; j<ny; j++) {
      x[j*nx+i] = x[j*nx+i] / 180. * 3.1415926535897932384626433832795028;
      y[j*nx+i] = y[j*nx+i] / 180. * 3.1415926535897932384626433832795028;
    }
  }
  if ( !( pj_latlon = pj_init_plus("+proj=latlong +ellps=WGS84"                                ) ) ) { printf("bad latlon\n"); exit(1); }
  if ( !( pj_bxy    = pj_init_plus("+proj=stere +ellps=WGS84 +datum=WGS84 +lat_ts=-71 +lat_0=-90 +lon_0=0") ) ) { printf("bad bxy   \n"); exit(1); }
  p = pj_transform( pj_latlon , pj_bxy , nx*ny , 1 , x , y , NULL );  if ( p != 0) { printf( "bad latlon2bxy   \n"); printf( "%s\n" , pj_strerrno(p) ); exit(1); }

  nc_check( nc_put_var_double( ncid , xvar , x ) );
  nc_check( nc_put_var_double( ncid , yvar , y ) );

  nc_check( nc_close( ncid ) );

  exit( EXIT_SUCCESS );

}


