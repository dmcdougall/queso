/* This is part of the netCDF package.
   Copyright 2016 University Corporation for Atmospheric Research/Unidata
   See COPYRIGHT file for conditions of use.

   Tests library ability to open files generated by a netcdf
   instance linked against libhdf5 1.10.0.  This is an issue at this point
   because, without specifying libver bounds, netcdf linked against
   libhdf5 1.8 cannot read those generated by libhdf5 1.10.0.

   This test will undoubtedly not age well, but right now it is
   fairly critical, and will help test if we've corrected the
   issue.

   Details: https://github.com/Unidata/netcdf-c/issues/250

   Files testing against (copied from nc_test4):

   * ref_hdf5_compat1.nc <- tst_vars.nc
   * ref_hdf5_compat2.nc <- tst_vars4.nc
   * ref_hdf5_compat3.nc <- tst_compounds.nc

   */

#include <nc_tests.h>
#include "err_macros.h"
#include "netcdf.h"


#define FILE_NAME1 "ref_hdf5_compat1.nc"
#define FILE_NAME2 "ref_hdf5_compat2.nc"
#define FILE_NAME3 "ref_hdf5_compat3.nc"


int main(int argc, char **argv) {

  int res = 0;
  int ncid = 0;

  printf("\n*** Testing libhdf5 file compatibility (open files generated by hdf5 1.10).\n");

  {
    printf("Testing %s\n",FILE_NAME1);
    if (nc_open(FILE_NAME1, NC_NOWRITE, &ncid)) ERR;
    if (nc_close(ncid)) ERR;
  }

  {
    printf("Testing %s\n",FILE_NAME2);
    if (nc_open(FILE_NAME2, NC_NOWRITE, &ncid)) ERR;
    if (nc_close(ncid)) ERR;
  }

  {
    printf("Testing %s\n",FILE_NAME3);
    if (nc_open(FILE_NAME3, NC_NOWRITE, &ncid)) ERR;
    if (nc_close(ncid)) ERR;
  }



  SUMMARIZE_ERR;
  FINAL_RESULTS;

}
