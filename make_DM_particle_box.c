#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
//CC make_DM_particle_box.c create_lognormal_halos.c integration.c util.c -lm -o dm_box -w
#include "header_DM.h"

char *pmfile = cvector(0,200);
char *fp_out = cvector(0,200);

int main(int argc, char *argv[]) {

  sprintf(pmfile, "/global/scratch2/sd/mwhite/QPM/new/pm_0.5882_0001.subsample");

  create_lognormal_halos(pmfile, fp_out);

}
  //Files.pmfile = 


