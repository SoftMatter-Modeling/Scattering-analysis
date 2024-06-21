#ifndef READ_DUMP_FILE
#define READ_DUMP_FILE
/*
CHANGE LOG
   July 3,2018: Started this code
OBJECTIVE
   allocate r array and
   read data from lammps dump file
   returns the number of atoms
    
*/
void alloc_r(double *& rr, float *& wrr, double *& correct, int nAtom, int *& atype);
int read_dump(char *dumpFile, double *& rr, float *& wrr,  double *& correct, float *pbcl, int *& atype,float *cm);
#endif
