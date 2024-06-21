/*
CHANGE LOG
   July 3,2018: Started this code
OBJECTIVE
   allocate r array and
   read data from lammps dump file
   returns the number of atoms 
   
   This sub file is especially for reading the MD trajectory of shear flow 
    
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

void alloc_r(double *& rr, float *& wrr, double *& correct, int nAtom, int *& atype){
    int i;
    rr= (double *) malloc (nAtom*3 * sizeof(double));
    for(i=0;i<nAtom*3;i++) rr[i]=0;
    
    wrr= (float *) malloc (nAtom*3 * sizeof(float));
    for(i=0;i<nAtom*3;i++) wrr[i]=0; //wrr is the wrapped coordinates
    
    correct= (double*) malloc (nAtom*3*sizeof(double));
    for(i=0;i<nAtom*3;i++) correct[i]=0; 
    
    
    atype=(int *) malloc (nAtom * sizeof(int));
}
int read_dump(char *dumpFile,double *& rr, float *& wrr, double *& correct, float *pbcl, int *& atype,float *cm){
    
    FILE *fpr0;
    int nAtom;
    int i,k;
    char buffer[256];

    fpr0=fopen(dumpFile,"r");
    for(i=0;i<3;i++) cm[i]=0.0;
    
    while(1){
        fgets(buffer,254,fpr0);
        if(feof(fpr0)!=0){
            break;
        }
        else{
            for(i=0;i<3;i++) fgets(buffer,254,fpr0);
            sscanf(buffer,"%i",&nAtom);
            
            alloc_r(rr,wrr,correct, nAtom,atype);   
            
            fgets(buffer,254,fpr0);
            for(i=0;i<3;i++){ 
                fgets(buffer,254,fpr0);
                sscanf(buffer,"%f %f %f",&pbcl[i*3+0],&pbcl[i*3+1],&pbcl[i*3+2]);
            }
            fgets(buffer,254,fpr0);
            
            
            for(i=0;i<nAtom;i++){
                fgets(buffer,254,fpr0);
                int ID,at;
                float p[3];
                float pinit0;
                
                double cc[3] ,r[3];
                sscanf(buffer,"%i %i %f %f %f %lf %lf %lf %lf %lf %lf %f",&ID,&at,&p[0],&p[1],&p[2],&r[0],&r[1],&r[2],&cc[0],&cc[1],&cc[2],&pinit0);

                for(k=0;k<3;k++){
                    wrr[(ID-1)*3+k]     = p[k];                 
                    rr[(ID-1)*3+k]      = r[k];
                    correct[(ID-1)*3+k] = cc[k];
                    cm[k] += r[k];
                }
                atype[(ID-1)]   = at;
                
                
            }
            for(k=0;k<3;k++) cm[k]=cm[k]/(1.0*nAtom);
            
        }
    }
    fclose(fpr0);
    return nAtom;
}

