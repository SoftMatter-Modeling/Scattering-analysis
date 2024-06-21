/*
 CHANGE LOG
  09.27.2021: based on the definition of G2 
 NOTES:
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "./lib/read_dump_file_CC_UniTenMD.h"
#include "./lib/declare_struct.h"
#include <mpi.h>
#include <openacc.h>


int main(int argc, char** argv){

    //mpi variables
    int me=0,nprocs=1;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&me);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    int pairs_per_mpi_thread = 0;    //Number of simulation frames per mpi thread
    
     // Set rank-gpu affinity
    int local_rank=0,local_size=1;
    MPI_Comm MPI_COMM_SHARE;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &MPI_COMM_SHARE);
    MPI_Comm_rank(MPI_COMM_SHARE,&local_rank);
    MPI_Comm_size(MPI_COMM_SHARE,&local_size);
    
    int ngpus;
    int gpunum;
    int nnodes = nprocs/local_size;

    ngpus = acc_get_num_devices (acc_device_nvidia);
    //printf ("# Number of GPUs per node: (%i) %i\n", me, ngpus);
    
    if (me==0){ 
        printf ("# Number of GPUs per node: %i\n", ngpus);
        printf ("# Number of nodes: %i\n", nnodes);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&ngpus,1,MPI_INT,0,MPI_COMM_WORLD);
    
   gpunum = (me/local_size)*ngpus+local_rank%ngpus;
   printf ("# me: %i , gpunum: %i\n", me, gpunum);
   acc_set_device_num (gpunum, acc_device_nvidia);
   //
    
    
    //timing variables
    clock_t start, end;
    if(me==0) start = clock();
    //end of timing variable declaration

    int natom_per_chain;             // degree-of-polymerization of the chain
    int dumpStep;       // how often are trajectories written to files
    int endStep;        // final trajectory frame that will be considered
    int breakStep;      // time where calculation breaks, useful only for time dependent calculations
    int initStep;       // initial trajectory frame that will be considered
    float dt;           // timestep used in integrating the equation of motion in MD (see https://lammps.sandia.gov/doc/timestep.html)
    
    int nStep;          // number of frames used in the calculation
    int nBreak;
    int nPair;
    
    int u,v,w,i,j,k;        
    
    int nChain=0;         // number of polymer chains in the simulation box
    
    double *rr0;          // !!!! attension: the rr is the unwrapped position in the flow frame. !!!!!!!!!!!!!!!!!!! 
    double *rr1;
    float  *wrr0;         // !!!! attension: the wrr is the unwrapped psoition in the upper triangle lammps frame  !!!!!!!!!!!!!!!!!!!
    float  *wrr1;
    double *correct0;     // the ccumulative correction computed on the fly of the simulation    
    double *correct1;   
    
    int *atype;          // pointer pointing to the ato types, atype={type_1,type_2,...,type_n}
    float pbcl[9];       // halfbox size of the simulation box. Origin or (0.0,0.0,0.0) is the location of the center of the simulation box
    float cm0[3];        // center of mass of all atoms in the simulation box
    float cm1[3];
    int nAtom;           // number of atoms in the simulation box    
    int flag_correct;    // the flag indicating whether takes the correction or not. The flag=1 if corrected required. 
    
    char directory[256];    //directory where the dump files are located  
    char buffer[256],TempC[30];
        
    // Get user input here
    if(me==0){
    while(1){
        fgets(buffer,254,stdin);
        if(feof(stdin)!=0){
            break;
        }
        else{
            if(strstr(buffer,"dumpStep")){
                sscanf(buffer, "%s %i", &TempC,&dumpStep);
                printf("# %s %i\n", TempC,dumpStep);
            }
            if(strstr(buffer,"endStep")){
                sscanf(buffer, "%s %i", &TempC,&endStep);
                printf("# %s %i\n",TempC,endStep);
            }
            if(strstr(buffer,"breakStep")){
                sscanf(buffer, "%s %i", &TempC,&breakStep);
                printf("# %s %i\n",TempC,breakStep);
            }
            if(strstr(buffer,"initStep")){
                sscanf(buffer, "%s %i", &TempC,&initStep);
                printf("# %s %i\n",TempC,initStep);
            }
            if(strstr(buffer,"dt")){
                sscanf(buffer, "%s %f", &TempC,&dt);
                printf("# %s %f\n",TempC,dt);
            }
            if(strstr(buffer,"natom_per_chain")){
                sscanf(buffer, "%s %i", &TempC,&natom_per_chain);
                printf("# %s %i\n",TempC,natom_per_chain);
            }          
            if(strstr(buffer,"flag_correct")){
                sscanf(buffer, "%s %i", &TempC,&flag_correct);
                printf("# %s %i\n",TempC,flag_correct);
            }            
            if(strstr(buffer,"directory")){
                sscanf(buffer, "%s %s", &TempC,&directory);
                printf("# %s %s\n",TempC,directory);
            }
        }
    }
    // end of "Get user input" section
    
    // Do preliminary calculations here
    nStep = (endStep-initStep)/dumpStep+1;
    nBreak = (breakStep-initStep)/dumpStep+1;
    nPair = nStep*(nStep+1)/2-(nStep-nBreak+1)*(nStep-nBreak)/2;
    pairs_per_mpi_thread = (int) ceilf(1.0*nPair/nprocs);
    printf("# nStep %i\n",nStep);
    printf("# nBreak %i\n",nBreak);
    printf("# nPair %i\n",nPair);
    printf("# MPI threads %i\n",nprocs);
    printf("# Pairs/mpi thread %i\n",pairs_per_mpi_thread);
    
    sprintf(buffer,"%sdump.%010i.txt",directory,initStep);
    nAtom=read_dump(buffer,rr0,wrr0,correct0,pbcl,atype,cm0);
    printf("# nAtom: %i\n",nAtom);
    nChain=nAtom/natom_per_chain;
    printf("# nChain: %i\n",nChain);
    } //end me==0 region
    
    
    //broad cast input variables
    MPI_Bcast(&dumpStep,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&endStep,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&breakStep,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&initStep,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&dt,1,MPI_FLOAT,0,MPI_COMM_WORLD);
    MPI_Bcast(&natom_per_chain,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&directory,256,MPI_CHAR,0,MPI_COMM_WORLD);
    
    MPI_Bcast(&nStep,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&nBreak,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&nPair,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&pairs_per_mpi_thread,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&nAtom,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&pbcl,9,MPI_FLOAT,0,MPI_COMM_WORLD);
    MPI_Bcast(&nChain,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&flag_correct,1,MPI_INT,0,MPI_COMM_WORLD);

    //end of broadcasting
    
    //allocate  Sq array for all threads
	int nQ=71;            
	double *q;	
	q=(double *) malloc (nQ*sizeof(double));
   
	for (i=0;i<nQ;i++) q[i]=0.01 * exp(1.0*i/10.0);
	printf("# Done with allocating q \n");
	MPI_Barrier(MPI_COMM_WORLD);

    int    nTheta=720;
    double dTheta=M_PI/nTheta;
    int    total_nq = nQ*nTheta;

    nPair=0;
    double *restrict Fave;
    Fave=(double *) malloc (total_nq*nBreak*sizeof(double));
    #pragma acc enter data create(Fave[0:total_nq*nBreak])
    #pragma acc parallel loop present(Fave)
    for (i=0; i<total_nq*nBreak; i++) Fave[i]=0.0;
   
    //determine exact number of frames per thread
	k=0;
    for(i=0;i<nBreak;i++){
        for(j=0;j<nStep-i;j++){
            int t_id=(int) floorf(1.0*k/pairs_per_mpi_thread);
            k++;
            if(me==t_id)             nPair++;
	    }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
     
    //allocate arrays
    double *restrict pair_traj; 
    int   *restrict pair_delta_t;
    
    if(nPair>0){
       pair_traj     = (double *) malloc (nPair*nAtom*6 * sizeof(double));
       pair_delta_t  = (int *) malloc (nPair*sizeof(int));
       for(i=0;i<nPair*nAtom*6;i++)    pair_traj[i]=0.0;
       for(i=0;i<nPair;i++) pair_delta_t[i]=0;
       }
     if(me==0) printf("# Done allocating trajectory arrays\n");
    
    //populate trajectory arrays
      int ii=0;
      int jj=0;
      int kk=0;
	  int ic=0; 
	  int ib=0; 
      double CCx, CCy, CCz;
	  
      for(i=0;i<nBreak;i++){
          for(j=0;j<nStep-i;j++){
              int t_id=(int) floorf(1.0*kk/pairs_per_mpi_thread);
              kk++;
              if(me==t_id){
                  sprintf(buffer,"%sdump.%010i.txt",directory,initStep+j*dumpStep);
                  read_dump(buffer,rr0,wrr0,correct0, pbcl,atype,cm0);
                  
                  sprintf(buffer,"%sdump.%010i.txt",directory,initStep+(j+i)*dumpStep);
                  read_dump(buffer,rr1,wrr1,correct1,pbcl,atype,cm1);
				  
				  for(ic=0; ic<nChain; ic++){
					  CCx=0.0; 
					  CCy=0.0; 
					  CCz=0.0; 
					  
					  for(ib=0; ib<natom_per_chain; ib++){ // caculate the COM total displacment 
						  ii=ic*natom_per_chain+ib;
                          CCx += (rr1[ii*3+0]-rr0[ii*3+0])/natom_per_chain;
						  CCy += (rr1[ii*3+1]-rr0[ii*3+1])/natom_per_chain;						  
						  CCz += (rr1[ii*3+2]-rr0[ii*3+2])/natom_per_chain; 						  
					  }
					  for(ib=0; ib<natom_per_chain; ib++){ // Do the correction based on the COM flow drift 
					      ii=ic*natom_per_chain+ib;
						  
						  pair_traj[jj*nAtom*6+ii*6+0]=rr0[ii*3+0];
                          pair_traj[jj*nAtom*6+ii*6+3+0]=rr1[ii*3+0]-CCx;
						  
						  pair_traj[jj*nAtom*6+ii*6+1]=rr0[ii*3+1];
                          pair_traj[jj*nAtom*6+ii*6+3+1]=rr1[ii*3+1]-CCy;
                          
                          pair_traj[jj*nAtom*6+ii*6+2]=rr0[ii*3+2];
                          pair_traj[jj*nAtom*6+ii*6+3+2]=rr1[ii*3+2]-CCz;					  
					  }									  
				  }
					  					  
	              pair_delta_t[jj]=i;
		 		  jj++;
                  free(rr0);                  
                  free(rr1);                  
                  free(wrr0);
                  free(wrr1);
                  free(atype);
                  free(correct0);
                  free(correct1);
                  }
              }
          }
        
        MPI_Barrier(MPI_COMM_WORLD);
        if(me==0) printf("# Done populating trajectory arrays\n");

  
 	double s_temp,s_cos,s_sin;
    double qr_0, qr_t, qx, qz;  
/*============================================================================    
calculate s(q,theta,t) here, consdiering the symetrix in the Uniaxial extension 
===============================================================================*/

    #pragma acc data present(Fave[0:total_nq*nBreak]) copyin(pair_traj[0:nPair*6*nAtom],pair_delta_t[0:nPair],q[0:nQ]) 
    {
        for(i=0;i<nPair;i++){
            #pragma acc parallel loop gang collapse(2) 
            for(ii=0;ii<nQ;ii++){
                for(jj=0;jj<nTheta;jj++){
                    
                        qx=q[ii]*sin(jj*dTheta); 
                        qz=q[ii]*cos(jj*dTheta);
                        
                        s_temp=0.0;
                        s_cos=0.0;
                        s_sin=0.0;   
                        #pragma acc loop vector                        
                        for(j=0;j<nAtom;j++){
                                qr_0=0;
                                qr_t=0;
                                                                        
                                qr_0 +=qx*pair_traj[i*nAtom*6+j*6];   
                                qr_0 +=qz*pair_traj[i*nAtom*6+j*6+2];
                                    
                                qr_t +=qx*pair_traj[i*nAtom*6+j*6+3];
                                qr_t +=qz*pair_traj[i*nAtom*6+j*6+3+2];

                                s_cos+=cos(qr_0)*cos(qr_t);
                                s_sin+=sin(qr_0)*sin(qr_t);                             
                            }                         
                            s_temp=s_cos+s_sin;  
                            int index=(ii*nTheta+jj)*nBreak+pair_delta_t[i];
                            #pragma acc atomic
                            Fave[index]+=s_temp/nChain/(nStep-pair_delta_t[i]);  

                    }
                }
            }
        #pragma acc wait            
    }      
    #pragma acc update self(Fave[0:total_nq*nBreak])  
    #pragma acc exit data delete(Fave[0:total_nq*nBreak])
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(me==0){
        printf("# Done calculating Fave(q)\n");
        MPI_Reduce(MPI_IN_PLACE, Fave, total_nq*nBreak ,MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  }
      else{
          MPI_Reduce(Fave, Fave, total_nq*nBreak, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  }
	 MPI_Barrier(MPI_COMM_WORLD);
    
    
    if(me==0){
        
        double  *s0, *s0m2, *s0m4; 
        s0             = (double *) malloc (nQ*nBreak*sizeof(double));        
        s0m2           = (double *) malloc (nQ*nBreak*sizeof(double));        
        s0m4           = (double *) malloc (nQ*nBreak*sizeof(double));
        
        for(i=0;i<nQ*nBreak;i++){
        s0[i]       = 0.0;       
        s0m2[i]     = 0.0;      
        s0m4[i]     = 0.0;
        }
               
        float sm20_const  = sqrt(5.0)/2.0;       
        float sm40_const  = 3.0/64.0;
        
        
        for(i=0;i<nQ;i++){// do the integration  or q vector average  
            for(j=0;j<nTheta;j++){
                     for(int v=0;v<nBreak;v++){
                         int index0  	 = (i*nTheta+j)*nBreak+v;
                         int index1   	 = i*nBreak+v;

                         //(1) Sl,m, = 0,0
                         s0[index1]     += Fave[index0]*sin(j*dTheta)*dTheta;
                        
                         //(3) Sl,m, = 2,0
                         s0m2[index1]   += Fave[index0]*sin(j*dTheta)*dTheta*\
                         sm20_const* (3.0*cos(j*dTheta)*cos(j*dTheta)-1.0);
                                             
                         //(7) Sl,m, = 4,0
                         s0m4[index1]  += Fave[index0]*sin(j*dTheta)*dTheta*\
                         sm40_const*(9.0+20.0*cos(2.0*j*dTheta) + 35.0*cos(4.0*j*dTheta));

                     }                    
            }
        }
        
        
        FILE *fpw0, *fpw2, *fpw6; // out put the file       
        fpw0=fopen("s0.txt","w");
        fpw2=fopen("s0m2.txt","w");
        fpw6=fopen("s0m4.txt","w");
        
        for(i=0;i<nQ;i++){
             fprintf(fpw0,"%f  ",q[i]);
             fprintf(fpw2,"%f  ",q[i]);
             fprintf(fpw6,"%f  ",q[i]);
            for(j=0;j<nBreak;j++){
                int index=i*nBreak+j;
                fprintf(fpw0,"   %0.12f",0.5*s0[index]/natom_per_chain); 
                fprintf(fpw2,"   %0.12f",0.5*s0m2[index]/natom_per_chain); 
                fprintf(fpw6,"   %0.12f",0.5*s0m4[index]/natom_per_chain); 

            } 
            fprintf(fpw0," \n"); 
            fprintf(fpw2," \n"); 
            fprintf(fpw6," \n"); 
            
            
        }        
        fclose(fpw0); 
        fclose(fpw2); 
        fclose(fpw6); 
        
        //====output the inter-mediate resluts just in case                
        fpw0=fopen("S_q_theta_phi_t.txt","w");
        for(i=0;i<total_nq*nBreak;i++) fprintf(fpw0,"%0.12f \n",Fave[i]);
        fclose(fpw0); 
                  
    }           

   if(me==0){
    end = clock();
    double total_time = ((double)(end-start))/CLOCKS_PER_SEC;
    printf("# total_time %f\n", total_time);
   }
   
   MPI_Finalize();
}
