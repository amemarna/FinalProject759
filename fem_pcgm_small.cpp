#include<mpi.h>
#include<math.h>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<stdio.h>

#include "MyVonMises.cpp"


using namespace std;



double invmatrix(double *A,double *Ainv,int n);
void invmatrix2(double *A,double *Ainv,int n);
void LUdecomp(double *A,double *L,double *U,double *P,double *I,int n);
void SolveSysOEqs(double* L,double* U,double* P,double* b,double* x, int n);
void checkmat(double* A, int n);
void matvecmult(double* a, double* b, double* c, int row, int col);
void matrixmult(double* a, double* b, double* c, int row, int inner, int col);
void Transpose(double *a,double *at,int row, int col);
void InitializeArray(double* A, int row, int col);
void Emodsix2three(double* De, double* Ce3);
void MyVonMises(double* sig_el,double* kap_el,double* eps_el,double* Etng_el,double* oldsig_el,double* oldkap_el,double* oldeps_el,double* Deps,double E, double nu,double H,double sy);
void AddArray(double* A, double* B, double* C, int row, int col, int op, int id);
void ApplyBCs(double* A,double* Af,int ndofs,int nfr,int* frdofs,int id);
double CompResNorm(double* A,int n);
void GlbAsmbly(double* Kel,double* fel,double* Kg,double* fg,int* dofs_el,int ndofs);
void dotproduct(double* a,double* b,double &c,int n);
void dotproduct_par(double* a,double* b,double &totalsum,int nfr,int comm_sz, int my_rank);



void matvecmult_par(double* A,double* x,double* y,int nfr, int comm_sz, int my_rank);  
void Mat_vect_mult(double* local_A, double* local_x, double* local_y,double* x,int n ,int local_n);





void readmatrix(FILE* fp, double* in, const int nrows, const int ncols)
{
  for (int i=0;i<nrows;i++){
    for (int j=0;j<ncols;j++){
      double r=fscanf(fp,"%lf",&in[i*ncols+j]);

    }
  }
  fclose(fp);
}


void Write2File(double* A, int r, int c, const char* fp)
{
  std::ofstream ofile(fp);
	for(unsigned int i = 0; i < r; i++){
	  for(unsigned int j = 0; j < c; j++){
	    ofile<<A[i*c+j]<<" ";

	  }
	  cout<<'\n';
	}
	ofile.close();

}


void hooke( const double E, const double nu)
{

  double fc = E/((1+nu)*(1-2*nu));
  double D_e[6][6] = {{fc*(1-nu),   fc*nu,       fc*nu,        0,            0,                0},
		     {fc*nu,       fc*(1-nu),   fc*nu,        0,            0,                0},
		     {fc*nu,       fc*nu,       fc*(1-nu),    0,            0,                0},
		     {   0,        0,           0,       fc*(1-2*nu)/2,     0,                0},
		     {   0,        0,           0,            0,      fc*(1-2*nu)/2,          0},
		     {   0,        0,           0,            0,            0,         fc*(1-2*nu)/2} };
  
}


void coordxtr(double* edof,double* coords, double* ex, double* ey, int nels)
{

  int nod;
  
  for (int i=0;i<nels;i++){
    for (int j=0;j<3;j++){

      nod = edof[i*7+(2*j+2)] / 2;
      ex[i*3+j]=coords[(nod-1)*2+0];
      ey[i*3+j]=coords[(nod-1)*2+1];
    }
  }

}


double prescribedu(double t, double maxdisp)
{
  return maxdisp*t;
}


void checkmat(double* A, int r, int c)
{
  for(int i=0;i<r;i++){
    for(int j=0;j<c;j++){


      A[i*c+j];

       }

  }

}






//int main()
int main(int argc,char** argv)  
{  


  int my_rank,comm_sz;

        MPI_Init(&argc, &argv);  
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  int bcload[11][2] ={{20, 1},
		      {22, 1},
		      {34, 1},
		      {36, 1},
		      {38, 1},
		      {40, 1},
		      {42, 1},
		      {44, 1},
		      {46, 1},
		      {48, 1},
		      {50, 1}};

  int bc0[24][2]={{23, 0},
             {17, 0}, 
	     {19, 0},
	     {25, 0},
	     {27, 0},
	     {29, 0},
	     {31, 0},
	     {16, 0},
	     {24, 0},
	     {70, 0},
	     {72, 0},
	     {74, 0},
	     {76, 0},
	     {20, 1},
	     {22, 1},
	     {34, 1},
	     {36, 1},
	     {38, 1},
	     {40, 1},
	     {42, 1},
	     {44, 1},
	     {46, 1},
	     {48, 1},
	     {50, 1}};




   
     
  FILE *fpe = fopen("edofdata.inp","r");
  FILE *fpc = fopen("coordata0.inp","r");
  FILE *fph = fopen("hooke0.inp","r");
  FILE *fpg = fopen("gamma.inp","r");
  
  const char* fpew = "edof.out";
  
  int nels=176;
  int neldofs=6;
  int nods=108;
  double* edofs  = (double*) malloc(nels*(neldofs+1)*sizeof(double));
  double* coords = (double*) malloc(2*nods*sizeof(double));
  double* De     = (double*) malloc(6*6*sizeof(double));
  double* gamma  = (double*) malloc(9*6*sizeof(double));
  readmatrix(fpe,edofs,nels,neldofs+1);
  
  readmatrix(fpc,coords,nods,2);
  checkmat(coords,nods,2);
  readmatrix(fph,De,6,6);
  
  readmatrix(fpg,gamma,9,6);
  
  
  double th    = 0.001;   // % plate thickness
  
  int ndofs = 216;
  
//% Analysis parameters
double NrOfTsteps = 1000.0; // % total number of timesteps
//double NrOfTsteps = 10.0; // % total number of timesteps 
double reltol     = 1e-6; //% tolerance for equilibrium iterations
int ptype      = 1;   // % 1 = plane stress 
int magn       = 10;  // % magnification factor for deformed mesh
int maxiter    = 200;  // % max number of iterations
int iter;
 
//% Material parameters for elastic plastic
double E          = 200*pow(10,9); 
double nu         = 0.3000000000000001;

double H          = 0.2*E;
double sy         = E/(1e2);


 double* ex = (double*) malloc(3*nels*sizeof(double));
 double* ey = (double*) malloc(3*nels*sizeof(double));

coordxtr(edofs,coords,ex,ey,nels);
checkmat(ex,nels,3);
checkmat(ey,nels,3); 
free(coords);

double u_el = 0.001;
double maxdisp = 3.0*u_el;



double* u   =(double*) malloc(ndofs*sizeof(double));

double* oldsig =(double*) malloc(9*nels*sizeof(double));
double* oldeps =(double*) malloc(9*nels*sizeof(double));

double* oldkap = (double*) malloc(nels*sizeof(double));
double* kap = (double*) malloc(nels*sizeof(double));
double* u_hist    =(double*) malloc(NrOfTsteps*sizeof(double)); 
double* FReaction =(double*) malloc(NrOfTsteps*sizeof(double));

double* eldisps=(double*)malloc(nels*(neldofs+1)*sizeof(double));

 InitializeArray(oldkap,1,nels);
 

//% Global timestepping
double olddisp = 0.0;
double* du = (double*) malloc(ndofs*sizeof(double));

 double dispcontrolled = 0.0;
 int nbcld=11;
 int nbc = 24;


      int* cstd =(int*) malloc(nbc*sizeof(int));
      int nfr0 = ndofs-nbc;
      int* frdofs =(int*) malloc(nfr0*sizeof(int));
      int ii = 0;
      int jj = 0;
      int nfr = 0;      
      for(int i=0;i<ndofs;i++){
	int cp = 0;
	for(int j=0;j<nbc;j++){
		if(i==bc0[j][0]-1){
		cstd[jj]=bc0[j][0]-1;
		jj++;
		cp=1;}
		
	  }
	  if(cp==0){
	    frdofs[ii]=i;
		  ii++;
		  nfr++;}
      }



	  
      	int ngp = 1;
	double gp1[1] = {0.333333333333333333};
	double gp2[1] = {0.333333333333333333};
	double gw[1]  = {0.50};

	double* fh = (double*) malloc(ndofs*sizeof(double));   
	InitializeArray(FReaction,1,NrOfTsteps);
	
	double err,prp;
	double Tstepc=-1.0;
	double t0 = MPI_Wtime();

 for (int Tstep=0;Tstep<NrOfTsteps;Tstep++){
   Tstepc++;
   if(Tstep>0){
     prp = (Tstepc)/NrOfTsteps + 1e-17;}
   else{
     prp = 0.0;}
    dispcontrolled = prescribedu(prp,maxdisp); 


    for (int m=0;m<ndofs;m++){
      for (int n=0;n<nbcld;n++){
	if(m==bcload[n][0]-1){
	  du[m]=dispcontrolled - olddisp;}
      }
    }

    
    for(int i=0;i<nels;i++){
	       eldisps[i*7+0]=edofs[i*7+0];
    }

    if(my_rank==0)
        printf("Timestep= %d Displacement= %lf \n",Tstep,dispcontrolled);

    

    double* sig =(double*) malloc(9*nels*sizeof(double));
    double* eps =(double*) malloc(9*nels*sizeof(double));
    double* kap = (double*) malloc(1*nels*sizeof(double));
    
    
    for(iter=0; iter<maxiter; iter++){

      
      for (int k=0;k<ndofs;k++){
	  for(int i=0;i<nels;i++){
	    for(int j=0;j<7;j++){
	      if( edofs[i*7+j]-1==k && j!=0){

		      eldisps[i*7+j]=du[k];}
	    }
	  }
      }
      


  double* fg = (double*) malloc(ndofs*sizeof(double));
  double* Kg = (double*) malloc(ndofs*ndofs*sizeof(double));          
      InitializeArray(Kg,ndofs,ndofs);
      InitializeArray(fg,1,ndofs);
      
      for(int iel=0; iel<nels; iel++){


	double* oldsig_el = (double*) malloc(9*sizeof(double));
	double* oldeps_el = (double*) malloc(9*sizeof(double));
	for(int j=0; j<9; j++){
	oldsig_el[j] = oldsig[iel*9+j];
	oldeps_el[j] = oldeps[iel*9+j];
	}


	double* oldkap_el = (double*) malloc(1*sizeof(double));
	oldkap_el[0]=oldkap[iel];
	
       	  
	double* disp_el = (double*) malloc(neldofs*sizeof(double));
	for(int j=1;j<7;j++){
	  disp_el[j-1]=eldisps[iel*7+j];}

	  
	double* elcrds =(double*) malloc(3*2*sizeof(double));
	for(int i=0;i<3;i++){
	 
	    elcrds[i*2+0]=ex[iel*3+i];
	    elcrds[i*2+1]=ey[iel*3+i];
	}


	double ksi,etha;
	double N1,N2,N3;
	double dN1dksi, dN1detha, dN2dksi, dN2detha, dN3dksi, dN3detha;
	double* dNdksi_v =(double*) malloc(6*sizeof(double));
	double* dNdksi =(double*) malloc(2*3*sizeof(double));
	double* Jmat  =(double*) malloc(2*2*sizeof(double));
	double* Jinv  =(double*) malloc(2*2*sizeof(double));
	double* dNdx  =(double*) malloc(2*3*sizeof(double));
	double* Bmat  =(double*) malloc(3*6*sizeof(double));
	double* Bmat_v  =(double*) malloc(18*sizeof(double));

	double* Kel = (double*) malloc(6*6*sizeof(double));
        double* fel = (double*) malloc(6*sizeof(double));

	InitializeArray(Kel,6,6);
	InitializeArray(fel,1,6);
	
	double* sig_el   = (double*) malloc(9*sizeof(double));
	double* eps_el   = (double*) malloc(9*sizeof(double));

	double* kap_el = (double*) malloc(1*sizeof(double));;

	for (int igp=0;igp<ngp;igp++){

	  
	  ksi=gp1[igp];
	  etha=gp2[igp];

	  N1=ksi;
	  N2=etha;
	  N3=1-ksi-etha;

	  dN1dksi=1;
	  dN1detha=0;
	  dN2dksi=0;
	  dN2detha=1;
	  dN3dksi=-1;
	  dN3detha=-1;


	  dNdksi_v[0]=dN1dksi;dNdksi_v[1]=dN2dksi;dNdksi_v[2]=dN3dksi;
	  dNdksi_v[3]=dN1detha;dNdksi_v[4]=dN2detha;dNdksi_v[5]=dN3detha;

	  int cnt=0;
	  for(int i=0;i<2;i++){
	    for(int j=0;j<3;j++){
	      dNdksi[i*3+j]=dNdksi_v[cnt];
	      cnt+=1;
	    }
	  }

	  matrixmult(dNdksi,elcrds,Jmat,2,3,2);


	  invmatrix2(Jmat,Jinv,2);
	  

	  matrixmult(Jinv,dNdksi,dNdx,2,2,3);

	  
 	for(int i=0;i<3;i++){
        for(int j=0;j<6;j++){
	Bmat[i*6+j]=0.0;
	}
	}


	Bmat[0]=dNdx[0];Bmat[2]=dNdx[1];Bmat[4]=dNdx[2];
	Bmat[7]=dNdx[3];Bmat[9]=dNdx[4];Bmat[11]=dNdx[5];
	Bmat[12]=dNdx[3];Bmat[14]=dNdx[4];Bmat[16]=dNdx[5];
	Bmat[13]=dNdx[0];Bmat[15]=dNdx[1];Bmat[17]=dNdx[2];
	  
	
	  double* Se = (double*) malloc(6*6*sizeof(double));
	
	  invmatrix(De,Se,6);
	  
	  
	  
	  double* Ce3 = (double*) malloc(3*3*sizeof(double));
	  
	  Emodsix2three(De,Ce3);
	  
	  
	  double* Deps3= (double*) malloc(3*sizeof(double));
	  matvecmult(Bmat,disp_el,Deps3,3,6);

	  
	  double* t_sig= (double*) malloc(6*sizeof(double));
          for(int j=0;j<6;j++){
	  t_sig[j]=0;
	  }

	  double* t_sig3= (double*) malloc(3*sizeof(double));
	  
	  
	  matvecmult(Ce3,Deps3,t_sig3,3,3);

	  
	  cnt=-1;
	  for(int j=0;j<6;j++){
		 if(j==0 || j==1 || j==3){
		   cnt+=1;
		   t_sig[j]=t_sig3[cnt];
		 }
	  }

	  
	  double* Depsv= (double*) malloc(6*sizeof(double));	  
	  matvecmult(Se,t_sig,Depsv,6,6);

	  
	  double* Deps= (double*) malloc(9*sizeof(double));	  
	  matvecmult(gamma,Depsv,Deps,9,6);

	  
	  double* Etng_el = (double*) malloc(9*9*sizeof(double));
		  
	  MyVonMises(sig_el,kap_el,eps_el,Etng_el,oldsig_el,oldkap_el,oldeps_el,Deps,E,nu,H,sy);

	  
	  double* Etngr = (double*) malloc(6*6*sizeof(double));
	  
	  for(int i=0;i<6;i++){
	    for(int j=0;j<6;j++){
	      if(i>2 && j>2){
	  	Etngr[i*6+j] = Etng_el[i*9+j]/2.0;}
	    else {
	      Etngr[i*6+j] = Etng_el[i*9+j];}
	    }
	  }

	  
	  double* BmatT = (double*) malloc(6*3*sizeof(double));
	  Transpose(Bmat,BmatT,3,6);

	  
	  double* CE3 = (double*) malloc(3*3*sizeof(double));
	  Emodsix2three(Etngr,CE3);


	  
	  double detJ=Jmat[0]*Jmat[3]-Jmat[1]*Jmat[2];

	  
	  double* BTC = (double*) malloc(6*3*sizeof(double));
	  matrixmult(BmatT,CE3,BTC,6,3,3);

	  
	  double* BTCB = (double*) malloc(6*6*sizeof(double));
	  matrixmult(BTC,Bmat,BTCB,6,3,6);

	  
	  for(int i=0;i<6;i++){
	    for(int j=0;j<6;j++){
	      Kel[i*6+j] = Kel[i*6+j] + BTCB[i*6+j]*detJ*gw[igp]*th;
	    }
	  }

	  const char* fpkw = "Kel.out";
	  
     	  double* sig_el3 = (double*) malloc(3*sizeof(double));
	  InitializeArray(sig_el3,1,3);
	  cnt=0;
	  for(int j=0;j<9;j++){
		 if(j==0 || j==1 || j==3){
		   sig_el3[cnt]=sig_el[j];
		   cnt+=1;
		 }
	  }


	  double* BTsig = (double*) malloc(6*sizeof(double));
	  InitializeArray(BTsig,1,6);
	  matvecmult(BmatT,sig_el3,BTsig,6,3);

	  
	  for(int i=0;i<6;i++){
	      fel[i] = fel[i] + BTsig[i]*detJ*gw[igp]*th;
	  }

 	free(BTsig);free(sig_el3);
	free(BTCB);free(BTC);free(CE3);free(BmatT);
	free(Etngr);free(Etng_el);free(Deps);free(Depsv);
	free(t_sig3);free(t_sig);free(Deps3);free(Ce3);free(Se);//free(Se3);
	free(oldsig_el);free(oldeps_el);free(disp_el);free(elcrds);
	free(dNdksi_v);free(dNdksi);free(Jmat);free(Jinv);free(dNdx);free(Bmat);free(Bmat_v); 
	}
	//}   //end of gp integration

	int* dofs_el = (int*) malloc(6*sizeof(int));
	int cnt=0;
	for (int j=1;j<7;j++){
	  dofs_el[cnt]=edofs[iel*7+j];
	  cnt++;
	}

	
	GlbAsmbly(Kel,fel,Kg,fg,dofs_el,ndofs);


	  const char* kgw = "Kg.dat";

	  
	

	for(int j=0; j<9; j++){
	sig[iel*9+j] = sig_el[j];
	eps[iel*9+j] = eps_el[j];
	}

	kap[iel]=kap_el[0];

	free(dofs_el);
	free(sig_el);free(eps_el);free(kap_el);
	free(Kel);free(fel);
      }  //end of element loop

      double* fext = (double*) malloc(ndofs*sizeof(double));
      InitializeArray(fext,1,ndofs);
      double* res  = (double*) malloc(ndofs*sizeof(double));
      InitializeArray(res,1,ndofs);
      double* resf = (double*) malloc(nfr*sizeof(double));
      InitializeArray(resf,1,nfr);
      
      
	  AddArray(fext,fg,res,ndofs,1,1,0);
      
      
      ApplyBCs(res,resf,ndofs,nfr,frdofs,0);

      
      err = CompResNorm(resf,nfr);
      double tol = 1e-6;
      double tv1 = *std::max_element(fg,fg+ndofs);
      double tv2 = *std::min_element(fg,fg+ndofs);

      if (err<=tol*max(fabs(tv1),fabs(tv2))){
	if(my_rank==0)
	printf("Tstep= %d iter= %d err= %lf \n",Tstep,iter,err);
	break;}
      if(my_rank==0)
      printf("Tstep= %d iter= %d err= %lf \n",Tstep,iter,err);


      double* Kgf = (double*) malloc(nfr*nfr*sizeof(double));
      ApplyBCs(Kg,Kgf,ndofs,nfr,frdofs,1);

      const char* kgfw = "Kgf.dat";


      double tcg0 = MPI_Wtime();
      //PCG method
      double* ro = (double*) malloc(nfr*sizeof(double));
      double* xv = (double*) malloc(nfr*sizeof(double));
      for(int i=0;i<nfr;i++){
	xv[i]=resf[i];}
      double* Axo = (double*) malloc(nfr*sizeof(double));

      matvecmult_par(Kgf,xv,Axo,nfr,comm_sz,my_rank);
      for(int i=0;i<nfr;i++){
	ro[i]=resf[i]-Axo[i];}

      double* LT = (double*) malloc(nfr*nfr*sizeof(double));
      for(int i=0;i<nfr;i++){
      for(int j=0;j<nfr;j++){
	if(i>=j){
	  LT[i*nfr+j]=Kgf[i*nfr+j];}
	else{
	  LT[i*nfr+j]=0.0;}
      }
      }

      double* LTT = (double*) malloc(nfr*nfr*sizeof(double));
      Transpose(LT,LTT,nfr,nfr);

      double* D = (double*) malloc(nfr*nfr*sizeof(double));
      double* Dinv = (double*) malloc(nfr*nfr*sizeof(double));
      double* Dinv2 = (double*) malloc(nfr*nfr*sizeof(double));
      double* In = (double*) malloc(nfr*nfr*sizeof(double));      
      for(int i=0;i<nfr;i++){
      for(int j=0;j<nfr;j++){
	if(i==j){
	  D[i*nfr+j]=Kgf[i*nfr+j];
	  Dinv[i*nfr+j]=1.0/Kgf[i*nfr+j];
	  Dinv2[i*nfr+j]=1.0/sqrt(Kgf[i*nfr+j]);
	  In[i*nfr+j]=1.0;}
	else{
	  D[i*nfr+j]=0.0;
	  Dinv[i*nfr+j]=0.0;
	  Dinv2[i*nfr+j]=0.0;
	  In[i*nfr+j]=0.0;}
      }
      }
      

      double* LTDinv = (double*) malloc(nfr*nfr*sizeof(double));
      matrixmult(LT,Dinv,LTDinv,nfr,nfr,nfr);
      double* ImLTDinv = (double*) malloc(nfr*nfr*sizeof(double));
      AddArray(In,LTDinv,ImLTDinv,nfr,nfr,1,1);
      double* Kb = (double*) malloc(nfr*nfr*sizeof(double));
      matrixmult(Dinv2,ImLTDinv,Kb,nfr,nfr,nfr);
      double* KbT = (double*) malloc(nfr*nfr*sizeof(double));
      Transpose(Kb,KbT,nfr,nfr);
      double* Mb = (double*) malloc(nfr*nfr*sizeof(double));


      for(int i=0;i<nfr;i++){
      for(int j=0;j<nfr;j++){
	  Mb[i*nfr+j]=Dinv[i*nfr+j];
      }
      }

      double* zo = (double*) malloc(nfr*sizeof(double));

      matvecmult_par(Mb,ro,zo,nfr,comm_sz,my_rank);
      
      double* pv = (double*) malloc(nfr*sizeof(double));
      double* ri = (double*) malloc(nfr*sizeof(double));
      for(int i=0;i<nfr;i++){
	pv[i]=zo[i];
        ri[i]=ro[i];}

      int cgit;
      double ierr_cg = CompResNorm(ri,nfr);
      double tolcg = 1e-10;
      double alpha,beta,err_cg;
      for(cgit=0;cgit<200;cgit++){
	double dot_rozo=0;

	dotproduct_par(ro,zo,dot_rozo,nfr,comm_sz,my_rank);
	double* Ap = (double*) malloc(nfr*sizeof(double));

	matvecmult_par(Kgf,pv,Ap,nfr,comm_sz,my_rank);
	
	double dot_App=0;

	dotproduct_par(Ap,pv,dot_App,nfr,comm_sz,my_rank);
	alpha = dot_rozo/dot_App;

       	double* rv = (double*) malloc(nfr*sizeof(double));
	for(int i=0;i<nfr;i++){
	  xv[i] = xv[i] + alpha*pv[i];
	  rv[i] = ro[i] - alpha*Ap[i];}

	double* zv = (double*) malloc(nfr*sizeof(double));

	matvecmult_par(Mb,rv,zv,nfr,comm_sz,my_rank);
	
	double dot_rz=0;

	dotproduct_par(rv,zv,dot_rz,nfr,comm_sz,my_rank);
       	beta = dot_rz/dot_rozo;

	for(int i=0;i<nfr;i++){
	  pv[i] = zv[i] + beta*pv[i];}
	
        err_cg = CompResNorm(rv,nfr);
	if(err_cg<=tolcg*ierr_cg){
	  break;}

	for(int i=0;i<nfr;i++){
	ro[i]=rv[i];
        zo[i]=zv[i];}

	free(rv);
	free(zv);
	free(Ap);
      }
      double tcg1 = MPI_Wtime();
      if(my_rank==0)
        printf("cg_it= %d ;cg_err= %lf ;time2conv= %lf \n",cgit,err_cg,tcg1-tcg0);

      free(ro);free(Axo);free(LT);free(LTT);
      free(D);free(Dinv);free(Dinv2);free(In);
      free(LTDinv);free(ImLTDinv);free(Kb);free(KbT);free(Mb);
      free(zo);free(pv);free(ri);


      double* ddu = (double*) malloc(nfr*sizeof(double));
	
	for(int i=0;i<nfr;i++){
        ddu[i]=xv[i];}

        free(xv);


	  int pos;
      for(int i=0;i<nfr;i++){
	pos  = frdofs[i];
	du[pos] = du[pos] + ddu[i];
      }

	  for (int i=0;i<ndofs;i++){
            fh[i]=fg[i];}

      free(Kgf);free(Kg);free(fg);
      free(ddu); free(res);free(resf);free(fext);
    }//end of iteration loop

    AddArray(u,du,u,ndofs,1,0,0);

    for(int i=0;i<nels;i++){
      for(int j=0; j<9; j++){
	oldsig[i*9+j] = sig[i*9+j];
	oldeps[i*9+j] = eps[i*9+j];
      }
    }

    for(int i=0;i<nels;i++){
      oldkap[i]=kap[i];
    }

    free(sig);
    free(eps);free(kap);
    olddisp = dispcontrolled;
    u_hist[Tstep]=dispcontrolled;
    int pos=0;
    for(int i=0;i<nbcld;i++){
      pos = bcload[i][0];
      FReaction[Tstep] = FReaction[Tstep] + 2*fh[pos];
    }
    
 }//end of loadsteps
 double t1 = MPI_Wtime();
 free(fh);
 //printf("%lf\n%lf\n%lf\n%lf\n",Tstepc,dispcontrolled,FReaction[99],t1-t0);
 if(my_rank==0)
 printf("total_disp= %lf\n total_time= %lf\n",dispcontrolled,t1-t0);

 MPI_Finalize();
 return 0;

}


void dotproduct(double* a,double* b,double &c,int n)
{
	
	for(int i=0;i<n;i++){
	  c  = c + a[i]*b[i];}

}


void dotproduct_par(double* a,double* b,double &totalsum,int nfr,int comm_sz, int my_rank)
{
  
   MPI_Bcast(&nfr, 1, MPI_INT, 0,MPI_COMM_WORLD);

   int loc_n = nfr/comm_sz;

   double* loc_a   = (double*) malloc(loc_n*sizeof(double));
   double* loc_b   = (double*) malloc(loc_n*sizeof(double));

   if(my_rank==0){
   MPI_Scatter(a, loc_n, MPI_DOUBLE, loc_a, loc_n, MPI_DOUBLE, 0,MPI_COMM_WORLD);}
   else{
     MPI_Scatter(a, loc_n, MPI_DOUBLE, loc_a, loc_n, MPI_DOUBLE, 0,MPI_COMM_WORLD);}

   if(my_rank==0){
     MPI_Scatter(b, loc_n, MPI_DOUBLE,loc_b, loc_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);}
   else{
     MPI_Scatter(b, loc_n, MPI_DOUBLE,loc_b, loc_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);}

   MPI_Barrier(MPI_COMM_WORLD);


  double loc_sum=0.0;
	for(int i=0;i<loc_n;i++){
	  loc_sum  = loc_sum + loc_a[i]*loc_b[i];}

	MPI_Allreduce(&loc_sum,&totalsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	free(loc_a);
	free(loc_b);
}

void GlbAsmbly(double* Kel,double* fel,double* Kg,double* fg,int* dofs_el,int ndofs)
{

  int ipos, jpos;
  for(int i=0;i<6;i++){
    for(int j=0;j<6;j++){
      ipos = dofs_el[i]-1;
      jpos = dofs_el[j]-1;
      Kg[ipos      *ndofs+      jpos] = Kg[ipos      *ndofs+      jpos] + Kel[i*6+j];

    }
  }


  for(int i=0;i<6;i++){
    ipos = dofs_el[i]-1;

	      fg[ipos] = fg[ipos] + fel[i];
  }
  
}



void InitializeArray(double* A, int row, int col)
{
  
 	for(int i=0;i<row;i++){
	  for(int j=0;j<col;j++){
	    A[i*col+j]=0.0;
	  }
	}

}


double CompResNorm(double* A,int n)
{

      double sum=0.0;
      for (int i=0;i<n;i++){
      sum+=A[i]*A[i];
      }

      return sqrt(sum);
}


void AddArray(double* A, double* B, double* C, int row, int col, int op, int id)
{


    if( id==0){ //id=0 --> vector
      	for(int i=0;i<row;i++){
	    if(op==0){
	      C[i]=A[i]+B[i];}
	    else{
	    C[i]=A[i]-B[i];}
	}
	
    }
    else{
 	for(int i=0;i<row;i++){
	  for(int j=0;j<col;j++){
	    if(op==0){
	      C[i*col+j]=A[i*col+j]+B[i*col+j];}
	    else{
	    C[i*col+j]=A[i*col+j]-B[i*col+j];}
	  }
	}
    }
}


void ApplyBCs(double* A,double* Af,int ndofs,int nfr,int* frdofs,int id)
{

  if( id==0){ //id=0 --> vector

    int jpos;
    for(int j=0;j<nfr;j++){
            jpos = frdofs[j]; 
	    Af[j] = A[jpos];
    }
    

  }
  else {
    int row = -1;
    int col;
    for(int i=0;i<ndofs;i++){
      for(int k=0;k<nfr;k++){
	if(i==frdofs[k]){
	  row++;
	  col=-1;
	for(int j=0;j<ndofs;j++){
	  for(int l=0;l<nfr;l++){
	    if(j==frdofs[l]){
	      col++;
	      Af[row*nfr+col]=A[i*ndofs+j];
	    }
	  }
	}
	}
      }
    }
  }
  
}



void Emodsix2three(double* A6, double* B3)
{  
	  int row=-1;
	  int col=-1;

	  double* A_inv6 = (double*) malloc(6*6*sizeof(double));
	  double* A_inv3 = (double*) malloc(3*3*sizeof(double));

	  invmatrix(A6,A_inv6,6);


          checkmat(A_inv6,6);
	  
	  for(int i=0;i<6;i++){
	    if(i<4 && i!=2){
	      row+=1;
	    for(int j=0;j<6;j++){
		 if(j<4 && j!=2){
		   col+=1;
		   A_inv3[row*3+col]=A_inv6[i*6+j];
		 }
	    }
	    }
	    col=-1;
	  }

          checkmat(A_inv3,3);
	  

	  invmatrix(A_inv3,B3,3);

	  free(A_inv6);
	  free(A_inv3);

}


void matrixmult(double* a, double* b, double* c, int row, int inner, int col){

  int i,j,k;
  double sum;

  for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			sum = 0;
			for (k = 0; k < inner; k++) {

				sum = sum + a[i*inner + k] * b[k*col + j];

			}
			c[i*col + j] = sum;
		}

  }

}


void matvecmult(double* a, double* b, double* c, int row, int col)
{


  int i,j;
  for (i = 0; i < row; i++) {
     double sum = 0;
                for (j = 0; j < col; j++) {
                                                                                                        
                                sum = sum + a[i*col + j] * b[j];
                                                                                                        
                }
c[i] = sum;
  }

}





void invmatrix2(double *A,double *Ainv2,int n)
{
  double detA = A[0]*A[3] - A[1]*A[2];

      Ainv2[0] = A[3]/detA;
      Ainv2[1] = A[1]*(-1.0)/detA;
      Ainv2[2] = A[2]*(-1.0)/detA;
      Ainv2[3] = A[0]/detA;

}


double invmatrix(double *A,double *Ainv,int n)
{


  double* LL = (double*) malloc(n*n*sizeof(double));
  double* UU = (double*) malloc(n*n*sizeof(double));
  double* II = (double*) malloc(n*n*sizeof(double));
  double* PP = (double*) malloc(n*n*sizeof(double));
  
  double* Ainv0 = (double*) malloc(n*n*sizeof(double));

  LUdecomp(A,LL,UU,PP,II,n);


          checkmat(II,n);

          
          double* I_col = (double*) malloc(n*sizeof(double));

	  
          for (int i=0;i<n;i++){

            int jj=0;
            for(int j=i;j<n*n;j+=n){
              I_col[jj]=II[j];

              jj+=1;
            }


	    double* A_col = (double*) malloc(n*sizeof(double));

            SolveSysOEqs(LL,UU,PP,I_col,A_col,n);


            jj=0;
            for(int j=i;j<n*n;j+=n){
              Ainv[j]=A_col[jj];
              jj++;
            }
            free(A_col);

          }

          checkmat(Ainv,n);

	  free(Ainv0);
	  free(I_col);
	  free(LL);
	  free(UU);
	  free(PP);
	  free(II);
	  
          return 0;
}



void LUdecomp(double *A,double *L,double *U,double *P,double *I,int n)
{



  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      if( i==j){
        L[i*n+j]=0;
        U[i*n+j]=A[i*n+j];
        P[i*n+j]=1;
        I[i*n+j]=1;}
      else{
      L[i*n+j]=0;
       U[i*n+j]=A[i*n+j];
        P[i*n+j]=0;
        I[i*n+j]=0;
      }
    }
  }





  //************************************************************                                        

for(int row = 0; row<n; row++){


    double* Acols = (double*) malloc((n-row)*sizeof(double));

    int ii=0;
    for(int i = row; i<n; i++){
       Acols[ii] = abs(A[i*n+row]);

     ii++;
    }



    int indx = distance(Acols,max_element(Acols,Acols+n-row));
    indx = indx + row;                                                                        

    free(Acols);


    double tmp;
  for (int j=0; j<n ; j++){


    if(j<n-1){
      tmp=L[row*n+j];
    L[row*n+j]=L[indx*n+j];
    L[indx*n+j]=tmp;}

    tmp=U[row*n+j];
    U[row*n+j]=U[indx*n+j];
    U[indx*n+j]=tmp;

    tmp=P[row*n+j];
    P[row*n+j]=P[indx*n+j];
    P[indx*n+j]=tmp;

  }


    for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      if( i==j){
        L[i*n+j]=1;}
      }
    }


  for (int col=row+1;col<n;col++){
    L[col*n+row] = U[col*n+row]/U[row*n+row];

      for (int k=0;k<n;k++){
        U[col*n+k]=U[col*n+k]-L[col*n+row]*U[row*n+k];
      }

  }


 }


checkmat(L,n);
 checkmat(U,n);
 checkmat(P,n);


}


                                                                                          

void SolveSysOEqs(double* L,double* U,double* P,double* b0,double* x, int n)
{


  double* y = (double*) malloc(n*sizeof(double));

  double* b = (double*) malloc(n*sizeof(double));
  matvecmult(P,b0,b,n,n);


  double sum;
  for(int i=0;i<n;i++){
    y[i]=b[i];
    for(int j=0;j<i;j++){
      y[i]-= (L[i*n+j]*(y[j]));
    }

  }



  for(int i=n-1;i>=0;i--){
    x[i]=y[i];
    for(int j=i+1;j<n;j++){
      x[i]-=U[i*n+j]*x[j];
    }
    x[i]/=(U[i*n+i]);
  }


  free(b);



}




void checkmat(double* A, int n)
{
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){


      A[i*n+j];

       }

  }

}













void Transpose(double *a,double *at,int row, int col)  
{
   int i,j;


   for (i=0;i<col;i++) {
      for (j=0;j<row;j++) {
	at[i*row+j]=a[j*col+i];
      }
   }
       
     
}








void matvecmult_par(double* A,double* x,double* y,int nfr, int comm_sz, int my_rank)  
{  



   MPI_Bcast(&nfr, 1, MPI_INT, 0,MPI_COMM_WORLD);
   int local_n = nfr/comm_sz;

   double* local_A   = (double*) malloc(local_n*nfr*sizeof(double));
   double* local_x   = (double*) malloc(local_n*sizeof(double));
   double* local_y   = (double*) malloc(local_n*sizeof(double));

   if(my_rank==0){
   MPI_Scatter(A, local_n*nfr, MPI_DOUBLE, local_A, local_n*nfr, MPI_DOUBLE, 0,MPI_COMM_WORLD);}
   else{
     MPI_Scatter(A, local_n*nfr, MPI_DOUBLE, local_A, local_n*nfr, MPI_DOUBLE, 0,MPI_COMM_WORLD);}

   if(my_rank==0){
     MPI_Scatter(x, local_n, MPI_DOUBLE,local_x, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);}
   else{
     MPI_Scatter(x, local_n, MPI_DOUBLE,local_x, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);}

   MPI_Barrier(MPI_COMM_WORLD);
   Mat_vect_mult(local_A, local_x, local_y, x, nfr, local_n);

   if(my_rank==0){
     MPI_Gather(local_y, local_n, MPI_DOUBLE, y, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);}
   else{
     MPI_Gather(local_y, local_n, MPI_DOUBLE, y, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);}


   free(local_A);
   free(local_x);
   free(local_y);
   

}



void Mat_vect_mult(double* loc_A, double* loc_x, double* loc_y,double* x, int n,int loc_n)
{

  int loc_i, j;

   MPI_Allgather(loc_x, loc_n, MPI_DOUBLE, x, loc_n, MPI_DOUBLE, MPI_COMM_WORLD);

   for (loc_i = 0; loc_i < loc_n; loc_i++) {
     loc_y[loc_i] = 0.0;
     for (j = 0; j < n; j++){
       loc_y[loc_i] += loc_A[loc_i*n+j]*x[j];
     }
   }
   
}  

