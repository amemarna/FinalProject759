#include<math.h>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<stdio.h>

#include "MyVonMises.cpp"


using namespace std;



float invmatrix(float *A,float *Ainv,int n);
void invmatrix2(float *A,float *Ainv,int n);
void LUdecomp(float *A,float *L,float *U,float *P,float *I,int n);
void SolveSysOEqs(float* L,float* U,float* b,float* x, int n);
void checkmat(float* A, int n);
void matvecmult(float* a, float* b, float* c, const int row, const int col);
void matrixmult(float* a, float* b, float* c, const int row, const int inner, const int col);
void Transpose(float *a,float *at,int row, int col);
void InitializeArray(float* A, int row, int col);
void Emodsix2three(float* De, float* Ce3);
//void mymises(float* sig_el,float kap_el,float* eps_el,float* Etng_el,float* oldsig_el,float oldkap_el,float* oldeps_el,float* Deps,float E,float nu,float H,float sy);
void AddArray(float* A, float* B, float* C, int row, int col, int op);
void ApplyBCs(float* A,float* Af,int ndofs,int nfr,float* frdofs,int id);
float CompResNorm(float* A,int n);
void GlbAsmbly(float* Kel,float* fel,float* Kg,float* fg,int* dofs_el,int ndofs);



void readmatrix(FILE* fp, float* in, const int nrows, const int ncols)
{
  for (int i=0;i<nrows;i++){
    for (int j=0;j<ncols;j++){
      int r=fscanf(fp,"%f",&in[i*nrows+j]);
      //fp << in[i*nrows+j];
    }
  }
}



//void hooke(float* De, const float E, const float nu)
void hooke( const float E, const float nu)
{

  // float* De=(float*)malloc(6*6*sizeof(float));
  float fc = E/((1+nu)*(1-2*nu));
  float D_e[6][6] = {{fc*(1-nu),   fc*nu,       fc*nu,        0,            0,                0},
		     {fc*nu,       fc*(1-nu),   fc*nu,        0,            0,                0},
		     {fc*nu,       fc*nu,       fc*(1-nu),    0,            0,                0},
		     {   0,        0,           0,       fc*(1-2*nu)/2,     0,                0},
		     {   0,        0,           0,            0,      fc*(1-2*nu)/2,          0},
		     {   0,        0,           0,            0,            0,         fc*(1-2*nu)/2} };
  

  //  return De;
}


void coordxtr(float* edof,float* coords, float* ex, float* ey, const int nels)
{
  int xi, yi;
  
  for (int i=0;i<nels;i++){
    for (int j=0;j<3;j++){
      xi = edof[i*7+(2*j+1)]+1;
      yi = edof[i*7+(2*j+2)]+1;
      ex[i*3+j]=coords[xi];
      ey[i*3+j]=coords[yi];
    }
  }

}


float prescribedu(const float t, const float maxdisp)
{
  return maxdisp*t;
}






int main()
{  


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

  int bc0[23][2]={{17, 0}, 
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



  
  //% Rectangelar plate with a centered hole given in meters
  //% width  = 0.2        
  //% height = 0.2
  //% radius = width/4
  //load metal_sheet
  FILE *fpe = fopen("edofdata.inp","r");
  FILE *fpc = fopen("coordata.inp","r");
  FILE *fph = fopen("hooke.inp","r");
  FILE *fpg = fopen("gamma.inp","r");
  
  int nels=176;
  int neldofs=6;
  int nods=108;
  float* edofs  = (float*) malloc(nels*(neldofs+1)*sizeof(float));
  float* coords = (float*) malloc(2*nods*sizeof(float));
  float* De     = (float*) malloc(6*6*sizeof(float));
  float* gamma  = (float*) malloc(9*6*sizeof(float));
  readmatrix(fpe,edofs,nels,neldofs+1);
  readmatrix(fpc,coords,nods,2);
  readmatrix(fph,De,6,6);
  readmatrix(fpg,gamma,9,6);
  
  float t      = 0.001;   // % plate thickness
  //  int ndofs  = max_element(edofs,edofs+nels*(neldofs+1));//numel(Dof);
// int nrelem = size(edof,1);
  int ndofs = 216;
  
//% Analysis parameters
int NrOfTsteps = 100; // % total number of timesteps
float reltol     = 1e-6; //% tolerance for equilibrium iterations
int ptype      = 1;   // % 1 = plane stress 
int magn       = 10;  // % magnification factor for deformed mesh
int maxiter    = 30;  // % max number of iterations

//% Material parameters for elastic plastic material
 float E          = 200*pow(10,9); 
float nu         = 0.3;
float H          = 0.2*E;
float sy         = E/(1e2);//%!!remember u changed it
//matmod     = 'plast'; % Change this to use "plast" instead of elast.
//float De[6][6];
//hooke(De,E,nu);
hooke(E,nu);
//% Converts from gamma notation to pure epsilon 9x1
//gam = [eye(3),zeros(3);zeros(3),eye(3)/2;zeros(3),eye(3)/2];

float gam[9][6]={ {1,0,0,0,0,0},
                  {0,1,0,0,0,0},
                  {0,0,1,0,0,0},
                  {0,0,0,0.5,0,0},
                  {0,0,0,0,0.5,0},
                  {0,0,0,0,0,0.5},
                  {0,0,0,0.5,0,0},
                  {0,0,0,0,0.5,0},
                  {0,0,0,0,0,0.5} };



//% sc = controlled stresses 
//sc = logical([0,0,1,0,1,1,0,1,1]); 
//ec = ~sc; %strain-controlled

 float* ex=(float*)malloc(3*nels*sizeof(float));
 float* ey=(float*)malloc(3*nels*sizeof(float));
//[ex,ey] = coordxtr(Edof,Coord,Dof,nen);
 coordxtr(edofs,coords,ex,ey,nels);

//% Solves for a linear material for reference.
//% As mentions in the project support, you may use this, or check during the main iterations.
//%Klin = spalloc(nrdof,nrdof,nen*2*nrelem);
//%for i=1:nrelem
//%    Ke = plante(ex(i,:),ey(i,:),[ptype t],De);
//%    Klin = assem(Edof(i,:),Klin,Ke);
//%end
//%u0 = solveq(Klin,zeros(nrdof,1),bc0); % Solves isotropic, elastic, problem for a unit displacement bc0.
//%u_el = ??? % Part of assignment 1.c to find this.
float u_el = 1e-3; //% You can still run with just any value. Might need to adjust in order to get decent time steps.
float maxdisp = 5*u_el;// % Maximum displacement 5 times the elastic limit.

  //modpar.E = E;
  //modpar.H = H;
  //modpar.nu = nu;
  //modpar.sy = sy;

float* u   =(float*) malloc(ndofs*sizeof(float));//zeros(nrdof,1);
float* sig =(float*) malloc(9*nels*sizeof(float));// zeros(9,nrelem);
float* eps =(float*) malloc(9*nels*sizeof(float));// zeros(9,nrelem);
float* oldsig =(float*) malloc(9*nels*sizeof(float));// zeros(9,nrelem);
float* oldeps =(float*) malloc(9*nels*sizeof(float));// zeros(9,nrelem);
float* oldeps_p =(float*) malloc(1*nels*sizeof(float));// zeros(1,nrelem);
float* eps_p = (float*) malloc(1*nels*sizeof(float));//=oldeps_p;
float* oldkap=(float*) malloc(1*nels*sizeof(float));//zeros(1,nrelem);
float* kap = (float*) malloc(1*nels*sizeof(float));//=oldkap; 
float* u_hist    =(float*) malloc(NrOfTsteps*sizeof(float));// zeros(NrOfTsteps,1); 
float* FReaction =(float*) malloc(NrOfTsteps*sizeof(float));// zeros(NrOfTsteps,1);

float* eldisps=(float*)malloc(nels*(neldofs+1)*sizeof(float));

//prescribedu = @(t) maxdisp*t;%sin(3*2*pi*t);

//% Global timestepping
float olddisp=0.f;
float* du =(float*) malloc(ndofs*sizeof(float));// zeros(nrdof,1);

 float dispcontrolled=0.f;
 int nbcld=11;
 int nbc = 23;


      float* cstd =(float*) malloc(nbc*sizeof(float));
      int nfr0 = ndofs-nbc;
      float* frdofs =(float*) malloc(nfr0*sizeof(float));
      int ii = 0;
      int jj = 0;
      int nfr = 0;
      for(int i=0;i<ndofs;i++){
	  for(int j=0;j<nbc;j++){
		if(i==bc0[j][0]){
		cstd[jj]=bc0[j][0];
		jj++;}
		else{
		  frdofs[ii]=i;
		  ii++;
		  nfr++;}
	  }
      }

 
 float gp1[3]={0.16666666, 0.66666666, 0.16666666};
 float gp2[3]={0.16666666, 0.16666666, 0.66666666};
 float gw[3]={0.16666666, 0.16666666, 0.16666666};


  float* fg = (float*) malloc(ndofs*sizeof(float));
 InitializeArray(FReaction,1,NrOfTsteps);
 
 for (int Tstep=1;Tstep<NrOfTsteps;Tstep++){
     
    dispcontrolled = prescribedu(Tstep/NrOfTsteps,maxdisp); //% Load scalar value of prescribed di\

    //du(bcload(:,1)) = dispcontrolled - olddisp;
    for (int m=0;m<ndofs;m++){
      for (int n=0;n<nbcld;n++){
	if(m==bcload[n][0]){
	  du[m]=dispcontrolled - olddisp;}
      }
    }

    
    //printf("Timestep %2.0f Displacement %.2f \n",Tstep,dispcontrolled);                                        

    for(int iter; iter<maxiter; iter++){


      	

      	  for(int i=0;i<nels;i++){
	       eldisps[i*7+0]=edofs[i*7+0];
	  }
	
      
      int pos;
      for (int k=0;k<nbcld;k++){
	  for(int i=0;i<nels;i++){
	    for(int j=1;j<7;j++){
	      if( edofs[i*7+j]==bcload[k][0]){
		      pos = edofs[i*7+j];
		      eldisps[i*7+j]=du[pos];}
	    }
	  }
      }


      float* Kg = (float*) malloc(ndofs*ndofs*sizeof(float));
      InitializeArray(Kg,ndofs,ndofs);
      InitializeArray(fg,1,ndofs);
      
      for(int iel=0; iel<nels-1; iel++){
      // [sig(:,i),kap(:,i),eps(:,i),E_tang]=mymises(oldsig(:,i),oldkap(:,i),oldeps(:,i),Deps,modpar);%,oldeps_p(:,i)


	float* oldsig_el = (float*) malloc(9*sizeof(float));
	float* oldeps_el = (float*) malloc(9*sizeof(float));
	for(int j=0; j<9; j++){
	oldsig_el[j] = oldsig[iel*9+j];
	oldeps_el[j] = oldeps[iel*9+j];
	}

	float oldkap_el=oldkap[iel];


	
	//Se = inv(De);
	//Ce = inv(Se[1 2 4],[1 2 4]);
	//trsig([1 2 4]) = Ce*;

	float* disp_el = (float*) malloc(neldofs*sizeof(float));
	for(int j=1;j<7;j++){
	  disp_el[j-1]=eldisps[iel*7+j];}


	//float elcrds[3][2];
	float* elcrds =(float*) malloc(3*2*sizeof(float));
	for(int i=0;i<3;i++){
	  // elcrds[i][0]=ex[iel*3+i];
	  // elcrds[i][1]=ey[iel*3+i];
	    elcrds[i*2+0]=ex[iel*3+i];
	    elcrds[i*2+1]=ey[iel*3+i];
	}


	int ngp=3;
	float ksi,etha;
	float N1,N2,N3;
	float dN1dksi, dN1detha, dN2dksi, dN2detha, dN3dksi, dN3detha;
	float* dNdksi_v =(float*) malloc(6*sizeof(float));
	float* dNdksi =(float*) malloc(2*3*sizeof(float));
	float* Jmat  =(float*) malloc(2*2*sizeof(float));
	float* Jinv  =(float*) malloc(2*2*sizeof(float));
	float* dNdx  =(float*) malloc(2*3*sizeof(float));
	float* Bmat  =(float*) malloc(3*6*sizeof(float));
	float* Bmat_v  =(float*) malloc(18*sizeof(float));

	float* Kel = (float*) malloc(6*6*sizeof(float));
        float* fel = (float*) malloc(6*sizeof(float));

	InitializeArray(Kel,6,6);
	InitializeArray(fel,1,6);
	
	float* sig_el   = (float*) malloc(9*sizeof(float));
	float* eps_el   = (float*) malloc(9*sizeof(float));
	float* Etang_el = (float*) malloc(9*9*sizeof(float));
	float kap_el;

	for (int igp=0;igp<ngp;igp++){
	  //for (int jgp=0;jgp<ngp;jgp++){

	  
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


	  // dNdksi={ {dN1dksi, dN2dksi, dN3dksi},
	  //	   {dN1detha, dN2detha, dN3detha} };

	  dNdksi_v[0]=dN1dksi;dNdksi_v[1]=dN2dksi;dNdksi_v[2]=dN3dksi;
	  dNdksi_v[3]=dN1detha;dNdksi_v[4]=dN2detha;dNdksi_v[5]=dN3detha;
	  //dNdksi_v={dN1dksi, dN2dksi, dN3dksi,dN1detha,dN2detha,dN3detha};

	  int cnt=0;
	  for(int i=0;i<2;i++){
	    for(int j=0;j<3;j++){
	      dNdksi[i*2+j]=dNdksi_v[cnt];
	      cnt+=1;
	    }
	  }
	  
	  // Jmat = dNdksi*[ex(i,:)^T, ey(i,:)^T];
	  matrixmult(dNdksi,elcrds,Jmat,2,3,2);
	  invmatrix2(Jmat,Jinv,2);
	  //dNdx=Jinv*dNdksi;
	  matrixmult(Jinv,dNdksi,dNdx,2,2,3);
	  
	  // Bmat = { {dN1dx, 0, dN2dx, 0, dN3dx, 0},
	  //         {0, dN1dy, 0, dN2dy, 0, dN3dy},
	  // 	   {dN1dy, dN1dx, dN2dy, dN2dx, dN3dy, dN3dx} };

	  // Bmat_v = {dN1dx, 0, dN2dx, 0, dN3dx, 0 , 0, dN1dy, 0, dN2dy, 0, dN3dy, dN1dy, dN1dx, dN2dy, dN2dx, dN3dy, dN3dx};


 	for(int i=0;i<3;i++){
        for(int j=0;j<6;j++){
	Bmat[i*6+j]=0;
	}
	}

	  
	  
	  int cnte=0;
	  int cnto=3;
	  for(int i=0;i<3;i++){
	    for(int j=0;j<6;j++){
	      if(i==0 && j%2==0){
		Bmat[i*6+j]=dNdx[cnt];
	        cnt+=1;}
	      else if (i==1 && j%2>0){
	      Bmat[i*6+j]=dNdx[cnt];
	      cnt+=1;}
	      else if (i==2 && j%2==0){
	      Bmat[i*6+j]=dNdx[cnto];
	      cnto+=1;}
	      else if (i==2 && j%2>0){
	      Bmat[i*6+j]=dNdx[cnte];
	      cnte+=1;}
	    }
	  }


	  cnt=0;
	  for(int i=0;i<3;i++){
	    for(int j=0;j<6;j++){
	      Bmat[i*6+j]=Bmat_v[cnt];
	      cnt+=1;
	    }
	  }


	  //int row=0;
	  //int col=0;

	  float* Se = (float*) malloc(6*6*sizeof(float));
	  float* Se3= (float*) malloc(3*3*sizeof(float));

	  invmatrix(De,Se,6);
	  
	  Emodsix2three(Se,Se3);
	  //for(int i=0;i<6;i++){
	  //  for(int j=0;j<6;j++){
	  //    if(i==0 || i==1 || i==3){
	  //	row+=1;
	  //	 if(j==0 || j==1 || j==3){
	  //	   col+=1;
	  //	   Se3[row*3+col]=Se[i*6+j];
	  //	 }
	  //	 col=0;
	  //    }
	  //  }
	  //}

	  
	  float* Ce3 = (float*) malloc(3*3*sizeof(float));
	  invmatrix(Se3,Ce3,3);
	  
	  
	  //Deps3 = Bmat*disp_el;
	  float* Deps3= (float*) malloc(3*sizeof(float));
	  matvecmult(Bmat,disp_el,Deps3,3,6);

	  float* t_sig= (float*) malloc(6*sizeof(float));
          for(int j=0;j<6;j++){
	  t_sig[j]=0;
	  }

	  float* t_sig3= (float*) malloc(3*sizeof(float));
	  
	  //t_sig3 = Ce*Deps;
	  matvecmult(Ce3,Deps3,t_sig3,3,3);

	  cnt=0;
	  for(int j=0;j<6;j++){
		 if(j==0 || j==1 || j==3){
		   cnt+=1;
		   t_sig[j]=t_sig3[cnt];
		 }
	  }
	
	  float* Depsv= (float*) malloc(6*sizeof(float));	  
	  matvecmult(Se,t_sig,Depsv,6,6);

	  float* Deps= (float*) malloc(9*sizeof(float));	  
	  matvecmult(gamma,Depsv,Deps,9,6);
	  
	  float* Etng_el = (float*) malloc(9*9*sizeof(float));
		  
	  MyVonMises(sig_el,kap_el,eps_el,Etng_el,oldsig_el,oldkap_el,oldeps_el,Deps,E,nu,H,sy);//%,oldeps_p(:,i);

	  float* Etngr = (float*) malloc(6*6*sizeof(float));
	  
	  for(int i=0;i<6;i++){
	    for(int j=0;j<6;j++){
	      if(i>2 && j>2){
		Etngr[i*6+j] =1/2*Etng_el[i*9+j];}
	      else {
		Etngr[i*6+j] = Etng_el[i*9+j];}
	    }
	  }

	  float* BmatT = (float*) malloc(6*3*sizeof(float));
	  Transpose(Bmat,BmatT,3,6);
	  float* CE3 = (float*) malloc(3*3*sizeof(float));
	  Emodsix2three(Etngr,CE3);

	  float detJ=Jmat[0]*Jmat[3]-Jmat[1]*Jmat[2];

	  float th = 1e-3;
	  
	  //Kel+ = BmatT*CE3*Bmat*detJ*gw[igp]*gw[jgp]*th;
	  float* BTC = (float*) malloc(6*3*sizeof(float));
	  matrixmult(BmatT,CE3,BTC,6,3,3);
	  float* BTCB = (float*) malloc(6*6*sizeof(float));
	  matrixmult(BTC,Bmat,BTCB,6,3,6);

	  for(int i=0;i<6;i++){
	    for(int j=0;j<6;j++){
	      Kel[i*6+j] = Kel[i*6+j] + BTCB[i*6+j]*detJ*gw[igp]*th;
	    }
	  }
	    //Kel+ = BTCB*detJ*gw[igp]*gw[jgp]*th;

     	  float* sig_el3 = (float*) malloc(3*sizeof(float));
	  cnt=0;
	  for(int j=0;j<9;j++){
		 if(j==0 || j==1 || j==3){
		   cnt+=1;
		   sig_el3[cnt]=sig_el[j];
		 }
	  }

	  float* BTsig = (float*) malloc(6*sizeof(float));
	  matvecmult(BmatT,sig_el3,BTsig,3,3);

	  for(int i=0;i<6;i++){
	      fel[i] = fel[i] + BTsig[i]*detJ*gw[igp]*th;
	  }
	  //fel+ = BTsig*detJ*gw[igp]*th;
	  //fel+ = BmatT*sig_el3*detJ*gw[igp]*th;
	  


	}
	//}   //end of gp integration

	int* dofs_el = (int*) malloc(6*sizeof(int));
	int cnt=0;
	for (int j=1;j<7;j++){
	  dofs_el[cnt]=edofs[iel*7+j];
	  cnt++;
	}

	
	GlbAsmbly(Kel,fel,Kg,fg,dofs_el,ndofs);
	//for(int i=0;i<6;i++){
	//for(int j=0;j<6;j++){
	//    Kg[dofs_el[i]*ndofs+dofs_el[j]] = Kg[dofs_el[i]*ndofs+dofs_el[j]] + Kel[i*6+j];
	//}
	//}

	for(int j=0; j<9; j++){
	oldsig[iel*9+j] = sig_el[j];
	oldeps[iel*9+j] = eps_el[j];
	}

	oldkap[iel]=kap_el;
	
      }  //end of element loop

      float* fext = (float*) malloc(ndofs*sizeof(float));
      InitializeArray(fext,ndofs,1);
      float* res  = (float*) malloc(ndofs*sizeof(float));
      float* resf = (float*) malloc(nfr*sizeof(float));
      
      AddArray(fext,fg,res,ndofs,1,1);
	//res = fext - fint;
      
      
      ApplyBCs(res,resf,ndofs,nfr,frdofs,0);

      float err = CompResNorm(resf,nfr);
      float tol = 1e-5;
      if (err<tol){
	break;}
      printf("%d %f \n",iter,err);

      float* Kgf = (float*) malloc(nfr*nfr*sizeof(float));
      ApplyBCs(Kg,Kgf,ndofs,nfr,frdofs,1);


      float* L = (float*) malloc(nfr*nfr*sizeof(float));
      float* U = (float*) malloc(nfr*nfr*sizeof(float));
      float* I = (float*) malloc(nfr*nfr*sizeof(float));
      float* P = (float*) malloc(nfr*nfr*sizeof(float));

      LUdecomp(Kgf,L,U,P,I,nfr);

      float* ddu = (float*) malloc(nfr*sizeof(float));
      SolveSysOEqs(L,U,resf,ddu,nfr);      

      
      for(int i=0;i<nfr;i++){
	pos  = frdofs[i];
	du[pos] = du[pos] + ddu[i];
      }
      //AddArray(du,ddu,du,ndofs,1,1);
      
    }

    AddArray(u,du,u,ndofs,1,1);
    olddisp = dispcontrolled;
    u_hist[Tstep]=dispcontrolled;
    int pos;
    for(int i=0;i<nbcld;i++){
      pos = bcload[i][0];
      FReaction[Tstep] = FReaction[Tstep] + 2*fg[pos];
    }
    
 }

 return 0;

}

//}


    
//return 0;

 
//}


void GlbAsmbly(float* Kel,float* fel,float* Kg,float* fg,int* dofs_el,int ndofs)
{


  for(int i=0;i<6;i++){
    for(int j=0;j<6;j++){
	      Kg[dofs_el[i]*ndofs+dofs_el[j]] = Kg[dofs_el[i]*ndofs+dofs_el[j]] + Kel[i*6+j];
    }
  }
	
  for(int i=0;i<6;i++){
	      fg[dofs_el[i]] = fg[dofs_el[i]] + fel[i];
  }
  
}



void InitializeArray(float* A, int row, int col)
{
  
 	for(int i=0;i<row;i++){
	  for(int j=0;j<col;j++){
	    A[i*col+j]=0;
	  }
	}

}


float CompResNorm(float* A,int n)
{

      float sum=0;
      for (int i=0;i<n;i++){
      sum+=sqrt(A[i]*A[i]);
      }

      return sum;
}


void AddArray(float* A, float* B, float* C, int row, int col, int op)
{

  //(if op=0 --> add)
  
 	for(int i=0;i<row;i++){
	  for(int j=0;j<col;j++){
	    if(op==0){
	      C[i*col+j]=A[i*col+j]+B[i*col+j];}
	    else{
	    C[i*col+j]=A[i*col+j]-B[i*col+j];}
	  }
	}

}


void ApplyBCs(float* A,float* Af,int ndofs,int nfr,float* frdofs,int id)
{

  if( id==0){ //id=0 --> vector
  
      int cnt=0;
      for(int i=0;i<ndofs;i++){
	for(int j=0;j<nfr;j++){
	  if(i==frdofs[j]){
	    Af[cnt]=A[i];
	    cnt++;}
	}
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



void Emodsix2three(float* De, float* Ce3)
{  
	  int row=0;
	  int col=0;

	  float* Se = (float*) malloc(6*6*sizeof(float));
	  float* Se3= (float*) malloc(3*3*sizeof(float));

	  invmatrix(De,Se,6);

	  for(int i=0;i<6;i++){
	    for(int j=0;j<6;j++){
	      if(i==0 || i==1 || i==3){
		row+=1;
		 if(j==0 || j==1 || j==3){
		   col+=1;
		   Se3[row*3+col]=Se[i*6+j];
		 }
		 col=0;
	      }
	    }
	  }
	  
	  //float* Ce3 = (float*) malloc(3*3*sizeof(float));
	  invmatrix(Se3,Ce3,3);

	  delete(Se);
	  delete(Se3);

	  //return Ce3;
}


void matrixmult(float* a, float* b, float* c, const int row, const int inner, const int col){

  int i,j,k;
  float sum;

  for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			sum = 0;
			for (k = 0; k < inner; k++) {

				sum = sum + a[i*row + k] * b[k*row + j];

			}
			c[i*row + j] = sum;
		}

  }

}


void matvecmult(float* a, float* b, float* c, const int row, const int col)
{


  int i,j;
  for (i = 0; i < row; i++) {
     float sum = 0;
                for (j = 0; j < col; j++) {
                  //sum = 0;                                                                           \
                                                                                                        
                        //for (k = 0; k < inner; k++) {                                                \
                                                                                                        

                  //sum = sum + a[i*row + k] * b[k*row + j];                                           \
                                                                                                        
                                sum = sum + a[i*row + j] * b[j];
                                //}                                                                    \
                                                                                                        
                                //c[i*row + j] = sum;                                                  \
                                                                                                        
                }
c[j] = sum;
  }

}





void invmatrix2(float *A,float *Ainv,int n)
{
  float det=A[0]*A[3]-A[1]*A[2];

      Ainv[0]=(1/det)*A[3];
      Ainv[1]=(-1/det)*A[1];
      Ainv[2]=(-1/det)*A[2];
      Ainv[3]=(1/det)*A[0];

}


float invmatrix(float *A,float *Ainv,int n)
{


  float* LL = (float*) malloc(n*n*sizeof(float));
  float* UU = (float*) malloc(n*n*sizeof(float));
  float* II = (float*) malloc(n*n*sizeof(float));
  float* PP = (float*) malloc(n*n*sizeof(float));

  LUdecomp(A,LL,UU,PP,II,n);

          cout<<"II ="<<'\n';
          checkmat(II,n);

          float* A_col = (float*) malloc(n*sizeof(float));
          float* I_col = (float*) malloc(n*sizeof(float));

          for (int i=0;i<n;i++){

            int jj=0;
            for(int j=i;j<n*n;j+=n){
              I_col[jj]=II[j];
              cout<<j<<"*"<<I_col[jj]<<'\n';
              jj+=1;
            }


            // matvecmult(P,I_col,I_col,3,3);                                                           
            SolveSysOEqs(LL,UU,I_col,A_col,n);

            //jj=0;                                                                                     
            //for(int j=0;j<n;j+=n){                                                                    
            //  cout<<j<<"#"<<A_col[jj]<<'\n';                                                          
            //}                                                                                         

            jj=0;
            for(int j=i;j<n*n;j+=n){
              Ainv[j]=A_col[jj];
              jj++;
            }
            delete(A_col);

          }

          cout<<"Ainv ="<<'\n';
          checkmat(Ainv,n);

	  delete(LL);
	  delete(UU);
	  delete(PP);
	  delete(II);
	  
          return 0;
}



void LUdecomp(float *A,float *L,float *U,float *P,float *I,int n)
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


    float* Acols = (float*) malloc((n-row)*sizeof(float));
    //int vsz = (n-row)*sizeof(float);                                                                  
    //vector<float> Acols(vsz);                                                                         
    int ii=0;
    for(int i = row; i<n; i++){
       Acols[ii] = abs(A[i*n+row]);
       //cout<<ii<<"  "<<Acols[ii]<<endl;                                                               
     ii++;
    }


    //float indx = 0;                                                                                   
    //int indx = max_element(Acols,Acols+n-row)-Acols;                                                  
    int indx = distance(Acols,max_element(Acols,Acols+n-row));
    indx = indx + row;//(row-1);                                                                        

    //cout<<row << indx <<'\n';                                                                         

    float tmp;
  for (int j=0; j<n ; j++){

    //if(row>0 && j<n-1){                                                                               
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

  cout<<"L ="<<'\n';
checkmat(L,n);
  cout<<"U ="<<'\n';
 checkmat(U,n);
  cout<<"P ="<<'\n';
 checkmat(P,n);


}





  //************************************************************                                       \
                                                                                                        

void SolveSysOEqs(float* L,float* U,float* b,float* x, int n)
{


  float* y = (float*) malloc(n*sizeof(float));
  //float* x = (float*) malloc(n*sizeof(float));                                                        


  //y[0] = b[0]/L[0];                                                                                   
  //float sum;                                                                                          
  //for(int i=1;i<n;i++){                                                                               
  //  sum = 0;                                                                                          
  //  for(int j=1;j<i-1;j++){                                                                           
  //    sum+=L[i*n+j]*y[j];                                                                             
  //  }                                                                                                 
  //  y[i]=(1/L[i*n+i])*(b[i]-sum);                                                                     
  // }                                                                                                  


  //x[n-1] = y[n-1]/U[n*n-1];                                                                           
  //for(int i=n-2;i>=0;i--){                                                                            
  //  sum = 0;                                                                                          
  //  for(int j=i+1;j<n;j++){                                                                           
  //    sum+=U[i*n+j]*x[j];                                                                             
  //  }                                                                                                 
  //  x[i]=(1/U[i*n+i])*(y[i]-sum);                                                                     
  // }                                                                                                  




  //y[0] = b[0]/L[0];
  float sum;
  for(int i=0;i<n;i++){
    y[i]=b[i];
    for(int j=0;j<i;j++){
      y[i]-= (L[i*n+j]*(y[j]));
    }
    //y[i]/=(L[i*n+i]);
  }


  //x[n-1] = y[n-1]/U[n*n-1];
  for(int i=n-1;i>=0;i--){
    x[i]=y[i];
    for(int j=i+1;j<n;j++){
      x[i]-=U[i*n+j]*x[j];
    }
    x[i]/=(U[i*n+i]);
  }


  //float det=Determinant(a,n);                                                                        \
                                                                                                        
  //CoFactor(a,n,b);                                                                                   \
                                                                                                        
  //return Transpose(cofmat,n);                                                                        \



}




void checkmat(float* A, int n)
{
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){

      cout<< A[i*n+j] << " ";

       }
    cout<<'\n';
  }

}













/*
   Transpose of a square matrix, do it in place
*/
//void Transpose(float **a,int n)
void Transpose(float *a,float *at,int row, int col)  
{
   int i,j;
   //float tmp;

   
   //   for (i=1;i<n;i++) {
   //   for (j=0;j<i;j++) {
   //	//tmp = a[i][j];
   //	 tmp = a[i*n+j];
   //      //a[i][j] = a[j][i];
   //	 a[i*n+j] = a[j*n+i];
   //         //a[j][i] = tmp;
   //	 a[j*n+i] = tmp;
   //      }
   //   }


   for (i=0;i<col;i++) {
      for (j=0;j<row;j++) {
	at[i*row+j]=a[j*col+i];
      }
   }
       
   //return at;
     
}



