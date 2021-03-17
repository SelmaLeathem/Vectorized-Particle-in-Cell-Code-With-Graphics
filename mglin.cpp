/*
* Source: Numerical Recipes, The Art of Scientific Computing 3rd edition by
*         William H. Press, Saul A. Teukolsky, William T. Vetterling, and Brian P. Flannery
*
*         Section 20.6: Multigrid Methods for Boundary Value Problems
*
*         Solves Poission equation on a two-dimensional grid using the multigrid method
*         which is O(n) as compared to O(n^3) for straight LU decomposition.
*          
*/
#include "mglin.h"
#include "parameters.h"

/*
On input u[0..n-1][0..n-1] contains the right-hand side rho, while on output it returns the
solution. The dimension n must be of the form 2^j + 1 for some integer j. (j is actually
the number of grid levels used in the solution, called ng below.) ncycle is the number of
V-cycles to be used at each level.
*/

Mglin::Mglin(MatDoub_IO &u, const int ncycle, double lengthIn) : n(u.nrows()), ng(0)
	{
		length = lengthIn;
		int nn=n;
		while (nn >>= 1) ng++;
		if ((n-1) != (1 << ng))
			throw("n-1 must be a power of 2 in mglin.");
		nn=n;
		int ngrid=ng-1;
		rho.resize(ng);
        /*
            Allocate storage for r.h.s. on grid ng - 1, and fill it with the input r.h.s.{ Similarly 
            allocate storage and fill r.h.s. on all coarse grids by restricting from finer grids.
        */
		rho[ngrid] = new MatDoub(nn,nn);  
		*rho[ngrid]=u;
		while (nn > 3) {
			nn=nn/2+1;
			rho[--ngrid]=new MatDoub(nn,nn);
			rstrct(*rho[ngrid],*rho[ngrid+1]);
		}
		nn=3;
		uj=new MatDoub(nn,nn);
		slvsml(*uj,*rho[0]);    //Initial solution on coarsest grid.
		for (int j=1;j<ng;j++) {  //Nested iteration loop.
			nn=2*nn-1;
			uj1=uj;
			uj=new MatDoub(nn,nn);
			interp(*uj,*uj1);   //Interpolate from grid j-1 to next finer grid j.
			delete uj1;
			for (int jcycle=0;jcycle<ncycle;jcycle++)   //V-cycle loop.
				mg(j,*uj,*rho[j]);
		}
		u = *uj;        //Return solution in u.
	}

    // Destructor deletes storage.
	Mglin::~Mglin()
	{
		if (uj != NULL) delete uj;
		for (int j=0;j<ng;j++)
			if (rho[j] != NULL) delete rho[j];
	}

/*
Coarse-to-fine prolongation by bilinear interpolation. If nf is the fine-grid dimension, the
coarse-grid solution is input as uc[0..nc-1][0..nc-1], where nc = nf/2 + 1. The finegrid
solution is returned in uf[0..nf-1][0..nf-1].
*/
	void Mglin::interp(MatDoub_O &uf, MatDoub_I &uc)
	{
		int nf=uf.nrows();
		int nc=nf/2+1;
		for (int jc=0;jc<nc;jc++)   //Do elements that are copies.
			for (int ic=0;ic<nc;ic++) uf[2*ic][2*jc]=uc[ic][jc];
		for (int jf=0;jf<nf;jf+=2)  //Do even-numbered columns, interpolating vertically.
			for (int iif=1;iif<nf-1;iif+=2)
				uf[iif][jf]=0.5*(uf[iif+1][jf]+uf[iif-1][jf]);
		for (int jf=1;jf<nf-1;jf+=2)    //Do odd-numbered columns, interpolating horizontally.
			for (int iif=0;iif<nf;iif++)
				uf[iif][jf]=0.5*(uf[iif][jf+1]+uf[iif][jf-1]);
	}

/*
Does coarse-to-fine interpolation and adds result to uf. If nf is the fine-grid dimension,
the coarse-grid solution is input as uc[0..nc-1][0..nc-1], where nc = nf/2 + 1. The
fine-grid solution is returned in uf[0..nf-1][0..nf-1]. res[0..nf-1][0..nf-1] is used
for temporary storage.
*/
	void Mglin::addint(MatDoub_O &uf, MatDoub_I &uc, MatDoub_O &res)
	{
		int nf=uf.nrows();
		interp(res,uc);
		for (int j=0;j<nf;j++)
			for (int i=0;i<nf;i++)
				uf[i][j] += res[i][j];
	}

/*
Solution of the model problem on the coarsest grid, where h = 1/2 . The right-hand side is
input in rhs[0..2][0..2] and the solution is returned in u[0..2][0..2].
*/
	void Mglin::slvsml(MatDoub_O &u, MatDoub_I &rhs)
	{
		double h=0.5*length; // h = 0.5;
		for (int i=0;i<3;i++)
			for (int j=0;j<3;j++)
				u[i][j]=0.0;
		u[1][1] = -h*h*rhs[1][1]/4.0;
	}

/*
Red-black Gauss-Seidel relaxation for model problem. Updates the current value of the
solution u[0..n-1][0..n-1], using the right-hand side function rhs[0..n-1][0..n-1].
*/
	void Mglin::relax(MatDoub_IO &u, MatDoub_I &rhs)
	{
		int n=u.nrows();
		double h=length/(n-1);  //1.0/(n-1);
		double h2=h*h;
		for (int ipass=0,jsw=1;ipass<2;ipass++,jsw=3-jsw) {  //Red and black sweeps.
			for (int j=1,isw=jsw;j<n-1;j++,isw=3-isw)
				for (int i=isw;i<n-1;i+=2)  //Gauss-Seidel formula.
					u[i][j]=0.25*(u[i+1][j]+u[i-1][j]+u[i][j+1]
						+u[i][j-1]-h2*rhs[i][j]);
		}
	}

/*
Returns minus the residual for the model problem. Input quantities are u[0..n-1][0..n-1]
and rhs[0..n-1][0..n-1], while res[0..n-1][0..n-1] is returned.
*/
	void Mglin::resid(MatDoub_O &res, MatDoub_I &u, MatDoub_I &rhs)
	{
		int n=u.nrows();
		double h=length/(n-1);  //1.0/(n-1);
		double h2i=1.0/(h*h);
		for (int j=1;j<n-1;j++) // Interior points.
			for (int i=1;i<n-1;i++)
				res[i][j] = -h2i*(u[i+1][j]+u[i-1][j]+u[i][j+1]
					+u[i][j-1]-4.0*u[i][j])+rhs[i][j];
		for (int i=0;i<n;i++)   // Boundary points.
			res[i][0]=res[i][n-1]=res[0][i]=res[n-1][i]=0.0;
	}

/*
Half-weighting restriction. If nc is the coarse-grid dimension, the fine-grid solution is input
in uf[0..2*nc-2][0..2*nc-2].The coarse-grid solution obtained by restriction is returned
in uc[0..nc-1][0..nc-1].
*/
	void Mglin::rstrct(MatDoub_O &uc, MatDoub_I &uf)
	{
		int nc=uc.nrows();
		int ncc=2*nc-2;
		for (int jf=2,jc=1;jc<nc-1;jc++,jf+=2) {  //Interior points.
			for (int iif=2,ic=1;ic<nc-1;ic++,iif+=2) {
				uc[ic][jc]=0.5*uf[iif][jf]+0.125*(uf[iif+1][jf]+uf[iif-1][jf]
					+uf[iif][jf+1]+uf[iif][jf-1]);
			}
		}
		for (int jc=0,ic=0;ic<nc;ic++,jc+=2) {  //Boundary points.
			uc[ic][0]=uf[jc][0];
			uc[ic][nc-1]=uf[jc][ncc];
		}
		for (int jc=0,ic=0;ic<nc;ic++,jc+=2) {
			uc[0][ic]=uf[0][jc];
			uc[nc-1][ic]=uf[ncc][jc];
		}
	}

/*
Recursive multigrid iteration. On input, j is the current level, u is the current value of the
solution, and rhs is the right-hand side. On output u contains the improved solution at the
current level.
*/
	void Mglin::mg(int j, MatDoub_IO &u, MatDoub_I &rhs)
	{
		//Number of relaxation sweeps before and after the coarse-grid correction is computed.
		const int NPRE=1,NPOST=1;
		int nf=u.nrows();
		int nc=(nf+1)/2;
		if (j == 0)		//Bottom of V: Solve on coarsest grid.
			slvsml(u,rhs);
		else {		//On downward stoke of the V.
			MatDoub res(nc,nc),v(nc,nc,0.0),temp(nf,nf);
			//v is zero for initial guess in each relaxation.
			for (int jpre=0;jpre<NPRE;jpre++)
				relax(u,rhs);	//Pre-smoothing.
			resid(temp,u,rhs);
			rstrct(res,temp);	//Restriction of the residual is the next r.h.s.
			mg(j-1,v,res);		//Recursive call for the coarse-grid correction.
			addint(u,v,temp);	//On upward stroke of V.
			for (int jpost=0;jpost<NPOST;jpost++)
				relax(u,rhs);	//Post-smoothing.
		}
	}