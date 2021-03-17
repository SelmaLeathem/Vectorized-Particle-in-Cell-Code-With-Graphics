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

#ifndef MGLIN_H
#define MGLIN_H

#define _CRT_SECURE_NO_WARNINGS   //Visual Studio Uses to Ignore VS related warnings

#include "nr3.h"


typedef const NRmatrix<double> MatDoub_I;
typedef NRmatrix<double> MatDoub, MatDoub_O, MatDoub_IO;

// Full multigrid algorithm for solution of linear elliptic equation
// on a square domain of side 1, so that delta = 1/(n-1).

class Mglin {
	private:
		int n,ng;
		double length;
		MatDoub *uj,*uj1;
		NRvector<NRmatrix<double> *> rho;  //Vector of pointers to rho on each level.

		void interp(MatDoub_O &uf, MatDoub_I &uc);
		void addint(MatDoub_O &uf, MatDoub_I &uc, MatDoub_O &res);
		void slvsml(MatDoub_O &u, MatDoub_I &rhs);
		void relax(MatDoub_IO &u, MatDoub_I &rhs);
		void resid(MatDoub_O &res, MatDoub_I &u, MatDoub_I &rhs);
		void rstrct(MatDoub_O &uc, MatDoub_I &uf);
		void mg(int j, MatDoub_IO &u, MatDoub_I &rhs);

	public:
	Mglin(MatDoub_IO &u, const int ncycle, double length);
	~Mglin();

};

#endif
