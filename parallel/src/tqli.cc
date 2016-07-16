/*
 * =====================================================================================
 *
 *       Filename:  tqli.cc
 *
 *    Description:  Tridiagonal QL Implicit (tqli) algorithm
 *        Created:  06/14/2016 21:19:24
 *
 *         Author:  Ken Hu, xnchnhu@gmail.com
 *
 * =====================================================================================
 */

#include <exception>
#include <iostream>
#include <cmath>
#include <utility>
#include <limits>
#include "tqli.h"
//#include "vt_user.h"

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a)) // Magnitude of a times sign of b

using namespace std;

double SQR(double a); 
double pythag(double a, double b);
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  pythag (P70, Numerical Recipes)
 *  Description:  Computes sqrt(a^2 + b^2) without destructive underflow or overflow.
 * =====================================================================================
 */

// Square a double value
static double sqrarg;
double SQR(double a) {
    return (sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg;
}

double pythag(double a, double b) {
    double absa,absb;
    absa=fabs(a);
    absb=fabs(b);
    if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
    else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  tqli (P480, Numerical Recipes)
 *  Description:  Calculate the eigenvalues and eigenvectors of a sysmetric triangular matrix
 * =====================================================================================
 */

/*-----------------------------------------------------------------------------
 *  Signiture of the funtion:
 * 	input:
 *		d[0..n-1] contains the diagonal elements of the tridiagonal matrix
 * 		e[0..n-1] inputs the subdiagonal elements of the tridiagonal matrix, with e[n-1] arbitrary
 *		z[0..n-1][0..n-1] is input as the identity matrix.
 *
 *	output:
 *		d[0..n-1] returns the eigenvalues
 *		e[0..n-1] is destroyed
 *		z[0..n-1][0..n-1] the kth column of z returns the normalized eigenvector corresponding to d[k].
 *-----------------------------------------------------------------------------*/

void tqli (vector<double>& d, vector<double>& e, unordered_map<int, vector<double>>& z) {

    //VT_TRACER("TQLI");
    int m,l,iter,i,k;
    double s,r,p,g,f,dd,c,b;
    const double EPS = numeric_limits<double>::epsilon();

    int n = d.size();

    vector<double> vinitial(n, 0);
    for(int i = 0; i < n; i++)	z[i] = vinitial;
    for(int i = 0; i < n; i++)	z[i][i] = 1;
    e.push_back(0.0);


    // Convenient to renumber the elements of e.
    //for (i = 2; i <= n; i++) 
    //	e[i-1] = e[i]; 
    //e[n] = 0.0;
    //
    //for (i = 0; i <= n; i++)
    //	cout << e[i] << " ";
    //cout << endl;

    for (l = 0; l < n; l++) {
        //cout << "tqli" << l << endl;
        iter = 0;
        do {
            // Look for a single small subdiagonal element to split the matrix.
            //for (m = l; m <= n-1; m++) 
            for (m = l; m < n-1; m++) 
            { 
                dd = std::abs(d[m])+std::abs(d[m+1]);
                //if ((double)(std::abs(e[m]) + dd) == dd) break;
                if (std::abs(e[m]) <= EPS * dd) break;
            }
            if (m != l) {
                if (iter++ == 30) throw std::runtime_error("Too many iterations in tqli.");
                g = (d[l+1] - d[l]) / (2.0 * e[l]); // Form shift.
                r = pythag(g,1.0);
                g = d[m] - d[l] + e[l] / (g + SIGN(r,g)); // This is dm − ks.
                s = c = 1.0;
                p = 0.0;

                // A plane rotation as in the original QL, 
                // followed by Givens rotations to restore tridiagonal form.
                for (i = m-1; i >= l; i--) { 
                    f = s * e[i];
                    b = c * e[i];
                    e[i+1] = (r = pythag(f,g));
                    // Recover from underflow.
                    if (r == 0.0) {
                        d[i+1] -= p;
                        e[m] = 0.0;
                        break;
                    }
                    s = f / r;
                    c = g / r;
                    g = d[i+1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    d[i+1] = g + (p = s * r);
                    g = c * r - b;
                    // Next loop can be omitted if eigenvectors not wanted
                    // Form eigenvectors.

                    //for (k = 1; k <= n; k++) { 
                    //	f = z[k][i+1];
                    //	z[k][i+1] = s * z[k][i] + c * f;
                    //	z[k][i] = c * z[k][i] - s * f;
                    //}
                    for (k = 0; k < n; k++) { 
                        f = z[k][i+1];
                        z[k][i+1] = s * z[k][i] + c * f;
                        z[k][i] = c * z[k][i] - s * f;
                    }
                }
                if (r == 0.0 && i >= l) continue;
                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            }
        } while (m != l);
    }
}
