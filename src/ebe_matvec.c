/* C Bind routine takes fortran assumed shape or allocated array as a two
 * dimensional array but it stores all the elements in second index.
 * Therefore the upperbound of the first index is always 0!
 * These C Bind arrays are implemented only in 201X compilers!*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* c routine to be called from fortran */
void ebe_matvec_c(int nelmt, int nedof, int neq, int **gdof, 
                      double **kmat, double **vin, double **vout){
int e,i,ig,j;
int en,enn,jn;
int im,iv;
double mv;
double km[nedof][nedof],g[nedof],v[nedof];

/* initialise vout */
for(j=0;j<neq+1;j++){
  vout[0][j]=0.0;
}

/* perform matrix-vector multiplication */
for(e=0;e<nelmt;e++){
  /* set gdof */
  en=e*nedof;
  enn=en*nedof;
  for(j=0;j<nedof;j++){
    iv=en+j;
    ig=gdof[0][iv]; /* C indices are fortran indices - 1 BUT fortran vector 
                         gin starts from 0 */
    g[j]=ig;
    v[j]=vin[0][ig];

    jn=j*nedof;
    for(i=0;i<nedof;i++){
      im=enn+jn+i;
      km[i][j]=kmat[0][im];
    }
  }

  /* local matrix vector product */
  for(i=0;i<nedof;i++){
    mv=0.0;
    for(j=0;j<nedof;j++){
      mv=mv+km[i][j]*v[j];
    }
    ig=g[i];
    vout[0][ig]=vout[0][ig]+mv;
  }
}

}
/*===========================================================================*/  

