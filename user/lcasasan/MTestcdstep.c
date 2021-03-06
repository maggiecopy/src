#include <stdio.h>

#include <rsf.h>
#include <rsfpwd.h>

int main (void) {
    float** a;
    float x[4], y[4], res[4];
    int i,j, iter;

    a = sf_floatalloc2(4,4);

    a[0][0] = 1.; a[0][1] = 1.;	a[0][2] = 1.; a[0][3] = 0.;
    a[1][0] = 1.; a[1][1] = 2.;	a[1][2] = 0.; a[1][3] = 0.;
    a[2][0] = 1.; a[2][1] = 3.;	a[2][2] = 1.; a[2][3] = 0.;
    a[3][0] = 1.; a[3][1] = 4.;	a[3][2] = 0.; a[3][3] = 1.;

    y[0]=3.; y[1]=3.; y[2]=5.; y[3]=7.;

    printf ("y = \n");
    for (i=0; i < 4; i ++) {
	printf (" %10.2f",y[i]);
    }
    printf ("\n");
    printf ("a = \n");
    for (j=0; j < 4; j ++) {
	for (i=0; i < 4; i ++) {
	    printf (" %10.2f",a[j][i]);
	}
	printf("\n");
    }
    printf("\n");

    matmult_init(a);

    printf ("cdstep\n------\n");
    for (iter =0; iter < 5; iter++) {
	sf_cdstep_init();
	sf_left_solver( matmult_lop, sf_cdstep, 4, x, y, iter,
			"res", res, "end");
	sf_cdstep_close();
	printf ("x = ");
	for (i=0; i < 4; i ++) {
	    printf (" %12.8f",x[i]);
	}
	printf ("\n");
	printf ("res = ");
	for (i=0; i < 4; i ++) {
	    printf (" %12.8f",res[i]);
	}
	printf ("\n");
    }

    free(a[0]);
    free(a);
    exit(0);
}
