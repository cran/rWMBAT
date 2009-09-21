#include <R.h>

void calcMI(int *v1, int *v2, int *n, double *MI)  
{  
	int *Nx, *Ny, *Nxy;
	Nx = (int *)R_alloc(3, sizeof(int));
	Ny = (int *)R_alloc(3, sizeof(int));
	Nxy = (int *)R_alloc(9, sizeof(int));
	int i,j,k;
	
	for (k=0; k<9; k++) Nxy[k] =0;
	for (k = 0; k < 3; k++) {
		Ny[k] = 0;
		for (j=0; j < *n; j++)  if (v2[j] == k+1) Ny[k]++;
	}
	for (k = 0; k < 3; k++) {
		Nx[k] =0;
		for (j=0; j< *n; j++) {
			if (v1[j] == k+1) {
				Nx[k]++;
				for (i=0; i< 3; i++) if (v2[j] == i+1) Nxy[ i*3 +k]++;
			}
		}
	}
	

	*MI = *n * log(*n);
	for (k = 0; k < 9; k++) if (Nxy[k]>0) *MI = *MI + (Nxy[k] * log(Nxy[k])) ;

	for (k = 0; k < 3; k++) {
		if (Ny[k]>0) *MI = *MI - (Ny[k] *log(Ny[k])) ;
		if (Nx[k]>0) *MI =*MI - (Nx[k] *log(Nx[k])) ;
	}

	*MI = *MI / (double) *n;
	*MI = *MI/ log(2);

}

