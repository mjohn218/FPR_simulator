#include "cell_neighbor_lists.h"

void cell_neighbor_listPBC(int Nx, int Ny, int Nz, int maxnbor, int *nbor, int *nborrev) {

	int i, j, k;
	int Nsofar = 0;
	int ix, iy, iz;
	int Ncell = Nx * Ny * Nz;
	int *Nrev = new int[Ncell]; //for keeping track of the other 13 neighbors that are ~back and down.
	int c = 0;
	int mybin;

	for (k = 0; k < Nz; k++) {
		for (j = 0; j < Ny; j++) {
			for (i = 0; i < Nx; i++) {
				Nsofar = 0;
				/*For each cell figure out its 13 neighbors that are ~forward and up*/
				if (i == (Nx - 1))
					ix = 0;
				else
					ix = i + 1;
				//count all the plus x boxes
				mybin = ix + j * Nx + k * (Nx * Ny);
				nbor[c * maxnbor + Nsofar] = mybin;
				nborrev[mybin * maxnbor + Nrev[mybin]] = c; //in case you need to know who also has your cell as a partner
				Nrev[mybin]++;
				Nsofar++;
				if (j == (Ny - 1))
					iy = 0;
				else
					iy = j + 1;
				mybin = ix + iy * Nx + k * (Nx * Ny);
				nbor[c * maxnbor + Nsofar] = mybin;
				nborrev[mybin * maxnbor + Nrev[mybin]] = c;
				Nrev[mybin]++;

				Nsofar++;

				if (k == (Nz - 1))
					iz = 0;
				else
					iz = k + 1;
				mybin = ix + iy * Nx + iz * (Nx * Ny);
				nbor[c * maxnbor + Nsofar] = mybin;
				nborrev[mybin * maxnbor + Nrev[mybin]] = c;
				Nrev[mybin]++;

				Nsofar++;

				if (k == (Nz - 1))
					iz = 0;
				else
					iz = k + 1;
				mybin = ix + j * Nx + iz * (Nx * Ny);
				nbor[c * maxnbor + Nsofar] = mybin;
				nborrev[mybin * maxnbor + Nrev[mybin]] = c;
				Nrev[mybin]++;

				Nsofar++;

				if (j == (Ny - 1))
					iy = 0;
				else
					iy = j + 1;
				mybin = i + iy * Nx + k * (Nx * Ny);
				nbor[c * maxnbor + Nsofar] = mybin;
				nborrev[mybin * maxnbor + Nrev[mybin]] = c;
				Nrev[mybin]++;

				Nsofar++;
				if (i > 0)
					ix = i - 1;
				else
					ix = Nx - 1;

				mybin = ix + iy * Nx + k * (Nx * Ny);
				nbor[c * maxnbor + Nsofar] = mybin;
				nborrev[mybin * maxnbor + Nrev[mybin]] = c;
				Nrev[mybin]++;

				Nsofar++;
				if (k == (Nz - 1))
					iz = 0;
				else
					iz = k + 1;
				mybin = ix + iy * Nx + iz * (Nx * Ny);
				nbor[c * maxnbor + Nsofar] = mybin;
				nborrev[mybin * maxnbor + Nrev[mybin]] = c;
				Nrev[mybin]++;

				Nsofar++;

				if (k == (Nz - 1))
					iz = 0;
				else
					iz = k + 1;
				mybin = i + iy * Nx + iz * (Nx * Ny);
				nbor[c * maxnbor + Nsofar] = mybin;
				nborrev[mybin * maxnbor + Nrev[mybin]] = c;
				Nrev[mybin]++;

				Nsofar++;

				if (k == (Nz - 1))
					iz = 0;
				else
					iz = k + 1;
				mybin = i + j * Nx + iz * (Nx * Ny);
				nbor[c * maxnbor + Nsofar] = mybin;
				nborrev[mybin * maxnbor + Nrev[mybin]] = c;
				Nrev[mybin]++;

				Nsofar++;

				if (i > 0)
					ix = i - 1;
				else
					ix = Nx - 1;
				mybin = ix + j * Nx + iz * (Nx * Ny);
				nbor[c * maxnbor + Nsofar] = mybin;
				nborrev[mybin * maxnbor + Nrev[mybin]] = c;
				Nrev[mybin]++;

				Nsofar++;

				if (j > 0)
					iy = j - 1;
				else
					iy = Ny - 1;

				mybin = ix + iy * Nx + iz * (Nx * Ny);
				nbor[c * maxnbor + Nsofar] = mybin;
				nborrev[mybin * maxnbor + Nrev[mybin]] = c;
				Nrev[mybin]++;

				Nsofar++;

				if (j > 0)
					iy = j - 1;
				else
					iy = Ny - 1;
				mybin = i + iy * Nx + iz * (Nx * Ny);
				nbor[c * maxnbor + Nsofar] = mybin;
				nborrev[mybin * maxnbor + Nrev[mybin]] = c;
				Nrev[mybin]++;

				Nsofar++;
				if (i == (Nx - 1))
					ix = 0;
				else
					ix = i + 1;
				mybin = ix + iy * Nx + iz * (Nx * Ny);
				nbor[c * maxnbor + Nsofar] = mybin;
				nborrev[mybin * maxnbor + Nrev[mybin]] = c;
				Nrev[mybin]++;

				Nsofar++;

				// cout <<"Cell: "<<c<<" Num neighbors: "<<Nsofar<<'\t';
// 	for(int sn=0;sn<Nsofar;sn++)
// 	  cout <<nbor[c*maxnbor+sn]<<" ";
// 	cout <<endl;

				c++;

			} //end looping over z cells
		} //end looping over y cells
	} //end looping over x cells

	delete[] Nrev;

}
