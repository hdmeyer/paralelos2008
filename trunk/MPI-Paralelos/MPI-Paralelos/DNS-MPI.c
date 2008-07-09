/*Se crea una grilla (mesh) lógica de procesos implementada
con la topología cartesiana de MPI (MPI_Cart_create) en 3 dimensiones.
Estos procesos son etiquetados desde P0,0,0 a Pn,n,n en un nuevo comunicador (comm_3D)
creado por la topología cartesiana */

#include <stdio.h>
#include "mpi.h"
#define n 10
int coords[3], dims[3], periods[3];
MPI_Comm comm_3d;
int id3D, tag =99;
int A[n][n], B[n][n], C[n][n];


void llenarMatriz(int m[n][n])
{
  static int k=0;
  int i, j;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      m[i][j] = k++;
}

void imprimirMatriz(int m[n][n])
{
  int i, j = 0;
  for (i=0; i<n; i++) {
    printf("\n\t| ");
    for (j=0; j<n; j++)
      printf("%2d ", m[i][j]);
    printf("|");
  }
}

int main(int argc, char** argv) {
    int mi_fila, mi_columna, mi_plano;
    int coords_envio[3];
    int rank_envio,size;
    int i,j;
    MPI_Status status;
    int buffA[1], buffB[1];


    MPI_Init(&argc, &argv);
    /* Cantidad de fil, col y plano = numero de procesos */
    dims[0]=dims[1]=dims[2]=n;
    periods[0]=periods[1]=periods[2]= 1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != n*n*n)
    {
        printf("Por favor corra con  %d procesos.\n", n*n*n);fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &comm_3d);


    printf(" cart create \n");
    /* Obtiene mi nuevo id en 3D */
    MPI_Comm_rank(comm_3d, &id3D);
    printf(" comm rank \n");
    /* Obtiene mis coordenadas */
    MPI_Cart_coords(comm_3d, id3D, 3, coords);
    printf(" CART COORDS\n");


    mi_fila = coords[0];
    mi_columna= coords[1];
    mi_plano = coords[2];

    printf(" mi fila %d \n", mi_fila);
    printf(" mi columna %d \n", mi_columna);
    printf(" mi plano %d \n", mi_plano);


    /*AHORA VERIFICAMOS SI SOMOS EL P[0][0][0] PARA CREAR LAS MATRICES Y EMPEZAR
    A DISTRIBUIR A LOS DEMAS PISOS*/
	if(mi_fila == 0 && mi_columna == 0 && mi_plano == 0){
	    printf("PROCESO PADRE \n");

        llenarMatriz(A);
        llenarMatriz(B);

        imprimirMatriz(A);
        printf("\n");
        imprimirMatriz(B);

		for (i=0; i<n; i++){
		    for (j=0; j<n; j++){
                coords_envio[0] = i;
                coords_envio[1] = j;
                coords_envio[2] = 0;
                MPI_Cart_rank(comm_3d, coords_envio, &rank_envio);
                //ACA ENVIAMOS A CADA PROCESO CORRESPONDIENTE DEL PLANO
                //0(DISTRIBUCIÓN N^2
                MPI_Send(&A[i][j], 1, MPI_INT, rank_envio, 99, comm_3d);
                MPI_Send(&B[i][j], 1, MPI_INT, rank_envio, 100, comm_3d);
		    }
		}



//		for (i=0;i<nproc-1;i++) {
//			first = i*chunksize;
//			MPI_Send(A,chunksize,MPI_INT, i, tag, MPI_COMM_WORLD);
//		}
//		first =my_id*chunksize;i++;
//		chunksize=size-1;
//		sum =0;
//		for (i=first;i<chunksize;i++) sum = sum +a[i];
//
//		for (i=0;i<nproc-1;i++) {
//			MPI_Recv(&psum,1,MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
//			printf("Respondio %d - envio %10.0f\n", status.MPI_SOURCE,psum);
//			sum =sum +psum;
//		}
//		printf("Resultado de la suma = %10.0f\n", sum);


	}else if(mi_plano == 0){
		MPI_Recv(buffA, 1, MPI_INT, 0, 99, comm_3d,&status);
		MPI_Recv(buffB, 1, MPI_INT, 0, 100, comm_3d,&status);

		printf("RECIBIDO EN A[%d][%d][%d] --> %d\n",mi_fila,mi_columna,mi_plano,buffA[0]);
		printf("RECIBIDO EN B[%d][%d][%d] --> %d\n",mi_fila,mi_columna,mi_plano,buffB[0]);
//		psum=0;
//		for (i=0;i<chunksize;i++) psum = psum +a[i];
//		MPI_Send(&psum, 1, MPI_DOUBLE, nproc-1, tag, MPI_COMM_WORLD);

        MPI_Send(&buffA, 1, MPI_INT, rank_envio, 99, comm_3d);
	}
	MPI_Finalize();
}
