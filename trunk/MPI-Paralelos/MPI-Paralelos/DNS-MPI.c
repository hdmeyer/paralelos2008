/*Se crea una grilla (mesh) l�gica de procesos implementada
con la topolog�a cartesiana de MPI (MPI_Cart_create) en 3 dimensiones.
Estos procesos son etiquetados desde P0,0,0 a Pn,n,n en un nuevo comunicador (comm_3D)
creado por la topolog�a cartesiana */

#include <stdio.h>
#include "mpi.h"
#define n 4
#include <math.h>

int coords[3], dims[3], periods[3];
MPI_Comm comm_3d;
int id3D, tag =99;
float A[n][n], B[n][n], C[n][n];


void llenarMatriz(float m[n][n])
{
  static float k=0;
  int i, j;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      m[i][j] = k++;
}

void imprimirMatriz(float m[n][n])
{
  int i, j = 0;
  for (i=0; i<n; i++) {
    printf("\n\t| ");
    for (j=0; j<n; j++)
      printf("%2f ", m[i][j]);
    printf("|");
  }
}

int main(int argc, char** argv) {
    int mi_fila, mi_columna, mi_plano;
    int coords_envio[3];
    int rank_envio,size;
    int i,j,k,l, cont_fila, cont_columna;
    MPI_Status statusA;
    MPI_Status statusB;
    //int buffA[1], buffB[1];


    MPI_Init(&argc, &argv);
    /* Cantidad de fil, col y plano = numero de procesos */
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //dimension del nro de procesos(cantidad de bloques)
    double m = pow((double)size,(double) 1/3);

    if (n % (int)m !=0 )
    {
        printf("Por favor corra con una cantidad de procesos multiplo de %d.\n", n);fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int tam_subM = n/m;
    printf("EL VALOR M ES: %d",(int) m);
    dims[0]=dims[1]=dims[2]=(int) m;
    periods[0]=periods[1]=periods[2]= 1;

    float subm_A[tam_subM][tam_subM];
    float subm_B[tam_subM][tam_subM];

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
    printf(" MI RANKING %d \n", id3D);


    /*AHORA VERIFICAMOS SI SOMOS EL P[0][0][0] PARA CREAR LAS MATRICES Y EMPEZAR
    A DISTRIBUIR A LOS DEMAS PISOS*/
	if(mi_fila == 0 && mi_columna == 0 && mi_plano == 0){
	    printf("PROCESO PADRE \n");

        llenarMatriz(A);
        llenarMatriz(B);

        imprimirMatriz(A);
        printf("\n");
        imprimirMatriz(B);
        /*Aqui basicamente lo que hacemos es enviar a cada proceso del plano cero
        una parte de A y B que es la que les corresponde, donde Pi,j,o tiene
        A[i][j] y B[i][j]*/

		for (i=0; i < m; i++){
		    for (j=0; j<m; j++){
                if(!(i==0 && j==0)){
                    coords_envio[0] = i;
                    coords_envio[1] = j;
                    coords_envio[2] = 0;
                    cont_fila=-1;
                    for (k=i*tam_subM; k < i*tam_subM+tam_subM; k++){
                        cont_fila++;
                        cont_columna=0;
                        for (l=j*tam_subM; l < j*tam_subM+tam_subM; l++){
                            subm_A[cont_fila][cont_columna]=A[k][l];
                            subm_B[cont_fila][cont_columna]=B[k][l];
                            cont_columna++;
                        }
                    }
                    MPI_Cart_rank(comm_3d, coords_envio, &rank_envio);
                    printf("RANKING AL Q ENVIO %d \n",rank_envio);
                    //ACA ENVIAMOS A CADA PROCESO CORRESPONDIENTE DEL PLANO
                    //0(DISTRIBUCI�N N^2
                    MPI_Send(subm_A, tam_subM*tam_subM, MPI_FLOAT, rank_envio, 1, comm_3d);
                    MPI_Send(subm_B, tam_subM*tam_subM, MPI_FLOAT, rank_envio, 2, comm_3d);
                }
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
		MPI_Recv(subm_A, tam_subM*tam_subM, MPI_FLOAT, 0, 1, comm_3d,&statusA);
		MPI_Recv(subm_B, tam_subM*tam_subM, MPI_FLOAT, 0, 2, comm_3d,&statusB);

		printf("RECIBIDO EN A[%d][%d][%d] --> %2f\n",mi_fila,mi_columna,mi_plano,subm_A[0][0]);
		printf("RECIBIDO EN A[%d][%d][%d] --> %2f\n",mi_fila,mi_columna,mi_plano,subm_A[(tam_subM-1)][(tam_subM-1)]);
		printf("RECIBIDO EN B[%d][%d][%d] --> %2f\n",mi_fila,mi_columna,mi_plano,subm_B[0][0]);
		printf("RECIBIDO EN B[%d][%d][%d] --> %2f\n",mi_fila,mi_columna,mi_plano,subm_B[(tam_subM-1)][(tam_subM-1)]);
//		psum=0;
//		for (i=0;i<chunksize;i++) psum = psum +a[i];
//		MPI_Send(&psum, 1, MPI_DOUBLE, nproc-1, tag, MPI_COMM_WORLD);
//
//        MPI_Send(&buffA, 1, MPI_INT, rank_envio, 99, comm_3d);
	}
	MPI_Finalize();
}
