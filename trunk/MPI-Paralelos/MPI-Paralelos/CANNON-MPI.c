#include <stdio.h>
#include "mpi.h"
#define n 9
#include <math.h>

int mi_id;
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
void imprimirSubMatriz(float m[3][3])
{
  int i, j = 0;
  for (i=0; i<3; i++) {
    printf("\n\t| ");
    for (j=0; j<3; j++)
      printf("%2f ", m[i][j]);
    printf("|");
  }
}
int main(int argc, char** argv) {
    int mi_fila, mi_columna;
    int rank_envio,size;
    int i,j,k,l;
    int m;              /*tamaño de cada submatriz*/
    MPI_Status statusA;
    MPI_Status statusB;

	MPI_Comm comm2d;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&mi_id);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	double m = pow((double)size,(double) 1/2);

	int tam_subM = n/m;
    printf("EL VALOR M ES: %d",(int) m);
    printf("EL VALOR TAMSUBM ES: %d",(int) tam_subM);
    periods[0]=periods[1]=periods[2]= 1;

    float subm_A[tam_subM][tam_subM];
    float subm_B[tam_subM][tam_subM];

}
