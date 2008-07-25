
#include <stdio.h>
#include "mpi.h"
#define n 10
#define tam_sum 5
#define imprimir 1
#include <math.h>

int coords[3], dims[3], periods[3];
MPI_Comm comm_3d, comm_col, comm_fil, comm_reduccion;
int id3D, tag =99;

static float subm_A[tam_sum][tam_sum];
static float subm_B[tam_sum][tam_sum];
static float subm_C[tam_sum][tam_sum];
static float subm_C_Plano0[tam_sum][tam_sum];
static float subm_C_aux[tam_sum][tam_sum];

//float A[n][n], B[n][n], C[n][n];


//void llenarMatriz(float m[n][n])
//{
//  static float k=0;
//  int i, j;
//  for (i=0; i<n; i++)
//    for (j=0; j<n; j++)
//      m[i][j] = k++;
//}
//
void imprimirMatriz(float m[n][n])
{
  int i, j = 0;
  for (i=0; i<n; i++) {
    printf("\n\t| ");
    for (j=0; j<n; j++)
      printf("%2f ", m[i][j]);
    printf("|");
  }
  printf("\n");
}
void imprimirSubMatriz(float m[tam_sum][tam_sum])
{
  int i, j = 0;
  for (i=0; i<tam_sum; i++) {
    printf("\n\t| ");
    for (j=0; j<tam_sum; j++)
      printf("%2f ", m[i][j]);
    printf("|");
  }
  printf("\n");
}
int main(int argc, char** argv) {
    int mi_fila, mi_columna, mi_plano, fila_recepcion, col_recepcion;
    int coords_envio[3], coords_recepcion[3], vector_logico[3];
    int rank_envio,size, valor_matriz;
    double timeIni, timeFin;
    int i,j,k,l, cont_fila, cont_columna;
    MPI_Status statusA;
    MPI_Status statusB;
    MPI_Status statusC;

	MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

	//dimension del nro de procesos(cantidad de bloques)
    double m = pow((double)size, (double)1/3);

    if (n % (int)m !=0 )
    {
        printf("Por favor corra con una cantidad de procesos multiplo de %d.\n", n);fflush(stdout);
		//printf(" mmmmmmmmmm %f\n", m);fflush(stdout);
		//printf(" mmmmmmmmmm %i\n", n);fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int tam_subM = n/m;


    //printf("EL VALOR M ES: %d",(int) m);
    //printf("EL VALOR TAMSUBM ES: %d",(int) tam_subM);
    dims[0]=dims[1]=dims[2]=(int) m;
    periods[0]=periods[1]=periods[2]= 1;


    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &comm_3d);



    //printf(" cart create \n");
    /* Obtiene mi nuevo id en 3D */
    MPI_Comm_rank(comm_3d, &id3D);
    //printf(" comm rank \n");
    /* Obtiene mis coordenadas */
    MPI_Cart_coords(comm_3d, id3D, 3, coords);
   // printf(" CART COORDS\n");


    mi_fila = coords[0];
    mi_columna= coords[1];
    mi_plano = coords[2];

   // printf(" mi fila %d \n", mi_fila);
    //printf(" mi columna %d \n", mi_columna);
   // printf(" mi plano %d \n", mi_plano);
   // printf(" MI RANKING %d \n", id3D);

    /*inicializamos submatrices C*/
    for (i=0; i<tam_subM; i++){
        for (j=0; j<tam_subM; j++){
            subm_C[i][j] =0;
            //subm_C_Plano0[i][j] =0;
        }
    }

    timeIni = MPI_Wtime();

    /*AHORA VERIFICAMOS SI SOMOS del PLANO 0 PARA CREAR LAS MATRICES Y EMPEZAR
    A DISTRIBUIR A LOS DEMAS PISOS*/

	/*A PARTIR DE AK VERIFICAMOS LOS PLANOS EN LOS QUE ESTAMOS Y HACEMOS LA OPERACION CORRECTA, ES DECIR
	SI ESTAMOS EN EL PLANO 0, ENVIAMOS A LOS DEMAS PLANOS LAS PARTES DE A Y B QUE CORRESPONDEN.*/

	if(mi_plano == 0){
	    valor_matriz=5;

        for (k=0; k<tam_subM; k++){
            for (l=0;l < tam_subM; l++){
                subm_A[k][l]=valor_matriz;
                subm_B[k][l]=valor_matriz;
                valor_matriz++;
            }
            valor_matriz =5;
        }
	    /*AQUI LO QUE HACEMOS ES ENVIAR LOS SUB-BLOQUES DE A, DE ACUERDO A LAS
	    COLUMNAS A DONDE CORRESPONDA Pijk = Aijj.*/
        if (mi_columna != 0) {
            coords_envio[0] = mi_fila;
            coords_envio[1] = mi_columna;
            coords_envio[2] = mi_columna;
            MPI_Cart_rank(comm_3d, coords_envio, &rank_envio);
            MPI_Send(subm_A, tam_subM*tam_subM, MPI_FLOAT, rank_envio, 99, comm_3d);
        }


        /*AQUI LO QUE HACEMOS ES ENVIAR LOS SUB-BLOQUES DE B, DE ACUERDO A LAS
	    FILAS A DONDE CORRESPONDA Pijk = Biji.*/
        if (mi_fila != 0) {
            coords_envio[0] = mi_fila;
            coords_envio[1] = mi_columna;
            coords_envio[2] = mi_fila;
            MPI_Cart_rank(comm_3d, coords_envio, &rank_envio);
            MPI_Send(subm_B, tam_subM*tam_subM, MPI_FLOAT, rank_envio, 100, comm_3d);
        }


    /*SI ESTOy EN OTRO PLANO QUIERE DECIR Q TENGO QUE RECIBIR EL BLOQUE Q SE ME ESTA ENVIANDO*/
	}else{
	    /*RECIBIMOS SIEMPRE DE UN PROCESO DEL PLANO 0*/
	    coords_recepcion[2] = 0;
	    coords_recepcion[0] = mi_fila;
	    coords_recepcion[1] = mi_columna;
	    MPI_Cart_rank(comm_3d, coords_recepcion, &rank_envio);
        if (mi_columna == mi_plano) {
            MPI_Recv(subm_A, tam_subM*tam_subM, MPI_FLOAT, rank_envio, 99, comm_3d, &statusA);
        }
        if (mi_fila == mi_plano) {
            MPI_Recv(subm_B, tam_subM*tam_subM, MPI_FLOAT, rank_envio, 100, comm_3d, &statusB);
        }
	}


	/*LUEGO DE HABER RECIBIDO NUESTRA PARTE, LO QUE DEBEMOS HACER EN EL CASO DE A ES
            REPARTIR MEDIANTE UN BROADCAST A TODA MI FILA*/
    vector_logico[0] = 0;
    vector_logico[1] = 1;
    vector_logico[2] = 0;
 //   printf("ENTRE PARA REPARTIR A \n");
    MPI_Cart_sub(comm_3d, vector_logico, &comm_col);
    MPI_Bcast(subm_A, tam_subM*tam_subM, MPI_FLOAT, mi_plano, comm_col);

    vector_logico[0] = 1;
    vector_logico[1] = 0;
    vector_logico[2] = 0;
//    printf("ENTRE PARA REPARTIR B \n");
    MPI_Cart_sub(comm_3d, vector_logico, &comm_fil);
    MPI_Bcast(subm_B, tam_subM*tam_subM, MPI_FLOAT, mi_plano, comm_fil);


    /*AHORA REALIZAMOS EL COMPUTO DE CADA SUB-MATRIZ QUE RECIBIMOS*/
    for (i = 0; i < tam_subM; i++) {
        for (j = 0; j < tam_subM; j++) {
            for (k = 0;  k < tam_subM; k++) {
                subm_C[i][j] += subm_A[i][k] * subm_B[k][j];
            }
        }
    }
printf("\n\n-------------------------------------------------- %i\n",tam_subM);fflush(stdout);
    /*UNA VEZ COMPUTADO EL CALCULO DE LAS SUBMATRICES PROCEDEMOS A REDUCIR
    TODO AL PLANO 0, CREANDO UNA SUBTOPOLOGÍA PARA HACER, Y REALIZAMOS LA SUMA DE
    LAS SECCIONES QUE CORRESPONDEN*/
    vector_logico[0] = 0;
    vector_logico[1] = 0;
    vector_logico[2] = 1;
    MPI_Cart_sub(comm_3d, vector_logico, &comm_reduccion);
    subm_C_Plano0[tam_subM][tam_subM];
    for (i=0; i<tam_subM; i++){
        for (j=0; j<tam_subM; j++){
            subm_C_Plano0[i][j] =0;
        }
    }
//printf("\n\n-------------------------------------------------- %i\n", id3D);fflush(stdout);

    MPI_Reduce(subm_C, subm_C_Plano0, tam_subM*tam_subM, MPI_FLOAT, MPI_SUM, 0, comm_reduccion);

//    printf("RECIBIDO EN A[%d][%d][%d] --> %2f\n",mi_fila,mi_columna,mi_plano,subm_A[0][0]);
//	printf("RECIBIDO EN A[%d][%d][%d] --> %2f\n",mi_fila,mi_columna,mi_plano,subm_A[(tam_subM-1)][(tam_subM-1)]);
//	printf("RECIBIDO EN B[%d][%d][%d] --> %2f\n",mi_fila,mi_columna,mi_plano,subm_B[0][0]);
//	printf("RECIBIDO EN B[%d][%d][%d] --> %2f\n",mi_fila,mi_columna,mi_plano,subm_B[(tam_subM-1)][(tam_subM-1)]);
    //MPI_Barrier(comm_3d);
    if(imprimir ==1 && mi_plano ==0){
        MPI_Send(subm_C_Plano0, tam_subM*tam_subM, MPI_FLOAT, 0, id3D, comm_3d);
        //printf("Matriz Resultado plano 0 P[%d][%d] \n",mi_fila,mi_columna);
        //imprimirSubMatriz(subm_C_Plano0);
        printf("\n");
    }
    MPI_Barrier(comm_3d);
	if(id3D == 0){
        //printf("RECIBIDO EN C[%d][%d][%d] --> %2f\n",mi_fila,mi_columna,mi_plano,subm_C_Plano0[0][0]);
        //printf("RECIBIDO EN C[%d][%d][%d] --> %2f\n",mi_fila,mi_columna,mi_plano,subm_C_Plano0[(tam_subM-1)][(tam_subM-1)]);
        /*GENERO DE VUELTA LAS MATRICES A Y B*/
        float A[n][n];
        float B[n][n];
        float C[n][n];


        int aux = tam_subM;
        if(imprimir == 1){
            for (i=0; i<n; i++){
                for (j=0; j<n; j++){
                    C[i][j] =0;
                }
            }
            int m_int = (int) m;
            for(i=0; i < m_int*m_int; i++){
                MPI_Recv(subm_C_aux,tam_subM*tam_subM, MPI_FLOAT,MPI_ANY_SOURCE,MPI_ANY_TAG,comm_3d,&statusC);
                MPI_Cart_coords(comm_3d, statusC.MPI_TAG ,3, coords_envio);
                //printf("RECIBIMOS DE [%d][%d] -->\n",coords_envio[0],coords_envio[1]);
                fila_recepcion = coords_envio[0];
                col_recepcion = coords_envio[1];
                cont_columna =0;
                cont_fila =-1;
                for(j=fila_recepcion*tam_subM; j<fila_recepcion*tam_subM+tam_subM; j++){
                    cont_fila++;
                    for(k=col_recepcion*tam_subM; k<col_recepcion*tam_subM+tam_subM; k++){
                        C[j][k]= subm_C_aux[cont_fila][cont_columna];
                        cont_columna++;
                    }
                    cont_columna=0;
                }
            }
        }
        valor_matriz=5;
        for (k=0; k<n; k++){
            for (l=0;l < n; l++){
                if(l == aux){
                    valor_matriz =5;
                    //aux = aux+aux;
                }
                A[k][l]=valor_matriz;
                B[k][l]=valor_matriz;
                valor_matriz++;
            }
            valor_matriz =5;
        }

        if(imprimir ==1){
            printf("MATRIZ A/n");
            imprimirMatriz(A);
            printf("MATRIZ B/n");
            imprimirMatriz(B);
            printf("MATRIZ C/n");
            imprimirMatriz(C);
        }
        timeFin = MPI_Wtime();
        printf("TIEMPO TARDADO---> %f segundos\n", timeFin-timeIni);

	}

	MPI_Finalize();
}
