
#include <stdio.h>
#include "mpi.h"
#define n 512
#include <math.h>

int coords[3], dims[3], periods[3];
MPI_Comm comm_3d;
int id3D, tag =99;
//float A[n][n], B[n][n], C[n][n];
void dimensionar(float**m, int tam){
    int i;
    m = (float **) malloc ( tam * sizeof(float) );
    for (i = 0; i <= tam; i++) {
		m[i] = (float *) malloc ( tam* sizeof(float) );
 	}
}

void liberar(float **m, int tam){
    int i;
    for (i = 0; i < tam; i++) {
        free(m[i]);
    }
    free(m);
}

int main(int argc, char** argv) {
    int mi_fila, mi_columna, mi_plano;
    int coords_envio[3], coords_recepcion[3];
    int size, valor_matriz;
    double timeIni, timeFin;
    int i,j,k,l, cont_fila, cont_columna;
    MPI_Status statusA;
    MPI_Status statusB;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id3D);

    //dimension del nro de procesos(cantidad de bloques)
    double m = pow((double)size,(double) 1/3);

    if (n % (int)m !=0 )
    {
        printf("Por favor corra con una cantidad de procesos multiplo de %d.\n", n);fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int tam_subM = n/m;
    //printf("EL VALOR M ES: %d",(int) m);
    //printf("EL VALOR TAMSUBM ES: %d",(int) tam_subM);
    dims[0]=dims[1]=dims[2]=(int) m;
    periods[0]=periods[1]=periods[2]= 1;

//    float subm_A[tam_subM][tam_subM];
//    float subm_B[tam_subM][tam_subM];
//    float subm_C[tam_subM][tam_subM];
    float **subm_A;
    float **subm_B;
    float **subm_C;
//    dimensionar(subm_A, tam_subM);
//    dimensionar(subm_B, tam_subM);
//    dimensionar(subm_C, tam_subM);
    //printf("DESPUES DE DIMENSIONAR");

    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &comm_3d);

    //printf(" cart create \n");
    /* Obtiene mi nuevo id en 3D */
    MPI_Comm_rank(comm_3d, &id3D);
    //printf(" comm rank \n");
    /* Obtiene mis coordenadas */
    MPI_Cart_coords(comm_3d, id3D, 3, coords);
    //printf(" CART COORDS\n");


    mi_fila = coords[0];
    mi_columna= coords[1];
    mi_plano = coords[2];

//    printf(" mi fila %d \n", mi_fila);
//    printf(" mi columna %d \n", mi_columna);
//    printf(" mi plano %d \n", mi_plano);
//    printf(" MI RANKING %d \n", id3D);

    /*inicializamos submatrices C*/
    subm_A = (float **) malloc ( tam_subM* sizeof(float) );
    subm_B = (float **) malloc ( tam_subM* sizeof(float) );
    subm_C = (float **) malloc ( tam_subM* sizeof(float) );
    for (i = 0; i < tam_subM; i++) {
		subm_A[i] = (float *) malloc ( tam_subM* sizeof(float) );
		subm_B[i] = (float *) malloc ( tam_subM* sizeof(float) );
		subm_C[i] = (float *) malloc ( tam_subM* sizeof(float) );
 	}
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
	    valor_matriz=1;
        for (k=0; k<tam_subM; k++){
            for (l=0;l < tam_subM; l++){
                subm_A[k][l]=valor_matriz*id3D;
                subm_B[k][l]=valor_matriz*id3D;
                valor_matriz++;
            }
        }
	}

	if(mi_plano == 0 && mi_columna >0){
	    /*AQUI LO QUE HACEMOS ES ENVIAR LOS SUB-BLOQUES DE A, DE ACUERDO A LAS
	    COLUMNAS A DONDE CORRESPONDA Pijk = Aijj.*/
	    int rank_envio;
        coords_envio[0] =coords[0];
        coords_envio[1] = coords[1];
        coords_envio[2] = coords[1];
        MPI_Cart_rank(comm_3d, coords_envio, &rank_envio);
        MPI_Send(subm_A, tam_subM*tam_subM, MPI_FLOAT, rank_envio, 99, comm_3d);
	}
	if (mi_plano > 0 && mi_columna == mi_plano) {
	    int rank_fuente;
	    coords_recepcion[0] = coords[0];
	    coords_recepcion[1] = coords[1];
	    coords_recepcion[2] = 0;
	    MPI_Cart_rank(comm_3d, coords_recepcion, &rank_fuente);
        MPI_Recv(subm_A, tam_subM*tam_subM, MPI_FLOAT, rank_fuente, 99, comm_3d, &statusA);
    }

        /*AQUI LO QUE HACEMOS ES ENVIAR LOS SUB-BLOQUES DE B, DE ACUERDO A LAS
	    FILAS A DONDE CORRESPONDA Pijk = Biji.*/
    if (mi_plano == 0 && mi_fila >0) {
        int rank_envio;
        coords_envio[0] = coords[0];
        coords_envio[1] = coords[1];
        coords_envio[2] = coords[0];
        MPI_Cart_rank(comm_3d, coords_envio, &rank_envio);
        MPI_Send(subm_B, tam_subM*tam_subM, MPI_FLOAT, rank_envio, 100, comm_3d);
    }

    if (mi_plano > 0 && mi_fila == mi_plano) {
	    int rank_fuente;
	    coords_recepcion[0] = coords[0];
	    coords_recepcion[1] = coords[1];
	    coords_recepcion[2] = 0;
	    MPI_Cart_rank(comm_3d, coords_recepcion, &rank_fuente);
        MPI_Recv(subm_B, tam_subM*tam_subM, MPI_FLOAT, rank_fuente, 100, comm_3d, &statusB);
    }


	MPI_Barrier(comm_3d);
	/*LUEGO DE HABER RECIBIDO NUESTRA PARTE, LO QUE DEBEMOS HACER EN EL CASO DE A ES
            REPARTIR MEDIANTE UN BROADCAST A TODA MI FILA*/
    MPI_Comm comm_col, comm_fil, comm_reduccion;
    int vector_logico_col[3];
    vector_logico_col[0] = 0;
    vector_logico_col[1] = 1;
    vector_logico_col[2] = 0;
    //printf("ENTRE PARA REPARTIR A \n");
    MPI_Cart_sub(comm_3d, vector_logico_col, &comm_col);

    MPI_Bcast(subm_A, tam_subM*tam_subM, MPI_FLOAT, mi_plano, comm_col);

    int vector_logico_filas[3];
    vector_logico_filas[0] = 1;
    vector_logico_filas[1] = 0;
    vector_logico_filas[2] = 0;
    //printf("ENTRE PARA REPARTIR B \n");
    MPI_Cart_sub(comm_3d, vector_logico_filas, &comm_fil);

    MPI_Bcast(subm_B, tam_subM*tam_subM, MPI_FLOAT, mi_plano, comm_fil);

    //MPI_Barrier(comm_fil);
    //MPI_Barrier(comm_col);
    /*AHORA REALIZAMOS EL COMPUTO DE CADA SUB-MATRIZ QUE RECIBIMOS*/
    for (i = 0; i < tam_subM; i++) {
        for (j = 0; j < tam_subM; j++) {
            for (k = 0;  k < tam_subM; k++) {
                subm_C[i][j] += subm_A[i][k] * subm_B[k][j];
            }
        }
    }

    /*UNA VEZ COMPUTADO EL CALCULO DE LAS SUBMATRICES PROCEDEMOS A REDUCIR
    TODO AL PLANO 0, CREANDO UNA SUBTOPOLOGÍA PARA HACER, Y REALIZAMOS LA SUMA DE
    LAS SECCIONES QUE CORRESPONDEN*/
    float **subm_C_Plano0;
    subm_C_Plano0 = (float **) malloc ( tam_subM * sizeof(float) );
    for (i = 0; i < tam_subM; i++) {
		subm_C_Plano0[i] = (float *) malloc ( tam_subM* sizeof(float) );
 	}

//    float subm_C_Plano0[tam_subM][tam_subM];
 	for (i=0; i<tam_subM; i++){
        for (j=0; j<tam_subM; j++){
            subm_C_Plano0[i][j] =0;
        }
    }
    int vector_logico_reduc[3];
    vector_logico_reduc[0] = 0;
    vector_logico_reduc[1] = 0;
    vector_logico_reduc[2] = 1;
    MPI_Cart_sub(comm_3d, vector_logico_reduc, &comm_reduccion);

    MPI_Reduce(subm_C, subm_C_Plano0, tam_subM*tam_subM, MPI_FLOAT, MPI_SUM, 0, comm_reduccion);

//    printf("RECIBIDO EN A[%d][%d][%d] --> %2f\n",mi_fila,mi_columna,mi_plano,subm_A[0][0]);
//	printf("RECIBIDO EN A[%d][%d][%d] --> %2f\n",mi_fila,mi_columna,mi_plano,subm_A[(tam_subM-1)][(tam_subM-1)]);
//	printf("RECIBIDO EN B[%d][%d][%d] --> %2f\n",mi_fila,mi_columna,mi_plano,subm_B[0][0]);
//	printf("RECIBIDO EN B[%d][%d][%d] --> %2f\n",mi_fila,mi_columna,mi_plano,subm_B[(tam_subM-1)][(tam_subM-1)]);
    //MPI_Barrier(comm_3d);
    MPI_Comm_free(&comm_col);
    MPI_Comm_free(&comm_fil);
    MPI_Comm_free(&comm_reduccion);

	if(id3D == 0){
        //printf("RECIBIDO EN C[%d][%d][%d] --> %2f\n",mi_fila,mi_columna,mi_plano,subm_C_Plano0[0][0]);
        //printf("RECIBIDO EN C[%d][%d][%d] --> %2f\n",mi_fila,mi_columna,mi_plano,subm_C_Plano0[(tam_subM-1)][(tam_subM-1)]);

        //printf("SUBMATRIZ A/n");
        //imprimirSubMatriz(subm_A);
        //printf("SUBMATRIZ B/n");
        //imprimirSubMatriz(subm_B);
        //printf("SUBMATRIZ C/n");
        //imprimirSubMatriz(subm_C_Plano0);
        timeFin = MPI_Wtime();
        printf("TIEMPO TARDADO---> %f segundos\n", timeFin-timeIni);

	}
//	if(mi_plano == 0){
//        for (i = 0; i < tam_subM; i++) {
//            free(subm_C_Plano0[i]);
//        }
//	}
//	for (i = 0; i < tam_subM; i++) {
//        free(subm_C[i]);
//    }
//    for (i = 0; i < tam_subM; i++) {
//        free(subm_A[i]);
//    }
//    for (i = 0; i < tam_subM; i++) {
//        free(subm_B[i]);
//    }

    MPI_Comm_free(&comm_3d);

	MPI_Finalize();
}
