
/*********************************************************************************************
 * Archivo: mm2d-mpi.c
 *
 * Multiplicacion de matrices utilizando la Interfaz de Paso de Mensajes MPI
 * y el particionamiento 2-D ciclico.
 *
 * Adaptacion del codigo: sam_mm.c
 * Disponible en: http://www.cs.ucsb.edu/~tyang/class/240b99/homework/sam_mm.c
 *
 * Para compilar:
 * cc mm2d-mpi.c -o mm2d-mpi -lmpi
 *
 * Adaptado por:
 * - Eduardo Rivas. erivas17@gmail.com
 * - Hugo Meyer. meyer.hugo@gmail.com
 *
 * Alumnos de la Universidad Nacional de Asuncion - Facultad Politecnica.
 * Carrera: Ingenieria Informatica.
 * Materia: Electiva V - Algoritmos Paralelos.
 * Profesor: Cristian Von Lucken.
 * Anho: 2008.
 **********************************************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <sys/time.h>
# include <mpi.h>

# define N 10		/* tamanho de la matriz */
# define SUBN 2		/* tamanho de la submatriz (o bloques) */


static int subnum;	 /* numero de sub-matrices en una fila/columna de la matriz */
static int nproc;	 /* numero de nodos MPI */
static int myid; 	 /* mi propio rank */

static float BUFFER_A[SUBN][N]; /* buffer utilizado para la distribucion de sub-matrices de A
                                 * necesarias para calcular un determinado bloque de la matriz C */

static float BUFFER_B[N][SUBN]; /* buffer utilizado para la distribucion de sub-matrices de B
                                 * necesarias para calcular un determinado bloque de la matriz C */

static float BUFFER_C[SUBN][SUBN]; /* buffer utilizado para almacenar el resulatdo de la multiplicacion
                                    * del bloque correspondiente a C */

struct indices    /* Tipo de dato para representar */
{	int i, j;     /* los ndices de una matriz */
};
typedef struct indices indice;



/**
 * Esta funcion mapea una tarea a un indice de la sub-matriz
 */
indice getIndice(int tarea)
{
	indice indi;
	indi.i = tarea / subnum;
	indi.j = tarea % subnum;
	return indi;
}

/**
 * Imprime el contenido de una matriz
 */
void imprimirMatriz(float m[N][N])
{
	int i, j = 0;
	for (i=0; i < N; i++) {
		printf("\n\t| ");
		for (j=0; j < N; j++)
			printf("%4f ", m[i][j]);
		printf("|");
	}
}

/**
 * Imprime el contenido del buffer A
 */
void imprimirBufferA()
{
	int i, j = 0;
	for (i=0; i < SUBN; i++) {
		printf("\n\t| ");
		for (j=0; j < N; j++)
			printf("%4f ", BUFFER_A[i][j]);
		printf("|");
	}
}

/**
 * Imprime el contenido del buffer B
 */
void imprimirBufferB()
{
	int i, j = 0;
	for (i=0; i < N; i++) {
		printf("\n\t| ");
		for (j=0; j < SUBN; j++)
			printf("%4f ", BUFFER_B[i][j]);
		printf("|");
	}
}

/**
 * Imprime el contenido del buffer C
 */
void imprimirBufferC()
{
	int i, j = 0;
	for (i=0; i < SUBN; i++) {
		printf("\n\t| ");
		for (j=0; j < SUBN; j++)
			printf("%4f ", BUFFER_C[i][j]);
		printf("|");
	}
}

/**
 * Inicializa las matrices en forma generica
 */
void inicializarMatrices(float A[N][N], float B[N][N], float C[N][N])
{
    int i, j;
	for (i=0; i < N; i++){
		for (j=0; j < N; j++) {
			A[i][j] = (j+3) * (i+1);
			B[i][j] = i+j;
			C[i][j]= 0;
		}
	}
}


/**
 * Rellena los buffers con los valores necesarios para calcular
 * el bloque correspondiente a la matriz resultado C
 * Parametros: A y B son las matrices de entrada
 *   indi: indica la primera posicion del bloque a calcular
 */
void cargarBuffers(float A[N][N], float B[N][N], indice indi){

    indice indiceA, indiceB;
	indiceA.i = indi.i * SUBN;
	indiceA.j = 0;
	indiceB.i = 0;
	indiceB.j = indi.j * SUBN;

    if ( N % SUBN != 0){
        /* Sobran bloques de tamanhos distintos a SUBN */
        if ( indiceA.i + SUBN > N )
            indiceA.i = N - SUBN;
        if ( indiceB.j + SUBN > N )
            indiceB.j = N - SUBN;
    }

	int i, j;
	/* Rellenamos el buufer A */
	for (i = 0; i < SUBN ; i++){
		for (j = 0; j < N; j++){
			BUFFER_A[i][j] = A[i+indiceA.i][j+indiceA.j];
		}
	}

	/* Rellenamos el buufer B */
	for (i = 0; i < N ; i++){
		for (j = 0; j < SUBN; j++){
			BUFFER_B[i][j] = B[i+indiceB.i][j+indiceB.j];
		}
	}
}

/**
 * Inicializa el buffer C
 */
void inicializarBufferC(float BUFFER_C[SUBN][SUBN]){

	int i, j;
	for (i = 0; i < SUBN ; i++)
		for (j = 0; j < SUBN; j++)
			BUFFER_C[i][j] = 0;
}

/**
 * Asigna los datos obtenidos en el BUFFER_C
 * a la matriz resultado C
 */
void asignarBUFFEResultado(float BUFFER_C[SUBN][SUBN], float C[N][N], indice indi){

     indice indiceC;
     indiceC.i = indi.i * SUBN;
	 indiceC.j = indi.j * SUBN;

     if ( N % SUBN != 0){
         /* Sobran bloques de tamanhos distintos a SUBN */
         if ( indiceC.i + SUBN > N )
             indiceC.i = N - SUBN;
         if ( indiceC.j + SUBN > N )
             indiceC.j = N - SUBN;
     }

     int i, j;
     for (i = 0; i < SUBN; i++)
         for (j = 0; j < SUBN; j++)
             C[i+indiceC.i][j+indiceC.j] = BUFFER_C[i][j];

}

/**
 * Realiza la multiplicacion correspondiente a un
 * determinado bloque de la matriz resultado C
 * Utiliza los BUFFER_A y BUFFER_B como parametros de entrada y
 * el resultado lo almacena en BUFFER_C
 */

void multiplicarSubMatriz(float BUFFER_A[SUBN][N], float BUFFER_B[N][SUBN], float BUFFER_C[SUBN][SUBN]){

    inicializarBufferC(BUFFER_C);

    int i, j, k;
    for (i = 0; i < SUBN; i++){
        for (j = 0; j < SUBN; j++){
            for (k = 0; k < N; k++){
                BUFFER_C[i] [j] += BUFFER_A[i] [k] * BUFFER_B[k] [j];
            }
        }
    }
}

int main(int argc, char *argv[])
{
	int i, j, k, token;
	double t1, t2;

    MPI_Status status;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	subnum = N / SUBN;
	if (subnum * SUBN != N ) {
		subnum++;
	}

    if (nproc > SUBN*SUBN)
        nproc = SUBN*SUBN; /* Maxima cantidad de nodos necesarios */

	if ( myid == 0){
		float A[N][N], B[N][N], C[N][N];
		int id_destino = 0;

        /*Se inicializan las matrices */
		printf("Inicializacion ...\n");
		inicializarMatrices(A, B, C);
		printf("Ok!!\n\n");

        /* Se realiza la multiplicacion de matrices */
		printf("Multiplicacion ...\n");
		t1 = MPI_Wtime();

		int tarea; /* Cantidad total de tareas disponibles */
		for (tarea = 0; tarea < subnum*subnum; tarea++ ){

            /* Se carga los buffers con los datos necesarios para cada terea */
			indice indi = getIndice(tarea);
			cargarBuffers(A, B, indi);
            id_destino = tarea % nproc;

            if ( id_destino != 0 ){
                /* Enviar las sub-matrices a multiplicar */
                MPI_Send(BUFFER_A, SUBN*N, MPI_FLOAT, id_destino, 1, MPI_COMM_WORLD);
                MPI_Send(BUFFER_B, N*SUBN, MPI_FLOAT, id_destino, 2, MPI_COMM_WORLD);
            }
            else{
                /* Se multiplican localmente las sub-matrices */
                multiplicarSubMatriz(BUFFER_A, BUFFER_B, BUFFER_C);
                asignarBUFFEResultado(BUFFER_C, C, indi);
            }

            if ( id_destino == nproc-1 ){
                /* Se reciben todas las tareas distribuidas hasta el momento */
                int tarea_recibida = tarea - (nproc-1);
                int proc;
                for (proc = 0; proc <= nproc-1; proc++){
                    id_destino = tarea_recibida % nproc;
                    if (id_destino != 0){
                        /* Si la tarea no corresponde al root (id=0) entonces se recibe */
                        MPI_Recv(BUFFER_C, SUBN*SUBN, MPI_FLOAT, id_destino,3, MPI_COMM_WORLD, &status);
                        indi = getIndice(tarea_recibida);
                        asignarBUFFEResultado(BUFFER_C, C, indi);
                    }
                    tarea_recibida++;
                }
            }
		}

        if ( (nproc < SUBN*SUBN) && (((SUBN*SUBN) % nproc) != 0) ){
            /* Se reciben las ultimas tareas pendientes */

            int tareas_pendientes = SUBN*SUBN % nproc; // cuantas tareas aun estan pendientes??
            int tarea_recibida = SUBN*SUBN - tareas_pendientes; // a partir de que tarea??
            int x;
            for (x = 0; x < tareas_pendientes; x++){
                id_destino = tarea_recibida % nproc;
                if (id_destino != 0){
                    /* Si la tarea no corresponde al root (id=0) entonces se recibe */
                    MPI_Recv(BUFFER_C, SUBN*SUBN, MPI_FLOAT, id_destino,3, MPI_COMM_WORLD, &status);
                    indice indi = getIndice(tarea_recibida);
                    asignarBUFFEResultado(BUFFER_C, C, indi);
                }
                tarea_recibida++;
            }
        }

		t2 = MPI_Wtime();
		/* Fin de la multiplicacion de matrices */
		printf("Ok!!\n\n");

		if ( N <= 10) {
		    /* Si el tamanho de la matriz es pequenha entonces se imprime el resultado*/
			printf("Imprimiendo Matriz A...\n");
			imprimirMatriz(A);
			printf("ok!!\n\n");
			printf("Imprimiendo Matriz B...\n");
			imprimirMatriz(B);
			printf("ok!!\n\n");
			printf("Imprimiendo Matriz C...\n");
			imprimirMatriz(C);
			printf("ok!!\n\n");
		}
		printf("\nDuracion total de la multilplicacion de matrices %4f segundos\n", t2-t1);
	}
	else{

		int tarea;
		int id_destino;
		for (tarea = 0; tarea < subnum*subnum; tarea++ ){

			indice indi = getIndice(tarea);
            id_destino = tarea % nproc;

            if ( id_destino == myid ){
                /* Se reciben los bloques de entrada, se multiplican y se envia el resultado */
                MPI_Recv(BUFFER_A, SUBN*N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(BUFFER_B, N*SUBN, MPI_FLOAT, 0, 2, MPI_COMM_WORLD, &status);
                multiplicarSubMatriz(BUFFER_A, BUFFER_B, BUFFER_C);
                MPI_Send(BUFFER_C, SUBN*SUBN, MPI_FLOAT, 0, 3, MPI_COMM_WORLD);
            }
		}
	}

    MPI_Finalize();
	return 0;
}

