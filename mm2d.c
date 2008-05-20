
/*********************************************************************************************
 * Archivo: mm2d.c
 *
 * Multiplicacion de matrices utilizando hilos y el particionamiento 2-D 
 *
 * Adaptacion del codigo: sam_mm.c  
 * Disponible en: http://www.cs.ucsb.edu/~tyang/class/240b99/homework/sam_mm.c
 * El cual es una multiplicacion de matrices cuadradas utilizando 
 * el particionamiento 1-D e implementado con la libreria MPI.
 *
 * Para compilar:
 * cc mm2d.c -o mm2d -lpthread
 *
 * Para ejecutar:
 * ./mm2d <cantidad_de_hilos>
 *
 * Adaptacion para su implementacion con hilos, particionamiento 2-D y matrices de tamanhos
 * arbitrarias (multiplicables) por:
 *
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
# include <pthread.h>
# include <sys/time.h>

# define M 10		/* Cantidad de filas de la matriz A */
# define S 5		/* Cantidad de columnas de la matriz A = Cantidad de filas de la matriz B */
# define N 8		/* Cantidad de columnas de la matriz B */
# define SUBN 3		/* Tamanho de la submatriz -- Obs. debe satisfacer la restriccion: SUBN > M && SUBN > S && SUBN > N  */

# define MAX(a,b) (((a) > (b)) ? (a) : (b))

static int subnumM;	 /* numero de sub-matrices en una fila de la matriz A */ 
static int subnumS;	 /* numero de sub-matrices en una columna/fila de las matrices A/B */ 
static int subnumN;	 /* numero de sub-matrices en una columna de la matriz B */ 
static int cant_hilos;	/* numbero total de hilos */

static int A[M][S], B[S][N], C[M][N];

static pthread_mutex_t lock;
static int tarea = -1;	/* indica el id de la tarea actual */


struct indices           /* Tipo de dato para representar */
{	int i, j;	/* los ndices de una matriz */
};
typedef struct indices indice;


/**
 * Obtiene el tiempo actual del sistema
 */
double getTime() {
	struct timeval t;
	gettimeofday(&t, NULL);
	return (double)t.tv_sec+t.tv_usec*0.000001;
}


/**
 * Esta funcion mapea una tarea a un indice de la sub-matriz 
 */
indice getIndice(int tarea)
{
	indice indi;
	indi.i = tarea / subnumN;
	indi.j = tarea % subnumN;
	return indi;
}

/**
 * Imprime el contenido de la matriz A
 */
void imprimirMatrizA()
{
	int i, j = 0;
	for (i=0; i < M; i++) {
		printf("\n\t| ");
		for (j=0; j < S; j++)
			printf("%2d ", A[i][j]);
		printf("|");
	}
}

/**
 * Imprime el contenido de la matriz B
 */
void imprimirMatrizB()
{
	int i, j = 0;
	for (i=0; i < S; i++) {
		printf("\n\t| ");
		for (j=0; j < N; j++)
			printf("%2d ", B[i][j]);
		printf("|");
	}
}

/**
 * Imprime el contenido de la matriz C
 */
void imprimirMatrizC()
{
	int i, j = 0;
	for (i=0; i < M; i++) {
		printf("\n\t| ");
		for (j=0; j < N; j++)
			printf("%2d ", C[i][j]);
		printf("|");
	}
}

/**
 * Inicializa las matrices en forma generica
 */
void inicializarMatrices()
{
	int i, j = 0;
	for (i=0; i < M; i++)
		for (j=0; j < S; j++)
			A[i][j] = (j+3) * (i+1);

	for (i=0; i < S; i++)
		for (j=0; j < N; j++)
			B[i][j] = i+j;

	for (i=0; i < M; i++)
		for (j=0; j < N; j++)
			C[i][j] = 0;
}

/**
 * Realiza la suma de una sola submatriz
 */
void sumarSubMatriz( indice indiceA, indice indiceB, indice indiceC )
{
	int i, j, k;

	/* Primero se calcula el tamanho de la submatrices */
	indice tamanho_submatriz_A, tamanho_submatriz_B;
	tamanho_submatriz_A.i = SUBN;
	tamanho_submatriz_A.j = SUBN;
	tamanho_submatriz_B.i = SUBN;
	tamanho_submatriz_B.j = SUBN;

	if ( indiceA.i + SUBN > M && subnumM * SUBN != M ){ // queda una porcion de bloque desigual en forma vertical (filas)?
		tamanho_submatriz_A.i  = M - (M / SUBN) * SUBN;
	}
	if ( indiceA.j + SUBN > S && subnumS * SUBN != S ){ // queda una porcion de bloque desigual en forma horizontal (columnas)?
		tamanho_submatriz_A.j  = S - (S / SUBN) * SUBN;
	}
	if ( indiceB.i + SUBN > S && subnumS * SUBN != S ){ // queda una porcion de bloque desigual en forma vertical (filas)?
		tamanho_submatriz_B.i  = S - (S / SUBN) * SUBN;
	}
	if ( indiceB.j + SUBN > N && subnumN * SUBN != N ){ // queda una porcion de bloque desigual en forma horizontal (columnas)?
		tamanho_submatriz_B.j  = N - (N / SUBN) * SUBN;
	}

	int tk = tamanho_submatriz_A.j;

	/* Luego se realiza la multiplicacion entre las submatrices */
	for (i=0; i < tamanho_submatriz_A.i; i++)
		for (j=0; j < tamanho_submatriz_B.j; j++)
			for (k=0; k < tk; k++){
				C[indiceC.i + i] [indiceC.j + j] += A[indiceA.i + i] [indiceA.j + k] * B[indiceB.i + k] [indiceB.j + j];
			};
}



/**
 * Realiza la suma de todas las submatrices de una tarea
 */
void realizarTarea(indice indi)
{
	int k;

	indice indiceA, indiceB, indiceC;
	indiceA.i = indi.i * SUBN;
	indiceA.j = 0;
	indiceB.i = 0;
	indiceB.j = indi.j * SUBN;
	indiceC.i = indi.i * SUBN;
	indiceC.j = indi.j * SUBN;
	
	/* Se calcula la submatriz resultante de C */
	for (k = 1; k <= subnumS; k++) {
		sumarSubMatriz(indiceA, indiceB, indiceC);
		indiceA.j += SUBN;
		indiceB.i += SUBN;
	}
}


/**
 * Asigna una tarea a un proceso con la 
 * submatriz corrspondiente
 */
void *mapearTarea(void *arg)
{
	/* Iteramos entre procesos hasta que se realicen todas las tareas */
	while (1) {

		pthread_mutex_lock(&lock);
		    tarea++;
		pthread_mutex_unlock(&lock);

		if ( tarea >= subnumM *subnumN)
			return NULL;

		indice indi = getIndice(tarea);
		realizarTarea(indi);
	}
}

int main(int argc, char *argv[])
{
	int i, j, k, token;
	pthread_t *hilos;
	double t1, t2;
	
	/* Inicializacion de parametros de la matriz/submatriz */
	if (argc!=2) {
		fprintf(stderr, "Modo de uso: ./mm2d <cantidad_de_hilos>\n");
		return -1;
	}

	cant_hilos=atoi(argv[1]);

	if (SUBN > M || SUBN > S || SUBN > N ||cant_hilos > M*N || SUBN <= 0) {
		fprintf(stderr, "Parametros no validos <M=%d> <S=%d> <N=%d> <SUBN=%d> <cantidad_de_hilos=%d>\n", M, S, N, SUBN, cant_hilos);
		return -1;
	}

	/* Se calcula el numero de submatrices por filas/columnas */
	subnumM = M / SUBN;
	subnumS = S / SUBN;
	subnumN = N / SUBN;
	if (subnumM * SUBN != M )
		subnumM++;
	if (subnumS * SUBN != S )
		subnumS++;
	if (subnumN * SUBN != N )
		subnumN++;

	
	/* Inicializacion de hilos */
	hilos = (pthread_t *)malloc(sizeof(pthread_t) * (cant_hilos-1));
	pthread_mutex_init(&lock, NULL);


	/* Se inicializan las matrices */
	printf("Inicializacion ...\n");
	inicializarMatrices();
	printf("Ok!!\n\n");

	/* Hora de multiplicar las matrices */
	printf("Multiplicacion ...\n");
	pthread_setconcurrency(cant_hilos);

	t1=getTime();
	/* Creamos y mapeamos todos los hilos */
	for (i=0; i < cant_hilos-1; i++)
		pthread_create(hilos+i, NULL, mapearTarea, NULL);
	mapearTarea(NULL);

	/* Esperamos a que finalicen todos los hilos */
	for (i=0; i < cant_hilos-1; i++)
		pthread_join(hilos[i], NULL);

	t2=getTime();

	printf("Ok!!\n\n");

	/* Se imprimien las matrices obtenidas si es que no son muy grandes */
	if ( M <= 10 && N <= 10) {
		printf("Imprimiendo Matriz A...\n");
		imprimirMatrizA();
		printf("ok!!\n\n");
		printf("Imprimiendo Matriz B...\n");
		imprimirMatrizB();
		printf("ok!!\n\n");
		printf("Imprimiendo Matriz C...\n");
		imprimirMatrizC();
		printf("ok!!\n\n");
	}

	printf("Duracion total de la multilplicacion de matrices %4f segundos\n", t2-t1);

	return 0;
}
