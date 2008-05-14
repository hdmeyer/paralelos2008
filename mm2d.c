/**
 * Archivo: mm2d.c
 *
 * Multiplicacion de matrices utilizando el particionamiento 2-D 
 *
 * Para compilar:
 * cc mm2d.c -o mm2d -lpthread
 *
 * Para ejecutar:
 * ./mm2d <cantidad_de_hilos>
 *
 */

# include <stdio.h>
# include <stdlib.h>
# include <pthread.h>
# include <sys/time.h>

# define N 9		/* tamanho de la matriz */
# define SUBN 3		/* tamanho de la submatriz  --  Obs: SUBN debe ser multiplo de N */


static int subnum;	 /* numero de sub-matrices en una fila/columna de la matriz */ 
static int cant_hilos;	/* numbero total de hilos */

static int A[N][N], B[N][N], C[N][N];

static pthread_mutex_t lock;
static int tarea = -1; 	/* indica el id de la tarea actual */


struct indices               /* Tipo de dato para representar */
{	int i, j;          		/* los ndices de una matriz */
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
	indi.i = tarea / subnum;
	indi.j = tarea % subnum;
	return indi;
}

/**
 * Imprime el contenido de una matriz
 */
void imprimirMatriz(int m[N][N])
{
	int i, j = 0;
	for (i=0; i < N; i++) {
		printf("\n\t| ");
		for (j=0; j < N; j++)
			printf("%2d ", m[i][j]);
		printf("|");
	}
}

/**
 * Realiza la suma de una sola submatriz
 */
void sumarSubMatriz( indice indiceA, indice indiceB, indice indiceC )
{
	int i, j, k;

	for (i=0; i < SUBN; i++)
		for (j=0; j < SUBN; j++)
			for (k=0; k < SUBN; k++)
				C[indiceC.i + i] [indiceC.j + j] += A[indiceA.i + i] [indiceA.j + k] * B[indiceB.i + k] [indiceB.j + j];
}


/**
 * Realiza la suma de todas las submatrices
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


	//printf("realizarTarea: %d...\n", tarea);
	//printf("i-j: %d.-.%d\n", indi.i, indi.j);

	for (k=0; k < subnum; k++) {
		sumarSubMatriz(indiceA, indiceB, indiceC);
		indiceA.j += SUBN;
		indiceB.i += SUBN;
	}
}


/**
 * Asigna una tarea a una submatriz para que luego 
 * se pueda realizar la suma en dicha submatriz
 */
void *mapearTarea(void *arg)
{
	while (1) {

		pthread_mutex_lock(&lock);
		    tarea++;
		pthread_mutex_unlock(&lock);

		if ( tarea >= subnum *subnum)
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
	subnum = N / SUBN;
	if (subnum * SUBN != N || cant_hilos > N*N) {
		fprintf(stderr, "Parametros no validos <N=%d> <SUBN=%d> <cantidad_de_hilos=%d>\n", N, SUBN, cant_hilos);
		return -1;
	}

	
	/* Inicializacion de hilos */
	hilos = (pthread_t *)malloc(sizeof(pthread_t) * (cant_hilos-1));
	pthread_mutex_init(&lock, NULL);


	printf("Inicializacion ...\n");

	/*Se inicializan las matrices */
	for (i=0; i < N; i++){
		for (j=0; j < N; j++) {
			A[i][j] = (j+3) * (i+1);
			B[i][j] = i+j;
			C[i][j]= 0;
		}
	}

	printf("Ok!!\n\n");
	printf("Multiplicacion ...\n");
	
	pthread_setconcurrency(cant_hilos);

	t1=getTime();
	for (i=0; i < cant_hilos-1; i++)
		pthread_create(hilos+i, NULL, mapearTarea, NULL);
	mapearTarea(NULL);

	for (i=0; i < cant_hilos-1; i++)
		pthread_join(hilos[i], NULL);

	t2=getTime();

	printf("Ok!!\n\n");

	if ( N <= 10) {
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

	printf("Duracion total de la multilplicacion de matrices %4f segundos\n", t2-t1);

	return 0;
}
