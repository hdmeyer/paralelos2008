
/*********************************************************************************************
 * Archivo: mm_serial.c
 *
 * Multiplicacion de matrices utilizando el algoritmo serial
 *
 *
 * Para compilar:
 * cc mm_serial.c -o mm_serial -lpthread
 *
 * Para ejecutar:
 * ./mm_serial 
 *
 * Autores:
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


static int A[M][S], B[S][N], C[M][N];


/**
 * Obtiene el tiempo actual del sistema
 */
double getTime() {
	struct timeval t;
	gettimeofday(&t, NULL);
	return (double)t.tv_sec+t.tv_usec*0.000001;
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
 * Realiza la multiplicacion de las matrices
 */
void multiplicarMatrices()
{
	int i, j, k;

	for (i=0; i < M; i++)
		for (j=0; j < N; j++)
			for (k=0; k < S; k++)
				C[i] [j] += A[i] [k] * B[k] [j];
}


int main(int argc, char *argv[])
{

	double t1, t2;

	/* Se inicializan las matrices */
	printf("Inicializacion ...\n");
	inicializarMatrices();
	printf("Ok!!\n\n");

	/* Hora de multiplicar las matrices */
	printf("Multiplicacion ...\n");
	t1=getTime();
	multiplicarMatrices();
	t2=getTime();

	printf("Ok!!\n\n");

	/* Se imprimien las matrices obtenidas si es que no son muy grandes */
	if ( M <= 10 && N <= 10 && S <= 10 ) {
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
