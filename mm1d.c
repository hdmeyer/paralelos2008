/**
 * File Name: sam_mm.c
 * Sample matrix multiplication code.
 * Algorithm: A and C are partitioned in 1D row blocks, B is partitioned in 
 * 1D column blocks. 
 *
 * The main multiplication involves bnum iterations, 
 * in iteration i, the ith sub-matrix column of B is broadcast to everybody 
 * then everybody uses this column and its own portion of A to update the ith 
 * column of its own portion of C.
 *
 * To compile:
 * cc sam_mm.c -o sam_mm -lmpi
 *
 * Comments: I've put lots of effort to make the code as easy to read as 
 * possible. So many opportunities for optimiaztions are possible, e.g.
 * all the dynmic checking of ownership can be eliminated, and even the
 * memcpy(...) can be removed.
 *
 * -Hong Tang, May. 1999
 */

# include <stdio.h>
# include <stdlib.h>
# include <pthread.h>
# include <sys/time.h>

static int n;	 	 // tamaño de la matriz
static int bsize;	 //tamaño de la submatriz
static int cant_hilos;	/* numbero total de hilos */
static int subnum;	 //numero de submatrices.
static int resto;
static pthread_mutex_t lock;
static int tarea = -1; 	// indica el id de la tarea actual 
int filasA;
int columnasB;
int comunAB;

double *A;
double *B;
double *C;

struct indices               /* Tipo de dato para representar */
{	int i, j;          		/* los ndices de una matriz */
};
typedef struct indices indice;
indice indiceGral;/* Aqui guardamos los indices de la matriz C para saber en que posición de ella estamos en un momento dado y nos
				sirve para movernos a traves de ella.*/

void cargarMatriz(double *, int, int);
void realizarTarea();
void *mapearTarea(void *);
void getIndice(int);
double getTime();
void imprimirMatriz(double *, int, int);

int main(int argc, char *argv[]){
	
	int i, j, k, token;
	pthread_t *hilos;
	double t1, t2;
	/*Para probar definimos arbitrariamente una matriz los tamaños de la matriz*/
	filasA=4;
	columnasB=5;
	comunAB = 3;
	
	/* Inicializacion de parametros de la matriz/submatriz */
	if (argc!=2) {
		fprintf(stderr, "Modo de uso: ./mm1d <cantidad_de_hilos>\n");
		return -1;
	}
	//leemos cantidad de procesos(hilos)
	cant_hilos=atoi(argv[1]);
	
	//dividimos la matriz de acuerdo a los procesos disponibles
	subnum=filasA;
	//resto=filasA%cant_hilos;
	
	/* Inicializacion de hilos */
	hilos = (pthread_t *)malloc(sizeof(pthread_t) * (cant_hilos-1));
	pthread_mutex_init(&lock, NULL);
	//Inicialización de las matrices.
	A = (double *)malloc (filasA*comunAB* sizeof(double));
	B = (double *)malloc (columnasB*comunAB* sizeof(double));
	C = (double *)malloc (filasA*columnasB* sizeof(double ));
	cargarMatriz(A, filasA ,comunAB);
	cargarMatriz(B,comunAB,columnasB);
	for (i=0; i < filasA; i++){
		for (j=0; j < columnasB; j++) {
			C[i*filasA+j] = 0;
		}
	}
	
	printf("Imprimimos matrices...\n");
	printf("Matriz A\n");
	imprimirMatriz(A,filasA,comunAB);
	printf("Matriz B\n");
	imprimirMatriz(B,comunAB,columnasB);
	
	
	printf("Ok!!\n\n");
	printf("Multiplicacion ...\n");
	pthread_setconcurrency(cant_hilos);
	
	t1=getTime();
	//creamos un hilo por cada proceso y le mandamos su parte para multiplicar
	for (i=0; i < cant_hilos-1; i++)
		pthread_create(hilos+i, NULL, mapearTarea, NULL);
	mapearTarea(NULL);

	for (i=0; i < cant_hilos-1; i++)
		pthread_join(hilos[i], NULL);

	t2=getTime();
	//AK HAY ERROR
	printf("Matriz C\n");
	imprimirMatriz(C,filasA,columnasB);
	printf("Duracion total de la multilplicacion de matrices %4f segundos\n", t2-t1);
	free(A);
	free(B);
	free(C);
	
}
void cargarMatriz(double X[], int fila, int columna){
	int i =0;
	int j = 0;
	printf(" Filas %d y columna %d \n", fila, columna);
	for (i=0; i < fila; i++){
		for (j=0; j < columna; j++) {
			X[i*fila+j]=(rand()%100)+1;
			printf(" Valor %2f \n",X[i*fila+j] );
			//printf("%2d",X[i*fila+j]);
		}
	}	
}
/*
void sumarSubMatriz( indice indiceA, indice indiceB, indice indiceC )
{
	int i, j, k;

	for (i=0; i < SUBN; i++)
		for (j=0; j < SUBN; j++)
			for (k=0; k < SUBN; k++)
				C[indiceC.i + i] [indiceC.j + j] += A[indiceA.i + i] [indiceA.j + k] * B[indiceB.i + k] [indiceB.j + j];
}
*/

/**
 * Realiza la suma de todas las submatrices
 */
void realizarTarea()
{
	int j, k;

	//printf("realizarTarea: %d...\n", tarea);
	//printf("i-j: %d.-.%d\n", indi.i, indi.j);
	// CORREGIMOS COLUMNASB POR COMUNAB
	for (k=0; k < comunAB; k++) {
		C[indiceGral.i*columnasB +indiceGral.j] += A[indiceGral.i*comunAB+k] * B[columnasB*k + indiceGral.j];
	}
	//sumarSubMatriz(indiceA, indiceB, indiceC);
	//indiceA.j += SUBN;
	//indiceB.i += SUBN;
	//Incrementamos aqui el valor para que la siguiente tarea, se mueva un lugar en la matriz C.
	printf("indice de la tarea J %d \n", indiceGral.j);
	indiceGral.j =indiceGral.j +1;

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

		if ( tarea >= filasA*columnasB)
			return NULL;

		getIndice(tarea);
		printf("indice de la tarea I %d \n", indiceGral.i);
		realizarTarea();
	}
}
void getIndice(int tarea)
{
	indiceGral.i = (int)tarea%subnum;
	//return indi;
}
double getTime() {
	struct timeval t;
	gettimeofday(&t, NULL);
	return (double)t.tv_sec+t.tv_usec*0.000001;
}

void imprimirMatriz(double m[], int fil, int col)
{
	int i, j = 0;
	for (i=0; i < fil; i++) {
		printf("\n\t| ");
		for (j=0; j < col; j++)
			printf("%2f ", m[i*fil + j]);
		printf("|");
	}
	printf("\n");
}

