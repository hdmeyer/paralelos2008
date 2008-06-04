
/*********************************************************************************************
 * Archivo: mm1d.c
 *
 * Multiplicacion de matrices utilizando hilos y el particionamiento 1-d
 *
 *
 * Para compilar:
 * cc mm1d.c -o mm1d -lpthread
 *
 * Para ejecutar:
 * ./mm1d <cantidad_de_hilos>
 *
 * Adaptacion para su implementacion con hilos, particionamiento 1-D y matrices de tamanhos
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
# define filasA 2500		/* Cantidad de filas de la matriz A */
# define columnasB 2500		/* Cantidad de columnas de la matriz A = Cantidad de filas de la matriz B */
# define comunAB 2500		/* Cantidad de columnas de la matriz B */

static int cant_hilos;	 /* numero total de hilos */
static int resto;
static pthread_mutex_t lock;
static int tarea = -1; 	// indica el id de la tarea actual 


static int A[filasA][comunAB], B[comunAB][columnasB], C[filasA][columnasB];
struct indices               /* Tipo de dato para representar */
{	int i, j;           /* los ndices de una matriz */
};
typedef struct indices indice;

void cargarMatrices();
void realizarTarea(indice indiceGral);
void *mapearTarea(void *);
indice getIndice(int);
double getTime();
void imprimirMatrices();

int main(int argc, char *argv[]){
	
	int i, j, k, token;
	pthread_t *hilos;
	double t1, t2;
	
	/* Inicializacion de parametros de la matriz/submatriz */
	if (argc!=2) {
		fprintf(stderr, "Modo de uso: ./mm1d <cantidad_de_hilos>\n");
		return -1;
	}
	//leemos cantidad de procesos(hilos)
	cant_hilos=atoi(argv[1]);
	
	
	/* Inicializacion de hilos */
	hilos = (pthread_t *)malloc(sizeof(pthread_t) * (cant_hilos-1));
	pthread_mutex_init(&lock, NULL);

	//Inicializaci� de las matrices.
	cargarMatrices();
	
	
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

	if (filasA <= 10 && columnasB <= 10 && comunAB <= 10){
		imprimirMatrices();
	}

	printf("\n\nDuracion total de la multilplicacion de matrices %4f segundos\n", t2-t1);
	
	
}
void cargarMatrices(){
	int i, j = 0;
	for (i=0; i < filasA; i++)
		for (j=0; j < comunAB; j++)
			A[i][j] = (j+3) * (i+1);

	for (i=0; i < comunAB; i++)
		for (j=0; j < columnasB; j++)
			B[i][j] = i+j;

	for (i=0; i < filasA; i++)
		for (j=0; j < columnasB; j++)
			C[i][j] = 0;	
}


/**
 * Realiza la suma de todas las submatrices
 */
void realizarTarea(indice indiceGral)
{
	int k;

	for (k=0; k < comunAB; k++)
		C[indiceGral.i][indiceGral.j] += A[indiceGral.i][k] * B[k][indiceGral.j];
	
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

		indice indiceGral = getIndice(tarea);
		realizarTarea(indiceGral);
	}
}

indice getIndice(int tarea)
{	
	indice indiceGral;
	indiceGral.i = tarea / columnasB;
	indiceGral.j = tarea % columnasB; 
	return indiceGral;
}

double getTime() {
	struct timeval t;
	gettimeofday(&t, NULL);
	return (double)t.tv_sec+t.tv_usec*0.000001;
}

/**
 * Imprime el contenido de las matrices
 */
void imprimirMatrices()
{
	int i, j = 0;
	printf("Imprimimos matrices...\n");
	printf("Matriz A\n");

	for (i=0; i < filasA; i++) {

		printf("\n\t| ");
		for (j=0; j < comunAB; j++)
			printf("%d ", A[i][j]);
		printf("|");
	}

	printf("\n");
	printf("Matriz B\n");

	for (i=0; i < comunAB; i++) {

		printf("\n\t| ");
		for (j=0; j < columnasB; j++)
			printf("%d ", B[i][j]);
		printf("|");
	}

	printf("\n");
	printf("Matriz C\n");

	for (i=0; i < filasA; i++) {

		printf("\n\t| ");
		for (j=0; j < columnasB; j++)
			printf("%d ", C[i][j]);
		printf("|");
	}

	printf("\n");
	
}

