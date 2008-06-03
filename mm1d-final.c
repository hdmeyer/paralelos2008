
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
# define filasA 10		/* Cantidad de filas de la matriz A */
# define columnasB 10		/* Cantidad de columnas de la matriz A = Cantidad de filas de la matriz B */
# define comunAB 10		/* Cantidad de columnas de la matriz B */

static int n;	 	 // tamaño de la matriz
static int bsize;	 //tamaño de la submatriz
static int cant_hilos;	/* numbero total de hilos */
static int subnum;	 //numero de submatrices.
static int resto;
static pthread_mutex_t lock;
static int tarea = -1; 	// indica el id de la tarea actual 
static int cont =0;


static int A[filasA][comunAB], B[comunAB][columnasB], C[filasA][columnasB];
struct indices               /* Tipo de dato para representar */
{	int i, j;          		/* los ndices de una matriz */
};
typedef struct indices indice;
indice indiceGral;/* Aqui guardamos los indices de la matriz C para saber en que posición de ella estamos en un momento dado y nos
				sirve para movernos a traves de ella.*/

void cargarMatrices();
void realizarTarea();
void *mapearTarea(void *);
void getIndice(int);
double getTime();
void imprimirMatrices();

int main(int argc, char *argv[]){
	
	int i, j, k, token;
	pthread_t *hilos;
	double t1, t2;
	/*Para probar definimos arbitrariamente una matriz los tamaños de la matriz*/
	//filasA=4;
	//columnasB=5;
	//comunAB = 3;
	
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
	//A = (double *)malloc (filasA*comunAB* sizeof(double));
	//B = (double *)malloc (columnasB*comunAB* sizeof(double));
	//C = (double *)malloc (filasA*columnasB* sizeof(double ));
	cargarMatrices();
	for (i=0; i < filasA; i++){
		for (j=0; j < columnasB; j++) {
			C[i][j] = 0;
		}
	}
	
	
	imprimirMatrices();
	
	printf("Ok!!\n\n");
	printf("Multiplicacion ...\n");
	pthread_setconcurrency(cant_hilos);
	/*ESTO LO HACEMOS PARA QUE LA PRIMERA VEZ QUE ENTRE A LLAMAR A GETINDICE, SE TOME
	EL INDICE 0*/
	indiceGral.j = columnasB -1;
	
	t1=getTime();
	//creamos un hilo por cada proceso y le mandamos su parte para multiplicar
	for (i=0; i < cant_hilos-1; i++)
		pthread_create(hilos+i, NULL, mapearTarea, NULL);
	mapearTarea(NULL);

	for (i=0; i < cant_hilos-1; i++)
		pthread_join(hilos[i], NULL);

	t2=getTime();
	//AK HAY ERROR
	imprimirMatrices();
	printf("Duracion total de la multilplicacion de matrices %4f segundos\n", t2-t1);
	printf("CANTIDAD DE VECES QUE LLAMAMOS A GETINDICE %d \n", cont);
	printf("VALOR FINAL DE INDICEGRAL.J %d \n", indiceGral.j);
	
	
}
void cargarMatrices(){
	int i =0;
	int j = 0;
	printf(" Filas %d y columna %d \n", filasA, comunAB);

	for (i=0; i < filasA; i++){
		for (j=0; j < comunAB; j++) {
			A[i][j]=(rand()%100)+1;
			printf(" Valor %d \n",A[i][j]);
			//printf("%2d",X[i*fila+j]);
		}
	}
	printf(" Filas %d y columna %d \n", comunAB, columnasB);

	for (i=0; i < comunAB; i++){
		for (j=0; j < columnasB; j++) {
			B[i][j]=(rand()%100)+1;
			printf(" Valor %d \n",B[i][j]);
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
	//pthread_mutex_lock(&lock);
	//indiceGral.j += 1;
	//if(indiceGral.j == columnasB){
	//	indiceGral.j = 0;
	//}
	//pthread_mutex_unlock(&lock);

	for (k=0; k < comunAB; k++) {
		pthread_mutex_lock(&lock);
		C[indiceGral.i][indiceGral.j] += A[indiceGral.i][k] * B[k][indiceGral.j];
		printf("VALORES MULTIPLICADOS %d * %d\n", A[indiceGral.i][k], B[k][indiceGral.j]);
		printf("VALOR DE LA MATRIZ PARCIAL C[%d][%d] -> %d \n", indiceGral.i, indiceGral.j, C[indiceGral.i][indiceGral.j]);
		pthread_mutex_unlock(&lock);
	}
	printf("VALOR FINAL DEL CAMPO C[%d][%d] -> %d \n", indiceGral.i,indiceGral.j, C[indiceGral.i][indiceGral.j]);
	//sumarSubMatriz(indiceA, indiceB, indiceC);
	//indiceA.j += SUBN;
	//indiceB.i += SUBN;
	//Incrementamos aqui el valor para que la siguiente tarea, se mueva un lugar en la matriz C.
	//printf("indice de la tarea J %d \n", indiceGral.j);
	//pthread_mutex_lock(&lock);
	//indiceGral.j =indiceGral.j +1;
	//pthread_mutex_unlock(&lock);
	

	printf("INDICE DE LA TAREA I %d \n", indiceGral.i);
	printf("INDICE DE LA TAREA J %d \n", indiceGral.j);
	
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
		//printf("indice de la tarea I %d \n", indiceGral.i);
		realizarTarea();
	}
}
void getIndice(int tarea)
{	
	cont++;
	//pthread_mutex_lock(&lock);
	indiceGral.i = (int)tarea / subnum;
	indiceGral.j = (int)tarea % subnum;
	//indiceGral.j += 1;
	//if(indiceGral.j == columnasB){
	//	indiceGral.j = 0;
	//}
	//pthread_mutex_unlock(&lock);
	
	//return indi;
}
double getTime() {
	struct timeval t;
	gettimeofday(&t, NULL);
	return (double)t.tv_sec+t.tv_usec*0.000001;
}

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

