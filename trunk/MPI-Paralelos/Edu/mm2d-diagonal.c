
/*********************************************************************************************
 * Archivo: mm2d-diagonal.c
 *
 * Implementacion del algoritmo 2-D Diagonal para la multiplicacion de matrices
 * densas utilizando la Interfaz de Paso de Mensajes MPI.
 *
 * Presentado en el paper:Communication Efficient Matrix Multiplication on Hipercubes (Himanshu Gupta)
 * Disponible en: http://portal.acm.org/citation.cfm?id=231664.231684&coll=GUIDE&dl=GUIDE&CFID=68476495&CFTOKEN=16169251
 *
 * Para compilar:
 * gcc mm2d-diagonal.c -o mm2d-diagonal -lmpi
 *
 * Autores:
 * - Eduardo Rivas. erivas17@gmail.com
 * - Hugo Meyer. meyer.hugo@gmail.com
 *
 * Alumnos de la Universidad Nacional de Asuncion - Facultad Politecnica.
 * Carrera: Ingenieria Informatica.
 * Materia: Electiva V - Algoritmos Paralelos.
 * Profesor: Cristian Von Lucken.
 * Anho: 2008.
 **********************************************************************************************/


 /*********************************************************************************************
 * Breve explicacion del algoritmo implementado:
 *
 * La implementacion propuesta se basa en el articulo mencionado anteriormente y de manera
 * a simplificar la explicacion lo hemos divido en 8 pasos basicos que se pueden apreciar
 * a lo largo de la implementacion. Los 8 pasos son los siguientes:
 *
 * Algoritmo 2-D Diagonal (Figura 5 del paper):
 *
 * Paso 1: Distribucion Inicial. Cada procesador diagonal P(i,i) almacena el i-esimo grupo
 *         de filas y columnas de las matrices de entrada A y B respectivamente.
 *
 * si (i = j)
 *    Paso 2: One To All Broadcast. Broadcast A(*,j) a todos los procesadores P(*,j).
 *          A(*,j) es el j-esimo grupo de columnas de A, inicialmente almacenado en P(j,j).
 *
 *    desde k=0 hasta q-1
 *       Paso 3: Enviar B(i,k) a todos los procesadores P(k,j).
 *         B(i,*) es el i-esimo grupo de columnas de B, inicialmente almacenado
 *         en P(i,i) y B(i,k) es el k-esimo grupo de columnas de B(i,*).
 *    fin desde
 * fin si
 *
 * Paso 4: Recibir A(*,j) y B(j,i) de P(j,j).
 *
 * Paso 5: Calcular I(*,i) = A(*,i) x B(j,i).
 *
 * Paso 6: Enviar I(*,i) a P(i,i).
 *
 * si (i = j)
 *    Pasos 7 y 8: All To One Reduction. Recibir I(*,i) de P(i,k) y Calcular C(*,i) = C(*,i) + I(*,i)
 *    desde k=0 hasta q-1
 *       Paso 7: Recibimos los resultados parciales. Recibir I(*,i) de P(i,k).
 *       Paso 8: Calculamos las sumas parciales. Calcular C(*,i) = C(*,i) + I(*,i).
 *    fin desde
 * fin si
 *
 **********************************************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <sys/time.h>
# include <mpi.h>
# include <math.h>


static int n;		 /* Tamanho de la matriz */
static int numgrup;	 /* Numero de grupos en una fila/columna de la matriz (q) */
static int tamgrup;	 /* Tamanho del bloque de cada grupo */
static int nproc;	 /* Numero de nodos MPI */
static int myid; 	 /* Mi propio rank */

struct indices    /* Tipo de dato para representar */
{	int i, j;     /* los ndices de una matriz */
};
typedef struct indices indice;



/**
 * Esta funcion mapea el id de un proceso a su grupo correspondiente (i,j)
 */
indice getIndice(int id_proceso)
{
	indice indi;
	indi.i = id_proceso / numgrup;
	indi.j = id_proceso % numgrup;
	return indi;
}

/**
 * Esta funcion, al contrario de getIndice(), devuelve el id de un proceso
 * a partir de su indice correspondiente (i,j)
 */
int getIdProceso(indice indi)
{
	return indi.i * numgrup + indi.j;
}

/**
 * Paso 1: Distribucion Inicial.
 * Inicializa las matrices en forma generica dejando las entradas
 * almacenadas en los procesadores diagonales del mesh 2D.
 */
void inicializarMatrices(float *a, float *b)
{
    int i, j;

    for (i=0; i < n; i++){
        for (j=0; j < n; j++) {
            a[i*n+j] = (j+3) * (i+1);
            b[i*n+j] = i+j;
        }
    }
}

/**
 * Imprime el contenido de un grupo de matrices de A
 */
void imprimirBloqueA(float *a)
{
	int i, j = 0;
	printf("\n");
	for (i=0; i < n; i++) {
		printf("\n\t| ");
		for (j=0; j < tamgrup; j++)
			printf("%4f ", a[i*tamgrup+j]);
		printf("|");
	}
}

/**
 * Imprime el contenido de un grupo de matrices de B
 */
void imprimirBloqueB(float *b)
{
	int i, j = 0;
	printf("\n");
	for (i=0; i < tamgrup; i++) {
		printf("\n\t| ");
		for (j=0; j < n; j++)
			printf("%4f ", b[i*n+j]);
		printf("|");
	}
}

/**
 * Imprime el contenido de un grupo de matrices de C
 */
void imprimirBloqueC(float *c)
{
	int i, j = 0;
	printf("\n");
	for (i=0; i < n; i++) {
		printf("\n\t| ");
		for (j=0; j < tamgrup; j++)
			printf("%4f ", c[i*tamgrup+j]);
		printf("|");
	}
}


/**
 * Imprime el contenido del BUFFER_B
 */
void imprimirBlufferB(float *b)
{
	int i, j = 0;
	printf("\n");
	for (i=0; i < tamgrup; i++) {
		printf("\n\t| ");
		for (j=0; j < tamgrup; j++)
			printf("%4f ", b[i*tamgrup+j]);
		printf("|");
	}
}

/**
 * Imprime el contenido de una matriz
 * matriz: es la matriz que se desea imprimir
 * filas: cantidad de filas de la matriz
 * columnas: catidad de columnas de la matriz
 */
void imprimirMatriz(float *matriz, int filas, int columnas)
{
	int i, j = 0;
	printf("\n");
	for (i=0; i < filas; i++) {
		printf("\n\t| ");
		for (j=0; j < columnas; j++)
			printf("%4f ", matriz[i*columnas+j]);
		printf("|");
	}
}


/**
 * Rellena los buffers con los valores necesarios para calcular
 * el bloque correspondiente a la matriz resultado C
 * Parametros: A y B son las matrices de entrada
 * BLOQUE_A y BLOQUE_B: son los bloques donde se almacenan los bloques
 * de entradas pertenecientes a A y a B respectivamente.
 * indi: indica la primera posicion del bloque a calcular
 */
void cargarBloques(float *A, float *B, float *BLOQUE_A, float *BLOQUE_B, int id_destino){

     indice indi = getIndice(id_destino);
     indice indiceA, indiceB;
	 indiceA.i = 0;
	 indiceA.j = indi.j * tamgrup;
	 indiceB.i = indi.i * tamgrup;
	 indiceB.j = 0;

	 if ( n % tamgrup != 0){
        /* Sobran bloques de tamanhos distintos a tamgrup */
        if ( indiceA.j + tamgrup > n )
            indiceA.j = n - tamgrup;
        if ( indiceB.i + tamgrup > n )
            indiceB.i = n - tamgrup;
     }

	 int i, j;
	 /* Rellenamos el bloque A */
	 for (i = 0; i < n ; i++)
		for (j = 0; j < tamgrup; j++)
			BLOQUE_A[i*tamgrup+j] = A[(i+indiceA.i)*n+j+indiceA.j];

	 /* Rellenamos el bloque B */
	 for (i = 0; i < tamgrup ; i++)
		for (j = 0; j < n; j++)
			BLOQUE_B[i*n+j] = B[(i+indiceB.i)*n+j+indiceB.j];

}


/**
 * Almacena los elementos, necesarios para el broadcast,
 * que se encuentran en BLOQUE_B al BUFFER_B
 */
void cargarBufferB(float *BLOQUE_B, float *BUFFER, indice indi){

     indice indice;
     indice.i = indi.i * tamgrup;
	 indice.j = indi.j * tamgrup;

     int i, j;
     for (i = 0; i < tamgrup; i++)
         for (j = 0; j < tamgrup; j++)
             BUFFER[i*tamgrup+j] = BLOQUE_B[i*n+(j+indice.j)];

}


/**
 * Inicializa el contenido de un bloque a cero
 * bloque: es el bloque que se desea inicializar
 * filas: cantidad de filas del bloque
 * columnas: catidad de columnas del bloque
 */
void inicializarBloque(float *bloque, int filas, int columnas)
{
	int i, j = 0;
	for (i=0; i < filas; i++)
		for (j=0; j < columnas; j++)
			bloque[i*columnas+j] = 0;

}


/**
 * Realiza la multiplicacion local y parcial correspondiente
 * a un determinado bloque de la matriz resultado C.
 * Utiliza los BLOQUE_A y BUFFER_B como parametros de entrada y
 * el resultado lo almacena en BLOQUE_C.
 */
void multiplicarBloques(float *BLOQUE_A, float *BUFFER_B, float *BLOQUE_C){

    inicializarBloque(BLOQUE_C, n, tamgrup);
    int i, j, k;
    for (i=0; i < n; i++){
        for (j=0; j < tamgrup; j++) {
            for (k = 0; k < tamgrup; k++){
                BLOQUE_C[i*tamgrup+j] += BLOQUE_A[i*tamgrup+k] * BUFFER_B[k*tamgrup+j];
            }
        }
    }
}


/**
 * Realiza la suma parcial por filas de procesos de manera a
 * obtener la columna correspondiente de la matriz resultado C.
 * Suma el contenido de BUFFER_C con BLOQUE_C y luego deja
 * el resultado en este ultimo.
 */
void sumarBloques(float *BLOQUE_C, float *BUFFER_C){

    int i, j;
    for (i=0; i < n; i++)
        for (j=0; j < tamgrup; j++)
            BLOQUE_C[i*tamgrup+j] += BUFFER_C[i*tamgrup+j];

}


/**
 * Asigna los datos obtenidos en el BLOQUE_C a la matriz resultado C.
 * id_origen: es el id del procesador quien obtuvo el resultado parcial
 */
void asignarBloqueResultado(float *C, float *BLOQUE_C, int id_origen){

     indice indi = getIndice(id_origen);
     indice indiceC;
     indiceC.i = indi.i * tamgrup;
	 indiceC.j = indi.j * tamgrup;

	 if ( n % tamgrup != 0){
         /* Sobran bloques de tamanhos distintos a SUBN */
         if ( indiceC.i + tamgrup > n )
             indiceC.i = n - tamgrup;
         if ( indiceC.j + tamgrup > n )
             indiceC.j = n - tamgrup;
     }

     int i, j;
     for (i=0; i < n; i++)
		for (j=0; j < tamgrup; j++)
			C[i*n+j+indiceC.j] = BLOQUE_C[i*tamgrup+j];

}

/**
 * Asigna los datos obtenidos en el BLOQUE_B a la matriz de entrada B
 * id_origen: es el id del procesador quien obtuvo el resultado parcial
 */
void asignarBloqueBResultado(float *B, float *BLOQUE_B, int id_origen){

     indice indi = getIndice(id_origen);
     indice indiceB;
     indiceB.i = indi.i * tamgrup;
	 indiceB.j = 0;

     if ( n % tamgrup != 0){
         /* Sobran bloques de tamanhos distintos a SUBN */
         if ( indiceB.i + tamgrup > n )
             indiceB.i = n - tamgrup;
         if ( indiceB.j + tamgrup > n )
             indiceB.j = n - tamgrup;
     }

     int i, j;
     for (i=0; i < tamgrup; i++)
		for (j=0; j < n; j++)
			B[(i+indiceB.i)*n + j] = BLOQUE_B[i*n+j];

}


int main(int argc, char *argv[])
{
	double t1, t2; /* Variables utilizadas para  medir los tiempos */
	static float *BLOQUE_A, *BLOQUE_B, *BLOQUE_C; /* Bloque de Matrices. De entrada: A y B. De resultado: C */

    static float *BUFFER_B; /* Buffer utilizado para almacenar solo algunos elementos de BLOQUE_B. */
                            /* Necesario para minimizar el tamanho de los mensajes durante el broadcast */

    MPI_Status status;
	MPI_Init(&argc, &argv);

    if (argc!=2) {
		if (myid==0)
			printf("Se debe pasar como parametro el tamanho de la matriz\n");
		MPI_Finalize();
		return -1;
	}

    n=atoi(argv[1]);

	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    if (nproc > n*n)
        nproc = n*n; /* Maxima cantidad de nodos necesarios */

    numgrup = sqrt(nproc); /* numero de grupos en una fila/columna de la matriz (q) */

    if ( numgrup != sqrt(nproc) || nproc < 3 ){
        if (myid==0){
            printf("\nNumero de procesos no es un numero cuadratico.");
            printf("\nSe debe contar con mesh de procesadores (q*q).");
            printf("\nEj: 4, 9, 16, 25, 36, ... procesadores.\n");
        }
        MPI_Finalize();
        return -1;
    }

    tamgrup = n / numgrup;

	if (n % tamgrup != 0){
		tamgrup++;
	}



    /* Cada proceso contiene solo un bloque del total de elementos de las matrices respectivas */
    BLOQUE_A=(float *)malloc(sizeof(float)*(n*tamgrup));
	BLOQUE_B=(float *)malloc(sizeof(float)*(tamgrup*n));
	BLOQUE_C=(float *)malloc(sizeof(float)*(n*tamgrup));
    BUFFER_B=(float *)malloc(sizeof(float)*(tamgrup*tamgrup));


    /************************** Generacion de matrices de entrada *******************************/
    /************************************** Paso 1 **********************************************/
    /*********************** Distribucion Inicial. Se inicializan las matrices ******************/

    indice indi = getIndice(myid); // Se obtiene la etiqueta del procesador en el mesh (i,j)
	t1 = 0;
    if ( myid == 0){
        /* Se reciben los resultados de los procesadores diagonales */

        float *A,*B;
        A=(float *)malloc(sizeof(float)*(n*n));
        B=(float *)malloc(sizeof(float)*(n*n));

        /* Se lee las matrices de entradas A y B */
		printf("Inicializacion ...\n");
        inicializarMatrices(A, B);
		printf("Ok!!\n\n");
		t1 = MPI_Wtime();
		printf("Multiplicando ...\n");

        /* Se envian los bloques iniciales a los procesadores diagonales */
        indice indi2;
        int id_destino, k;
        for (k = 1; k < numgrup; k++){

            /* Se calcula el id del proceso destino */
            indi2.i = k;
            indi2.j = k;
            id_destino = getIdProceso(indi2);

            /* Se calcula los bloques que le pertenece */
            cargarBloques(A, B, BLOQUE_A, BLOQUE_B, id_destino);

            /* Se envian los bloques iniciales a los procesadores diagonales */
            MPI_Send(BLOQUE_A, n*tamgrup, MPI_FLOAT, id_destino, 10, MPI_COMM_WORLD);
            MPI_Send(BLOQUE_B, tamgrup*n, MPI_FLOAT, id_destino, 11, MPI_COMM_WORLD);

        }

        /* Se carga el ultimo bloque que pertenece al root */
        cargarBloques(A, B, BLOQUE_A, BLOQUE_B, 0);
        free(A);
        free(B);
    }
    else if ( indi.i == indi.j ){

        /* Se reciben los bloques iniciales */
        /* Solo los procesadores diagonales almacenan los datos de entrada */
        /* Cada procesador diagonal P(i,i) almacena el i-esimo grupo de filas y columnas */
        /* de las matrices de entrada A y B respectivamente */

        MPI_Recv(BLOQUE_A, n*tamgrup, MPI_FLOAT, 0, 10, MPI_COMM_WORLD, &status);
        MPI_Recv(BLOQUE_B, tamgrup*n, MPI_FLOAT, 0, 11, MPI_COMM_WORLD, &status);
    }


    if ( indi.i == indi.j ){ // Es un proceso diagonal

        /************************************** Paso 2 **********************************************/
        /********* One To All Broadcast. Broadcast A(*,j) a todos los procesadores P(*,j) ***********/
        int id_destino;
        for (id_destino = 0; id_destino < nproc; id_destino++){
            indice indi2 = getIndice(id_destino);
            if ( (indi.j == indi2.j) && (indi.i != indi2.i) ){ // Si es diagonal y no soy yo mismo
                MPI_Send(BLOQUE_A, n*tamgrup, MPI_FLOAT, id_destino, 1, MPI_COMM_WORLD);
            }
        }


        /************************************** Paso 3 **********************************************/
        /******************** Enviar B(i,k) a todos los procesadores P(k,j) *************************/
        /* B(i,*) es el i-esimo grupo de columnas de B, inicialmente almacenado */
        /* en P(i,i) y B(i,k) es el k-esimo grupo de columnas de B(i,*) */

        int k, mi_columna = 0;
        indice indi2;
        for (k = 0; k < numgrup; k++){

            /* Se calcula el id del proceso destino */
            indi2.i = k;
            indi2.j = indi.j;
            id_destino = getIdProceso(indi2);

            if ( indi.i != k ){ // Si no soy yo mismo

                /* Se calcula el bloque de columnas a enviar */
                indi2.i = indi.i;
                indi2.j = k;
                cargarBufferB(BLOQUE_B, BUFFER_B, indi2);

                /* Se envia el buffer */
                MPI_Send(BUFFER_B, tamgrup*tamgrup, MPI_FLOAT, id_destino, 2, MPI_COMM_WORLD);
            }
            else{
                /* Guardo mi columna para luego cargar mi buffer */
                mi_columna = k;
            }
        }

        /* Cargo mi buffer que me corresponde */
        indi2.i = indi.i;
        indi2.j = mi_columna;
        cargarBufferB(BLOQUE_B, BUFFER_B, indi2);

    }
    else{ /* No es un proceso diagonal; se reciben los bloques enviados */

        /************************************** Paso 4 **********************************************/
        /*********************** Recibir A(*,j) y B(j,i) de P(j,j) **********************************/
        indi.i = indi.j; // Se obtiene el id del proceso emisor (diagonal)
        int id_origen = getIdProceso(indi);

        /* Recibimos el bloque de matriz A */
        MPI_Recv(BLOQUE_A, n*tamgrup, MPI_FLOAT, id_origen,1, MPI_COMM_WORLD, &status);

        /* Recibimos el bloque de matriz B */
        MPI_Recv(BUFFER_B, tamgrup*tamgrup, MPI_FLOAT, id_origen,2, MPI_COMM_WORLD, &status);
    }


    /************************************** Paso 5 **********************************************/
    /*********************** Calcular I(*,i) = A(*,i) x B(j,i) **********************************/
    multiplicarBloques(BLOQUE_A, BUFFER_B, BLOQUE_C);

    indi = getIndice(myid);
    if ( indi.i != indi.j ){ /* No es un proceso diagonal */

        /************************************** Paso 6 **********************************************/
        /******************************* Enviar I(*,i) a P(i,i) *************************************/
        indi.j = indi.i; // Se obtiene el id del proceso receptor (diagonal)
        int id_destino = getIdProceso(indi);

        /* Se envia el resultado parcial al proceso diagonal */
        MPI_Send(BLOQUE_C, n*tamgrup, MPI_FLOAT, id_destino, 3, MPI_COMM_WORLD);
    }
    else{ /* Es un proceso diagonal */

        /********************************** Pasos 7 y 8 *********************************************/
        /**** All To One Reduction. Recibir I(*,i) de P(i,k) y Calcular C(*,i) = C(*,i) + I(*,i) ****/
        float *BUFFER_C; /* Creamos un buffer temporal para recibir los resultados parciales */
        BUFFER_C=(float *)malloc(sizeof(float)*(n*tamgrup));

        indice indi2;
        int id_origen,k;
        for (k = 0; k < numgrup; k++){

            /* Se calcula el id del proceso origen */
            indi2.i = indi.i;
            indi2.j = k;

            if ( indi.j != k ){ // Si no soy yo mismo

                /* Se recive los resultados parciales */
                id_origen = getIdProceso(indi2);

                /* Paso 7: Recibimos los resultados parciales. Recibir I(*,i) de P(i,k) */
                MPI_Recv(BUFFER_C, n*tamgrup, MPI_FLOAT, id_origen,3, MPI_COMM_WORLD, &status);

                /* Paso 8: Calculamos las sumas parciales. Calcular C(*,i) = C(*,i) + I(*,i) */
                sumarBloques(BLOQUE_C, BUFFER_C);
            }
        }

        free(BUFFER_C);
    }


    /********************************************************************************************/
    /************************ Fin de la Multiplicacion de Matrices ******************************/
    /********************************************************************************************/
    /* El paper explica hasta este punto, y finaliza diciendo que los bloques resultados se quedan */
    /* en los procesadores diagonales del mesh. A partir de aqui lo que se realiza es enviar dichos bloques */
    /* al proceso root de manera a realizar la recoleccion de datos (asumiendo que el procesador root */
    /* cuenta con memoria suficiente para almacenar toda la matriz). Esto nos permite imprimir o almacenar */
    /* el resultado de la matriz C resultante, asi como las matrices de entrada, en pantalla o en un archivo */
    /* para comprobar el correcto funcionamiento del algoritmo implementado */


    indi = getIndice(myid);
    if ( myid == 0){
        /* Se reciben los resultados de los procesadores diagonales */

        float *A,*B,*C; /* Matriz donde se almacena el resultado final */
        A=(float *)malloc(sizeof(float)*(n*n));
        B=(float *)malloc(sizeof(float)*(n*n));
        C=(float *)malloc(sizeof(float)*(n*n));

        /* Se almacena el resultado parcial en la matriz de resultado final, C */
        asignarBloqueResultado(C, BLOQUE_C, 0);
        asignarBloqueResultado(A, BLOQUE_A, 0);
        asignarBloqueBResultado(B, BLOQUE_B, 0);

        /* Se reciben los demas resultados parciales de los demas procesadores */
        indice indi2;
        int id_origen,k;
        for (k = 1; k < numgrup; k++){

            /* Se calcula el id del proceso origen */
            indi2.i = k;
            indi2.j = k;

            /* Se reciben los resultados parciales */
            id_origen = getIdProceso(indi2);
            MPI_Recv(BLOQUE_C, n*tamgrup, MPI_FLOAT, id_origen, 4, MPI_COMM_WORLD, &status);
            MPI_Recv(BLOQUE_A, n*tamgrup, MPI_FLOAT, id_origen, 5, MPI_COMM_WORLD, &status);
            MPI_Recv(BLOQUE_B, tamgrup*n, MPI_FLOAT, id_origen, 6, MPI_COMM_WORLD, &status);

            /* Se almacena el resultado parcial en la matriz de resultado final, C */
            asignarBloqueResultado(C, BLOQUE_C, id_origen);
            asignarBloqueResultado(A, BLOQUE_A, id_origen);
            asignarBloqueBResultado(B, BLOQUE_B, id_origen);
        }
				
		t2 = MPI_Wtime();

		if ( n <= 10){
            printf("\n\n------------------------------------------------------------------------");
            printf("\n*********************** Matrices de Entrada ******************************");
            printf("\n------------------------------------------------------------------------");
            printf("\n\nLa matriz de entrada A fue: ");
            imprimirMatriz(A,n,n);
            printf("\n\nLa matriz de entrada B fue: ");
            imprimirMatriz(B,n,n);
            printf("\n\n\n------------------------------------------------------------------------");
            printf("\n***************  La matriz resultado final es: **************************");
            printf("\n------------------------------------------------------------------------");
            imprimirMatriz(C,n,n);
		}
		printf("\nDuracion total de la multilplicacion de matrices %4f segundos\n", t2-t1);

        free(A);
        free(B);
        free(C);
    }
    else if ( indi.i == indi.j ){

        /* Se envia el resultado parcial al proceso root */
        MPI_Send(BLOQUE_C, n*tamgrup, MPI_FLOAT, 0, 4, MPI_COMM_WORLD);
        MPI_Send(BLOQUE_A, n*tamgrup, MPI_FLOAT, 0, 5, MPI_COMM_WORLD);
        MPI_Send(BLOQUE_B, tamgrup*n, MPI_FLOAT, 0, 6, MPI_COMM_WORLD);
    }

    free(BLOQUE_A);
    free(BLOQUE_B);
    free(BLOQUE_C);
    free(BUFFER_B);
    MPI_Finalize();
	return 0;
}

