// Celem tego programu jest prezentacja pomiaru i analizy
//efektywnosci programu za pomocÄ?CodeAnalyst(tm).
// Implementacja mnoÅŸenia macierzy jest realizowana za pomoca typowego
// algorytmu podrÄ?cznikowego.
#include <stdio.h>
#include <time.h>
#include <windows.h>
#include "omp.h"

#define USE_MULTIPLE_THREADS true
#define MAXTHREADS 128
int NumThreads;
double start;

static const int MAT_SIZE = 3024;     // liczba wierszy macierzy
static const int R = 50;     // liczba wierszy macierzy

float matrix_a[MAT_SIZE + R][MAT_SIZE + R];    // lewy operand
float matrix_b[MAT_SIZE + R][MAT_SIZE + R];    // prawy operand
float matrix_r[MAT_SIZE + R][MAT_SIZE + R];    // wynik

FILE* result_file;

void initialize_matrices()
{
	// zdefiniowanie zawarosci poczatkowej macierzy
	//#pragma omp parallel for
	for (int i = 0; i < MAT_SIZE; i++) {
		for (int j = 0; j < MAT_SIZE; j++) {
			matrix_a[i][j] = (float)rand() / RAND_MAX;
			matrix_b[i][j] = (float)rand() / RAND_MAX;
			matrix_r[i][j] = 0.0;
		}
	}
}

void initialize_matricesZ()
{
	// zdefiniowanie zawarosci poczatkowej macierzy
#pragma omp parallel for
	for (int i = 0; i < MAT_SIZE; i++) {
		for (int j = 0; j < MAT_SIZE; j++) {
			matrix_r[i][j] = 0.0;
		}
	}
}
void print_result()
{
	// wydruk wyniku
	for (int i = 0; i < MAT_SIZE; i++) {
		for (int j = 0; j < MAT_SIZE; j++) {
			fprintf(result_file, "%6.4f ", matrix_r[i][j]);
		}
		fprintf(result_file, "\n");
	}
}

void multiply_matrices_KIJ()
{
	// mnozenie macierzy KIJ
	int i, j, k;
#pragma omp parallel for private(i,j)
	for (k = 0; k < MAT_SIZE; k++)
		for (i = 0; i < MAT_SIZE; i++)
			for (j = 0; j < MAT_SIZE; j++)
				matrix_r[i][j] += matrix_a[i][k] * matrix_b[k][j];
}

void multiply_matrices_KIJ_SEQ()
{
	// mnozenie macierzy KIJ
	int i, j, k;
	for (k = 0; k < MAT_SIZE; k++)
		for (i = 0; i < MAT_SIZE; i++)
			for (j = 0; j < MAT_SIZE; j++)
				matrix_r[i][j] += matrix_a[i][k] * matrix_b[k][j];
}

void multiply_matrices_IKJ_sequential()
{
	// mnozenie macierzy sekwencyjnie
	for (int i = 0; i < MAT_SIZE; i++)
		for (int k = 0; k < MAT_SIZE; k++)
			for (int j = 0; j < MAT_SIZE; j++)
				matrix_r[i][j] += matrix_a[i][k] * matrix_b[k][j];
}

void multiply_matrices_6() {
	for (int i = 0; i < MAT_SIZE; i += R) {
		for (int j = 0; j < MAT_SIZE; j += R) {
			for (int k = 0; k < MAT_SIZE; k += R) {
#pragma omp parallel for
				for (int ii = i; ii < i + R; ++ii) {
					for (int kk = k; kk < k + R; ++kk) {
						for (int jj = j; jj < j + R; ++jj) {
							matrix_r[ii][jj] += matrix_a[ii][kk] * matrix_b[kk][jj];
						}
					}
				}
			}
		}
	}
}



void print_elapsed_time()
{
	double elapsed;
	double resolution;

	// wyznaczenie i zapisanie czasu przetwarzania
	elapsed = (double)clock() / CLK_TCK;
	resolution = 1.0 / CLK_TCK;
	printf("Czas: %8.4f sec \n",
		elapsed - start);

	fprintf(result_file,
		"Czas wykonania programu: %8.4f sec (%6.4f sec rozdzielczosc pomiaru)\n",
		elapsed - start, resolution);
}

int main(int argc, char* argv[])
{
	//   start = (double) clock() / CLK_TCK ;
	if ((result_file = fopen("classic.txt", "a")) == NULL) {
		fprintf(stderr, "nie mozna otworzyc pliku wyniku \n");
		perror("classic");
		return(EXIT_FAILURE);
	}


	//Determine the number of threads to use
	if (USE_MULTIPLE_THREADS) {
		SYSTEM_INFO SysInfo;
		GetSystemInfo(&SysInfo);
		NumThreads = SysInfo.dwNumberOfProcessors;
		if (NumThreads > MAXTHREADS)
			NumThreads = MAXTHREADS;
	}
	else
		NumThreads = 1;
	fprintf(result_file, "Klasyczny algorytm mnozenia macierzy, liczba watkow %d \n", NumThreads);
	printf("liczba watkow  = %d\n\n", NumThreads);

	initialize_matrices();
	start = (double)clock() / CLK_TCK;
	multiply_matrices_KIJ();
	printf("KIJ ");
	print_elapsed_time();

	initialize_matricesZ();
	start = (double)clock() / CLK_TCK;
	multiply_matrices_KIJ_SEQ();
	printf("KIJ seq ");
	print_elapsed_time();


	//start = (double)clock() / CLK_TCK;
	//multiply_matrices_IKJ_sequential();
	//printf("IKJ_sequential ");
	//print_elapsed_time();

	//initialize_matricesZ();


	//start = (double)clock() / CLK_TCK;
	//multiply_matrices_6();
	//printf("ijk ii_kk_jj");
	//print_elapsed_time();


	fclose(result_file);

	return(0);
}