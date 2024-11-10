#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

#include "fft.h"

#define SIZE 9

void print(double complex* data, size_t size) {
	for (size_t i = 0; i < size; i++) {
		printf("%lf + %lfI\n", creal(data[i]), cimag(data[i]));
	}
}

int main() {
	printf("\nTime:\n");
	double complex data[SIZE] = {
		1 + 0 * I,
		2 + 0 * I,
		3 + 0 * I,
		4 + 0 * I,
		5 + 0 * I,
		6 + 0 * I,
		7 + 0 * I,
		8 + 0 * I,
		9 + 0 * I
	};

	print(data, SIZE);

	if (fft_fft(data, SIZE, 1)) {
		return 1;
	}

	printf("\nFFT:\n");
	print(data, SIZE);

	if (fft_ifft(data, SIZE, 1)) {
		return 1;
	}

	printf("\nInverse Transform:\n");
	print(data, SIZE);

	return 0;
}
