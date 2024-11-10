#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062
#endif

void swap(double complex* a, double complex* b) {
	double complex temp = *a;
	*a = *b;
	*b = temp;
}

void shuffle(double complex* data, size_t size, size_t stride) {
	size_t target = 0;
	for (size_t position = 0; position < size; position++) {
		if (target > position) {
			swap(&data[position * stride], &data[target * stride]);
		}

		size_t mask = size >> 1;

		while (target & mask) {
			target &= ~mask;
			mask >>= 1;
		}
		target |= mask;
	}
}

int radixfft(double complex* data, size_t size, size_t stride) {
	size_t halfsize = size / 2;
	shuffle(data, size, stride);

	double complex* twiddles = calloc(halfsize, sizeof(double complex));
	if (twiddles == NULL) {
		fprintf(stderr, "[ERROR]: Radix FFT could not allocate auxillary space\n");
		return 1;
	}

	for (size_t i = 0; i < halfsize; i++) {
		twiddles[i] = cexp(-I * M_PI * i / ((double)halfsize));
	}

	for (size_t step = 1; step < size; step <<= 1) {
		const size_t jump = step << 1;
		double complex twiddle = 1.0 + 0.0 * I;
		for (size_t group = 0; group < step; group++) {
			for (size_t pair = group; pair < size; pair += jump) {
				const size_t match = pair + step;
				const double complex product = twiddle * data[match * stride];
				data[match * stride] = data[pair * stride] - product;
				data[pair * stride] += product;
			}

			if (group + 1 == step) {
				continue;
			}

			twiddle = twiddles[(group + 1) * halfsize / step];
		}
	}

	free(twiddles);

	return 0;
}

int radixifft(double complex* data, size_t size, size_t stride) {
	size_t halfsize = size / 2;
	shuffle(data, size, stride);

	double complex* twiddles = calloc(halfsize, sizeof(double complex));
	if (twiddles == NULL) {
		fprintf(stderr, "[ERROR]: Radix FFT could not allocate auxillary space\n");
		return 1;
	}

	for (size_t i = 0; i < halfsize; i++) {
		twiddles[i] = cexp(I * M_PI * i / ((double)halfsize));
	}

	for (size_t step = 1; step < size; step <<= 1) {
		const size_t jump = step << 1;
		double complex twiddle = 1.0 + 0.0 * I;
		for (size_t group = 0; group < step; group++) {
			for (size_t pair = group; pair < size; pair += jump) {
				const size_t match = pair + step;
				const double complex product = twiddle * data[match * stride];
				data[match * stride] = data[pair * stride] - product;
				data[pair * stride] += product;
			}

			if (group + 1 == step) {
				continue;
			}

			twiddle = twiddles[(group + 1) * halfsize / step];
		}
	}

	free(twiddles);

	return 0;
}

int fft(double complex* data, size_t size, size_t stride) {
	size_t larger_size = 1;
	while (larger_size >> 1 < size) {
		larger_size <<= 1;
	}

	double complex* chirp = calloc(size, sizeof(double complex));
	if (chirp == NULL) {
		fprintf(stderr, "[ERROR]: could not allocate auxillary space for chirp\n");
		return 1;
	}

	for (size_t i = 0; i < size; i++) {
		chirp[i] = cexp(-I * M_PI * i * i / ((double)size));
	}

	double complex* a = calloc(larger_size, sizeof(double complex));
	if (a == NULL) {
		free(chirp);
		fprintf(stderr, "[ERROR]: could not allocate auxillary space for a\n");
		return 1;
	}

	for (size_t i = 0; i < size; i++) {
		a[i] = chirp[i] * data[i * stride];
	}

	double complex* b = calloc(larger_size, sizeof(double complex));

	if (b == NULL) {
		free(chirp);
		free(a);
		fprintf(stderr, "[ERROR]: could not allocate auxillary space for b\n");
		return 1;
	}

	b[0] = chirp[0];
	for (size_t i = 1; i < size; i++) {
		b[i]               = conj(chirp[i]);
		b[larger_size - i] = conj(chirp[i]);
	}

	if (radixfft(a, larger_size, 1)) {
		free(chirp);
		free(a);
		free(b);
		fprintf(stderr, "[ERROR]: Failed during internal RADIX FFT call\n");
		return 1;
	}
	if (radixfft(b, larger_size, 1)) {
		free(chirp);
		free(a);
		free(b);
		fprintf(stderr, "[ERROR]: Failed during internal RADIX FFT call\n");
		return 1;
	}

	for (size_t i = 0; i < larger_size; i++) {
		a[i] *= b[i];
	}

	if (radixifft(a, larger_size, 1)) {
		free(chirp);
		free(a);
		free(b);
		fprintf(stderr, "[ERROR]: Failed during internal RADIX IFFT call\n");
		return 1;
	}

	for (size_t i = 0; i < size; i++) {
		data[i * stride] = a[i] * chirp[i] / larger_size;
	}

	free(chirp);
	free(a);
	free(b);

	return 0;
}

int ifft(double complex* data, size_t size, size_t stride) {
	size_t larger_size = 1;
	while (larger_size >> 1 < size) {
		larger_size <<= 1;
	}

	double complex* chirp = calloc(size, sizeof(double complex));
	if (chirp == NULL) {
		fprintf(stderr, "[ERROR]: could not allocate auxillary space from chirp\n");
		return 1;
	}

	for (size_t i = 0; i < size; i++) {
		chirp[i] = cexp(I * M_PI * i * i / ((double)size));
	}

	double complex* a = calloc(larger_size, sizeof(double complex));
	if (a == NULL) {
		free(chirp);
		fprintf(stderr, "[ERROR]: could not allocate auxillary space from a\n");
		return 1;
	}

	for (size_t i = 0; i < size; i++) {
		a[i] = chirp[i] * data[i * stride];
	}

	double complex* b = calloc(larger_size, sizeof(double complex));
	if (b == NULL) {
		free(chirp);
		free(a);
		fprintf(stderr, "[ERROR]: could not allocate auxillary space from b\n");
		return 1;
	}

	b[0] = chirp[0];
	for (size_t i = 1; i < size; i++) {
		b[i]               = conj(chirp[i]);
		b[larger_size - i] = conj(chirp[i]);
	}

	if (radixfft(a, larger_size, 1)) {
		free(chirp);
		free(a);
		free(b);
		fprintf(stderr, "[ERROR]: Failed during internal RADIX FFT call\n");
		return 1;
	}
	if (radixfft(b, larger_size, 1)) {
		free(chirp);
		free(a);
		free(b);
		fprintf(stderr, "[ERROR]: Failed during internal RADIX FFT call\n");
		return 1;
	}

	for (size_t i = 0; i < larger_size; i++) {
		a[i] *= b[i];
	}

	if (radixifft(a, larger_size, 1)) {
		free(chirp);
		free(a);
		free(b);
		fprintf(stderr, "[ERROR]: Failed during internal RADIX IFFT call\n");
		return 1;
	}

	for (size_t i = 0; i < size; i++) {
		data[i * stride] = a[i] * chirp[i] / larger_size;
	}

	free(chirp);
	free(a);
	free(b);

	return 0;
}

void print(double complex* data, size_t size) {
	for (size_t i = 0; i < size; i++) {
		printf("%lf + %lfI\n", creal(data[i]), cimag(data[i]));
	}
}

#define SIZE 9

int main() {
	printf("\nTime:\n");
	double complex data[SIZE] = {
		1 + 9 * I,
		2 + 8 * I,
		3 + 7 * I,
		4 + 6 * I,
		5 + 5 * I,
		6 + 6 * I,
		7 + 3 * I,
		8 + 2 * I,
		9 + 1 * I
	};

	print(data, SIZE);

	if (fft(data, SIZE, 1)) {
		return 1;
	}

	printf("\nFFT:\n");
	print(data, SIZE);

	if (ifft(data, SIZE, 1)) {
		return 1;
	}

	printf("\nInverse Transform:\n");
	print(data, SIZE);

	return 0;
}
