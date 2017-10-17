#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>

void parse_args(int argc, char *argv[], int *n, int *plots);
void reset(double *matrix, int n, int seed);
void print_matrix(double *matrix, int n);
void print_sol(double *matrix, int n);
void swap_rows(int r1, int r2, double *matrix, int n);
double serial(double *matrix, int n);
double mp(double *matrix, int n, int procs);


int main(int argc, char *argv[]){
	//set default arguments
	int n=500;
	int plots=0;
	parse_args(argc, argv, &n, &plots);

	double stime, mptime;
	int seed=time(NULL);



	// sets up augmented matrix
	double *matrix=malloc(n*(n+1)*sizeof(double));
	reset(matrix, n, seed);

	stime=serial(matrix, n);
	print_sol(matrix,n);

	reset(matrix, n, seed);

	mptime=mp(matrix, n, omp_get_num_procs());
	print_sol(matrix,n);


	printf("Matrix Size:\t%d\n",n);
        printf("Serial time:\t%.2es\n",stime);
	printf("OpenMP time:\t%.2es\n",mptime);
	printf("Speed-up:\t%.2lfx\n",stime/mptime);
	printf("Procs Utilised:\t%d\n",omp_get_num_procs());
	free(matrix);

	// creates some data files if requested
	if(plots==1){
		int i;
		int max=1000;
		int max_procs=omp_get_num_procs();

		// sets up new primes array
		matrix=malloc(max*(max+1)*sizeof(double));
		reset(matrix, n, seed);
	
		// Plot by size
		FILE *plot;
		if((plot=fopen("gauss_size.dat", "w"))==NULL){
			printf("*** ERROR 1 ***\nCouldn't open gauss_size.dat\n");
			exit(1);	
		}

		// loops through increasing size by 1000
		for(i=20;i<=max;i+=20){
			reset(matrix, i, seed);
			stime=serial(matrix, i);
			reset(matrix, i, seed);
			mptime=mp(matrix, i, max_procs);
			fprintf(plot, "%d %lf\n",i,stime/mptime);
		}
		fclose(plot);

		// Plot by number of processes
		FILE *plot2;
		if((plot2=fopen("gauss_procs.dat", "w"))==NULL){
			printf("*** ERROR 1 ***\nCouldn't open gauss_procs.dat\n");
			exit(1);	
		}

		reset(matrix, max, seed);
		stime=serial(matrix, max);
		// loops through increasing number of processes
		for(i=1;i<=max_procs;i++){
			reset(matrix, max, seed);
			mptime=mp(matrix, max, i);
			fprintf(plot2, "%d %lf\n",i,stime/mptime);
		}
		fclose(plot2);

		free(matrix);
	}


	return 0;
}


void parse_args(int argc, char *argv[], int *n, int *plots){
	//parse command line arguments
	int opt;
	while((opt=getopt(argc,argv,"n:p"))!=-1){
		switch(opt){
			case 'n':
				*n=atoi(optarg);
				break;
			case 'p':
				*plots=1;
				break;
			default:
				fprintf(stderr,"Usage: %s -n\n",argv[0]);
				exit(EXIT_FAILURE);
		}
	}
	return;
}

void reset(double *matrix, int n, int seed){
	srand48(seed);
	int i;
	for(i=0;i<n*(n+1);i++){
		matrix[i]=4*(drand48()-0.5);
	}
}

double serial(double *matrix, int n){
	int i,j,k;
	double mult;
	struct timeval t1, t2;
	double stime;

	gettimeofday(&t1, NULL);
	for(i=0;i<n-1;i++){

		// checks if the row has a zero in the (i,i) spot
		if(matrix[i*(n+1)+i]==0){
			for(j=i+i;j<n;j++){
				// finds a suitable row to swap with
				if(matrix[j*(n+1)+i]!=0){
					swap_rows(i, j, matrix, n);
					break;
				}

			}
		}

		// if all the rows have a zero in the ith column then either infinite or no solutions
		if(matrix[i*(n+1)+i]==0){
			printf("Unique solution does not exist\n");
			exit(1);
		}
		
		// runs gaussian elimination
		for(j=i+1;j<n;j++){
			mult=matrix[j*(n+1)+i];
			if(mult!=0){
				mult/=matrix[i*(n+1)+i];
				for(k=i;k<n+1;k++){
					matrix[j*(n+1)+k]-=mult*matrix[i*(n+1)+k];
				}
			}
		}
	}

	// if (n-1, n-1) entry is zero then either infinite or no solution
	if(matrix[(n-1)*(n+1)+n-1]==0){
		printf("Unique solution does not exist\n");
		exit(1);
	}

	// BACK SUBSTITUTION
	for(i=n-1;i>=0;i--){
		// solves backwards for each value of xi
		for(j=n-1;j>i;j--){
			matrix[i*(n+1)+n]-=matrix[i*(n+1)+j]*matrix[j*(n+1)+n];
		}
		matrix[i*(n+1)+n]/=matrix[i*(n+1)+i];
	}

	gettimeofday(&t2, NULL);
	stime=(t2.tv_sec-t1.tv_sec)+(t2.tv_usec-t1.tv_usec)/1000000.0;
	return stime;
}



double mp(double *matrix, int n, int procs){
	int i,j,k;
	double mult;
	struct timeval t1, t2;
	double mptime;

	gettimeofday(&t1, NULL);
	for(i=0;i<n-1;i++){

		// checks if the row has a zero in the (i,i) spot
		if(matrix[i*(n+1)+i]==0){
			for(j=i+i;j<n;j++){
				// finds a suitable row to swap with
				if(matrix[j*(n+1)+i]!=0){
					swap_rows(i, j, matrix, n);
					break;
				}

			}
		}

		// if all the rows have a zero in the ith column then either infinite or no solutions
		if(matrix[i*(n+1)+i]==0){
			printf("Unique solution does not exist\n");
			exit(1);
		}
		
		// runs gaussian elimination
		#pragma omp parallel for private(k,mult) num_threads(procs) 
		for(j=i+1;j<n;j++){
			mult=matrix[j*(n+1)+i];
			if(mult!=0){
				mult/=matrix[i*(n+1)+i];
				for(k=i;k<n+1;k++){
					matrix[j*(n+1)+k]-=mult*matrix[i*(n+1)+k];
				}
			}
		}
	}

	// if (n-1, n-1) entry is zero then either infinite or no solution
	if(matrix[(n-1)*(n+1)+n-1]==0){
		printf("Unique solution does not exist\n");
		exit(1);
	}
		

	// BACK SUBSTITUTION
	for(i=n-1;i>=0;i--){
		// solves backwards for each value of xi
		double b=matrix[i*(n+1)+n];
		#pragma omp parallel for reduction(+:b) num_threads(omp_get_num_procs())
		for(j=n-1;j>i;j--){
			b-=matrix[i*(n+1)+j]*matrix[j*(n+1)+n];
		}
		matrix[i*(n+1)+n]=b/matrix[i*(n+1)+i];
	}
	gettimeofday(&t2, NULL);
	mptime=(t2.tv_sec-t1.tv_sec)+(t2.tv_usec-t1.tv_usec)/1000000.0;
	return mptime;

}


// prints augmented matrix
void print_matrix(double *matrix, int n){
	if(n<=10){
		int i,j;
		for(i=0;i<n;i++){
			for(j=0;j<n+1;j++){
				printf("%lf\t",matrix[i*(n+1)+j]);
			}
			printf("\n");
		}
		printf("\n");
	}
}

// prints last column ie solution
void print_sol(double *matrix, int n){
	if(n<=10){
		int i;
		printf("The solution is:\n");
		for(i=0;i<n;i++){
			printf("x%d\t%lf\n",i+1,matrix[i*(n+1)+n]);
		}
		printf("\n");
	}
}

// swaps row r1 with r2 in the matrix
void swap_rows(int r1, int r2, double *matrix, int n){
	double temp;
	int i;
	for(i=0;i<n+1;i++){
		temp=matrix[r1*(n+1)+i];
		matrix[r1*(n+1)+i]=matrix[r2*(n+1)+i];
		matrix[r2*(n+1)+i]=temp;
	}
}

