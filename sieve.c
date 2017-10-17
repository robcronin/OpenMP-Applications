#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include <sys/time.h>

void parse_args(int argc, char *argv[], int *n, int *plots);
void print_list(char *list, int n);
double serial(char *primes, int n);
double mp(char *primes, int n, int procs);
void reset(char *primes, int n);




int main(int argc, char *argv[]){
	//set default arguments
	int n=1000000;
	int plots=0;
	parse_args(argc, argv, &n, &plots);
	double stime, mptime;

	// sets up a char array indicating prime or not
	char *primes=malloc((n+1)*sizeof(char));
	primes[0]='n';
	primes[1]='n';
	reset(primes, n);	// initialises 2-n as y

	stime=serial(primes, n);

	reset(primes,n);

	mptime=mp(primes, n, omp_get_num_procs());




	print_list(primes, n);
	printf("List Size:\t%d\n",n);
	printf("Serial time:\t%.2es\n",stime);
	printf("OpenMP time:\t%.2es\n",mptime);
	printf("Speed-up:\t%.2lfx\n",stime/mptime);
	printf("Procs Utilised:\t%d\n",omp_get_num_procs());
	free(primes);

	// creates some data files if requested
	if(plots==1){
		int i;
		int max=1000000;
		int max_procs=omp_get_num_procs();

		// sets up new primes array
		primes=malloc((max+1)*sizeof(char));
		primes[0]='n';
		primes[1]='n';
		reset(primes, max);
		
		// Plot by size
		FILE *plot;
		if((plot=fopen("sieve_size.dat", "w"))==NULL){
			printf("*** ERROR 1 ***\nCouldn't open sieve_size.dat\n");
			exit(1);	
		}

		// loops through increasing size by 1000
		for(i=1000;i<=max;i+=1000){
			stime=serial(primes, i);
			reset(primes, i);
			mptime=mp(primes, i, max_procs);
			reset(primes, i);
			fprintf(plot, "%d %lf\n",i,stime/mptime);
		}
		fclose(plot);

		// Plot by number of processes
		FILE *plot2;
		if((plot2=fopen("sieve_procs.dat", "w"))==NULL){
			printf("*** ERROR 1 ***\nCouldn't open sieve_procs.dat\n");
			exit(1);	
		}

		stime=serial(primes, max);
		// loops through increasing number of processes
		for(i=1;i<=max_procs;i++){
			reset(primes, max);
			mptime=mp(primes, max, i);
			reset(primes, max);
			fprintf(plot2, "%d %lf\n",i,stime/mptime);
		}
		fclose(plot2);

		free(primes);
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

double serial(char *primes, int n){
	int i,j;
	struct timeval t1, t2;
	// sift through current primes and mark multiples of it up to n as not prime
	//	can stop sifting once i^2>n
	gettimeofday(&t1, NULL);
	for(i=2;i*i<=n;i++){
		if(primes[i]=='y'){
			// can begin sifting from i^2
			for(j=i*i;j<=n;j+=i){
				primes[j]='n';
			}
		}
	}
	gettimeofday(&t2, NULL);
	double stime=(t2.tv_sec-t1.tv_sec)+(t2.tv_usec-t1.tv_usec)/1000000.0;
	return stime;
}

double mp(char *primes, int n, int procs){
	int i,j;
	struct timeval t1, t2;
	// sift through current primes and mark multiples of it up to n as not prime
	//	can stop sifting once i^2>n
	gettimeofday(&t1, NULL);
	for(i=2;i*i<=n;i++){
		if(primes[i]=='y'){
			// can begin sifting from i^2
			#pragma omp parallel for num_threads(procs)
			for(j=i*i;j<=n;j+=i){
				primes[j]='n';
			}
		}
	}
	gettimeofday(&t2, NULL);
	double mptime=(t2.tv_sec-t1.tv_sec)+(t2.tv_usec-t1.tv_usec)/1000000.0;
	return mptime;
}

void reset(char *primes, int n){
	int i;
	for(i=2;i<=n;i++){
		primes[i]='y';
	}
}

// prints list of numbers with a y
void print_list(char *list, int n){
	if(n<=100){
		int i;
		for(i=0;i<=n;i++){
			if(list[i]=='y'){
				printf("%d\n",i);
			}
		}
	}
}
