#include<stdio.h>
#include<stdlib.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <omp.h>
#include <immintrin.h>
#include <avx2intrin.h>
int main(int argc, char * argv)
{
	
	double gflop;
	double *A;
	double *B;
	double *BT;
	double*C_golden;
	float*C_goldenf;
	double*C1;
	double*C2;
	double *C3;
	double *C4;
	double*C5;
	double*C6;
	float*C7;
	float *AF;
	float *BF;
	FILE *fp,*ff;
		
	if((fp=fopen("matrixfile.txt","wb"))==NULL)
	{
		//printf("\nopen file Fail,close!");
		getchar();
		exit(1);
	}
	if((ff=fopen("plot.txt","wb"))==NULL)
	{
		//printf("\nopen file Fail,close!");
		getchar();
		exit(1);
	}
	int g=1600;
	for(int m=g;m<=2000;m+=100)
	{
		
	fprintf(ff,"a%d=[",m);
	fprintf(fp,"\\multirow{5}{*}{%d}",m);
	for(int nthreads=1;nthreads<=16;nthreads*=2)
	{
		
		int k=m;
		int n=m;
		omp_set_num_threads(nthreads);
		/*int nthreads=1;//atoi(argv[1]);
		
		int m=200;//atoi(argv[2]);
		int k=200;///atoi(argv[3]);
		int n=200;//atoi(argv[4]);*/
		//printf("%d******\n",nthreads);
		
		//printf("Matrix A is %i x %i, matrix B is %i x %i\n", m, k, k, n);
		gflop =2.0* (double)m *(double)n *(double)k/ 1000000000.0;
		A=(double*) malloc(sizeof(double)*m*k);
		B=(double*) malloc( sizeof(double)*k*n);
		BT=(double*) malloc(sizeof(double)*k*n);
		C_golden=( double*)malloc(sizeof(double)*m*n);
		C_goldenf=( float*)malloc(sizeof(float)*m*n);
		C1=(double*) malloc(sizeof(double)*m*n);
		C2=(double*) malloc(sizeof(double)*m*n);
		C3 =(double*) malloc(sizeof(double)*m *n);
		C4 =(double*) malloc(sizeof(double)*m*n);
		C5=( double*) malloc(sizeof(double)*m*n);
		C6=(double*) malloc(sizeof(double)*m*n);
		C7=(float*) malloc(sizeof(double)*m*n);
		AF=(float*) malloc(sizeof(float)*m*k);
		BF=(float*) malloc( sizeof(float)*k*n);
		
		srand(0);
		for (int i= 0; i<m; i++)
			for (int j= 0; j <k; j++)
				{
					AF[i+k+j]=rand()%m;
					A[i+k+j]=rand()%m;
				}
		for (int i= 0; i<k; i++)
			for(int j=0;j<n;j++)
				{
					BF[i+n+j]=rand()%k;
					B[i+n+j]=rand()%k;
				}
		//warmup
		memset(C_golden,0,sizeof(double)*m*n);
		memset(C_goldenf,0,sizeof(float)*m*n);
		#pragma omp parallel for
		for(int mi=0;mi<m;mi++)
		{
			for(int ni=0;ni<n;ni++)
			{
				for(int ki=0;ki<k;ki++)
				{
					C_golden[mi*n+ni]+=A[mi*k+ki]*B[ki*n+ni];
					C_goldenf[mi*n+ni]+=AF[mi*k+ki]*BF[ki*n+ni];
					}	
			}
		}
		struct timeval t1,t2;
		//fprintf(ff,"%d\r\n",nthreads);
		fprintf(fp, "&threads=%d",nthreads);
		// row-col method, A in row-major, B in row-major :m-n-k
		memset(C1, 0, sizeof(double)* m*n); 
		gettimeofday(&t1,NULL);
		 #pragma omp parallel for
		 for (int mi =0; mi<m; mi++)
		 	for (int ni =0; ni<n; ni++)
		 		for (int ki =0; ki<k; ki++)
		 			C1[mi*n+ni] +=A[mi*k+ ki] *B[ki * n+ ni];
		gettimeofday(&t2, NULL);
		double time_use=(t2.tv_sec-t1.tv_sec)*1000+(t2.tv_usec-t1.tv_usec)/1000;
		//printf("\nGEMM 1 (row-row) used %4.5f s, %4.2f GFlop/s\n", time_use, gflop/time_use);
		fprintf(ff, "[%4.5f,",time_use);
		fprintf(fp, "&%4.5f",time_use);
		/*int count =0;
		for(int i=0;i<m*n;i++)
			if(C_golden[i]!=C1[i])
		 		count++;
		 if (count ==0)
		 	printf("GEMM 1 (row-row) PASS! \n\n"); 
		 else
		 	printf("GEMM 1 row-row) NOT PASS! \n\n");*/
		
		//method 2 row-col method, A in row-major.B in col-major
		for(int i=0;i<k;i++)
			for (int j=0;j<n;j++)
				BT[j*k+1]=B[i*n+1];
		memset(C2,0,sizeof(double)*m*n);
		gettimeofday(&t1,NULL);
		#pragma omp parallel for
		 for (int mi =0; mi<m; mi++)
		 	for (int ni =0; ni<n; ni++)
		 		for (int ki =0; ki<k; ki++)
		 			C2[mi *n+ni] +=A[mi *k+ ki] *BT[ni * k+ ki];
		gettimeofday(&t2, NULL);
		time_use=(t2.tv_sec-t1.tv_sec)*1000+(t2.tv_usec-t1.tv_usec)/1000;
		//printf("GEMM 2 (row-col) used %4.5f s, %4.2f GFlop/s\n", time_use, gflop/time_use);
		fprintf(ff, "%4.5f,",time_use);
		fprintf(fp, "&%4.5f",time_use);
	/*	count =0;
		for(int i=0;i<m*n;i++)
			if(C_golden[i]!=C1[i])
		 		count++;
		 if (count ==0)
		 	printf("GEMM 2 (row-col) PASS! \n\n"); 
		 else
		 	printf("GEMM 2 row-col) NOT PASS! \n\n");*/
		
		//method 3. col-col method .A in row-majoe,B in row-major:n-k-m
		memset(C3,0,sizeof(double)*m*n);
		gettimeofday(&t1, NULL);
		 #pragma omp parallel for
	 	for (int ni =0; ni<n; ni++)
	 		for (int ki =0; ki<k; ki++)
	 			for (int mi =0; mi<m; mi++)
	 				C3[mi *n+ni] +=A[mi *k+ ki] *B[ki * n+ ni];
		gettimeofday(&t2, NULL);
		time_use=(t2.tv_sec-t1.tv_sec)*1000+(t2.tv_usec-t1.tv_usec)/1000;
		//printf("\nGEMM 3 (col-col) used %4.5f s, %4.2f GFlop/s\n", time_use, gflop/time_use);
		fprintf(fp, "&%4.5f",time_use);
		fprintf(ff, "%4.5f,",time_use);
	/*	count =0;
		for(int i=0;i<m*n;i++)
			if(C_golden[i]!=C3[i])
		 		count++;
		 if (count ==0)
		 	printf("GEMM 3 (col-col) PASS! \n\n"); 
		 else
		 	printf("GEMM 3 (col-col) NOT PASS! \n\n");*/
		
		//method 4 row-row method A in row -major,B in row-major:m-k-n
		memset(C4,0,sizeof(double)*m*n);
		gettimeofday(&t1, NULL);
		 #pragma omp parallel for
	 	for (int mi =0; mi<m; mi++)
	 		for (int ki =0; ki<k; ki++)
	 			for (int ni =0; ni<n; ni++)
	 				C4[mi *n+ni] +=A[mi *k+ ki] *B[ki * n+ ni];
		gettimeofday(&t2, NULL);
		time_use=(t2.tv_sec-t1.tv_sec)*1000+(t2.tv_usec-t1.tv_usec)/1000;
		fprintf(fp, "&%4.5f",time_use);
		fprintf(ff, "%4.5f,",time_use);
		//printf("\nGEMM 4 (row-row) used %4.5f s, %4.2f GFlop/s\n", time_use, gflop/time_use);
		/*count =0;
		for(int i=0;i<m*n;i++)
			if(C_golden[i]!=C4[i])
		 		count++;
		 if (count ==0)
		 	printf("GEMM 4 (row-row) PASS! \n\n"); 
		 else
		 	printf("GEMM 4 (row-row) NOT PASS! \n\n");*/
		
		//method 5 row-row avx2 method A in row -major,B in row-major:m-k-n
		memset(C5,0,sizeof(double)*m*n);
		gettimeofday(&t1, NULL);
		int nloop=n/4;
		int tail=n-4*nloop;
		 #pragma omp parallel for
	 	for (int mi =0; mi<m; mi++)
	 		for (int ki =0; ki<k; ki++)
	 		{
	 			double aval=A[mi*k+ki];
	 			__m256d aavx2=_mm256_set1_pd(aval);
	 			//avx2 part
	 			for (int ni =0; ni<nloop; ni++)
	 			{
	 				__m256d bavx2=_mm256_loadu_pd(&B[ki*n+ni*4]);
	 				bavx2=_mm256_mul_pd(aavx2,bavx2);
	 				__m256d cavx2=_mm256_loadu_pd(&C5[mi*n+ni*4]);
	 				cavx2=_mm256_add_pd(cavx2,bavx2);
	 				_mm256_storeu_pd(&C5[mi*n+ni*4],cavx2);
				 }
				 //tail part
				 for (int ni=4*nloop;ni<n;ni++)
	 				C5[mi *n+ni] +=aval*B[ki*n+ni];
			 }	
		gettimeofday(&t2, NULL);
		time_use=(t2.tv_sec-t1.tv_sec)*1000+(t2.tv_usec-t1.tv_usec)/1000;
		//printf("\nGEMM 5 (row-row avx2) used %4.5f s, %4.2f GFlop/s\n", time_use, gflop/time_use);
		fprintf(fp, "&%4.5f",time_use);
		fprintf(ff, "%4.5f,",time_use);
		/*count =0;
		for(int i=0;i<m*n;i++)
			if(C_golden[i]!=C5[i])
		 		count++;
		 if (count ==0)
		 	printf("GEMM 5 (row-row avx2) PASS! \n\n"); 
		 else
		 	printf("GEMM 5 (row-row avx2) NOT PASS! \n\n");*/
		int BIT=4;
		//method 6 row -row avx2 block method.A in row-major,B in row -major
		if(m%BIT==0&&k%BIT==0&&n%BIT==0)
		{
			memset(C6,0,sizeof(double)*m*n);
			gettimeofday(&t1, NULL);
			int mb=m/BIT;
			int kb=k/BIT;
			int nb=n/BIT;
			 #pragma omp parallel for
		 	for (int mbi =0; mbi<mb; mbi++)
		 		for (int nbi =0; nbi<nb; nbi++)
		 		{
		 			int cbase=(mbi*4+0)*n+nbi*4;
		 			__m256d cr1=_mm256_loadu_pd(&C6[cbase]);
		 			__m256d cr2=_mm256_loadu_pd(&C6[cbase+n]);
		 			__m256d cr3=_mm256_loadu_pd(&C6[cbase+2*n]);
		 			__m256d cr4=_mm256_loadu_pd(&C6[cbase+3*n]);
		 			
		 			for (int kbi =0; kbi<kb; kbi++)
		 			{
		 				int bbase=(kbi*4+0)*n+nbi*4;
			 			__m256d br1=_mm256_loadu_pd(&B[bbase]);
			 			__m256d br2=_mm256_loadu_pd(&B[bbase+n]);
			 			__m256d br3=_mm256_loadu_pd(&B[bbase+2*n]);
			 			__m256d br4=_mm256_loadu_pd(&B[bbase+3*n]);
			 			__m256d a;
			 			//row 1
			 			int abase=(mbi*4+0)*k+kbi*4+0;
			 			a=_mm256_set1_pd(A[abase]);
			 			a=_mm256_mul_pd(a,br1);
			 			cr1=_mm256_add_pd(cr1,a);
			 			
			 			a=_mm256_set1_pd(A[abase+1]);
			 			a=_mm256_mul_pd(a,br2);
			 			cr1=_mm256_add_pd(cr1,a);
			 			
			 			a=_mm256_set1_pd(A[abase+2]);
			 			a=_mm256_mul_pd(a,br3);
			 			cr1=_mm256_add_pd(cr1,a);
			 			
			 			a=_mm256_set1_pd(A[abase+3]);
			 			a=_mm256_mul_pd(a,br4);
			 			cr1=_mm256_add_pd(cr1,a);
			 			
			 			//row 2
			 			abase+=k;
			 			a=_mm256_set1_pd(A[abase]);
			 			a=_mm256_mul_pd(a,br1);
			 			cr2=_mm256_add_pd(cr2,a);
			 			
			 			a=_mm256_set1_pd(A[abase+1]);
			 			a=_mm256_mul_pd(a,br2);
			 			cr2=_mm256_add_pd(cr2,a);
			 			
			 			a=_mm256_set1_pd(A[abase+2]);
			 			a=_mm256_mul_pd(a,br3);
			 			cr2=_mm256_add_pd(cr2,a);
			 			
			 			a=_mm256_set1_pd(A[abase+3]);
			 			a=_mm256_mul_pd(a,br4);
			 			cr2=_mm256_add_pd(cr2,a);
			 			
			 			//row 3
			 			abase+=k;
			 			a=_mm256_set1_pd(A[abase]);
			 			a=_mm256_mul_pd(a,br1);
			 			cr3=_mm256_add_pd(cr3,a);
			 			
			 			a=_mm256_set1_pd(A[abase+1]);
			 			a=_mm256_mul_pd(a,br2);
			 			cr3=_mm256_add_pd(cr3,a);
			 			
			 			a=_mm256_set1_pd(A[abase+2]);
			 			a=_mm256_mul_pd(a,br3);
			 			cr3=_mm256_add_pd(cr3,a);
			 			
			 			a=_mm256_set1_pd(A[abase+3]);
			 			a=_mm256_mul_pd(a,br4);
			 			cr3=_mm256_add_pd(cr3,a);
			 			
			 			//row 4
			 			abase+=k;
			 			a=_mm256_set1_pd(A[abase]);
			 			a=_mm256_mul_pd(a,br1);
			 			cr4=_mm256_add_pd(cr4,a);
			 			
			 			a=_mm256_set1_pd(A[abase+1]);
			 			a=_mm256_mul_pd(a,br2);
			 			cr4=_mm256_add_pd(cr4,a);
			 			
			 			a=_mm256_set1_pd(A[abase+2]);
			 			a=_mm256_mul_pd(a,br3);
			 			cr4=_mm256_add_pd(cr4,a);
			 			
			 			a=_mm256_set1_pd(A[abase+3]);
			 			a=_mm256_mul_pd(a,br4);
			 			cr4=_mm256_add_pd(cr4,a);
					 }
					 
					 _mm256_storeu_pd(&C6[cbase],cr1);
					 _mm256_storeu_pd(&C6[cbase+n],cr2);
					 _mm256_storeu_pd(&C6[cbase+2*n],cr3);
					 _mm256_storeu_pd(&C6[cbase+3*n],cr4);
				 }	
		
			gettimeofday(&t2, NULL);
			time_use=(t2.tv_sec-t1.tv_sec)*1000+(t2.tv_usec-t1.tv_usec)/1000;
			//printf("\nGEMM 6 (row-row avx2-block) used %4.5f s, %4.2f GFlop/s\n", time_use, gflop/time_use);
			fprintf(fp, "&%4.5f",time_use);
			fprintf(ff, "%4.5f,",time_use);
			/*count =0;
			for(int i=0;i<m*n;i++)
				if(C_golden[i]!=C6[i])
			 		count++;
			 if (count ==0)
			 	printf("GEMM 6 (row-row avx2=block) PASS! \n\n"); 
			 else
			 	printf("GEMM 6 (row-row avx2-block) NOT PASS! \n\n");*/
		}
		//method 7 row -row avx2 block=8 method.A in row-major,B in row -major
		
		
		BIT=8;
		if(m%BIT==0&&k%BIT==0&&n%BIT==0)
		{
			memset(C7,0,sizeof(float)*m*n);
			gettimeofday(&t1, NULL);
			int mb=m/BIT;
			int kb=k/BIT;
			int nb=n/BIT;
			 #pragma omp parallel for
		 	for (int mbi =0; mbi<mb; mbi++)
		 		for (int nbi =0; nbi<nb; nbi++)
		 		{
		 			int cbase=(mbi*BIT+0)*n+nbi*BIT;
		 			__m256 cr[BIT];
		 			for(int ii=0;ii<BIT;ii++)
		 				cr[ii]=_mm256_load_ps(&C7[cbase+ii*n]);
		 			
		 			for (int kbi =0; kbi<kb; kbi++)
		 			{
		 				int bbase=(kbi*BIT+0)*n+nbi*BIT;
		 				__m256 br[BIT];
		 				for(int jj=0;jj<BIT;jj++)
		 					br[jj]=_mm256_loadu_ps(&BF[bbase+jj*n]);
			 			__m256 a;
			 			
			 			int abase=(mbi*BIT+0)*k+kbi*BIT+0;
			 			for(int row=0;row<BIT;row++)
			 			{
			 				for(int col=0;col<BIT;col++)
			 				{
			 					a=_mm256_set1_ps(AF[abase+col]);
					 			a=_mm256_mul_ps(a,br[col]);
					 			cr[row]=_mm256_add_ps(cr[row],a);
							 }
							 abase+=k;
						 }
			 			
					 }
					 for(int f=0;f<BIT;f++)
					 	_mm256_storeu_ps(&C7[cbase+f*n],cr[f]);
				 }	
		
			gettimeofday(&t2, NULL);
			double time_use=(t2.tv_sec-t1.tv_sec)*1000+(t2.tv_usec-t1.tv_usec)/1000;
			//printf("\nGEMM 7 (row-row avx2-block-8) used %4.5f s, %4.2f GFlop/s\n", time_use, gflop/time_use);
			fprintf(fp, "&%4.5f \\\\\r\n",time_use);
			if(nthreads<16)
			{
				fprintf(fp,"\\cline{2-9}");
				fprintf(ff, "%4.5f],",time_use);
			}
			
		else
		fprintf(ff, "%4.5f]]",time_use);
			/*int count =0;
			for(int i=0;i<m*n;i++)
				if(C_goldenf[i]!=C7[i])
			 		count++;
			 if (count ==0)
			 	printf("GEMM 7 (row-row avx2-block-8) PASS! \n\n"); 
			 else
			 	printf("GEMM 7 (row-row avx2-block-8) NOT PASS! \n\n");*/
		//fprintf(fp, "\r\n");
		}
		else
		{
			fprintf(fp,"&\\\\\r\n");
			
			if(nthreads<16)
			{
				fprintf(fp,"\\cline{2-9}");
				fprintf(ff, "%4.5f],",time_use);
			}
			
		else
		fprintf(ff, "%4.5f]]",time_use);
		}
		
	}
	fprintf(fp,"\\hline");
	fprintf(ff,"\r\n");
	
}
	fclose(fp);	
 	fclose(ff);
	
	
	free(A);
	free(B);
	free(AF);
	free(BF);
	free(BT);
	free(C1);
	free(C2);
	free(C3);
	free(C4);
	free(C5);
	free(C6);
	free(C7);
}
 
