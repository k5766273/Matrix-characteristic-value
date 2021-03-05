#include <iostream>
#include <map>
#include<math.h>
#include <iomanip>
using namespace std;
#define  N 8  //方阵的宽度
#define pi 3.1415926
//幂法求特征值
//求向量的无穷范数
double getMax(double* A){
	double max = 0;
	for (int i = 0; i < N; ++i){
		if (max < fabs(*(A + i))){
			max = fabs(*(A + i));
		}
	}
	return max;
}
void multi(double array[N], double* a0){
	int i, j;
	double* result = new double[N];
	memset(result, 0, N * sizeof(double));
	for (i = 0; i < N; ++i){
		for (j = 0; j < N; ++j){
			result[i] += array[i*N+j] * a0[j];
		}
		cout << "a0[" << i << "]=" << result[i] << endl;
	}
	for (i = 0; i < N; ++i){
		*(a0 + i) = result[i];
	}
	delete[] result;
}
bool power_png(double &pld, double* env, double* a, int n) {
	int count = 0;		//迭代次数
	double maxElemt1, maxElemt2;
	do{
		++count;
		maxElemt1 = getMax(env);
		for (int j = 0; j < n; ++j){
			*(env + j) = *(env + j) / maxElemt1;	
		}
		multi(a, env);
		maxElemt2 = getMax(env);
		cout << "最大特征值 = " << maxElemt2 << endl << endl;
		if (count > 500)
			break;//迭代次数检验
	} while (fabs((maxElemt2 - maxElemt1) / maxElemt1) > 0.0001);//精度
	pld = maxElemt2;
	cout << "迭代次数为 " << count << endl;
	cout << "最大特征值为 " << pld << endl;
	return false;
}
//雅克比求特征值
bool Jacobi(double* ev,double* a, int n  ){
	int nCount = 0;		
	while (1) {
		double prec = a[1];
		int nRow = 0;
		int nCol = 1;
		for (int i = 0; i < n; i++) {		
			for (int j = 0; j < n; j++) {		
				double d = fabs(a[i * n + j]);
				if ((i != j) && (d > prec)) {
					prec = d;
					nRow = i;
					nCol = j;
				}
			}
		}
		if (prec < 1e-10)     //精度
			break;
		if (nCount > 10000)     //迭代次数
			break;
		nCount++;
		double Aii = a[nRow * n + nRow];
		double Aij = a[nRow * n + nCol];
		double Ajj = a[nCol * n + nCol];
		double Angle = 0.5 * atan2(-2 *Aij,Ajj -Aii);
		double SinTheta = sin(Angle);
		double CosTheta = cos(Angle);
		double Sin2Theta = sin(2 * Angle);
		double Cos2Theta = cos(2 * Angle);
		a[nRow * n + nRow] =Aii * CosTheta * CosTheta + Ajj * SinTheta * SinTheta + 2 *Aij * CosTheta * SinTheta;
		a[nCol * n + nCol] =Aii * SinTheta * SinTheta + Ajj * CosTheta * CosTheta - 2 *Aij * CosTheta * SinTheta;
		a[nRow * n + nCol] = 0.5 * (Ajj -Aii) * Sin2Theta +Aij * Cos2Theta;
		a[nCol * n + nRow] = a[nRow * n + nCol];
		for (int i = 0; i < n; i++) {
			if ((i != nCol) && (i != nRow)) {
				int u = i * n + nRow;	 
				int w = i * n + nCol;	
				prec = a[u];
				a[u] = a[w] * SinTheta + prec * CosTheta;
				a[w] = a[w] * CosTheta - prec * SinTheta;
			}
		}
		for (int j = 0; j < n; j++) {
			if ((j != nCol) && (j != nRow)) {
				int u = nRow * n + j;	
				int w = nCol * n + j;	
				prec = a[u];
				a[u] = a[w] * SinTheta + prec * CosTheta;
				a[w] = a[w] * CosTheta - prec * SinTheta;
			}
		}
		for (int i = 0; i < n; i++) {
			int u = i * n + nRow;		 
			int w = i * n + nCol;		
		}
	}
	map<double, int> mapEigen;
	for (int i = 0; i < n; i++) {
		ev[i] = a[i * n + i];
		mapEigen.insert(make_pair(ev[i], i));
	}
	double* pdbTmpVec = new double[n * n];
	map<double, int>::reverse_iterator iter = mapEigen.rbegin();
	for (int j = 0; iter != mapEigen.rend(), j < n; ++iter, ++j) {
		ev[j] = iter->first;
	}

	for (int i = 0; i < n; i++) {
		double dSumVec = 0;
		for (int j = 0; j < n; j++)
			dSumVec += pdbTmpVec[j * n + i];
		if (dSumVec < 0) {
			for (int j = 0; j < n; j++)
				pdbTmpVec[j * n + i] *= -1;
		}
	}
	delete[]pdbTmpVec;
	return false;
}
//Hensenberg矩阵
struct arrays {
	double a[N * N] = { 0 };
	double c[N * N];//[N * N] = { 0 };
	double b[N * N];//[N * N] = { 0 };
} arr[N - 1];
void getInverse(double* a) {
	for (int i = 0; i < N - 1; i++)
		for (int j = i + 1; j < N; j++)
			a[j * N + i] = -a[j * N + i];
}
void getTemp(double* a, double* b, int n, int p) {//p为迭代次数 n为矩阵宽度
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {
			if (i == j)
				a[i * N + j] = 1;
			else
				a[i * N + j] = 0;
		}
	for (int i = p + 1, t = 1; i < n; i++, t++)
		if (b[p * N + p - 1] != 0)
			a[i * N + p] = -b[(p + t) * n + p - 1] / b[p * N + p - 1];
}
void multiply(double* c, double* a, double* b) {
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j) {
			c[i * N + j] = 0;
			for (int s = 0; s < N; ++s) {
				c[i * N + j] += a[i * N + s] * b[s * N + j];
			}
			if (c[i * N + j] < 2.84217e-5 && c[i * N + j]>-2.84217e-5)//精度
				c[i * N + j] = 0;
		}
}
double* gauss_hessen2(double* b, int n) {
	for (int i = 0; i < N * N; i++)
		arr[0].b[i] = b[i];
	for (int i = 0; i < n - 2; i++) {
		getTemp(arr[i].a, arr[i].b, n, i + 1);
		multiply(arr[i].c, arr[i].a, arr[i].b);
		getInverse(arr[i].a);
		multiply(arr[i + 1].b, arr[i].c, arr[i].a);
	}
	return  arr[N - 2].b;
}
//QR分解求特征值
bool qr_eng(double* A, double* Q, double* R) {
	int i, j, k, m;
	double temp;
	double a[N], b[N];
	for (j = 0; j < N; j++) {
		for (i = 0; i < N; i++) {
			a[i] = A[i * N + j];
			b[i] = A[i * N + j];
		}
		for (k = 0; k < j; k++) {
			R[k * N + j] = 0;
			for (m = 0; m < N; m++)
				R[k * N + j] += a[m] * Q[m * N + k];
			for (m = 0; m < N; m++)
				b[m] -= R[k * N + j] * Q[m * N + k];
		}
		temp = 0;
		for (i = 0; i < N; i++)
			temp += b[i] * b[i];
		R[j * N + j] = sqrt(temp);
		for (i = 0; i < N; i++)
			Q[i * N + j] = b[i] / sqrt(temp);
	}
	return false;
}
void Multiplicate(double A[N * N], double R[N * N], double Q[N * N]) {
	int i, j, k;
	double temp;
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++) {
			temp = 0;
			for (k = 0; k < N; k++)
				temp += R[i * N + k] * Q[k * N + j];
			A[i * N + j] = temp;
		}
}
int main(){
	double A[64] = {
		611,196,-192,407,-8,-52,-49,29,
		196,899,113,-192,-71,-43,-8,-44,
		-192,113,899,196,61,49,8,52,
		407,-192,196,611,8,44,59,-23,
		-8,-71,61,8,411,-599,208,208,
		-52,-43,49,44,-599,411,208,208,
		-49,-8,8,59,208,208,99,-911,
		29,-44,52,-23,208,208,-911,99
	};
	double B[100] = {
	1.0, 1.0 / 2.0,1.0 / 3.0,1.0 / 4.0,1.0 / 5.0, 1.0 / 6.0,  1.0 / 7.0,  1.0 / 8.0,  1.0 / 9.0,  1.0 / 10.0,
	1.0 / 2.0,  1.0 / 3.0, 1.0 / 4.0,   1.0 / 5.0,1.0 / 6.0, 1.0 / 7.0,  1.0 / 8.0,  1.0 / 9.0,  1.0 / 10.0, 1.0 / 11.0,
	1.0 / 3.0,  1.0 / 4.0,  1.0 / 5.0,   1.0 / 6.0, 1.0 / 7.0, 1.0 / 8.0,  1.0 / 9.0,  1.0 / 10.0, 1.0 / 11.0, 1.0 / 12.0,
	 1.0 / 4.0,  1.0 / 5.0,1.0 / 6.0,1.0 / 7.0,1.0 / 8.0, 1.0 / 9.0,  1.0 / 10.0, 1.0 / 11.0, 1.0 / 12.0, 1.0 / 13.0,
	 1.0 / 5.0,  1.0 / 6.0,1.0 / 7.0,1.0 / 8.0,1.0 / 9.0, 1.0 / 10.0, 1.0 / 11.0, 1.0 / 12.0, 1.0 / 13.0, 1.0 / 14.0,
	 1.0 / 6.0,  1.0 / 7.0, 1.0 / 8.0,1.0 / 9.0,1.0 / 10.0, 1.0 / 11.0, 1.0 / 12.0, 1.0 / 13.0, 1.0 / 14.0, 1.0 / 15.0,
	 1.0 / 7.0,   1.0 / 8.0,1.0 / 9.0,   1.0 / 10.0,1.0 / 11.0,  1.0 / 12.0,   1.0 / 13.0, 1.0 / 14.0, 1.0 / 15.0, 1.0 / 16.0,
	 1.0 / 8.0,   1.0 / 9.0, 1.0 / 10.0,  1.0 / 11.0,1.0 / 12.0,  1.0 / 13.0,  1.0 / 14.0, 1.0 / 15.0, 1.0 / 16.0, 1.0 / 17.0,
	 1.0 / 9.0,   1.0 / 10.0,1.0 / 11.0,  1.0 / 12.0, 1.0 / 13.0, 1.0 / 14.0,   1.0 / 15.0,1.0 / 16.0, 1.0 / 17.0, 1.0 / 18.0,
	 1.0 / 10.0,  1.0 / 11.0 ,1.0 / 12.0,  1.0 / 13.0,1.0 / 14.0, 1.0 / 15.0,   1.0 / 16.0,1.0 / 17.0, 1.0 / 18.0, 1.0 / 19.0
	};
	double C[144] = {
		12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		11, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		10, 10, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,
		9,  9,  9,   9,  8,  7,  6,  5,  4,  3,  2,  1,
		8,  8,  8,   8,  8,  7,  6,  5,  4,  3,  2,  1,
		7,  7,  7,   7,  7,  7,  6,  5,  4,  3,  2,  1,
		6,  6,  6,   6,  6,  6,  6,  5,  4,  3,  2,  1,
		5,  5,  5,   5,  5,  5,  5,  5,  4,  3,  2,  1,
		4,  4,  4,   4,  4,  4,  4,  4,  4,  3,  2,  1,
		3,  3,  3,   3,  3,  3,  3,  3,  3,  3,  2,  1,
		2,  2,  2,   2,  2,  2,  2,  2,  2,  2,  2,  1,
		1,  1,  1,   1,  1,  1,  1,  1,  1,  1,  1,  1,
	};
	double* D = new double[400];
	double* E = new double[2500];
	//D、E数组初始化
	for (int i = 0; i < 20; i++) {
		for (int j = 0; j < 20; j++) {
			double num = sqrt(2.0 / 21.0) * sin((((i + 1.0) * (j + 1.0) * pi) / 21.0));
			D[i * 20 + j] = num;
		}
	}
	for (int i = 0; i < 50; i++) {
		for (int j = 0; j < 50; j++) {
			if (j < i)
				E[i * 50 + j] = -1;
			else if (i == j)
				E[i * 50 + j] = 1;
			else if (j == 49)
				E[i * 50 + j] = 1;
			else
				E[i * 50 + j] = 0;
		}
	}
	//幂法求特征值
/*	double copy[N*N];//复制数组用于求误差
	memcpy(copy, A, N*N * sizeof(double));
	double pld=0;
	double env[N] = { 1,0.6,0.3,0.2};
	clock_t start, finish;
	double  duration;
	start = clock();
	power_png(pld, env, A, N);
	finish = clock();
	duration = (float)(finish - start) / CLOCKS_PER_SEC;
	printf("运行时间%f seconds\n", duration);
	cout << pld;
	//求幂法的误差
	double c[N * N] = { 0 };
	double a[N * N] = { 0 };
	double b[N * N] = { 0 };
	for (int i = 0; i < N; i++) {
		for (int s = 0; s < N; s++) {
			c[i] += copy[N * i + s] * env[s];//原先未操作过的矩阵A乘向量
			
		}
	}
	double fenzi = 0;
	double fenmu = 0;
	for (int i = 0; i < N; i++) {
		c[i] -= pld * env[i];
		fenzi += c[i]*c[i];
		fenmu += env[i] * env[i];
	}
	fenzi = sqrt(fenzi);
	fenmu= sqrt(fenmu);
    double offset = fenzi / fenmu;
	cout << "误差=" << offset << endl;*/
    //雅克比求特征值
	clock_t start, finish;
	double  duration;
	start = clock();
	double v[N] = { 0 };
	bool re = Jacobi(v ,&A[0], N );
	if (!re) {
		cout << "矩阵特征值:" << endl;
		for (int i = 0; i < N; i++) {
			cout  << v[i]<<",";
	}
	}
	else
		cout << "false" << endl;
	finish = clock();
	duration = (float)(finish - start) / CLOCKS_PER_SEC;
	printf("运行时间%f seconds\n", duration); 
	//QR分解求特征值
	/*clock_t start, finish;
	double  duration;
	start = clock();
	int i, j;
	double Q[N * N] = { 0 };
	double R[N * N] = { 0 };
	gauss_hessen2(A, N);
	for (i = 1; i <= 2000; i++) {
		qr_eng(E, Q, R);
		Multiplicate(E, R, Q);
	}
	cout << "矩阵特征值:\n";
	for (i = 0; i < N; i++) //输出特征值
		cout <<  R[i * N + i] <<",";
	finish = clock();
	duration = (float)(finish - start) / CLOCKS_PER_SEC;
	printf("运行时间%f seconds\n", duration);*/
	system("pause");
}