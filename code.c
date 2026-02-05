#include <stdio.h>

//2.
void GJE(int n, double A[n][n], double num[n])
{
    int r, c, k;
    double temp, pivot;
    for (c = 0; c < n; c++)
    {
        pivot = A[c][c];
        for (k = 0; k < n; k++)
        {
            A[c][k] = A[c][k] / pivot;
        }
        num[c] = num[c] / pivot;
        
        for (r = 0; r < n; r++)
        {
            if (c != r && A[r][c] != 0)
            {
                temp = A[r][c];
                num[r] = num[r] - (temp * num[c]);

                for (int k = 0; k < n; k++)
                {  
                    A[r][k] = A[r][k] - (temp * A[c][k]);
                }
            }
        }
    }
    
    for (int x = 0; x < n; x++)
    {
        printf("X[%d] = %.2f\n", x + 1, x[x]);
    }
}

//3.
void LU(int n, double U[n][n], double num[n]){
    int r, c;
    double L[n][n];
    double temp;
    for (c = 0; c < n; c++)
    {
        for (r = 0; r < n; r++)
        {
            temp = U[r][c] / U[c][c];
            
            if (c < r && U[r][c] != 0)
            {
                L[r][c] = temp;
                for (int k = 0; k < n; k++)
                {  
                    U[r][k] = U[r][k] - (temp * U[c][k]);
                }
            } else if (c > r && U[r][c] != 0){
                L[r][c] = 0;
            } else if (c == r){
                L[r][c] = 1;
            }
        }
    }
    
    printf("L = \n");
    for (int r = 0; r < n; r++)
    {
        for (int c = 0; c < n; c++)
        {
            printf("%.2f ", L[r][c]);
        }
        printf("\n");
    }
    printf("\n");
    
    printf("U =\n");
    for (int r = 0; r < n; r++)
    {
        for (int c = 0; c < n; c++)
        {
            printf("%.2f ", U[r][c]);
        }
        printf("\n");
    }
    printf("\n");
    
    double x[n], y[n];
    for (int i = 0; i < n; i++) {
        double sum_L = 0.0;
        for (int j = 0; j < i; j++) {
            sum_L += L[i][j] * y[j];
        } 
        y[i] = num[i] - sum_L;
    }
    
    for (int i = n-1; i >= 0; i--) {
        double sum_U = 0.0;
        for (int j = i+1; j < n; j++) {
            sum_U += U[i][j] * x[j];
        } 
        x[i] = (y[i] - sum_U) / U[i][i];
    }
    
    for (int i = 0; i < n; i++) {
        printf("x[%d] = %.2f\n", i+1, x[i]);
    }
}


int main()
{
    int n;
    printf("Enter matrix size n: ");
    scanf("%d", &n);

    double A[n][n], num[n];

    printf("Enter matrix A (%d x %d):\n", n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            scanf("%lf", &A[i][j]);
        }
    }

    printf("Enter num (%d values):\n", n);
    for (int i = 0; i < n; i++)
    {
        scanf("%lf", &num[i]);
    }
    printf("\n");

    printf("GJE is:\n");
    GJE(n, A, num);

    printf("LU is:\n");
    LU(n, A, num);
    return 0;
}