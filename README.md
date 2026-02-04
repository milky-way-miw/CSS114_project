    #include <stdio.h>
# 2. GJE
    void GJE(int n, float A[n][n], float num[n]){
        int r, c;
        float temp;
        for (c = 0; c < n; c++)
        {
            for (r = 0; r < n; r++)
            {
                if (c != r && A[r][c] != 0)
                {
                    temp = A[r][c] / A[c][c];
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
            printf("X[%d] = %.2f\n", x + 1, num[x] / A[x][x]);
        }
    }

# Main
 
    int main()
    {
        int n;
        printf("Enter matrix size n: ");
        scanf("%d", &n);

        float A[n][n], num[n];

        printf("Enter matrix A (%d x %d):\n", n, n);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                scanf("%f", &A[i][j]);
            }
        }

        printf("Enter num (%d values):\n", n);
        for (int i = 0; i < n; i++)
        {
            scanf("%f", &num[i]);
        }
        printf("\n");

        printf("GJE is:\n");
        GJE(n, A, num);
        return 0;
    }
