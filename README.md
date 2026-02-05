    #include <stdio.h>
# 1. Gauss Elimination with pivoting
    น้ำฟ้าเพิ่มโค้ดตงนี้น้าค่าบ
    
# 2. Gauss Jordan Elimination
    void GJE(int n, double A[n][n], double num[n]){
        int r, c;    //ใช้วนแถวกับหลัก r = row, c = column
        double temp;    
        for (c = 0; c < n; c++)
        {
            for (r = 0; r < n; r++)    //เพราะว่าเราต้องดูจากหลักก่อน
            {
                if (c != r && A[r][c] != 0)    //เราจะไม่เปลี่ยนค่าของตัวเองแนวทแยงให้เป็นค่า 0 และหน้าตัวเลขที่เราจะเปลี่ยนเป็นค่า 0 อยู่แล้วก็ไม่ต้องทำอะไร
                {
                    temp = A[r][c] / A[c][c];    //ใช้เก็บค่าส่วนต่างของตัวที่เราเลือกกับตัวที่เราจะให้มันเป็นค่า 0
                    num[r] = num[r] - (temp * num[c]);

                    for (int k = 0; k < n; k++)    //วนเผื่อลบค่าทั้งแถว
                    {  
                        A[r][k] = A[r][k] - (temp * A[c][k]);    //เอาค่าเมทริกแถวที่เลือกมาคูณกับtemp แล้วลบด้วยเราจะทำให้เป็นค่า 0
                    }
                }
            }
         }
    
        for (int x = 0; x < n; x++)    //เป็นการวนซ้ำเผื่อแสดงคำตอบ
        {
            printf("X[%d] = %.2f\n", x + 1, num[x] / A[x][x]);
        }
    }

# 3. LU Factorization
    void LU(int n, double U[n][n], double num[n]){
        int r, c;
        double L[n][n];
        double temp;
        for (c = 0; c < n; c++)
        {
            for (r = 0; r < n; r++)
            {
                temp = U[r][c] / U[c][c];
            
                if (r > c && U[r][c] != 0)
                {
                    L[r][c] = temp;
                    for (int k = 0; k < n; k++)    //วนเผื่อลบค่าทั้งแถว
                    {  
                        U[r][k] = U[r][k] - (temp * U[c][k]);
                    }
                } else if r < c && U[r][c] != 0){
                    L[r][c] = 0;
                } else if (r == c){
                    L[r][c] = 1;
                }
            }
        }
    
        printf("L = \n");
        for (int r = 0; r < n; r++)    //วนซ้ำเผื่อ print เมทริก L
        {
            for (int c = 0; c < n; c++)
            {
                printf("%.2f ", L[r][c]);
            }
            printf("\n");
        }
        printf("\n");
    
        printf("U =\n");    //วนซ้ำเผื่อ print เมทริก U
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
            double sum_L = 0;
            for (int j = 0; j < i; j++) {
                sum_L = sum_L + L[i][j] * y[j];
            } 
            y[i] = num[i] - sum_L;
        }
    
        for (int i = n-1; i >= 0; i--) {
            double sum_U = 0;
            for (int j = i+1; j < n; j++) {
                sum_U = sum_U + U[i][j] * x[j];
            }
            x[i] = (y[i] - sum_U) / U[i][i];
        }
    
        for (int i = 0; i < n; i++) {
            printf("x[%d] = %.2f\n", i+1, x[i]);
        }
    }
    
# 4. Inverse Matrix
    น้ำฟ้าเพิ่มโค้ดตงนี้น้าค่าบ

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
        GJE(n, A, num);    เรียกใฃ้ฟังก์ชัน GJE

        printf("LU is:\n");
        LU(n, A, num);    เรียกใช้ฟิงก์ชัน LU
        return 0;
    }
