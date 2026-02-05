    #include <stdio.h>
# 1. Gauss Elimination with pivoting
    น้ำฟ้าเพิ่มโค้ดตงนี้น้าค่าบ
    
# 2. Gauss Jordan Elimination
    void GJE(int n, float A[n][n], float num[n]){
        int r, c;    ใช้วนแถวกับหลัก r = row, c = column
        float temp;    
        for (c = 0; c < n; c++)
        {
            for (r = 0; r < n; r++)    เพราะว่าเราต้องดูจากหลักก่อน
            {
                if (c != r && A[r][c] != 0)    เราจะไม่เปลี่ยนค่าของตัวเองแนวทแยงให้เป็นค่า 0 และหน้าตัวเลขที่เราจะเปลี่ยนเป็นค่า 0 อยู่แล้วก็ไม่ต้องทำอะไร
                {
                    temp = A[r][c] / A[c][c];    ใช้เก็บค่าส่วนต่างของตัวที่เราเลือกกับตัวที่เราจะให้มันเป็นค่า 0
                    num[r] = num[r] - (temp * num[c]);

                    for (int k = 0; k < n; k++)    วนเผื่อลบค่าทั้งแถว
                    {  
                        A[r][k] = A[r][k] - (temp * A[c][k]);    เอาค่าเมทริกแถวที่เลือกมาคูณกับtemp แล้วลบด้วยเราจะทำให้เป็นค่า 0
                    }
                }
            }
         }
    
        for (int x = 0; x < n; x++)    เป็นการวนซ้ำเผื่อแสดงผล คำตอบ
        {
            printf("X[%d] = %.2f\n", x + 1, num[x] / A[x][x]);
        }
    }

# 3. LU Factorization

    
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
        return 0;
    }
