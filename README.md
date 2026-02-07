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

# 4. Inverse Matrix
#include <stdio.h>

int main() {
    int b_matrix = 0;
    double A[3][3] ,b[3];
    double det, inv[3][3], x[3];
    
    // รับค่าเมทริก A
    printf("Enter matrix A (3 x 3):\n");
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            scanf("%lf", &A[i][j]);
        }
    }
    
    // เช็คว่ามีเมทริก b มั้ย
    printf("Enter 1 if have b matrix or 0 if dont have b matrix :");
    scanf("%d",&b_matrix);

    // ถ้ามีให้รับค่าเมทริก b
    if(b_matrix == 1){
        printf("Enter b matrix (3 values):\n");
        for (int i = 0; i < 3; i++){
            scanf("%lf", &b[i]);
        }
    printf("\n");
    }
    
    // คำนวณ determinant ของ A
    det = A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
        - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
        + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);

    if(det == 0) {
        printf("Dont have Inverse\n");
        return 0;
    }

    // หา cofactor matrix
    double cof[3][3];
    cof[0][0] =  (A[1][1]*A[2][2] - A[1][2]*A[2][1]);
    cof[0][1] = -(A[1][0]*A[2][2] - A[1][2]*A[2][0]);
    cof[0][2] =  (A[1][0]*A[2][1] - A[1][1]*A[2][0]);

    cof[1][0] = -(A[0][1]*A[2][2] - A[0][2]*A[2][1]);
    cof[1][1] =  (A[0][0]*A[2][2] - A[0][2]*A[2][0]);
    cof[1][2] = -(A[0][0]*A[2][1] - A[0][1]*A[2][0]);

    cof[2][0] =  (A[0][1]*A[1][2] - A[0][2]*A[1][1]);
    cof[2][1] = -(A[0][0]*A[1][2] - A[0][2]*A[1][0]);
    cof[2][2] =  (A[0][0]*A[1][1] - A[0][1]*A[1][0]);

    // adjoint
    double adj[3][3];
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            adj[i][j] = cof[j][i];

    // inverse = adj / det
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            inv[i][j] = adj[i][j] / det;

    // คูณ A^-1 กับ b
    if(b_matrix == 1){
        for(int i=0;i<3;i++) {
        x[i] = 0;
            for(int j=0;j<3;j++) {
                x[i] += inv[i][j] * b[j];
            }
        } 
        printf("Answer x:\n");
        for(int i=0;i<3;i++) {
        printf("x[%d] = %.4f\n", i, x[i]);
        }
    }

    // แสดงผล
    printf("Inverse Matrix =\n");
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            printf("%.2f ", inv[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    return 0;
}