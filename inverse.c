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