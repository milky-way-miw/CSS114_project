#include <stdio.h>

//1.
void gaussElimination(int n, double A[n][n], double b[n]){
double x[n];
for (int i = 0; i < n - 1; i++){ // ไล่ทีละคอลัมน์ (pivot column)
    int maxRow = i; // สมมติว่าแถว i เป็น pivot ก่อน
    // หาแถวที่มีค่ามากที่สุดในคอลัมน์ i
    for (int r = i + 1; r < n; r++){
        if (A[r][i] * A[r][i] > A[maxRow][i] * A[maxRow][i]) {
            maxRow = r; // ถ้าใหญ่กว่า เปลี่ยน pivot
        }
    }
    // สลับแถวถ้า pivot ไม่ได้อยู่แถวบน
    if (maxRow != i){
        for (int c = 0; c < n; c++){
            double temp = A[i][c];
            A[i][c] = A[maxRow][c];
            A[maxRow][c] = temp;
        }
        double tempb = b[i];
        b[i] = b[maxRow];
        b[maxRow] = tempb;        // สลับค่าฝั่งขวา b ให้ตรงกัน
    }
    // ตรวจสอบว่า pivot เป็น 0 หรือไม่
    if (A[i][i] == 0){
        printf("pivot is 0\n");    // ไม่สามารถคำนวณต่อได้
        return;                   // ออกจากฟังก์ชัน
    }

    // กำจัดค่าด้านล่าง pivot ให้เป็นศูนย์
    for (int j = i + 1; j < n; j++){
        double factor = A[j][i] / A[i][i];  // ตัวคูณในการกำจัด
        for (int k = i; k < n; k++){
            A[j][k] = A[j][k] - (factor * A[i][k]); //แถวล่าง = แถวล่าง − factor × แถว pivot
        }
        b[j] = b[j] - (factor * b[i]);  // ปรับฝั่งขวา
    }
}
for (int i = n - 1; i >= 0; i--){
    x[i] = b[i];// เริ่มจากค่าฝั่งขวา
    for (int j = i + 1; j < n; j++){
        x[i] = x[i] - (A[i][j] * x[j]);
    }
    x[i] = x[i] / A[i][i]; // หารด้วยสัมประสิทธิ์ของตัวแปร
}
for (int i = 0; i < n; i++){
    printf("x%d = %.4f\n", i + 1, x[i]);
}
}

//2.
void GJE(int n, double A[n][n], double num[n]) {
    int r, c;
    double temp;

    // เปลี่ยนค่านอกจากแนวทแยงให้เป็น 0
    for (c = 0; c < n; c++) {
        for (r = 0; r < n; r++) {
            if (c != r && A[r][c] != 0) {
                temp = A[r][c] / A[c][c];
                num[r] = num[r] - (temp * num[c]);

                for (int k = 0; k < n; k++) {
                    A[r][k] = A[r][k] - (temp * A[c][k]);
                }
            }
        }
    }

    // แสดงผล
    for (int x = 0; x < n; x++) {
        printf("X[%d] = %.2f\n", x + 1, num[x] / A[x][x]);
    }
}

//3.
void LU(int n, double U[n][n], double num[n]){
    int r, c;
    double L[n][n];
    double temp;
    
    //เปลี่ยนค่าให้เป็น 0
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
    
    //แสดงผล ของเมทริก L
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

    //แสดงผล ของเมทริก U
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
    //จากสูตร Ly = b
    for (int i = 0; i < n; i++) {
        double sum_L = 0.0;
        for (int j = 0; j < i; j++) {
            sum_L = sum_L + L[i][j] * y[j];
        } 
        y[i] = num[i] - sum_L;
    }
    
    //จากสูตร Ux = y
    for (int i = n-1; i >= 0; i--) {
        double sum_U = 0.0;
        for (int j = i+1; j < n; j++) {
            sum_U = sum_U + U[i][j] * x[j];
        } 
        x[i] = (y[i] - sum_U) / U[i][i];
    }
    
    //แสดงผล
    for (int i = 0; i < n; i++) {
        printf("x[%d] = %.2f\n", i+1, x[i]);
    }
}


int main()
{
    int n;
    //รับค่าขนาดของเมทริก A
    printf("Enter matrix size n: ");
    scanf("%d", &n);

    double A[n][n], num[n];
    //รับค่าเมทริก A
    printf("Enter matrix A (%d x %d):\n", n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            scanf("%lf", &A[i][j]);
        }
    }

    //รับค่าเมริก num
    printf("Enter num (%d values):\n", n);
    for (int i = 0; i < n; i++)
    {
        scanf("%lf", &num[i]);
    }
    printf("\n");

    //เรียกใช้ฟังก์ชัน GE
    printf("GE is:\n");
    gaussElimination(n, A, num);
    printf("\n");

    //เรียกใช้ฟังก์ชัน GJE
    printf("GJE is:\n");
    GJE(n, A, num);
    printf("\n");

    //เรียกใช้ฟังก์ชัน LU
    printf("LU is:\n");
    LU(n, A, num);
    return 0;
}