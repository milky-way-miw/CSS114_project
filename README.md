# file >> code.c
    #include <stdio.h>
## 1. Gauss Elimination with pivoting
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
    
## 2. Gauss Jordan Elimination
    void GJE(int n, double A[n][n], double num[n]) {
        int r, c;
        double temp;

        // เปลี่ยนค่านอกจากแนวทแยงให้เป็น 0
        for (c = 0; c < n; c++) {
            for (r = 0; r < n; r++) {
                if (c != r && A[r][c] != 0) {
                    temp = A[r][c] / A[c][c];   //ค่าส่วนต่าง
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

## 3. LU Factorization
    void LU(int n, double U[n][n], double num[n]){
        int r, c;
        double L[n][n];
        double temp;
    
        //เปลี่ยนค่าให้เป็น 0
        for (c = 0; c < n; c++)
        {
            for (r = 0; r < n; r++)
            {
                temp = U[r][c] / U[c][c];   //ค่าส่วนต่าง
            
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

## Main
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

# file >> inverse.c
## 4. Inverse Matrix
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