    #include <stdio.h>
# 1. Gauss Elimination with pivoting
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

        printf("GE is:\n");
        gaussElimination(n, A, num); //เรียกใฃ้ฟังก์ชัน GE
        
        printf("GJE is:\n");
        GJE(n, A, num);    //เรียกใฃ้ฟังก์ชัน GJE

        printf("LU is:\n");
        LU(n, A, num);    //เรียกใช้ฟิงก์ชัน LU
        return 0;
    }
