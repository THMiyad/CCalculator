#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#define PI 3.14159265358979323846
#define E  2.71828182845904523536

#define C_LIGHT 299792458.0          
#define C_PLANCK 6.62607015e-34      
#define C_GRAV 6.67430e-11           
#define C_ME 9.1093837e-31           
#define C_MP 1.6726219e-27           
#define C_MN 1.6749275e-27           
#define C_ECHARGE 1.602176634e-19    
#define C_NA 6.02214076e23           
#define C_GASR 8.314462618           
#define C_BOLTZ 1.380649e-23         
#define C_FARADAY 96485.33212        
#define C_EPS0 8.854187817e-12       
#define C_MU0 1.25663706212e-6       
#define C_GOLDEN 1.618033988749895   

const char *expr;
double lastAns = 0.0;
double memoryM = 0.0;
int angleMode = 0;   
int fixDecimals = 4;

double parseExpression();
double parseTerm();
double parseFactor();

void skipSpaces() {
    while (*expr == ' ' || *expr == '\t' || *expr == ',') expr++;
}

double factorial(double n) {
    if (n < 0 || n > 170) return 0;
    long long fact = 1;
    for (int i = 1; i <= (int)n; i++) fact *= i;
    return (double)fact;
}

double nPr(double n, double r) {
    if (r > n || n < 0 || r < 0) return 0;
    return factorial(n) / factorial(n - r);
}

double nCr(double n, double r) {
    if (r > n || n < 0 || r < 0) return 0;
    return factorial(n) / (factorial(r) * factorial(n - r));
}

long long binaryToDecimal(char *binStr) {
    long long dec = 0;
    while (*binStr) {
        if (*binStr == '0' || *binStr == '1') {
            dec = dec * 2 + (*binStr - '0');
        }
        binStr++;
    }
    return dec;
}

void baseNConverter() {
    long long num = 0;
    int choice;
    printf("\n[BASE-N CONVERTER]\n");
    printf("1. Decimal (DEC)\n2. Hexadecimal (HEX)\n3. Octal (OCT)\n4. Binary (BIN)\nChoice: ");
    if (scanf("%d", &choice) != 1) { getchar(); return; }
    getchar();

    if (choice == 1) {
        printf("Enter Decimal Num: ");
        scanf("%lld", &num); getchar();
    } else if (choice == 2) {
        printf("Enter Hexadecimal Num: ");
        scanf("%llx", &num); getchar();
    } else if (choice == 3) {
        printf("Enter Octal Num: ");
        scanf("%llo", &num); getchar();
    } else if (choice == 4) {
        char binStr[65];
        printf("Enter Binary Num: ");
        if (fgets(binStr, sizeof(binStr), stdin) != NULL) {
            binStr[strcspn(binStr, "\n")] = 0;
            num = binaryToDecimal(binStr);
        }
    } else {
        printf("Invalid Choice.\n\n");
        return;
    }

    printf("\n--- RESULTS ---\n");
    printf(" DEC: %lld\n HEX: %llX\n OCT: %llo\n", num, num, num);
    printf(" BIN: ");
    for (int i = 15; i >= 0; i--) {
        int k = (num >> i) & 1;
        printf("%d", k);
        if (i % 4 == 0) printf(" ");
    }
    printf("\n\n");
}

void solveQuadratic() {
    double a, b, c;
    printf("\n[QUADRATIC: ax^2+bx+c=0]\nEnter a, b, c (space or comma separated): ");
    if (scanf("%lf %*[, ] %lf %*[, ] %lf", &a, &b, &c) != 3) {
        getchar();
        return;
    }
    getchar();
    if (a == 0) { printf("Invalid (a can't be 0).\n\n"); return; }
    double d = b*b - 4*a*c;
    if (d >= 0) {
        printf("x1 = %.*lf\nx2 = %.*lf\n\n", fixDecimals, (-b + sqrt(d)) / (2*a), fixDecimals, (-b - sqrt(d)) / (2*a));
    } else {
        printf("x1 = %.*lf + %.*lfi\nx2 = %.*lf - %.*lfi\n\n", fixDecimals, -b / (2*a), fixDecimals, sqrt(-d) / (2*a), fixDecimals, -b / (2*a), fixDecimals, sqrt(-d) / (2*a));
    }
}

void solveSimultaneous() {
    double a, b, c, d, e, f;
    printf("\n[SIMULTANEOUS EQ]\nEnter a, b, c, d, e, f: ");
    if (scanf("%lf %*[, ] %lf %*[, ] %lf %*[, ] %lf %*[, ] %lf %*[, ] %lf", &a, &b, &c, &d, &e, &f) != 6) {
        getchar();
        return;
    }
    getchar();
    double det = a * e - b * d;
    if (det == 0) { printf("No unique solution.\n\n"); return; }
    printf("x = %.*lf, y = %.*lf\n\n", fixDecimals, (c * e - b * f) / det, fixDecimals, (a * f - c * d) / det);
}

void statisticsMode() {
    int n;
    printf("\n[STATISTICS]\nEnter number of data points: ");
    if (scanf("%d", &n) != 1 || n <= 0) { getchar(); return; }
    double *arr = (double *)malloc(n * sizeof(double)), sum = 0;
    printf("Enter %d numbers (space or comma separated):\n", n);
    for (int i = 0; i < n; i++) { 
        if (i == 0) scanf("%lf", &arr[i]);
        else scanf("%*[, ] %lf", &arr[i]);
        sum += arr[i]; 
    }
    getchar();
    double mean = sum / n, varianceSum = 0;
    for (int i = 0; i < n; i++) varianceSum += pow(arr[i] - mean, 2);
    printf("Count = %d\nSum = %.*lf\nMean = %.*lf\nStd Dev = %.*lf\n\n", n, fixDecimals, sum, fixDecimals, mean, fixDecimals, sqrt(varianceSum / n));
    free(arr);
}

void matrixMode() {
    int r1, c1, r2, c2;
    printf("\n[MATRIX MODE (Up to 5x5)]\n");
    printf("Enter rows and columns for Matrix A (e.g., 3 3 or 3,3): ");
    if (scanf("%d %*[, ] %d", &r1, &c1) != 2 || r1 > 5 || c1 > 5 || r1 <= 0 || c1 <= 0) {
        printf("Invalid size (Max 5x5).\n\n");
        getchar();
        return;
    }
    
    double A[5][5], B[5][5], Res[5][5];
    printf("Enter elements of Matrix A (%dx%d):\n", r1, c1);
    for(int i = 0; i < r1; i++) {
        for(int j = 0; j < c1; j++) {
            scanf("%lf", &A[i][j]);
        }
    }

    printf("Enter rows and columns for Matrix B (e.g., 3 3 or 3,3): ");
    if (scanf("%d %*[, ] %d", &r2, &c2) != 2 || r2 > 5 || c2 > 5 || r2 <= 0 || c2 <= 0) {
        printf("Invalid size (Max 5x5).\n\n");
        getchar();
        return;
    }
    
    printf("Enter elements of Matrix B (%dx%d):\n", r2, c2);
    for(int i = 0; i < r2; i++) {
        for(int j = 0; j < c2; j++) {
            scanf("%lf", &B[i][j]);
        }
    }
    getchar();

    printf("1. Add (A + B)\n2. Multiply (A * B)\nChoice: ");
    int ch; 
    scanf("%d", &ch); 
    getchar();

    if (ch == 1) {
        if (r1 != r2 || c1 != c2) {
            printf("Error: Order mismatch for addition!\n\n");
            return;
        }
        printf("Result Matrix (%dx%d):\n", r1, c1);
        for(int i = 0; i < r1; i++) {
            printf("[ ");
            for(int j = 0; j < c1; j++) {
                Res[i][j] = A[i][j] + B[i][j];
                printf("%.*lf ", fixDecimals, Res[i][j]);
            }
            printf("]\n");
        }
        printf("\n");
    } 
    else if (ch == 2) {
        if (c1 != r2) {
            printf("Error: Column of A must match Row of B for multiplication!\n\n");
            return;
        }
        printf("Result Matrix (%dx%d):\n", r1, c2);
        for(int i = 0; i < r1; i++) {
            printf("[ ");
            for(int j = 0; j < c2; j++) {
                Res[i][j] = 0;
                for(int k = 0; k < c1; k++) {
                    Res[i][j] += A[i][k] * B[k][j];
                }
                printf("%.*lf ", fixDecimals, Res[i][j]);
            }
            printf("]\n");
        }
        printf("\n");
    } else {
        printf("Invalid choice.\n\n");
    }
}

void complexMode() {
    double r1, i1, r2, i2;
    printf("\n[COMPLEX NUMBER]\nEnter C1 (Real Imag): ");
    scanf("%lf %*[, ] %lf", &r1, &i1);
    printf("Enter C2 (Real Imag): ");
    scanf("%lf %*[, ] %lf", &r2, &i2);
    getchar();
    printf("1. Add\n2. Multiply\nChoice: ");
    int ch; scanf("%d", &ch); getchar();
    if(ch == 1) {
        printf("Result: %.*lf + %.*lfi\n\n", fixDecimals, r1 + r2, fixDecimals, i1 + i2);
    } else {
        printf("Result: %.*lf + %.*lfi\n\n", fixDecimals, (r1*r2 - i1*i2), fixDecimals, (r1*i2 + r2*i1));
    }
}

void vectorMode() {
    double a[3], b[3];
    printf("\n[VECTOR 3D]\nEnter Vector A (x y z): ");
    scanf("%lf %*[, ] %lf %*[, ] %lf", &a[0], &a[1], &a[2]);
    printf("Enter Vector B (x y z): ");
    scanf("%lf %*[, ] %lf %*[, ] %lf", &b[0], &b[1], &b[2]);
    getchar();
    double dot = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    double cx = a[1]*b[2] - a[2]*b[1];
    double cy = a[2]*b[0] - a[0]*b[2];
    double cz = a[0]*b[1] - a[1]*b[0];
    printf("Dot = %.*lf\nCross = (%.*lf, %.*lf, %.*lf)\nMagnitude A = %.*lf\n\n", 
           fixDecimals, dot, fixDecimals, cx, fixDecimals, cy, fixDecimals, cz, fixDecimals, sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]));
}

void unitConverter() {
    printf("\n[UNIT CONVERTER]\n1. Celsius to Fahrenheit\n2. Fahrenheit to Celsius\nChoice: ");
    int ch; scanf("%d", &ch); getchar();
    double val;
    if(ch == 1) {
        printf("Enter Celsius: "); scanf("%lf", &val); getchar();
        printf("Fahrenheit = %.*lf F\n\n", fixDecimals, (val * 9.0/5.0) + 32);
    } else {
        printf("Enter Fahrenheit: "); scanf("%lf", &val); getchar();
        printf("Celsius = %.*lf C\n\n", fixDecimals, (val - 32) * 5.0/9.0);
    }
}

void financialMode() {
    printf("\n[FINANCIAL MODE]\n1. Simple Interest\n2. Compound Interest\nChoice: ");
    int ch; scanf("%d", &ch); getchar();
    double p, r, t;
    printf("Enter Principal (P): "); scanf("%lf", &p);
    printf("Enter Rate %% (R): "); scanf("%lf", &r);
    printf("Enter Time in years (T): "); scanf("%lf", &t);
    getchar();
    if(ch == 1) {
        double si = (p * r * t) / 100.0;
        printf("Simple Interest = %.*lf\nTotal Amount = %.*lf\n\n", fixDecimals, si, fixDecimals, p + si);
    } else {
        double ci = p * pow((1 + r / 100.0), t) - p;
        printf("Compound Interest = %.*lf\nTotal Amount = %.*lf\n\n", fixDecimals, ci, fixDecimals, p + ci);
    }
}

void constantsMode() {
    printf("\n[PHYSICAL & MATH CONSTANTS]\n");
    printf("c    = Speed of Light (%.5e m/s)\n", C_LIGHT);
    printf("h    = Planck's Constant (%.5e J·s)\n", C_PLANCK);
    printf("G    = Gravitational Constant (%.5e)\n", C_GRAV);
    printf("me   = Electron Mass (%.5e kg)\n", C_ME);
    printf("mp   = Proton Mass (%.5e kg)\n", C_MP);
    printf("mn   = Neutron Mass (%.5e kg)\n", C_MN);
    printf("e    = Elementary Charge (%.5e C)\n", C_ECHARGE);
    printf("Na   = Avogadro Number (%.5e)\n", C_NA);
    printf("R    = Gas Constant (%.5lf)\n", C_GASR);
    printf("k    = Boltzmann Constant (%.5e)\n", C_BOLTZ);
    printf("F    = Faraday Constant (%.5lf)\n", C_FARADAY);
    printf("eps0 = Vacuum Permittivity (%.5e)\n", C_EPS0);
    printf("mu0  = Vacuum Permeability (%.5e)\n", C_MU0);
    printf("phi  = Golden Ratio (%.5lf)\n\n", C_GOLDEN);
}

double parseFactor() {
    skipSpaces();
    if (*expr == '\0') return 0;
    if (*expr == '-') { expr++; return -parseFactor(); }
    if (*expr == '+') { expr++; return parseFactor(); }
    if (*expr == '(') {
        expr++;
        double val = parseExpression();
        if (*expr == ')') expr++;
        return val;
    }
    if ((unsigned char)expr[0] == 0xE2 && (unsigned char)expr[1] == 0x88 && (unsigned char)expr[2] == 0x9A) {
        expr += 3; return sqrt(parseFactor());
    }
    if ((unsigned char)expr[0] == 0xCF && (unsigned char)expr[1] == 0x80) {
        expr += 2; return PI;
    }

    if (isalpha((unsigned char)*expr)) {
        char func[15];
        int i = 0;
        while ((isalpha((unsigned char)*expr) || isdigit((unsigned char)*expr) || *expr == '^' || *expr == '-') && i < 14) {
            func[i++] = *expr++;
        }
        func[i] = '\0';

        if (strcmp(func, "pi") == 0) return PI;
        if (strcmp(func, "e") == 0) return E;
        if (strcmp(func, "c") == 0) return C_LIGHT;
        if (strcmp(func, "h") == 0) return C_PLANCK;
        if (strcmp(func, "G") == 0) return C_GRAV;
        if (strcmp(func, "me") == 0) return C_ME;
        if (strcmp(func, "mp") == 0) return C_MP;
        if (strcmp(func, "mn") == 0) return C_MN;
        if (strcmp(func, "Na") == 0) return C_NA;
        if (strcmp(func, "R") == 0) return C_GASR;
        if (strcmp(func, "k") == 0) return C_BOLTZ;
        if (strcmp(func, "F") == 0) return C_FARADAY;
        if (strcmp(func, "eps0") == 0) return C_EPS0;
        if (strcmp(func, "mu0") == 0) return C_MU0;
        if (strcmp(func, "phi") == 0) return C_GOLDEN;

        if (strcmp(func, "Ans") == 0 || strcmp(func, "ans") == 0) return lastAns;
        if (strcmp(func, "M") == 0 || strcmp(func, "m") == 0) return memoryM;

        if (strcmp(func, "sqrt") == 0) return sqrt(parseFactor());
        if (strcmp(func, "cbrt") == 0) return cbrt(parseFactor());
        if (strcmp(func, "ln") == 0) return log(parseFactor());
        if (strcmp(func, "abs") == 0) return fabs(parseFactor());

        if (strcmp(func, "log") == 0) {
            skipSpaces();
            if (*expr == '(') {
                expr++;
                double b = parseExpression();
                if (*expr == ',') { expr++; double v = parseExpression(); if (*expr == ')') expr++; return log(v) / log(b); }
                if (*expr == ')') expr++;
                return log10(b);
            }
            return log10(parseFactor());
        }

        if (strcmp(func, "sinh") == 0) return sinh(parseFactor());
        if (strcmp(func, "cosh") == 0) return cosh(parseFactor());
        if (strcmp(func, "tanh") == 0) return tanh(parseFactor());

        if (strcmp(func, "asin") == 0 || strcmp(func, "sin-1") == 0 || strcmp(func, "sin^-1") == 0) {
            double res = asin(parseFactor());
            return (angleMode == 0) ? res * (180.0 / PI) : res;
        }
        if (strcmp(func, "acos") == 0 || strcmp(func, "cos-1") == 0 || strcmp(func, "cos^-1") == 0) {
            double res = acos(parseFactor());
            return (angleMode == 0) ? res * (180.0 / PI) : res;
        }
        if (strcmp(func, "atan") == 0 || strcmp(func, "tan-1") == 0 || strcmp(func, "tan^-1") == 0) {
            double res = atan(parseFactor());
            return (angleMode == 0) ? res * (180.0 / PI) : res;
        }

        double val = parseFactor();
        int localDeg = angleMode;
        skipSpaces();
        if ((unsigned char)expr[0] == 0xC2 && (unsigned char)expr[1] == 0xB0) { expr += 2; localDeg = 0; }
        double angleRad = (localDeg == 0) ? val * (PI / 180.0) : val;

        if (strcmp(func, "sin") == 0) return sin(angleRad);
        if (strcmp(func, "cos") == 0) return cos(angleRad);
        if (strcmp(func, "tan") == 0) return tan(angleRad);
        if (strcmp(func, "cot") == 0) return 1.0 / tan(angleRad);
        if (strcmp(func, "sec") == 0) return 1.0 / cos(angleRad);
        if (strcmp(func, "csc") == 0) return 1.0 / sin(angleRad);
    }

    char *endPtr;
    double val = strtod(expr, &endPtr);
    if (expr != endPtr) {
        expr = endPtr;
        skipSpaces();
        if ((unsigned char)expr[0] == 0xC2 && (unsigned char)expr[1] == 0xB0) expr += 2;
        else if (*expr == '%') { expr++; val /= 100.0; }
        else if (*expr == '!') { expr++; val = factorial(val); }
        else if (*expr == 'P' || *expr == 'p') { expr++; val = nPr(val, parseFactor()); }
        else if (*expr == 'C' || *expr == 'c') { expr++; val = nCr(val, parseFactor()); }
        return val;
    }
    return 0;
}

double parseTerm() {
    double val = parseFactor();
    while (1) {
        skipSpaces();
        if (isalpha((unsigned char)*expr) || *expr == '(' || 
            ((unsigned char)expr[0] == 0xE2 && (unsigned char)expr[1] == 0x88 && (unsigned char)expr[2] == 0x9A) ||
            ((unsigned char)expr[0] == 0xCF && (unsigned char)expr[1] == 0x80)) {
            val *= parseFactor();
            continue;
        }
        if (*expr == '*' || *expr == 'x' || ((unsigned char)expr[0] == 0xC3 && (unsigned char)expr[1] == 0x97)) {
            if ((unsigned char)expr[0] == 0xC3 && (unsigned char)expr[1] == 0x97) expr += 2; else expr++;
            val *= parseFactor();
        } 
        else if (*expr == '/' || ((unsigned char)expr[0] == 0xC3 && (unsigned char)expr[1] == 0xB7)) {
            if ((unsigned char)expr[0] == 0xC3 && (unsigned char)expr[1] == 0xB7) expr += 2; else expr++;
            double denom = parseFactor();
            if (denom != 0) val /= denom; else printf("\nMath Error\n");
        } 
        else if (*expr == '^' || ((unsigned char)expr[0] == 0xCB && (unsigned char)expr[1] == 0x86)) {
            if ((unsigned char)expr[0] == 0xCB && (unsigned char)expr[1] == 0x86) expr += 2; else expr++;
            val = pow(val, parseFactor());
        }
        else break;
    }
    return val;
}

double parseExpression() {
    double val = parseTerm();
    while (1) {
        skipSpaces();
        if (*expr == '+') { 
            expr++; 
            double nextVal = parseTerm();
            skipSpaces();
            if (*expr == '%') { expr++; nextVal = val * (nextVal / 100.0); }
            val += nextVal; 
        }
        else if (*expr == '-') { 
            expr++; 
            double nextVal = parseTerm();
            skipSpaces();
            if (*expr == '%') { expr++; nextVal = val * (nextVal / 100.0); }
            val -= nextVal; 
        }
        else break;
    }
    return val;
}

void showModeMenu() {
    int ch;
    do {
        printf("\n=== MENU ===\n");
        printf("1. Angle Mode (Deg/Rad)\n");
        printf("2. Decimal Fix (Decimals: %d)\n", fixDecimals);
        printf("3. Memory Control (M: %.*lf)\n", fixDecimals, memoryM);
        printf("4. Base-N (Dec/Hex/Oct/Bin)\n");
        printf("5. Quadratic Eq\n");
        printf("6. Simultaneous Eq\n");
        printf("7. Statistics\n");
        printf("8. Matrix Mode (Up to 5x5)\n");
        printf("9. Complex Number\n");
        printf("10. Vector (3D)\n");
        printf("11. Unit Converter\n");
        printf("12. Financial (Interest)\n");
        printf("13. Physical Constants\n");
        printf("0. Return to Main Interface\n");
        printf("Choice: ");
        
        if (scanf("%d", &ch) != 1) { getchar(); break; }
        getchar();

        if (ch == 1) {
            printf("1: Deg  2: Rad -> ");
            int a; scanf("%d", &a); getchar();
            angleMode = (a == 2) ? 1 : 0;
            printf("Mode: %s\n\n", angleMode == 0 ? "DEG" : "RAD");
        } else if (ch == 2) {
            printf("Enter decimal places (0 to 6): ");
            scanf("%d", &fixDecimals); getchar();
            if (fixDecimals < 0) fixDecimals = 0;
            if (fixDecimals > 6) fixDecimals = 6;
            printf("Decimals set to %d\n\n", fixDecimals);
        } else if (ch == 3) {
            printf("\n[MEMORY TOOLS]\n");
            printf("1. Store Ans to M\n2. Add Ans to M (M+)\n3. Clear Memory (MC)\nChoice: ");
            int mch; scanf("%d", &mch); getchar();
            if (mch == 1) { memoryM = lastAns; printf("Stored. M = %.*lf\n\n", fixDecimals, memoryM); }
            else if (mch == 2) { memoryM += lastAns; printf("Added. M = %.*lf\n\n", fixDecimals, memoryM); }
            else if (mch == 3) { memoryM = 0.0; printf("Memory Cleared (M = 0)\n\n"); }
        } else if (ch == 4) baseNConverter();
        else if (ch == 5) solveQuadratic();
        else if (ch == 6) solveSimultaneous();
        else if (ch == 7) statisticsMode();
        else if (ch == 8) matrixMode();
        else if (ch == 9) complexMode();
        else if (ch == 10) vectorMode();
        else if (ch == 11) unitConverter();
        else if (ch == 12) financialMode();
        else if (ch == 13) constantsMode();
        else if (ch == 0) {
            printf("Returning to Main Interface...\n\n");
            break;
        }
    } while (ch != 0);
}

int main() {
    char input[200];
    printf("=========================================\n");
    printf("     UNIVERSAL SCIENTIFIC CALCULATOR     \n");
    printf("=========================================\n");
    printf("  * Math  : +, -, *, /, x, ÷, ^, ˆ, %, ()\n");
    printf("  * Trig  : sin, cos, tan, cot, sec, csc\n");
    printf("  * Inv   : asin, acos, atan, sin^-1\n");
    printf("  * Root  : sqrt, √, cbrt\n");
    printf("  * Sym   : pi, π, e, Ans, M, °\n");
    printf("  * Menu  : 'm' (Setup/Tools)\n");
    printf("=========================================\n\n");

    while (1) {
        printf(angleMode == 0 ? "D " : "R ");
        if (fgets(input, sizeof(input), stdin) == NULL) break;
        input[strcspn(input, "\n")] = 0;

        if (strcmp(input, "off") == 0) break;
        if (strcmp(input, "m") == 0) {
            showModeMenu();
            continue;
        }
        if (strlen(input) == 0) continue;

        expr = input;
        lastAns = parseExpression();
        printf("Ans = %.*lf\n\n", fixDecimals, lastAns);
    }
    return 0;
}
    
