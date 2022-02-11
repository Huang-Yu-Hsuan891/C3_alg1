#include <stdio.h>
#include <stdlib.h>

int main() {
    // for H1
    int n1, rc1;
    int e1;
    int N1, RC1;
    // for H2
    int n2, rc2;
    int e2;
    int N2, RC2;
    //for H3
    int n3, rc3;
    int e3;
    int N3, RC3;

    int i, j, m, k;

    int **H1Cmask;
    int H1Cmaskrow, H1Cmaskcolumn;
    int **H2Cmask;
    int H2Cmaskrow, H2Cmaskcolumn;
    int **H3Cmask;
    int H3Cmaskrow, H3Cmaskcolumn;
    int **H1C;
    int H1Crow, H1Ccolumn;
    int **H2C;
    int H2Crow, H2Ccolumn;
    int **H3C;
    int H3Crow, H3Ccolumn;


    int **H;
    int Hrow, Hcolumn;
    int **H1;
    int H1row, H1column;


    FILE *fpr;
    fpr = fopen("H1.txt", "r");
    fscanf(fpr,"%d",&e1);
    fscanf(fpr, "%d",&n1);
    fscanf(fpr, "%d",&rc1);
    fscanf(fpr, "%d",&RC1);
    fscanf(fpr, "%d",&N1);
    N1 = N1 * 2;
    H1Cmaskrow = RC1;
    H1Cmaskcolumn = N1;
    H1Cmask = (int **)malloc(H1Cmaskrow * sizeof(int *));
    for (i = 0; i < H1Cmaskrow; i++) H1Cmask[i] = (int *)malloc(H1Cmaskcolumn * sizeof(int));
    for (i = 0; i < RC1; i++) {
        for (j = 0; j < N1; j++) {
            fscanf(fpr, "%d", &H1Cmask[i][j]);
        }
    }
    fclose(fpr);

    FILE *fpr2;
    fpr2 = fopen("H2.txt", "r");
    fscanf(fpr2,"%d",&e2);
    fscanf(fpr2, "%d",&n2);
    fscanf(fpr2, "%d",&rc2);
    fscanf(fpr2, "%d",&RC2);
    fscanf(fpr2, "%d",&N2);
    N2 = N2 * 2;
    H2Cmaskrow = RC2;
    H2Cmaskcolumn = N2;
    H2Cmask = (int **)malloc(H2Cmaskrow * sizeof(int *));
    for (i = 0; i < H2Cmaskrow; i++) H2Cmask[i] = (int *)malloc(H2Cmaskcolumn * sizeof(int));
    for (i = 0; i < RC2; i++) {
        for (j =0; j < N2; j++) {
            fscanf(fpr2, "%d", &H2Cmask[i][j]);
        }
    }
    fclose(fpr2);

    FILE *fpr3;
    fpr3 = fopen("H3.txt", "r");
    fscanf(fpr3,"%d",&e3);
    fscanf(fpr3, "%d",&n3);
    fscanf(fpr3, "%d",&rc3);
    fscanf(fpr3, "%d",&RC3);
    fscanf(fpr3, "%d",&N3);
    N3 = N3 * 2;
    H3Cmaskrow = RC3;
    H3Cmaskcolumn = N3;
    H3Cmask = (int **)malloc(H3Cmaskrow * sizeof(int *));
    for (i = 0; i < H3Cmaskrow; i++) H3Cmask[i] = (int *)malloc(H3Cmaskcolumn * sizeof(int));
    for (i = 0; i < RC3; i++) {
        for (j =0; j < N3; j++) {
            fscanf(fpr3, "%d", &H3Cmask[i][j]);
        }
    }
    fclose(fpr3);
   
    H1Crow = rc1;
    H1Ccolumn = n1;
    H1C = (int **)malloc(H1Crow * sizeof(int *));
    for (i = 0; i < H1Crow; i++) H1C[i] = (int *)malloc(H1Ccolumn * sizeof(int));
    N1 = N1 / 2;
    for (i = 0; i < H1Crow; i++) {
        for (j = 0; j < H1Ccolumn; j++) {
            H1C[i][j] = 0;
        }
    }
    int temp;
    for (i = 0; i < RC1; i++) {
        for (j = 0; j < N1; j++) {
            if (H1Cmask[i][2 * j] == 0) {
                printf("no\n");
                continue;
            }
            if (H1Cmask[i][2 * j] == 1) {
                printf("yes\n");
                for (m = 0; m < e1; m++) {
                   H1C[e1 * i + m][e1 * j + (m + H1Cmask[i][2 * j + 1]) % e1] = 1; 
                }    
            }
               
        }
    }
    H2Crow = rc2;
    H2Ccolumn = n2;
    H2C = (int **)malloc(H2Crow * sizeof(int *));
    for (i = 0; i < H2Crow; i++) H2C[i] = (int *)malloc(H2Ccolumn * sizeof(int));
    N2 = N2 / 2;
    for (i = 0; i < H2Crow; i++) {
        for (j = 0; j < H2Ccolumn; j++) {
            H2C[i][j] = 0;
        }
    }
    for (i = 0; i < RC2; i++) {
        for (j = 0; j < N2; j++) {
            if (H2Cmask[i][2 * j] == 0) {
                printf("no\n");
                continue;
            }
            if (H2Cmask[i][2 * j] == 1) {
                printf("yes\n");
                for (m = 0; m < e2; m++) {
                   H2C[e2 * i + m][e2 * j + (m + H2Cmask[i][2 * j + 1]) % e2] = 1; 
                }    
            }
               
        }
    }
    H3Crow = rc3;
    H3Ccolumn = n3;
    H3C = (int **)malloc(H3Crow * sizeof(int *));
    for (i = 0; i < H3Crow; i++) H3C[i] = (int *)malloc(H3Ccolumn * sizeof(int));
    N3 = N3 / 2;
    for (i = 0; i < H3Crow; i++) {
        for (j = 0; j < H3Ccolumn; j++) {
            H3C[i][j] = 0;
        }
    }
    for (i = 0; i < RC3; i++) {
        for (j = 0; j < N3; j++) {
            if (H3Cmask[i][2 * j] == 0) {
                printf("no\n");
                continue;
            }
            if (H3Cmask[i][2 * j] == 1) {
                printf("yes\n");
                for (m = 0; m < e3; m++) {
                   H3C[e3 * i + m][e3 * j + (m + H3Cmask[i][2 * j + 1]) % e3] = 1; 
                }    
            }
               
        }
    }

    Hrow = e1 * RC1 + e2 * RC2/*+e3 * RC3*/;
    if (Hrow == (rc3 + rc2)) printf("yes row = %d\n", Hrow);
    Hcolumn = e1 * N1 + e2 * N2;
    if (Hcolumn == (n1 + n2)) printf("yes column = %d\n", Hcolumn);
    H = (int **)malloc(Hrow * sizeof(int *));
    for (i = 0;i < Hrow; i++) H[i] = (int *)malloc(Hcolumn * sizeof(int));

    H1row = Hrow;
    H1column = Hcolumn;
    H1 = (int **)malloc(H1row * sizeof(int *));
    for (i = 0; i < H1row; i++) H1[i] = (int *)malloc(H1column  * sizeof(int));

    for (i = 0; i < Hrow; i++) {
        for (j = 0; j < Hcolumn; j++) {
            H[i][j] = 0;
        }
    }
    printf("H1Crow = %d\n", H1Crow);
    for (i = 0; i < H1Crow; i++) {
        for (j = 0; j < H1Ccolumn; j++) {
            if (H1C[i][j] == 1){
                H[i][j] = H1C[i][j];
                //H[rc + i][n + j] = H1C[i][j];
            }
        }
    }
    for (i = 0; i < H2Crow; i++) {
        for (j = 0; j < H2Ccolumn; j++) {
            if (H2C[i][j] == 1) H[i][n1 + j] = H2C[i][j];
        }
    }
    for (i = 0; i < H3Crow; i++) {
        for (j = 0; j < H3Ccolumn; j++) {
            if (H3C[i][j] == 1){
                H[rc1 + i][n1 + j] = H3C[i][j];
            }
        }
    }
    int dc1 = 3;
    int dc2 = 6;
    int dv1 = 12;
    int dv2 = 5;
    int dv3 = 6;
    for (i = 0; i < Hrow; i++) {
        for (j = 0; j < Hcolumn; j++) 
            H1[i][j] = H[i][j];
    }
    FILE *outfp1;
    outfp1 = fopen("C3parchematrix.txt","w");
    fprintf(outfp1,"%d ",Hcolumn);
    fprintf(outfp1,"%d ",Hrow);
    fprintf(outfp1,"\n");
    fprintf(outfp1,"12 5 6 3 6 ");
    //fprintf(outfp1,"%d ", num1/(Hcolumn/2));
    //fprintf(outfp1,"%d ", num2/(Hcolumn/2));
    fprintf(outfp1,"\n");
    for (i = 0; i < Hcolumn; i++) {
        for (j = 0; j < Hrow; j++) {
            if (H[j][i] == 1) fprintf(outfp1,"%d ",j+1);
        }
        fprintf(outfp1,"\n");
    }
    for (i = 0; i < Hrow; i++) {
        for (j = 0; j < Hcolumn; j++) {
            if (H[i][j] == 1) fprintf(outfp1,"%d ", j+1);
        }
        fprintf(outfp1, "\n");
    }
    fclose(outfp1);
    // GAUSS JORDAN METHOD
    int temprow = 0;    // for store pivot row 
    int stop = 0;
    int temptrans;
    int s = 0;
    printf("gauss jordan form!\n");
    for (i = 0; i < 4936 + s && temprow < 4936/*rc - s*//*n*/; i++) {                        // find pivot
        stop = 0;
        //printf("i = %d\n", i);
        for (j = temprow/*0*/; j < 4936 && stop == 0; j++) {      // where the pivot in the row
            if (H[j][i]  ==  1 && j == temprow) {
                temprow = j;
                //printf("j = %d temprow = %d \n", j, temprow);
                stop = 1;
            }
            else if (H[j][i]  ==  1 && j != temprow) {
                stop = 1;
                //printf("j = %d temprow = %d\n", j, temprow);
                for (k = 0; k < 9872; k++) {
                    temptrans = H[j][k];
                    H[j][k] = H[temprow][k];
                    H[temprow][k] = temptrans;
                }
            }
        }
        if(j == 4936 && stop == 0) {
            //i++;
            s++;
            printf("s = %d\n",s);
            continue;
        }
        // Eliminate the column have 1
        printf(" i = %d\n", i);
        for (j = 0; j < 4936; j++) {
            if (j == temprow) continue;                    // skip the original row because the row no need Elimination
            else {
                //printf("h\n");
                if (H[j][i] == 1) {
                    //printf("h\n");
                    for (k = 0; k < 9872; k++) {
                        H[j][k] = (H[j][k] + H[temprow][k]) % 2;
                    }
                }
            }
        }
        temprow++;
        //printf("temprow = %d\n", temprow);
    }
    printf("temprow = %d\n", temprow);
    int a = 0;
    int a1 = 0;
    int c[4936];
    int c1[4936];
    m = 0;
    int **Hsyst;
    int Hsystrow = 4936;
    int Hsystcolumn = 4936;
    int **Hsyst1;
    int Hsyst1row = 4936;
    int Hsyst1column = 4936;
    int **Gsys;
    int Gsysrow = 4936;
    int Gsyscolumn = 9872;
    int **G;
    int Grow = 4936;
    int Gcolumn = 9872;

    Hsyst = (int **)malloc(Hsystrow * sizeof(int *));
    for (i = 0; i < Hsystrow; i++) Hsyst[i] = (int *)malloc(Hsystcolumn * sizeof(int));
    Hsyst1 = (int **)malloc(Hsyst1row * sizeof(int *));
    for (i = 0; i < Hsyst1row; i++) Hsyst1[i] = (int *)malloc(Hsyst1column * sizeof(int));
    Gsys = (int **) malloc (Gsysrow * sizeof(int *));
    for (i = 0; i < Gsysrow; i++) Gsys[i] = (int *) malloc(Gsyscolumn *sizeof(int));
    G = (int **)malloc(Grow * sizeof(int *));
    for (i = 0; i < Grow; i++) G[i] = (int *)malloc(Gcolumn * sizeof(int));
    int q = 0;
    int q1 = 0;
    //int check[4936];
    for (j = 0; j < 9872; j++) {
        m = 0;
        q1 = 0;
        printf("j = %d ", j);
        for (i = 0; i <4936; i++) {
            if (H[i][j] == 0) m = m + 1;
            if (H[i][j] == 1) q1 = i;
        }
        printf("m = %d; a = %d; a1 = %d ;", m, a, a1);
        if (m == 4936 - 1) {
            if (q1 >= q) {
                for (i = 0; i < 4936; i++) {
                    Hsyst[i][a] = H[i][j]; 
                }
                c[a] = j;
                a++;
                q++;
            } else {
                for (i = 0; i < 4936; i++) {
                    Hsyst1[i][a1] = H[i][j];
                }
                c1[a1] = j;
                a1++;
            }
            //c[a] = j;
            //a++;
            //q++;
        }
        else {
            for (i = 0; i < 4936; i++) {
                Hsyst1[i][a1] = H[i][j];
            }
            c1[a1] = j;
            a1++;
        }
    }
    printf("yes!");
    for (i = 0; i < 4936; i++) printf("%d ", c[i]);
    printf("\n");
    for (i = 0; i < 4936; i++) printf("%d ",c1[i]);
    printf("\n");
    m = 0;
    for (j = 0; j < 9872; j++) {
        if (j < 4936) {
            for (i = 0; i < 4936; i++) {
                Gsys[i][j] = Hsyst[i][j];
            }
        } else {
                for (k = 4936; k < 9872; k++) {
                    Gsys[m][k] = Hsyst1[k-4936][j-4936];
                }  
                m++;
        }
        
    }
    for (i = 0; i < 4936; i++) {
        for (j = 0; j < 9872; j++) {
            printf("%d ", Gsys[i][j]);
        }
    }
    a = 0;
    a1 = 0;
    for (j = 0; j < 9872; j++) {
        if (c[a] == j) {
            printf("nna\n");
            for (i = 0; i < 4936; i++) {
                G[i][j] = Gsys[i][4936 + a];
            }
            a++;
        }
        if (c1[a1] == j) {
            printf("nna1\n");
            for (i = 0; i < 4936; i++) {
                G[i][j] = Gsys[i][a1];
            }
            a1++;
        }

    }
    printf("y\n");
    for (i = 0; i < 4936; i++) {
        for (j = 0; j < 9872; j++) {
            if (G[i][j] == 1) printf("%d ", j + 1);
        }
        printf("\n");
    }

    FILE *outfp;
    outfp = fopen("generator1.txt","w");
    for (int i = 0; i < 4936; i++) {
        for (int j = 0; j < 9872; j++) {
            fprintf(outfp,"%d ",G[i][j]);
        }
        fprintf(outfp,"\n");
    }
    fclose(outfp);
    printf("\n\n");
    int temp1;
    for (i = 0; i < 4936; i++) {
        temp1 = 0;
        for (j = 0; j < 9872; j++) {
            temp1 += (G[i][j] * H1[i][j]);
            temp1 = temp1 % 2;   
        }
        printf("%d ", temp1);
        if(temp1 == 1) printf("yes\n");
    }

    for (i = 0; i < H1Cmaskrow; i++) free(H1Cmask[i]);
    free(H1Cmask);
    printf("a\n");
    for (i = 0; i < H2Cmaskrow; i++) free(H2Cmask[i]);
    free(H2Cmask);
    for (i = 0; i < H3Cmaskrow; i++) free(H3Cmask[i]);
    free(H3Cmask);
    printf("a\n");
    for (i = 0; i < Hrow; i++) free(H[i]);
    free(H);
    printf("a\n");
    for (i = 0; i < H1row; i++) free(H1[i]);
    free(H1);
    printf("a\n");
    for (i = 0; i < H1Crow; i++) free(H1C[i]);
    printf("a\n");
    free(H1C);
    printf("a\n");
    for (i = 0; i < H2Crow; i++) free(H2C[i]);
    free(H2C);
    for (i = 0; i < H3Crow; i++) free(H3C[i]);
    free(H3C);
    printf("a\n");
    for (i = 0; i < Hsystrow; i++) free(Hsyst[i]);
    free(Hsyst);
    for (i = 0; i < Hsyst1row; i++) free(Hsyst1[i]);
    free(Hsyst1);
    for (i = 0; i < Gsysrow; i++) free(Gsys[i]);
    free(Gsys);
    for (i = 0; i < Grow; i++) free(G[i]);
    free(G);


    return 0;

}