#include "matrixBasics.c"

#define MAXREC 33

//On recopie du contenu d'une grande matrice dans une plus petite a partir d'une certaine position dans la grande
inline void copieTo(struct matrix * G, struct matrix * P, size_t line, size_t col){
    if((G->nbCols < (P->nbCols+col)) || (G->nbLines < (P->nbLines+line))){//faire un mod debeug avec le prépocesseur qui active ce genre de vérification pour gagner en vitesse si besoin
        fprintf(stderr,"-Unvalid matrix sizes in CopieTo-\n");
        exit(-1);
        }
    for(unsigned long int i=0; i<P->nbLines; i++){
        for(unsigned long int j=0; j<P->nbCols; j++){
            //printf("\n %ld", i);
            //printf(" %ld \n", j);
            *get(P,i,j) = *get(G, i+line, j+col);
        }
    }
}

struct matrix multRec(struct matrix * m1, struct matrix * m2){
    if(m1->nbCols != m2->nbLines){
        fprintf(stderr,"-Uncompatible matrix multiply intent on multRec-\n");
        exit(-1);
        }
    if(m1->nbLines < MAXREC || m1->nbCols< MAXREC || m2->nbCols< MAXREC) return multiply(m1, m2);
    //cas général
    //calcul du découpage
    size_t nbLines1 = m1->nbLines/2;
    size_t nbLinesCols = m2 -> nbLines/2;//m1Cols = m2Lines
    size_t nbCols2 = m2->nbCols/2;
    //initialisation des sous-matrices
    struct matrix A1 = newMatrix(nbLines1, nbLinesCols);
    struct matrix A2 = newMatrix(nbLines1, m1->nbCols-nbLinesCols);
    struct matrix A3 = newMatrix(m1->nbLines-nbLines1, nbLinesCols);
    struct matrix A4 = newMatrix(m1->nbLines-nbLines1, m1->nbCols-nbLinesCols);
    struct matrix B1 = newMatrix(nbLinesCols, nbCols2);
    struct matrix B2 = newMatrix(nbLinesCols, m2->nbCols-nbCols2);
    struct matrix B3 = newMatrix(m2->nbLines-nbLinesCols, nbCols2);
    struct matrix B4 = newMatrix(m2->nbLines-nbLinesCols, m2->nbCols-nbCols2);
    copieTo(m1,&A1,0,0);
    copieTo(m1,&A2,0,nbLinesCols);
    copieTo(m1,&A3,nbLines1,0);
    copieTo(m1,&A4,nbLines1,nbLinesCols);
    copieTo(m2,&B1,0,0);
    copieTo(m2,&B2,0,nbCols2);
    copieTo(m2,&B3,nbLinesCols,0);
    copieTo(m2,&B4,nbLinesCols,nbCols2);
    //cas général
    //calcul des quatre sous matrice
    struct matrix subProd1 = multRec(&A1,&B1);
    struct matrix subProd2 = multRec(&A2,&B3);
    struct matrix m11 = add(&subProd1, &subProd2);
    deleteMatrix(&subProd1);
    deleteMatrix(&subProd2);
    subProd1 = multRec(&A1,&B2);
    subProd2 = multRec(&A2,&B4);
    struct matrix m12 = add(&subProd1, &subProd2);
    deleteMatrix(&subProd1);
    deleteMatrix(&subProd2);
    deleteMatrix(&A1);
    deleteMatrix(&A2);
    subProd1 = multRec(&A3,&B1);
    subProd2 = multRec(&A4,&B3);
    struct matrix m21 = add(&subProd1, &subProd2);
    deleteMatrix(&subProd1);
    deleteMatrix(&subProd2);
    deleteMatrix(&B1);
    deleteMatrix(&B3);
    subProd1 = multRec(&A3,&B2);
    subProd2 = multRec(&A4,&B4);
    struct matrix m22 = add(&subProd1, &subProd2);
    deleteMatrix(&subProd1);
    deleteMatrix(&subProd2);
    deleteMatrix(&A3);
    deleteMatrix(&A4);
    deleteMatrix(&B2);
    deleteMatrix(&B4);
    struct matrix result = newMatrix(m1->nbLines,m2->nbCols);
    for(size_t i = 0; i<nbLines1; i++){
        for(size_t j = 0; j< nbCols2; j++){
            *get(&result, i, j) = *get(&m11, i, j);
            }
        for(size_t j = nbCols2; j<result.nbCols; j++){
            *get(&result, i, j) = *get(&m12, i, j-nbCols2);
            }
    }
    for(size_t i = nbLines1; i<result.nbLines; i++){
        for(size_t j = 0; j < nbCols2; j++){
            *get(&result, i, j) = *get(&m21, i-nbLines1, j);
            }
        for(size_t j = nbCols2; j<result.nbCols; j++){
            *get(&result, i, j) = *get(&m22, i-nbLines1, j-nbCols2);
            }
    }
    return result;
}

//Utilitaires pour l'algo strassen
//simple calcul de max et de min
inline size_t min(size_t i, size_t j){
    if(j>i)return i;
    return j;
}
inline size_t max(size_t i, size_t j){
    if(j<i)return i;
    return j;
}

//calcul de la plus grande puissance de 2 inférieure en O(1);
inline size_t powerOf2(size_t l){
    size_t max = 9223372036854775808;//2^63 pour des meilleurs performance sur la multiplication matricielle on utilisera 2^31
    for(size_t result = 0; max > 0; max = max >> 1){
        result = l & max;
        if(result != 0) return result;
    }
    return 0;
}

//maltriplication par l'agorithme de strassen
struct matrix strassenCore(struct matrix * m1, struct matrix* m2){
    size_t n = m1->nbLines;
    size_t split = m1->nbLines/2;
    if(n< MAXREC) return multiply(m1,m2);
    //découpage des matrices initialles
    struct matrix A1 = newMatrix(split, split);
    struct matrix A2 = newMatrix(split, split);
    struct matrix A3 = newMatrix(split, split);
    struct matrix A4 = newMatrix(split, split);
    struct matrix B1 = newMatrix(split, split);
    struct matrix B2 = newMatrix(split, split);
    struct matrix B3 = newMatrix(split, split);
    struct matrix B4 = newMatrix(split, split);
    copieTo(m1,&A1,0,0);
    copieTo(m1,&A2,0,split);
    copieTo(m1,&A3,split,0);
    copieTo(m1,&A4,split,split);
    copieTo(m2,&B1,0,0);
    copieTo(m2,&B2,0,split);
    copieTo(m2,&B3,split,0);
    copieTo(m2,&B4,split,split);
    //matrices de calcul intermédiaire
    struct matrix S1 = add(&A4,&A2);
    struct matrix S2 = subtract(&A4,&A3);
    struct matrix S3 = add(&S2,&A2);
    struct matrix S4 = subtract(&S3,&A1);
    struct matrix T1 = add(&B4,&B2);
    struct matrix T2 = subtract(&B4,&B3);
    struct matrix T3 = add(&T2,&B2);
    struct matrix T4 = subtract(&T3,&B1);
    //partie récursivité on libère toute la mémoire inutile avant chaque appel récursif
    deleteMatrix(&A4);
    deleteMatrix(&B4);
    struct matrix P1 = strassenCore(&S1,&T1);
    deleteMatrix(&S1);
    deleteMatrix(&T1);
    struct matrix P2 = strassenCore(&S2,&T2);
    deleteMatrix(&S2);
    deleteMatrix(&T2);
    struct matrix P3 = strassenCore(&S3,&T3);
    deleteMatrix(&S3);
    deleteMatrix(&T3);
    struct matrix P4 = strassenCore(&A1,&B1);
    deleteMatrix(&A1);
    deleteMatrix(&B1);
    struct matrix P5 = strassenCore(&A2,&B3);
    deleteMatrix(&A2);
    deleteMatrix(&B3);
    struct matrix P6 = strassenCore(&S4, &B2);
    deleteMatrix(&S4);
    deleteMatrix(&B2);
    struct matrix P7 = strassenCore(&A3,&T4);
    deleteMatrix(&T4);
    deleteMatrix(&A3);
    //calcul de la matrice finale
    struct matrix U1 = add(&P3, &P5);
    deleteMatrix(&P3);
    struct matrix U2 = subtract(&P1, &U1);
    deleteMatrix(&P1);
    struct matrix U3 = subtract(&U1, &P2);
    deleteMatrix(&U1);
    struct matrix C1 = add(&P4, &P5);
    deleteMatrix(&P4);
    deleteMatrix(&P5);
    struct matrix C2 = subtract(&U3, &P6);
    deleteMatrix(&U3);
    deleteMatrix(&P6);
    struct matrix C3 = subtract(&U2, &P7);
    deleteMatrix(&P7);
    struct matrix C4 = add(&P2, &U2);
    deleteMatrix(&P2);
    deleteMatrix(&U2);
    //copie les données dans la matrice de résultat finale
    struct matrix result = newMatrix(n, n);
    for(size_t i = 0; i<split; i++){
        for(size_t j = 0; j<split; j++){
            *get(&result, i, j) = *get(&C1, i, j);
            }
        for(size_t j = split; j<n; j++){
            *get(&result, i, j) = *get(&C2, i, j-split);
            }
    }
    for(size_t i = split; i<n; i++){
        for(size_t j = 0; j < split; j++){
            *get(&result, i, j) = *get(&C3, i-split, j);
            }
        for(size_t j = split; j<n; j++){
            *get(&result, i, j) = *get(&C4, i-split, j-split);
            }
    }
    deleteMatrix(&C1);
    deleteMatrix(&C2);
    deleteMatrix(&C3);
    deleteMatrix(&C4);
    return result;
}




struct matrix multStrassen(struct matrix * m1, struct matrix * m2){
    if(m1->nbCols != m2->nbLines){
        fprintf(stderr,"-Uncompatible matrix multiply intent on multRec-\n");
        exit(-1);
        }
    //cas de base on a plus qu'une seule ligne ou une seule colone
    if(m1->nbLines < MAXREC || m1->nbCols< MAXREC || m2->nbCols< MAXREC) return multiply(m1, m2);
    //cas général
    //calcul du découpage
    size_t split = min(m1->nbLines, m2->nbCols);
    split = min(split, m1->nbCols);
    split = powerOf2(split-1);//Le -1 sert a éviter de devoir géré le cas où une longueur est une puissance de 2.
                              //Dans un ago 100% optimiser il ne faut pas faire comme ça.
    //initialisation des sous-matrices
    struct matrix A1 = newMatrix(split, split);
    struct matrix A2 = newMatrix(split, m1->nbCols-split);
    struct matrix A3 = newMatrix(m1->nbLines-split, split);
    struct matrix A4 = newMatrix(m1->nbLines-split, m1->nbCols-split);
    struct matrix B1 = newMatrix(split, split);
    struct matrix B2 = newMatrix(split, m2->nbCols-split);
    struct matrix B3 = newMatrix(m2->nbLines-split, split);
    struct matrix B4 = newMatrix(m2->nbLines-split, m2->nbCols-split);   //printf("\n Check point 4 \n");
    copieTo(m1,&A1,0,0);
    copieTo(m1,&A2,0,split);
    copieTo(m1,&A3,split ,0);
    copieTo(m1,&A4,split ,split );
    copieTo(m2,&B1,0,0);
    copieTo(m2,&B2,0,split );
    copieTo(m2,&B3,split ,0);
    copieTo(m2,&B4,split ,split );
    //cas général
    //calcul des quatre sous matrice
    struct matrix subProd1 = strassenCore(&A1,&B1);
    struct matrix subProd2 = multRec(&A2,&B3);
    struct matrix m11 = add(&subProd1, &subProd2);
    deleteMatrix(&subProd1);
    deleteMatrix(&subProd2);
    subProd1 = multRec(&A1,&B2);
    subProd2 = multRec(&A2,&B4);
    struct matrix m12 = add(&subProd1, &subProd2);
    deleteMatrix(&subProd1);
    deleteMatrix(&subProd2);
    deleteMatrix(&A1);
    deleteMatrix(&A2);
    subProd1 = multRec(&A3,&B1);
    subProd2 = multRec(&A4,&B3);
    struct matrix m21 = add(&subProd1, &subProd2);
    deleteMatrix(&subProd1);
    deleteMatrix(&subProd2);
    deleteMatrix(&B1);
    deleteMatrix(&B3);
    subProd1 = multRec(&A3,&B2);
    subProd2 = multRec(&A4,&B4);
    struct matrix m22 = add(&subProd1, &subProd2);
    deleteMatrix(&subProd1);
    deleteMatrix(&subProd2);
    deleteMatrix(&A3);
    deleteMatrix(&A4);
    deleteMatrix(&B2);
    deleteMatrix(&B4);
    struct matrix result = newMatrix(m1->nbLines,m2->nbCols);
    for(size_t i = 0; i<split ; i++){
        for(size_t j = 0; j< split ; j++){
            *get(&result, i, j) = *get(&m11, i, j);
            }
        for(size_t j = split ; j<result.nbCols; j++){
            *get(&result, i, j) = *get(&m12, i, j-split);
            }
    }
    for(size_t i = split ; i<result.nbLines; i++){
        for(size_t j = 0; j < split ; j++){
            *get(&result, i, j) = *get(&m21, i-split, j);
            }
        for(size_t j = split ; j<result.nbCols; j++){
            *get(&result, i, j) = *get(&m22, i-split, j-split);
            }
    }
    return result;
}
