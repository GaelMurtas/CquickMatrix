# include <stdio.h>
# include <stdlib.h>
#ifndef _BASICS
#define _BASICS
#define _DEBUG 1

struct matrix{
    double * tab;
    unsigned long int nbLines;
    unsigned long int nbCols;
};

struct matrix newMatrix(const unsigned long int nbl,const unsigned long int nbc){
    struct matrix M;
    M.tab = (double *) malloc(sizeof(double)*nbl*nbc);
    M.nbLines = nbl;
    M.nbCols = nbc;
    return M;
}

void deleteMatrix(struct matrix * M){
    free(M->tab);
}

inline double* get(struct matrix * M, const unsigned long int i, const unsigned long int j){
    if((i<M->nbLines) && (j<M->nbCols)){
        return &M->tab[(i)*M->nbCols+j];
    }
    fprintf(stderr,"-out of bound access intant-\n");
    exit(-1);
}

void affiche(struct matrix M){
    for(unsigned long int i = 0; i < M.nbLines; i++){
            for(unsigned long int j = 0; j < M.nbCols; j++){
                printf("%f ",*get(&M,i,j));}
            printf("\n");
        }
}

inline struct matrix add(struct matrix * m1, struct matrix * m2){
    if ((m1->nbLines != m2->nbLines)||(m1->nbCols != m2->nbCols)){
    fprintf(stderr,"-Differents sizes matrix add intent-\n");
    exit(-1);
    }
    struct matrix m = newMatrix(m1->nbLines, m1->nbCols);
    for(unsigned long int i = 0; i<(m.nbLines*m.nbCols); i++){
        m.tab[i] = m1->tab[i]+m2->tab[i];
    }
    return m;
}

inline struct matrix subtract(struct matrix * m1, struct matrix * m2){
    if ((m1->nbLines != m2->nbLines)||(m1->nbCols != m2->nbCols)){
    fprintf(stderr,"-Differents sizes matrix add intent-\n");
    exit(-1);
    }
    struct matrix m = newMatrix(m1->nbLines, m1->nbCols);
    for(unsigned long int i = 0; i<(m.nbLines*m.nbCols); i++){
        m.tab[i] = m1->tab[i]-m2->tab[i];
    }
    return m;
}

inline struct matrix multiply(struct matrix * m1, struct matrix * m2){
    if (m2->nbLines != m1->nbCols){
        fprintf(stderr,"-Uncompatible matrix multiply intent-\n");
        exit(-1);
        }
    struct matrix M = newMatrix(m1->nbLines, m2->nbCols);
    for(unsigned long int i = 0; i < m1->nbLines; i++){
            for(unsigned long int j = 0; j < m2->nbCols; j++){
                double e = 0;
                for(unsigned long int k = 0; k < m2->nbLines; k++){
                    e += *get(m1,i,k)*(*get(m2,k,j));
                }
                *get(&M,i,j) = e;
            }
    }
    return M;
}

#endif
