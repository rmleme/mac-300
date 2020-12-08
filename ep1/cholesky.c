/******************************************************************************
 *** Regis de Abreu Barbosa                 Numero USP: 3135701             ***
 *** Rodrigo Mendes Leme                    Numero USP: 3151151             ***
 *** Exercicio-Programa 1                                                   ***
 *** Descricao: recebe um sistema linear do tipo Gx = b e, usando armazena- ***
 ***            mento envelope, encontra sua solucao, atravez de fatoracao  ***
 ***            Cholesky, forward substitution e back substitution.         ***
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ERRO -1
#define MAX 1010

double DIAG[MAX],
       ENV[MAX * (MAX - 1) / 2];
int IENV[MAX + 1];

// Funcao: le_matriz
// Entrada: o nome de um arquivo.
// Saida: ERRO se nao conseguiu ler o arquivo; 1 caso contrario. Alem disso,
//        devolve uma matriz e sua dimensao.
// Descricao: le de um arquivo os elemento de uma matriz.

int le_matriz(char *nome_arq, double G[][MAX], int *n)
{
  FILE *arquivo;
  int i,
      j,
      cont;
  double valor;

  if ((arquivo = fopen(nome_arq,"r")) == NULL)
    return ERRO;
  else
  {
    fscanf(arquivo,"%d",n);
    for (i = 1; i <= *n; i++)
      for (j = 1; j <= *n; j++)
        G[i][j] = 0;
    for (cont = 1; cont <= *n * *n; cont++)     // Le a matriz G
    {
      fscanf(arquivo,"%d %d %lf",&i,&j,&valor);
      G[i][j] = valor;
    }
    fclose(arquivo);
    return 1;
  }
}

// Funcao: le_vetor
// Entrada: a dimensao de um vetor.
// Saida: ERRO se nao conseguiu ler o arquivo; 1 caso contrario. Alem disso,
//        devolve um vetor.
// Descricao: le de um arquivo os elementos de um vetor.

int le_vetor(double b[MAX], int *n)
{
  FILE *arquivo;
  char nome_arq[50];
  int i,
      cont;
  double valor;

  printf ("Digite o nome do arquivo que contem b: ");
  scanf ("%s",nome_arq);
  if ((arquivo = fopen(nome_arq,"r")) == NULL)
    return ERRO;
  else
  {
    fscanf(arquivo,"%d",n);
    for (cont = 1; cont <= *n; cont++)          // Le o vetor b
    {
      fscanf(arquivo,"%d %lf",&i,&valor);
      b[i] = valor;
    }
    fclose(arquivo);
    return 1;
  }
}

// Funcao: monta_envelope
// Entrada: uma matriz e sua dimensao.
// Saida: os tres vetores que compoem o envelope da matriz dada.
// Descricao: monta um esquema de armazenamento do tipo envelope para uma dada
//            matriz.

void monta_envelope(double G[][MAX], int n)
{
  int qtos_env,
      i,
      j,
      k;

  for (i = 1; i <= n; i++)       // Guarda a diagonal principal
    DIAG[i] = G[i][i];

  IENV[1]  = 1;
  qtos_env = 0;
  for (i = 2; i <= n; i++)       // "Percorre" a matriz linha a linha
  {
    j = 1;
    while ((G[i][j] == 0) && (j <= i - 1))      // Procura o primeiro elemento
      j++;                                      // nao nulo da linha i
    if (j <= i - 1)         // A linha i contem elementos nao todos nulos
    {
      IENV[i] = qtos_env + 1;
      for (k = j; k <= i - 1; k++)     // Guarda em ENV os elementos
      {                                // nao nulos da linha i
        qtos_env++;
        ENV[qtos_env] = G[i][k];
      }
    }
    else                 // A linha i so contem elementos nulos
      IENV[i] = -1;      // "Marca" a linha i para posterior manipulacao
  }
  IENV[n + 1] = qtos_env + 1;
  for (i = n; i >= 2; i--)      // Faz as linhas nulas da matriz "apontarem"
    if (IENV[i] == -1)          // para a linha seguinte
      IENV[i] = IENV[i + 1];
}

// Funcao: A
// Entrada: a linha e a coluna de um elemento de uma matriz.
// Saida: um elemento do envelope da matriz.
// Descricao: recebe as coordenadas de um elemento da matriz e as converte para
//            coordenadas do envelope.

double A(int i, int j)
{
  if (i == j)          // E um elemento da diagonal principal
    return DIAG[i];
  else                 // A linha i (ou a parte inicial dela) e nula
    if ((IENV[i] == IENV[i + 1]) || (j < i - IENV[i + 1] + IENV[i]))
      return 0;
    else               // E um elemento do envelope
      return ENV[IENV[i + 1] - i + j];
}

// Funcao: forward_subst
// Entrada: um vetor b, o numero da linha e da coluna iniciais e finais e o ta-
//          manho da semibanda.
// Saida: ERRO se a matriz e singular; 1 caso contrario. Alem disso, devolve o
//        numero de flops.
// Descricao: resolve sistemas triangulares inferiores pelo metodo de forward
//            substitution.

int forward_subst(double b[MAX], int inic, int fim, int s, long int *flops)
{
  int i,
      j;

  for (i = inic; i <= fim; i++)
  {
    if (A(i,i) == 0)       // A matriz e singular
      return ERRO;
    for (j = (i - s > 0 ? i - s : inic); j <= i - 1; j++)
    {
      b[i] -= A(i,j) * b[j];
      (*flops)++;
    }
    b[i] /= A(i,i);
  }
  return 1;
}

// Funcao: back_subst
// Entrada: um vetor b, sua dimensao e o tamanho da semibanda.
// Saida: ERRO se a matriz e singular; 1 caso contrario. Alem disso, devolve o
//        numero de flops.
// Descricao: resolve sistemas triangulares superiores pelo metodo de back
//            substitution.

int back_subst(double b[MAX], int n, int s, long int *flops)
{
  int i,
      j;

  for (j = n; j >= 1; j--)
  {
    if (A(j,j) == 0)       // A matriz e singular
      return ERRO;
    b[j] /= A(j,j);
    for (i = (j - s > 0 ? j - s : 1); i <= j - 1; i++)
    {
      b[i] -= A(j,i) * b[j];
      (*flops)++;
    }
  }
  return 1;
}

// Funcao: cholesky
// Entrada: o envelope de uma matriz e sua dimensao.
// Saida: ERRO se a matriz nao e positiva definida; 0 caso contrario. Alem dis-
//        so, devolve o numero de flops.
// Descricao: calcula o fator de Cholesky de uma matriz positiva definada, pela
//            forma de borda.

int cholesky(int n, int s, long int *flops)
{
  int j,
      k,
      semibanda,
      zeros_inic;
  double *h,
         prod_int_h;
  long int flops_forward;

  for (j = 1; j <= n; j++)
  {
    h = (double *) malloc((j + 1) * sizeof(double));
    for (k = 1; k <= j - 1; k++)         // Copia da matriz A o vetor coluna c
      h[k] = A(j,k);

    flops_forward = 0;
    if (j - s <= 1)
      semibanda = 1;
    else
      semibanda = j - s;
    forward_subst(h,semibanda,j - 1,s,&flops_forward);
    *flops += flops_forward;

    prod_int_h = 0;
    for (k = 1; k <= j - 1; k++)
      prod_int_h += h[k] * h[k];
    if (A(j,j) - prod_int_h <= 0)        // A matriz nao eh positiva definida
      return ERRO;
    DIAG[j] = sqrt(A(j,j) - prod_int_h);

    for (zeros_inic = 0; h[zeros_inic + 1] == 0; zeros_inic++);
    for (k = 1 + zeros_inic; k <= j - 1; k++)     // Joga na matriz A o vetor h
      ENV[IENV[j] + k - 1 - zeros_inic] = h[k];   // Guarda na submatriz linha
    free(h);
  }
  return 1;
}

// Funcao: calcula_banda
// Entrada: o nome de um arquivo.
// Saida: o tamanho da banda da matriz.
// Descricao: calcula o tamanho da banda de uma dada matriz.

void calcula_banda(char *nome_arq, int *s)
{
  FILE *arquivo;
  int i,
      j,
      x;
  double valor;

  arquivo = fopen(nome_arq,"r");
  fscanf(arquivo,"%d",&i);
  *s = 0;    
  while (!feof(arquivo))
  {
    fscanf(arquivo,"%d %d %lf",&i,&j,&valor);
    if ((i > j) && (valor != 0))
    {
      x = i - j;
      if (x > *s)
        *s = x;
    } 
    else
      x = j - i;
  }
  fclose(arquivo);
}

// Funcao: main
// Descricao: recebe e resolve um sistema linear do tipo Gx = b lido de um ar-
//            quivo, chamando as rotinas apropriadas.

void main(int argc, char *argv[])
{
  double G[MAX][MAX],
         b[MAX];
  int n,          // Dimensao da matriz
      s,          // Tamanho da semibanda da matriz
      k;
  long int flops_cholesky = 0,     // Numero de flops da funcao cholesky()
           flops_forward  = 0,     // Numero de flops da funcao forward_subst()
           flops_back     = 0;     // Numero de flops da funcao back_subst()

  if (le_matriz(argv[1],G,&n) == ERRO)
    printf("Erro: nao conseguiu ler a matriz.\n");
  else
    if (le_vetor(b,&n) == ERRO)
      printf("Erro: nao conseguiu ler o vetor b.\n");
    else
    {
      calcula_banda(argv[1],&s);
      monta_envelope(G,n);
      if (cholesky(n,s,&flops_cholesky) == ERRO)
        printf("Erro: a matriz nao e positiva definida.\n");
      else
        if (forward_subst(b,1,n,s,&flops_forward) == ERRO)
          printf("Erro: a matriz nao e singular.\n");
        else
          if (back_subst(b,n,s,&flops_back) == ERRO)
            printf("Erro: a matriz nao e singular.\n");
          else
          {
            printf("Solucao do sistema:\n");
            for (k = 1; k <= n; k++)
              printf("x[%d] = %g\n",k,b[k]);
            printf("Numero de flops da fatoracao Cholesky: %ld.\n",
                   flops_cholesky);
            printf("Numero de flops do forward substitution: %ld.\n",
                   flops_forward);
            printf("Numero de flops do back substitution: %ld.\n",flops_back);
          }
    }
}
