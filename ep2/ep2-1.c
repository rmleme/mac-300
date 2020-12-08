/******************************************************************************
 *** Regis de Abreu Barbosa                 Numero USP: 3135701             ***
 *** Rodrigo Mendes Leme                    Numero USP: 3151151             ***
 *** Exercicio-Programa 2                                                   ***
 *** Descricao: resolve sistemas lineares quadrados por decomposicao QR,    ***
 ***            utilizando rotacoes, rotacoes rapidas ou reflexoes.         ***
 *****************************************************************************/

#include <stdio.h>
#include <math.h>

#define ERRO -1
#define MAX 100

// Funcao: back_subst
// Entrada: uma matriz A, um vetor b e a dimensao de ambos.
// Saida: ERRO se a matriz e singular; 1 caso contrario.
// Descricao: resolve sistemas triangulares superiores pelo metodo de back
//            substitution.

int back_subst(double A[][MAX], double b[MAX], int n)
{
  int i,
      j;

  for (j = n; j >= 1; j--)
  {
    if (A[j][j] == 0)       // A matriz e singular
      return ERRO;
    b[j] /= A[j][j];
    for (i = 1; i <= j - 1; i++)
      b[i] -= A[i][j] * b[j];
  }
  return 1;
}

// ROTACAO

// Funcao: rotac
// Entrada: A matriz A, o vetor b e a dimensao n da matriz.
// Saida: o vetor x (em b) tal que Ax = b.
// Descricao: Resolve o sistema linear Ax = b por decomposicao QR
//            utilizando o algoritmo de rotacoes.

int rotac(double A[][MAX], double b[], int n)
{
  double sen_theta,             // Seno do angulo de rotacao em Q_t
         cos_theta,             // Cosseno do angulo de rotacao em Q_t
         aux;
  int i,
      j,
      col,
      qtos_prod  = 0,
      qtas_somas = 0;

  // Faz o produto de rotacoes
  for (j = 1; j <= n - 1; j++)
  {
    for (i = j + 1; i <= n; i++)
    {
      cos_theta = A[j][j] / sqrt(A[j][j] * A[j][j] + A[i][j] * A[i][j]);
      sen_theta = A[i][j] / sqrt(A[j][j] * A[j][j] + A[i][j] * A[i][j]);

      // Calculando R = Q_t * A; guarda o resultado na propria matriz A
      for (col = 1; col <= n; col++)
      {
        aux       = A[j][col];
        A[j][col] = cos_theta * aux + sen_theta * A[i][col];
        A[i][col] = -sen_theta * aux + cos_theta * A[i][col];
      }
      qtos_prod  += 4 * n;
      qtas_somas += 2 * n;

      // Calculando c = Q_t * b; guarda o resultado no proprio vetor b
      aux  = b[j];
      b[j] = cos_theta * aux + sen_theta * b[i];
      b[i] = -sen_theta * aux + cos_theta * b[i];
    }
  }

  // Resolvendo Rx = c (Ax = b)
  back_subst(A,b,n);

  for (i = 1; i <= n; i++)
    printf("x[%d] = %g\n",i,b[i]);
  printf("Numero de somas: %d\n",qtas_somas);
  printf("Numero de produtos: %d\n",qtos_prod);
  return 0;
}

// ROTACAO RAPIDA

// Funcao: interchange
// Entrada:  A matriz A, o vetor b, as linhas i e j e  a dimensao n da matriz.
// Descricao: Troca as linhas i e j da matriz A e do vetor b.

void interchange (double A[][MAX], double b[MAX], int i, int j, int n)
{
  int k;
  double temp;

  for (k = 1; k <= n; k++) {
    temp = A[i][k];
    A[i][k] = A[j][k];
    A[j][k] = temp;
  }
  temp = b[i]; b[i] = b[j]; b[j] = temp;
}

// Funcao: updB_tg
// Entrada: A matriz A, o vetor b, os indices i e j das linhas, a dimensao n
//          do vetor e os valores de p e q do algoritmo.
// Descricao: Atualiza os valores de A e b para o algoritmo de rotacoes rapidas
//            pela tangente do angulo theta.

void updB_tg (double A[][MAX], double b[MAX], int i, int j, int n, double p,
              double q)
{
  int k;
  double temp;

  for (k = 1; k <= n; k++) {
    temp = A[i][k];
    A[i][k] += (p*A[j][k]);
    A[j][k] -= (q*temp);
  }
  temp = b[i];
  b[i] += (p*b[j]);
  b[j] -= (q*temp);
}

// Funcao: updB_ctg
// Entrada: A matriz A, o vetor b, os indices i e j das linhas, a dimensao n
//          do vetor e os valores de p e q do algoritmo.
// Descricao: Atualiza os valores de A e b para o algoritmo de rotacoes rapidas
//            pela cotangente do angulo theta.

void updB_ctg (double A[][MAX], double b[MAX], int i, int j, int n, double p,
               double q)
{
  int k;
  double temp;

  for (k = 1; k <= n; k++) {
    temp = A[i][k];
    A[i][k] = (p*A[i][k]) + A[j][k];
    A[j][k] = (q*A[j][k]) - temp;
  }
  temp = b[i];
  b[i] = (p*b[i]) + b[j];
  b[j] = (q*b[j]) - temp;
}

// Funcao: fastrotator
// Entrada: A matriz A, o vetor b e a dimensao n da matriz.
// Saida: O vetor x (em b) tal que Ax = b.
// Descricao: Resolve o sistema linear Ax = b por decomposicao QR
//            utilizando o algoritmo de rotacoes rapidas.

void fastrotator (double A[][MAX], double b[], int n) 
{
  int i, j, m, flops = 0;
  double p, q, r, temp;
  double g[n+1];

  for (i = 0; i <= n; i++)
    g[i] = 1;

  for (i = 1; i < n; i++) {
    for (j = i+1; j <= n; j++) {
      m = i;
      if (A[j][m] != 0) {
	if (A[i][m] == 0) {
	  interchange (A, b, i, j, n);
	  continue;
	}
	q = A[j][m]/A[i][m];
	p = q*g[j]/g[i];
	r = p*q;
	if (r <= 1) {
	  r++; r = 1/r;
	  g[i] = r*g[i];
	  g[j] = r*g[j];
	  updB_tg (A, b, i, j, n, p, q);
          flops += 2 * n;
	} else {
	  p = 1/p;
	  q = 1/q;
	  r = p*q;
	  r++; r = 1/r;
	  temp = g[i];
	  g[i] = r*g[j];
	  g[j] = r*temp;
	  updB_ctg (A, b, i, j, n, p, q);
	  flops += 2 * n;
	}
      }
    }
  }
  back_subst (A, b, n);
  for (i = 1; i <= n; i++)
    printf ("x[%d] = %g\n", i, b[i]);
  printf("Numero de flops: %d\n", flops);
}

// REFLEXAO

// Funcao: obtem_Q
// Entrada: uma matriz A, uma coluna da mesma e sua dimensao.
// Saida: os vetores sigma e gama.
// Descricao: obtem os vetores sigma e gama usados no calculo do refletor Q.

void obtem_Q(double A[][MAX], double sigma[MAX], double gama[MAX], int k,
             int n)
{
  double m;
  int i;

  m = 0;
  for (i = k; i <= n; i++)      // Calcula a norma infinita da coluna k de A
    if (fabs(A[i][k]) > m)
      m = A[i][k];

  if (m == 0)      // A coluna k inteira vale 0
    gama[k] = 0;
  else
  {
    for (i = k; i <= n; i++)
      A[i][k] /= m;
    sigma[k] = 0;

    for (i = k; i <= n; i++)             // Calcula a norma 2 da coluna k de A
      sigma[k] += A[i][k] * A[i][k];
    sigma[k] = -sqrt(sigma[k]);
    if (A[k][k] < 0)
      sigma[k] = -sigma[k];

    A[k][k]  += sigma[k];
    gama[k  ] = 1 / (sigma[k] * A[k][k]);
    sigma[k] *= m;
  }
}

// Funcao: back_subst_reflex
// Entrada: uma matriz A, um vetor b, o vetor sigma e a dimensao de todos.
// Descricao: resolve sistemas triangulares superiores pelo metodo de back
//            substitution.

void back_subst_reflex(double A[][MAX], double b[MAX], double sigma[MAX],
                       int n)
{
  int i,
      j;

 for (j = n; j >= 1; j--)
 {
   b[j] /= -sigma[j];
   for (i = 1; i <= j - 1; i++)
     b[i] -= A[i][j] * b[j];
 }
}

// Funcao: reflex
// Entrada: A matriz A, o vetor b e a dimensao n da matriz.
// Saida: O vetor x (em b) tal que Ax = b.
// Descricao: resolve sistemas lineares do tipo Ax = b, atraves de decomposicao
//            QR com reflexoes.

int reflex(double A[][MAX], double b[], int n)
{
  double sigma[MAX],        // Usado na obtencao de Q
         gama[MAX],         // Usado na obtencao de Q
         pi;                // Usado na multiplicacao de Q por um vetor
  int i,
      j,
      k,
      qtos_flops = 0;

  // Calculando a decomposicao QR
  for (k = 1; k <= n - 1; k++)
  {
    obtem_Q(A,sigma,gama,k,n);
    for (j = k + 1; j <= n; j++)      // Multiplica Q pelas outras
    {                                 // colunas da matriz
      pi = 0;
      for (i = k; i <= n; i++)
        pi += A[i][k] * A[i][j];
      pi *= gama[k];
      for (i = k; i <= n; i++)
      A[i][j] -= pi * A[i][k];
      qtos_flops += 2 * n;
    }
  }
  sigma[n] = -A[n][n];
  gama[n]  = A[n][n];

  // Testa se a matriz e singular
  for (i = 1; i <= n; i++)
  {
    if (gama[i] == 0)
      return ERRO;
  }

  // Calculando c = Q_t * b; guarda o resultado no proprio vetor b
  for (i = 1; i <= n - 1; i++)
  {
    pi = 0;
    for (k = i; k <= n; k++)
      pi += A[k][i] * b[k];
    pi *= gama[i];
    for (k = i; k <= n; k++)
    b[k] -= pi * A[k][i];
  }

  // Resolvendo Rx = c (Ax = b)
  back_subst_reflex (A,b,sigma,n);

  for (i = 1; i <= n; i++)
    printf("x[%d] = %g\n",i,b[i]);
  printf("Numero de flops: %d\n",qtos_flops);
  return 0;
}

// Funcao: main
// Descricao: permite que o usuario escolha um metodo de resolucao de sistemas
//            lineares e chama as funcoes que executam o mesmo.

int main (void)
{
  FILE *arquivoA,
       *arquivob;
  int opcao,
      i, 
      j,
      n,
      tam_n = 0;
  double A[MAX][MAX],
         b[MAX],
         temp;
  char nomeA[20], 
       nomeb[20];

  printf ("Resolucao de sistemas lineares por decomposicao QR utilizando:\n");
  printf ("1 - rotacoes\n2 - rotacoes rapidas\n3 - reflexoes\n");
  printf ("Escolha o metodo de resolucao: ");
  scanf ("%d", &opcao);
  while (opcao < 1 || opcao > 3) {
    printf ("Opcao invalida.\n");
    return ERRO;
  }
  switch (opcao) {
    case 1:  
      printf ("Rotacoes.\n");
      break;
    case 2:  
      printf ("Rotacoes rapidas.\n");
      break;
    case 3: 
      printf ("Reflexoes.\n");
      break;
  }

  printf ("\nDigite o nome do arquivo que contem a matriz A: ");
  scanf ("%s", nomeA);
  if (!(arquivoA = fopen (nomeA, "r"))) {
      printf ("Arquivo nao encontrado.\n");
      return ERRO;
  } 
  fscanf (arquivoA, "%d", &tam_n);
  n = tam_n;
  while (fscanf (arquivoA, "%d %d %lf", &i, &j, &temp) == 3) {
    if (i <= n && j <= n)
      A[i][j] = temp;
    else printf ("Erro\n");
  }

  printf ("\nDigite o nome do arquivo que contem o vetor b: ");
  scanf ("%s", nomeb);
  if (!(arquivob = fopen (nomeb, "r"))) {
      printf ("Arquivo nao encontrado.\n");
      return ERRO;
  } 
  fscanf (arquivob, "%d", &tam_n);
  if (tam_n != n) {
    printf ("Erro no tamanho do vetor.\n");
    return ERRO;
  }
  while (fscanf (arquivob, "%d %lf", &i, &temp) == 2) {
    if (i <= n)
      b[i] = temp;
  }

  switch (opcao) {
    case 1:
      rotac (A, b, n);
      break;
    case 2:
      fastrotator (A, b, n);
      break;
    case 3:
      reflex (A, b, n);
      break;
    }
  return 0;
}
