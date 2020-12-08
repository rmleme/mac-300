/******************************************************************************
 *** Regis de Abreu Barbosa                 Numero USP: 3135701             ***
 *** Rodrigo Mendes Leme                    Numero USP: 3151151             ***
 *** Exercicio-Programa 2                                                   ***
 *** Descricao: resolve o problema de minimos quadrados (casos posto comple-***
 ***            to e incompleto), via decomposicao QR por reflexoes.        ***
 *****************************************************************************/

#include <stdio.h>
#include <math.h>

#define MAX 100
#define EPSILON 1e-10        // Aproximacao numerica de zero

// Funcao: le_matriz
// Saida: 1 se conseguiu ler a matriz do arquivo; 0 caso contrario. Alem disso,
//        retorna a matriz lida e a dimensao da mesma.
// Descricao: le um arquivo contendo numeros e, a partir dos mesmos, monta uma
//            matriz.

int le_matriz(double A[][MAX], int *n, int *m)
{
  FILE *arquivoA;
  int i, 
      j;
  char nomeA[30];
  double temp;

  printf ("Digite o nome do arquivo que contem a matriz: ");
  scanf ("%s", nomeA);
  if (!(arquivoA = fopen (nomeA, "r")))
  {
    printf ("Arquivo nao encontrado.\n");
    return 0;
  } 
  fscanf (arquivoA, "%d", n);
  fscanf (arquivoA, "%d", m);
  while (fscanf (arquivoA, "%d %d %lf", &i, &j, &temp) == 3)
    if (i <= *n && j <= *m)       // O elemento lido encontra-se na matriz
      A[i][j] = temp;
    else
    {
      printf ("Erro no arquivo de entrada.\n");
      return 0;
    }
  fclose(arquivoA);
  return 1;
}

// Funcao: le_vetor
// Entrada: a dimensao que o vetor deve possuir.
// Saida: 1 se conseguiu ler o vetor do arquivo; 0 caso contrario. Alem disso,
//        retorna o vetor lido.
// Descricao: le um arquivo contendo numeros e, a partir dos mesmos, monta um
//            vetor.

int le_vetor(double b[MAX], int n)
{
  FILE *arquivob;
  int i, 
      tam_b;
  double temp;
  char nomeb[30];

  printf ("Digite o nome do arquivo que contem o vetor: ");
  scanf ("%s", nomeb);
  if (!(arquivob = fopen (nomeb, "r")))
  {
    printf ("Arquivo nao encontrado.\n");
    return 0;
  } 
  fscanf (arquivob, "%d", &tam_b);
  if (tam_b != n)
  {
    printf ("Erro no tamanho do vetor.\n");
    return 0;
  }
  while (fscanf (arquivob, "%d %lf", &i, &temp) == 2)
    if (i <= n)
      b[i] = temp;
  fclose (arquivob);
  return 1;
}

// Funcao: obtem_Q
// Entrada: uma matriz A, uma linha e uma coluna da mesma e sua dimensao.
// Saida: os vetores sigma e gama.
// Descricao: obtem os vetores sigma e gama usados no calculo do refletor Q.

void obtem_Q(double A[][MAX], double sigma[MAX], double gama[MAX], int lin, 
	     int col, int n)
{
  double m;
  int i;

  m = 0;
  for (i = lin; i <= n; i++)      // Calcula a norma infinita da coluna k de A
    if (fabs(A[i][col]) > m)
      m = A[i][col];

  if (m == 0 || fabs(m) < EPSILON)     // A coluna k inteira vale zero
    gama[col] = 0;
  else
  {
    for (i = lin; i <= n; i++)
      A[i][col] /= m;
    sigma[col] = 0;

    for (i = lin; i <= n; i++)           // Calcula a norma 2 da coluna k de A
      sigma[col] += A[i][col] * A[i][col];
    sigma[col] = -sqrt(sigma[col]);
    if (A[lin][col] < 0)
      sigma[col] = -sigma[col];

    A[lin][col] += sigma[col];
    gama[col]    = 1 / (sigma[col] * A[lin][col]);
    sigma[col]  *= m;
  }
}

// Funcao: back_subst
// Entrada: uma matriz A, um vetor b, o vetor sigma e a dimensao de todos.
// Descricao: resolve sistemas triangulares superiores pelo metodo de back
//            substitution.

void back_subst(double A[][MAX], double b[MAX], double sigma[MAX], int n)
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

// Funcao: troca
// Entrada: uma matriz A, os vetores sigma e gama, as colunas a serem permuta-
//          das e a dimensao de todos.
// Descricao: troca as colunas da matriz A e os valores de sigma e gama das
//            colunas coli e colj.

void troca (double A[][MAX], double sigma[MAX], double gama[MAX], int coli,
            int colj, int n)
{
  int i;
  double temp;
  
  for (i = 1; i <= n; i++)
  {
    temp       = A[i][coli];
    A[i][coli] = A[i][colj];
    A[i][colj] = temp;
  }

  temp        = sigma[coli];           
  sigma[coli] = sigma[colj];
  sigma[colj] = temp;

  temp       = gama[coli];
  gama[coli] = gama[colj];
  gama[colj] = temp;
}

// Funcao: quad_min
// Entrada: uma matriz A, um vetor b e a dimensao de ambos.
// Descricao: resolve sistemas lineares do tipo Ax = b, atraves de decomposicao
//            QR com reflexoes (casos posto completo e incompleto).

void quad_min(double A[][MAX], double b[MAX], int n, int m)
{
  double sigma[MAX],        // Usado na obtencao de Q
         gama[MAX],         // Usado na obtencao de Q
         pi,                // Usado na multiplicacao de Q por um vetor
         temp,
         norma,
         maior_norma;
  int i,
      j,
      k, 
      col,
      rr = m,               // Iteracao a partir da qual a matriz A e singular
      c_troca[MAX];         // Guarda as permutacoes feitas pela decomposicao

  for (i = 1; i <= MAX; i++)    // No inicio, nenhuma coluna de A foi permutada
    c_troca[i] = i;

  // Calculando a decomposicao QR
  for (k = 1; k <= m; k++)
  {
    obtem_Q(A,sigma,gama,k,k,n);

    if (gama[k] == 0)         // Caso posto incompleto; deve permutar colunas
    {
      maior_norma = 0;
      for (j = k + 1; j <= m; j++)       // Calcula coluna de maior norma
      {
        norma = 0;
        col   = j;
        for (i = k; i <= n; i++)
          norma += fabs(A[i][j]);
        if (norma > maior_norma)
          maior_norma = norma;
      }
      if (maior_norma != 0 && fabs(maior_norma) > EPSILON)   // Ainda nao che-
      {                                                      // gou numa matriz
        obtem_Q(A,sigma,gama,k,col,n);                       // singular
        troca(A,sigma,gama,k,col,n);	
        c_troca[k] = col;
        continue;
      }
      else            // Termina o calculo de R[1][1]
      {
        rr = k - 1;
        break;
      }
    }
    else           // Caso posto completo
      rr = m;

    for (j = k + 1; j <= m; j++)      // Multiplica Q pelas outras
    {                                 // colunas da matriz
      pi = 0;
      for (i = k; i <= n; i++)
        pi += A[i][k] * A[i][j];
      pi *= gama[k];
      for (i = k; i <= n; i++)
        A[i][j] -= pi * A[i][k];
    }
  }

  // Calculando c = Q_t * b; guarda o resultado no proprio vetor b
  for (i = 1; i <= m; i++)
  {
    pi = 0;
    for (k = i; k <= n; k++)
      pi += A[k][i] * b[k];
    pi *= gama[i];
    for (k = i; k <= n; k++)
      b[k] -= pi * A[k][i];
  }

  // Resolvendo Rx = c (Ax = b)
  back_subst(A,b,sigma,rr);

  // Atribui zero ao complemento do vetor x
  for (i = rr + 1; i <= m; i++)
    b[i] = 0.0;

  // Transforma x^ em x, desfazendo o pivoteamento do caso posto incompleto
  for (i = rr; i >= 1; i--)
  {
    if (c_troca[i] != i)     // A coluna i foi permutada; deve desfazer
    {
      temp          = b[i];
      b[i]          = b[c_troca[i]];
      b[c_troca[i]] = temp;
    }
  }
}

// Funcao: main
// Descricao: carrega arquivos contendo um sistema linear e, atraves de decom-
//            posicao QR, acha a solucao do problema de quadrados minimos.

int main (void)
{
  int i, 
      n,
      m;
  double A[MAX][MAX],
         b[MAX];

  printf ("Resolucao do problema dos quadrados minimos.\n");
  if (le_matriz(A,&n,&m))
    if (le_vetor(b,n))
    {
      quad_min(A,b,n,m);
      for (i = 1; i <= m; i++)
        printf("x[%d] = %g\n",i,b[i]);
      return 0;
    }
  return -1;
}
