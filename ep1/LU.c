/******************************************************************************
 *** Regis de Abreu Barbosa                 Numero USP: 3135701             ***
 *** Rodrigo Mendes Leme                    Numero USP: 3151151             ***
 *** Exercicio-Programa 1                                                   ***
 *** Descricao: Recebe uma matriz do tipo banda, faz a decomposicao da      ***
 *** matriz em LU, L triangular inferior e U triangular superior, e resolve ***
 *** o sistema Ax = b da seguinte forma: primeiro resolve o sistema Ly = b  ***
 *** e depois Ux = y.                                                       ***
 *****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <time.h>

static int n,    // Numero de linhas da matriz de entrada
           tam,  // Numero de colunas da matriz orientada a linha
           n_A;  // Numero de elementos da matriz A

// Funcao: val
// Entrada: a matriz, a linha i e coluna j.
// Saida: o valor (i,j) da matriz

double val (double matriz[], int i, int j) {
  int t = j + ((i-1)*tam);
  if (t > 0 && t < n_A) 
    return matriz[t];
  else return 0;
}

// Funcao: entrada
// Entrada: A matriz, a linha i, a coluna j, e o valor de (i,j).
// Descricao: Atribui a posicao (i,j) da matriz o valor da variavel valor

void entrada (double matriz[], int i, int j, double valor) {
  int t = j + ((i-1)*tam);
  if (t > 0 && t < n_A)
    matriz[t] = valor;
}

// Funcao: banda
// Entrada: o arquivo que contem a matriz.
// Saida: as bandas inferior s e superior t.

void banda (FILE *arq, int *s, int *t)
{
  int i, j, x;
  double valor;

  *s = *t = 0;
  fscanf (arq, "%d", &i);
  n = i;

  while (!feof(arq)) {
    fscanf (arq, "%d %d %lf", &i, &j, &valor);
    if (valor != 0) {
      if (i > j) {
        x = i - j;
        if (x > *s) *s = x;
      } else {
        x = j - i;
        if (x > *t) *t = x;
      }
    }
  }
  printf ("Matriz A:\nDimensao: %d\n", n);
  printf ("Banda inferior: %d\nBanda superior: %d\n\n", *s, *t);
}

// Funcao: preenche_matriz
// Entrada: O arquivo que contem a matriz, o vetor A que guarda a matriz,
//          as bandas s e t.
// Descricao: Preenche a matriz A pela forma esquematizada para matrizes banda

void preenche_matriz (FILE *arq, double *A, int s, int t)
{
  int i, j;
  double valor;

  fscanf (arq, "%d", &i);
  if (n != i) {
    printf ("Erro de leitura!\n");
    return ;
  }
  s++;
  while (!feof(arq)) {
    if (fscanf (arq, "%d %d %lf", &i, &j, &valor) == 3) {
      if (valor != 0) entrada (A, i, j-i+s, valor);
    }
  }
}

// Funcao: preenche_b
// Saida: O vetor b na variavel de entrada *b.
//        Retorna 0 se o arquivo for invalido, 1 caso contrario.


int preenche_b (double *b)
{
  int i = 0;
  double temp;
  char nome[20];
  FILE *arquivo;

  printf ("Digite o nome do arquivo que contem b: ");
  scanf ("%s", nome);
  if (!(arquivo = fopen (nome, "r"))) {
      printf ("Arquivo nao encontrado.\n");
      return 0;
  } else {
    fscanf (arquivo, "%d", &i);
    if (i != n) {
      printf ("Erro no tamanho do vetor b.\n");
      return 0;
    }
    while (fscanf (arquivo, "%d %lf", &i, &temp) > 1)
      b[i] = temp;
  }
  return 1;
}


// Funcao: decompLU
// Entrada: o arquivo que contem a matriz e as bandas s e t.
// Saida: o tempo da decomposicao LU
// Descricao: Decompoe a matriz A contida no arquivo de entrada em A = LU,
//            L triangular inferior e U triangular inferior.
//            Apos a decomposicao, resolve o sistema Ax = b, enquanto o
//            usuario fornecer o b na entrada.

int decompLU (FILE *arq, int s, int t)
{
  double A[(n+1)*(tam+1)+1], b[n+2];
  int i, j, k, temp_int, tempo_LU, cont;
  double temp, temp1;
  char c[5];

  n_A = (n+1)*(tam+1);
  for (i = 1; i <= n_A; i++)
    A[i] = 0.0;

  preenche_matriz (arq, A, s, t);

  cont = 0;
  for (k = 1; k < n; k++) {        // Decomposicao LU sem pivoteamento.
    for (i = k+1; i <= n; i++) {
      j = k-i+s+1;
      if (j > 0) {
        if (val(A, k, s+1) == 0) { // A matriz e singular
          printf ("A matriz e singular.\n");
          return 0;
        }
        temp = val(A, i, j)/val(A, k, s+1);
        entrada (A, i, j, temp);
        for (j = s+2; j <= tam; j++) {
          temp1 = val (A, i, j+k-i) - (temp*val(A, k, j));
          entrada (A, i, j+k-i, temp1);
	  cont++;
        }
      }
    }
  }

  printf ("Decomposicao LU: %d flops\n\n", cont );

  while (temp_int != 0) {

    if (!preenche_b(b)) return tempo_LU;
    
    cont = 0;
    for (i = 2; i <= n; i++) {     // foward-substitution
      for (j = 1; j < i; j++) {
        if ((temp_int = j-i+s+1) > 0) {
          b[i] -= (val(A, i, temp_int)*b[j]);
	  cont++;
        }
      }
    }
    printf ("Foward substitution: %d flops.\n", cont);

    cont = 0;
    for (i = n; i > 0; i--) {      // back-substitution
      for (j = i+1; j <= n; j++) { 
        if ((temp_int = j-i+s+1) <= tam) {
          b[i] -= (val(A, i, temp_int)*b[j]);
	  cont++;
        } else break;
      }
      b[i] /= val (A, i, s+1);
    }
    printf ("Back substitution: %d flops.\n\n", cont);

    printf ("Solucao do sistema:\n");
    for (i = 1; i <= n; i++) {
      printf ("x[%d] = %f\n", i, b[i]);
    }

    printf ("\nDeseja continuar (s/n)? ");
    scanf ("%s[sn]", c);
    if (c[0] == 'n') temp_int = 0;
    else {
      temp_int = 1;
      printf ("\n"); 
    }
  }  
  return tempo_LU;
}

// Funcao: main

int main (int argc, char *argv[])
{
  FILE *arquivo;
  int t, s, tempo;

  if ((arquivo = fopen (argv[1], "r")) == 0) {
    printf ("ERRO ao abrir o arquivo.\n");
    printf ("Digite o nome do arquivo na linha de entrada.\n");
    return 0;
  }

  banda (arquivo, &s, &t);
  tam = s + t + 1;
  fclose (arquivo);
  arquivo = fopen (argv[1], "r");
  tempo = decompLU (arquivo, s, t);
  return 1;
}
