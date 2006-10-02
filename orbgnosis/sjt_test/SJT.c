/* $Id: SJT.c,v 1.2 2006/10/02 03:40:45 trs137 Exp $ */
/*===================================================================*/
/* C program for distribution from the Combinatorial Object Server.  */
/* Generate permutations by transposing adjacent elements            */
/* via the Steinhaus-Johnson-Trotter algorithm.  This is             */
/* the same version used in the book "Combinatorial Generation."     */
/* Both the permutation (in one-line notation) and the positions     */
/* being transposed (as a 2-cycle) are output.                       */
/* The program can be modified, translated to other languages, etc., */
/* so long as proper acknowledgement is given (author and source).   */  
/* Programmer: Frank Ruskey, 1995.                                   */
/* The latest version of this program may be found at the site       */
/* http://theory.cs.uvic.ca/inf/perm/PermInfo.html                    */
/*===================================================================*/

int  NN, i, count=0 ;
int  p[100], pi[100];   /* The permutation and its inverse */
int  dir[100];          /* The directions of each element  */

void PrintPerm() {
  int i;
  /* uncomment if you want to print the index of each perm */
  /*
     count = count + 1;
     printf( "[%8d] ", count ); 
  */
  printf("0");
  for (i=1; i <= NN; ++i) printf( "%d", p[i] );
} /* PrintPerm */;

void PrintTrans( int x, int y ) {
  /*printf( "    (%d %d)", x, y );*/
  printf( "\n" );
} /* PrintTrans */;

void Move( int x, int d ) {
  int z;
  PrintTrans( pi[x], pi[x]+d ); 
  z = p[pi[x]+d];
  p[pi[x]] = z;
  p[pi[x]+d] = x;
  pi[z] = pi[x];
  pi[x] = pi[x]+d;  
} /* Move */;

void Perm ( int n ) { 
  int i;
  if (n > NN) PrintPerm();
  else {
    Perm( n+1 );
    for (i=1; i<=n-1; ++i) {
       Move( n, dir[n] );  Perm( n+1 ); 
    }
    dir[n] = -dir[n];
  }
} /* of Perm */;

int main (int argc, char **argv) {
  if (argc < 2)
  {
    printf ("\nUsage: enter integer 1-9.\n");
    exit(1);
  }
  NN = (int)atoi(argv[1]);
  for (i=1; i<=NN; ++i) {
    dir[i] = -1;  p[i] = i;  pi[i] = i;
  }
  Perm ( 1 );
  printf( "\n" );
}
