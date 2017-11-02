/* countFstPi.c - Takes a .count file (of # of reads supporting ref vs. alt 
allele for each individual, without chr or pos information) and outputs pi 
and Fst.  Calculations choose a random read from each individual at each 
site.  A seed file called "seedms" must be present for the program to work.  
Usage:

countFstPi k pop1 pop2

k = # of individuals in input file  (each line should have 2k fields)

Output:

piW1 piW2 piB Fst

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define SMAX 15000000

int main (int argc, char *argv[])
{
  int a, b, k, n1, n2, S=0, *poplist1, *poplist2, temp1, allele, sites=0;
  double piW1=0., piW2=0., piB=0.;
  FILE *pfseed, *FILE1, *FILE2;
  unsigned short seedv[3], *pseed;
  double ran1();
  int count1, count2, **count, flag, ref1, alt1, ref2, alt2;

  if (argc != 4) {
    printf("Usage: countFstPi k pop1 pop2\n");
    exit(0);
  }
  k = atoi (argv[1]);
  count = (int **) malloc (SMAX * sizeof (int *));
  for (a=0; a<SMAX; ++a)
    count[a] = (int *) malloc (2*k * sizeof (int));
  poplist1 = (int *) malloc ((k+1) * sizeof (int));
  for (a=0; a<=k; ++a)
    poplist1[a] = 0;
  poplist2 = (int *) malloc ((k+1) * sizeof (int));
  for (a=0; a<=k; ++a)
    poplist2[a]= 0;
 
  FILE1 = fopen (argv[2], "r");
  fscanf (FILE1, "%d ", &n1);
  for (a=0; a<n1; ++a) {
    fscanf (FILE1, "%d ", &temp1);
    poplist1[temp1] = 1;
  }
  fclose (FILE1);
  FILE2 = fopen (argv[3], "r");
  fscanf (FILE2, "%d ", &n2);
  for (a=0; a<n2; ++a) {
    fscanf (FILE2, "%d ", &temp1);
    poplist2[temp1] = 1;
  }
  fclose (FILE2);

  pfseed = fopen("seedms","r");
  for(a=0;a<3;a++) fscanf(pfseed," %hd",seedv+a);
  fclose( pfseed);
  seed48( seedv );

  while (scanf ("%d\t%d", &count1, &count2) != EOF) {
    count[S][0] = count1;
    count[S][1] = count2;
    for (a=1; a<k; ++a) 
      scanf ("\t%d\t%d", &count[S][2*a], &count[S][2*a+1]);
    scanf ("\n");
    ++S;
  }
      
  for (a=0; a<S; ++a) {
    ref1 = ref2 = alt1 = alt2 = 0;
    for (b=0; b<k; ++b) {
      if (count[a][2*b]+count[a][2*b+1]>0 && poplist1[b+1]>0)
	{
	  flag = count[a][2*b]+count[a][2*b+1];
	  allele = (int) (ran1() * flag);
	  if (allele < count[a][2*b])
	    ++ref1;
	  else
	    ++alt1;
	}
      else if (count[a][2*b]+count[a][2*b+1]>0 && poplist2[b+1]>0)
	{
	  flag = count[a][2*b]+count[a][2*b+1];
	  allele = (int) (ran1() * flag);
	  if (allele < count[a][2*b])
	    ++ref2;
	  else
	    ++alt2;
	}
    }
    /*
    printf("%d %d %d %d\n", ref1, alt1, ref2, alt2);
    */
    if ((ref1+alt1)>1 && (ref2+alt2)>1) {
      piW1 += (ref1*alt1) / ((ref1+alt1) * (ref1+alt1-1) * .5);
      piW2 += (ref2*alt2) / ((ref2+alt2) * (ref2+alt2-1) * .5);
      piB += (double) (ref1*alt2 + ref2*alt1) / ((ref1+alt1) * (ref2+alt2));
      ++sites;
    }
  }

  printf("piW1 = %.3f\npiW2 = %.3f\npiB = %.3f\nFst = %.3f\n#sites = %d\n",
	 piW1, piW2, piB, 1. - (piW1+piW2) / (2. * piB), sites);

  pfseed = fopen("seedms","w");
  pseed = seed48(seedv);
  fprintf(pfseed,"%d %d %d\n",pseed[0], pseed[1],pseed[2]);
  fclose( pfseed);
}
