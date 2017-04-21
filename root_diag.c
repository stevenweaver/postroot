/* Program to create data for the root diagnostic  */
/* described in Boykin et al. 2006. Comparison of  */
/* methods for rooting phylogenetic trees: A case  */
/* study using Orcuttieae (Poaceae: Choloridoideae)*/

/* Written by L. Salter Kubatko, February 2006     */

#include <stdio.h>
#include <stdlib.h>

int main() {

  int i, incr=100, total_it, skip, ntax, burnin, nlisttax;
  char *taxlist[200];
  FILE *testfile, *awkfile, *tf;

  tf = fopen("rootpp.dat","r");
  fclose(tf);
  
  printf("\n\nThis program will create a data set that can be used to  \n");
  printf("construct a diagnostic plot for convergence of MrBayes \n");
  printf("with respect to root position. You must have already run \n");
  printf("the post_root program on your output from MrBayes,and the\n");
  printf("resulting rootpp.dat file must be stored in this directory.\n");
  /*
  printf("directory, and you must copy your .nex.t file to a\n");
  printf("file called posttrees in this directory and remove\n");
  printf("the header information.\n");

  printf("\n  Please enter the number of taxa: ");
  scanf("%d",&ntax);

  printf("\n  Please enter the total number of iterations (ngen) you collected (no commas): ");
  scanf("%d",&total_it);

  printf("\n  Please enter the number of iterations discarded as burnin (burnin) (no commas): ");
  scanf("%d",&burnin);

  printf("\n  Please enter the skip value (samplefreq) (no commas): ");
  scanf("%d",&skip);

  printf("\n\nIt is inefficient to tabulate root proportions for every\n");
  printf("saved iteration of the MCMC, unless your run was very short.\n");
  printf("It is more convenient to compute posterior root probabilities \n");
  printf("for, say, every 100th saved iteration.\n");

  printf("\n  Please enter the increment at which you'd like to compute \n");
  printf("  posterior root probabilities (default is 100): ");
  scanf("%d",&incr);
  
  printf("\n\nTabulating posterior root probabilities for \n");
  printf("every %dth saved iteration (this will take a few minutes) ....\n\n",incr);

  for (i=incr; i<=(total_it-burnin)/skip; i=i+incr) {
    testfile = fopen("testfile","w");
    fprintf(testfile,"%d %d\n",ntax,i);
    fclose(testfile);
    system("post_root < testfile > treefile");
    system("cat rootpp.dat treefile > rootpp2.dat");
    system("/bin/cp -f rootpp2.dat rootpp.dat");
    system("/bin/rm -f rootpp2.dat");
  }

  printf("\n\nRoot positions have been tabulated. Data for each \n");
  printf("iteration are saved in the file rootpp.dat.\n\n");
  */


  printf("\nTo create a data set that can be easily plotted you must\n");
  printf("specify the numbers corresponding to the taxa on one side of\n");
  printf("the root of interest. Use the header at the top of the \n");
  printf(".nex.t file to find these numbers.  Enter them in a list \n");
  printf("separated by spaces, and use the side of the root that\n"); 
  printf("contains taxon 1. Enter a space and then a # at the end of \n");
  printf("your list.\n");

  awkfile = fopen("awkfile","w");
  fprintf(awkfile,"/");

  printf("\n  Enter the number of taxa in your list: ");
  scanf("%d",&nlisttax);
  printf("\n  Enter the taxon list here: ");
  for (i=0; i<nlisttax-1; i++) {
    scanf("%s",&taxlist);
    fprintf(awkfile,"%s ",taxlist);
  }
  
  scanf("%s",&taxlist);
  fprintf(awkfile,"%s",taxlist);
  fprintf(awkfile,"/ && NF==%d {print;}\n",nlisttax+2);
  fclose(awkfile);
  system("awk -f awkfile rootpp.dat > rootdiag.dat");
  
  printf("\n\nData to create the graph the root diagnostic have been written\n");
  printf("to the file rootdiag.dat. To create the plot, graph the first column\n");
  printf("of data (generation number) on the x-axis, and the second column \n");
  printf("divided by the first column on the y-axis.  This forms a plot of \n");
  printf("posterior probability of the root of interest versus generation number.\n");
  printf("If several roots are of interest, run this program multiple times on \n");
  printf("the same rootpp.dat file.\n\n");

  printf("Bye-bye!\n\n");

  return(1);

}
