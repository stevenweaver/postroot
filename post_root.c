/* This program computes the frequency of each rooting  */
/* in the posterior distribution of trees obtained from */
/* MrBayes.  The input file is the *.out.t file from    */
/* MrBayes, with the top lines removed and a sinlge     */
/* line displaying the number of taxa and the number of */
/* trees in the *.out.t file.                           */

/* Written by L. Salter, July 3, 2003                   */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <string.h>
#include <assert.h>
#include "post_root.h"


int *ppTwoRow[2], **ppMatrix, *filled_ind, *parents, nseq;
char *aa;
int ***pppTR, **rootlist, *rootlist_curr;
char **matnumlist;


int FindParentI(int target) {

  int i, j, parent=0;

  for (j=0; j<2; j++) {
    for (i=0; i<nseq; i++) {
      if (ppTwoRow[j][i]==target) {
	parent = i+nseq+1;
	break;
      }
    }
  }

  /*  if (parent==0) { printf("Could not find parent for node %d\n",target); }*/
  return(parent);

}


int NextNode(int next) {

  int stop=0;

  while (stop!=1) {

    if (filled_ind[next+1-(nseq+1)] != 1) {

      next+=1;
      stop=1;

    }

    else next=+1;

  }

  return next;

}

double ReadLength(FILE* fp) {

  float length=0.0;

  fscanf(fp,"%8f",&length);
  return (double)length;

}

int AddToNode(int current,int next,int last_tip,FILE* fp,double read_length) {

  int sss, junk, i, ss2, ss3, toadd;
  double temp;

  sss = fgetc(fp); 

  if (ppTwoRow[0][current-(nseq+1)] == 0) i = 0;
  else i=1;

  if (i==0) {

    ppTwoRow[i][current-(nseq+1)] = last_tip-48;
    /*ppLengthMat[current][last_tip-48] = ppLengthMat[last_tip-48][current] = read_length;*/
    if (sss!='(') {

      ss2 = fgetc(fp);
      if (ss2 == ':') {
            ungetc(ss2,fp);
            toadd = sss;
      }
      else {
            ss3 = fgetc(fp);
            if (ss3 == ':') {
              ungetc(ss3,fp);
              toadd = 10*(sss-48) + ss2;
            }
            else toadd = 100*(sss-48) + 10*(ss2-48) + ss3;
      }
      junk = fgetc(fp);
      temp = ReadLength(fp);
      ppTwoRow[1][current-(nseq+1)] = toadd-48;
      /*ppLengthMat[current][toadd-48] = ppLengthMat[toadd-48][current] = temp;*/
      sss = fgetc(fp); 
      filled_ind[current-(nseq+1)] = 1;
      
    }

  }

  if (i==1 && sss!='(') {
    
    ss2 = fgetc(fp);
    if (ss2 == ':') {
      ungetc(ss2,fp);
      toadd = sss;
    }
    else {
      ss3 = fgetc(fp);
      if (ss3 == ':') {
	ungetc(ss3,fp);
	toadd = 10*(sss-48) + ss2;
      }
      else toadd = 100*(sss-48) + 10*(ss2-48) + ss3;
    }
    /*sss = toadd;*/
    ppTwoRow[1][current-(nseq+1)] = toadd-48;
    filled_ind[current-(nseq+1)] = 1;
  }

  return sss;

}

int FindLastOpen(int current) {

  int node, stop=0;

  node = current-1;

  while(node>=nseq+1 && stop==0) {

    if (filled_ind[node-(nseq+1)]==1) node = node-1;
    else stop=1;

  }

  return node;

}

int CloseBack(int open, int current, int last_char, int last_tipp, FILE* fp, double read_length) {

  int sss;
  int i, j;

  sss=getc(fp);

  if (sss!=';') {
  
    if (ppTwoRow[0][open-(nseq+1)] == 0) j=0;
    else {
      
      j=1;
      filled_ind[open-(nseq+1)] = 1;
      
    }
    
    if (sss == ':') {

      read_length = ReadLength(fp);
      sss = fgetc(fp); 

    }
      
    ppTwoRow[j][open-(nseq+1)] = current;  
    /*ppLengthMat[open][current] = ppLengthMat[current][open] = read_length;*/
   
  }
    
  return sss;
 
}

void ReadTreePrior(FILE* fp) {

  int ss, ss1, ss2, ss3, last_ss, lastlast_ss, last_tip;
  int next_avail = nseq+1, last_open=nseq+1, curr_node=nseq+1;
  int parent;
  double read_length;

  ss = fgetc(fp); 
  last_ss = ss;
  ss = fgetc(fp); 

  while (!feof(fp) && ss!=';') {

    if (ss == '(') {

      next_avail = NextNode(next_avail);
      curr_node = next_avail;
      last_ss = ss;
      ss = fgetc(fp); 

    }

    else {

      if (ss == ',') {

	last_ss = ss;
	ss = AddToNode(curr_node,next_avail,last_tip,fp,read_length);

      }

      else {

	if (ss == ')') {

	  last_open = FindLastOpen(curr_node);
	  lastlast_ss = last_ss;
	  last_ss = ss;
	  ss = CloseBack(last_open,curr_node,lastlast_ss,last_tip,fp,read_length);
	  if (ss!=';' && filled_ind[last_open-(nseq+1)]==1) curr_node = FindParentI(curr_node);
	  else curr_node = last_open;

	}
	
	else {

	  ss2 = fgetc(fp);
	  if (ss2 == ':') {
	    ungetc(ss2,fp);
	    last_tip = ss;
	  }
	  else {
	    ss3 = fgetc(fp);
	    if (ss3 == ':') {
	      ungetc(ss3,fp);
	      last_tip = 10*(ss-48) + ss2;
	    }
	    else last_tip = 100*(ss-48) + 10*(ss2-48) + ss3;
	  }
	  last_ss = last_tip;
	  ss = fgetc(fp); 
	  if (ss == ':') read_length = ReadLength(fp);
	  parent = FindParentI(last_tip-48);
	  /*if (parent != 0) ppLengthMat[parent][last_tip-48] = ppLengthMat[last_tip-48][parent] = read_length; */
	  ss = fgetc(fp); 
	  
	}

      }

    }
      
  }

}



/***  Function to search through ppTwoRow to find parent ***/
/***  of a specified node.                               ***/

int find_parent(int target) {

int i, j, parent=0;

 for (j=0; j<2; j++) {
   for (i=0; i<nseq; i++) {
     if (ppTwoRow[j][i]==target) {
      parent = i+nseq+1;
      break;
     }
   }
 }

 /*if (parent==0) { printf("Could not find parent in find_parent for target = %d\n",target); }*/
 return(parent);

}



/** Make the parents vector **/

void make_parents () {

  int j, par;
  
    for (j=1; j<nseq+1; j++) parents[j] = find_parent(j);
    for (j=nseq+2; j<2*nseq+1; j++) parents[j] = find_parent(j);

}



/***  Function to find the generation of the input node   ***/

int find_gen(int current_node) {

  int gen_counter=0, parent;

    parent = current_node;

    while (parent != nseq+1) {

      parent = parents[parent];
      gen_counter++;

    }

    return gen_counter;

}




/***  Function to construct the IND matrix  ***/

void make_indmat() {

  int q, r;
  int parent1, parent2, check=0;
  int dist1=0, dist2=1;

  for (q=1; q<nseq+1; q++) {

    for (r=q+1; r<nseq+1; r++) {

      parent1 = parents[q];
      parent2 = parents[r];

      while (parent1 != nseq+1 && check != 1) {

        while (parent2 != nseq+1 && check != 1) {

          if (parent1 == parent2) {

            check = 1;
            break;

          }

          parent2 = parents[parent2];
          if (check != 1) dist2++;

        }

        if (check == 1) break;

        else {

          parent1 = parents[parent1];
          dist1++;
          dist2 = 1;
          parent2 = parents[r];
          
        }

      }

      ppMatrix[q][r] = dist1 + dist2;

      if ((parent1 == nseq+1 || parent2 == nseq+1) && check == 0) {

        ppMatrix[q][r] = find_gen(q) + find_gen(r) - 1;
	if (q==1) rootlist_curr[r]=0;
        check = 1;

      }
 
      check = 0;
      dist1 = 0;
      dist2 = 1;

    }

    check = 0;
    dist1 = 0;
    dist2 = 1;

  }

}




/***  Convert IND matrix to a number representing the tree.  This   ***/
/***  number is a unique representation of the tree and is used to  ***/
/***  compare trees                                                 ***/

void make_matnum() {

  int q, r, temp;
  int counter=0;

  for (q=1; q<nseq+1; q++) {

    for (r=q+1; r<nseq+1; r++) {

      aa[counter] = (char)((ppMatrix[q][r]/100)+48);
      temp = ppMatrix[q][r]%100;
      aa[counter+1] = (char)(temp/100+48);
      aa[counter+2] = (char)((temp%100)+48);

      counter = counter + 3;

    }

  }

  aa[counter] = '\0';

}


int main() {
  
  FILE *ptrees,*rd;
  int i, j, k, l, ii, kk, ll, numtreesprior;
  int find_next, check, rootcheck1, rootcheck2, check_count;
  int num_unique = 0, num_roots = 0;
  long numtreestotal;

  ptrees=fopen("posttrees","r");
  fscanf(ptrees,"%d %d",&nseq,&numtreesprior);

  printf("\n\nReading trees from file posttrees .......\n");
  printf("The number of sequences is %d, and the number of trees is %d\n\n",nseq, numtreesprior);

  rd = fopen("rootpp.dat","w");
  
  /* Allocate memory for trees and counts */
  filled_ind = (int*)malloc((2*nseq)*sizeof(int));
  for (i=0; i<2; i++)
    {
      ppTwoRow[i] = (int*)malloc((2*nseq)*sizeof(int));
      if (ppTwoRow[i]==NULL)
	{
	  printf("Can't memalloc ppTwoRow[%d]\n",i);
	  break;
	}
    }
  parents = (int*)malloc((2*nseq+1)*sizeof(int));
  if (parents==NULL) 
    {
      printf("Can't memalloc parents\n");
    }
  ppMatrix = (int**)malloc(nseq*sizeof(int*));
  if (ppMatrix==NULL)
    {
      printf("Can't memalloc ppMatrix\n");
    }

  for (i=0; i<nseq; i++)
    {
      ppMatrix[i] = (int*)malloc(nseq*sizeof(int));
      if (ppMatrix[i]==NULL)
        {
          printf("Can't memalloc ppMatrix[%d]\n",i);
          break;
        }
    }
  aa = (char*)malloc(((nseq*(nseq-1)/2)*3+1)*sizeof(char));
  if (aa==NULL) {
    printf("can't memalloc aa\n");
    exit(1);
  }
  pppTR = (int***)malloc((2*numtreesprior+1)*sizeof(int**));
   if (pppTR==NULL) {
    printf("Can't memalloc pppTR\n");
   }
   for (i=0; i<2*numtreesprior+1; i++) {
    pppTR[i] = (int**)malloc(2*sizeof(int*));
    for (j=0; j<2; j++) {
      pppTR[i][j] = (int*)malloc((nseq)*sizeof(int));
    }
   }
   matnumlist = (char**)malloc((2*numtreesprior+1)*sizeof(char*));
   for (i=0; i<2*numtreesprior; i++) {     
     matnumlist[i] = (char*)malloc(((nseq*(nseq-1)/2)*3+1)*sizeof(char));
   }
   rootlist = (int**)malloc((2*numtreesprior+1)*sizeof(int*));
   for (i=0; i<2*numtreesprior; i++) {
     rootlist[i] = (int*)malloc((nseq+2)*sizeof(int));
   }
   rootlist_curr = (int*)malloc((nseq+2)*sizeof(int));


   /* printf("Reading tree number ");
      fflush(0);*/

   for (k=0; k<numtreesprior; k++) {
   
     for (i=0; i<nseq+2; i++) {

       ppTwoRow[0][i]=0;
       ppTwoRow[1][i]=0;
       filled_ind[i]=0;
       filled_ind[i+nseq]=0;
       parents[i]=0;
       parents[i+nseq]=0;
       
     }


     find_next = fgetc(ptrees);

     while (find_next!=61) {
       find_next = fgetc(ptrees);
     }
     
     ReadTreePrior(ptrees);


     for (i=0; i<2; i++) {    
       for (j=0; j<nseq; j++) {
	 
	 ppTwoRow[i][j] = ppTwoRow[i][j+1];
	 if (ppTwoRow[i][j]>nseq+1) ppTwoRow[i][j] = ppTwoRow[i][j]-1;
	 /*printf("%d ",ppTwoRow[i][j]);*/

       }       
       
       /*    printf("%d\n",ppTwoRow[i][j]);
	     fflush(0);*/

     }


     for (ii=0;ii<nseq+1;ii++) rootlist_curr[ii]=1;

     make_parents();
     make_indmat();
     make_matnum();

     /*for (ii=2; ii<nseq+1; ii++) printf("%d %d\n",ii,rootlist_curr[ii]);*/
     
     check = 1;
     
     /*  printf(" %d",k+1);
	 fflush(0);*/


     /* Add to list of new roots if it hasn't come up before */
     /* Not the same as unique trees */
     
     if (k==0) {
       for (ii=2; ii<nseq+1; ii++) rootlist[0][ii] = rootlist_curr[ii];
       rootlist[0][nseq+1] = 1;
       num_roots+=1;
     }

     else {

       rootcheck1 = 1;
       if (k>0) ll=0;
       else ll=1;
       
       /* Check to see if the rooting is unique */
       
       while (rootcheck1==1 && ll<num_roots) {
	 rootcheck2=1;
	 for (ii=2; ii<nseq+1; ii++) {
	   if (rootlist[ll][ii] != rootlist_curr[ii]) rootcheck2 = 99;
	 }
	 if (rootcheck2==1) {
	   rootlist[ll][nseq+1] += 1;
	   rootcheck1=99;
	 }
	 ll += 1;
       }
       
       if (ll==num_roots && rootcheck1!=99) {
	 
	 for (ii=2; ii<nseq+1; ii++) rootlist[num_roots][ii] = rootlist_curr[ii];
	 rootlist[num_roots][nseq+1] = 1;
	 num_roots += 1;
	 
       }

     }

     /* Try printing each time through loop for root_diag */
      
     check_count = 0;
     for (kk=0; kk<num_roots; kk++) {
       
       fprintf(rd,"    %d     %d         1 ",k+1,rootlist[kk][nseq+1]);
       for (ii=2; ii<nseq+1; ii++) {
	 if (rootlist[kk][ii]==1) fprintf(rd,"%d ",ii);
       }
       fprintf(rd,"\n");
       
       check_count += rootlist[kk][nseq+1];
       
     }
     fprintf(rd,"\n");

   }

    printf("The number of unique rootings was %d\n\n\n",num_roots);

   check_count = 0;

   printf("Frequency   Taxa on one side\n");

   for (k=0; k<num_roots; k++) {

     /*printf("Rooting %d occurred %d times - the taxa to the left of the root are:\n",k+1,rootlist[k][nseq+1]);*/
     printf("    %d         1 ",rootlist[k][nseq+1]);
     for (ii=2; ii<nseq+1; ii++) {
       if (rootlist[k][ii]==1) printf("%d ",ii);
     }
     printf("\n");
     
     check_count += rootlist[k][nseq+1];

     }


   printf("\n\n%d trees have been tallied.\n",check_count);
   printf("\nRoot probabilities at each generation have been written to the file\nrootpp.dat.  These can be used with the root_diag program to check \nconvergence with respect to root position."); 

   printf("\n\n");

   printf("Bye-bye!\n\n");

   
}
















