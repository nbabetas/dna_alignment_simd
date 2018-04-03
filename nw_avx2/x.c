/***********************************
  Optimal Global Alignment Program 
     written by Benjamin Corum
       email: ben@corum.org	 
             3/14/01			 
***********************************/

/*
THIS SOFTWARE WAS DESIGNED FOR AN UNDERGRADUATE COMPUTATIONAL BIOLOGY TUTORIAL.
IN NO EVENT SHALL THE AUTHOR BE LIABLE TO ANY PARTY FOR DAMAGES
ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION,
EVEN IF THE AUTHOR HAS BEEN ADVISED OF THE POSSIBILITY DAMAGE.
*/

/*
This program creates a global alignment of two sequences of symbols.
Using the Needleman-Wunsch Algorithm all possible alignments are scored
on a 2 dimensional matrix, and an optimal sequence is returned to the user.

To run this program correctly 2 FASTA formatted test files named 'str1.fa' and 'str2.fa' 
also must be downloaded or created and placed into the same folder as the executable.


This program was designed to display a 20*20 NW matrix on the screen at run time. This feature
seriously limits its ability to align longer sequences. By removing the displayed matrix and 
expanding certain character arrays, longer sequences can be aligned.
*/


//Library Files
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*File pointers*/
FILE *ptr_file_1, *ptr_file_2;	//pointer to test files

/*Definitions*/
#define TRUE 1
#define FALSE 0

/*Global Variables*/
char inputC[5];		//user input character
int inputI;			//integer value of user input
int intcheck = TRUE;	//used to check validity of 1st user input

char holder;			//used to measure file lengths
char ch;
int filelen1 = 0;		//measured file1 length
int filelen2 = 0;		//measured file2 length
int i,j;				//counters
int lenA,lenB;			//string lengths
char FASTA1[200];		//test file1 --> string
char FASTA2[200];		//test file2 --> string 								
											
int HiScore;			//HiScore holder
int HiScorePos[2];		//position of ^ holder

int MaxAdd;			//maximum score of column and row 1 below & to the right 

char dash = '-';		//character used in additionof dashes to strings

char strA[20];				//holds 1st string to be aligned in character array		
char strB[20];				//holds 2nd string to be aligned in character array 
char MaxAlignmentA[50];		//holds final alignment	of string A	
char MaxAlignmentB[50];		//holds final alignment of string B

int NWArray[21][21];		//N-W matrix

int PositionA;
int PositionB;

/*Matrix Rescoring Function*/
int MaxLookUp(int Row, int Column) {
	
	int MaxVal;			//integer to hold maximum value found
	static int rc;		//internal row counter				
	static int cc;		//internal column counter

	if((Row == lenB) && (Column == lenA)) {	//if cell is the bottom corner of the matrix 
		MaxVal = 0;							//no column or rows below to the right (no score)
	}
												
	if((Row == lenB) || (Column == lenA)) {	//if cell is along bottom or right side
		MaxVal = 0;							//no column or rows also (no score)
	}

	else {

		for(rc = Row+1; rc <= lenB; ++rc) {			//scan row +1 
		
			if(NWArray[rc][Column+1] > MaxVal) {
				MaxVal = NWArray[rc][Column+1];		//& return greatest value in that row
			}
		}
		
		for(cc = Column+1; cc <= lenA; ++cc) {		//scan column +1
		
			if(NWArray[Row+1][cc] > MaxVal) {		//& if column score > maxVal reset maxVal
				MaxVal = NWArray[Row+1][cc];
	
			}
		}											
	}	
	return MaxVal;
}

/*Make Alignment Function*/
int MakeAlignment(PositionA,PositionB) {

	/*Function Variables*/
	int strAcounter = 0;
	int strBcounter = 0;
	int MaxAcounter = 0;
	int MaxBcounter = 0;

	int relmax; 
	int relmaxpos[2];


	while((PositionA != lenA) && (PositionB != lenB)) {	//until the bottom or far right of the matrix is reached

		relmax = -1;	//reset relmax for every new position
		
		//FIND HIGHEST SCORING POSITION IN THE NEXT COLUMN OR ROW
		for(i = PositionA + 1; i < lenA + 1; ++i) {		//scan +1 column
			
			if(NWArray[i][PositionB + 1] > relmax) {	//row changes column stays the same
				relmax = NWArray[i][PositionB + 1];
				relmaxpos[0] = i;
				relmaxpos[1] = (PositionB + 1);
			}
		}
		
		for(i = PositionB + 1; i < lenB + 1; ++i) {		//scan +1 row
		
			if(NWArray[PositionA + 1][i] > relmax) {	//column changes row stays the same
				relmax = NWArray[PositionA + 1][i];
				relmaxpos[0] = (PositionA + 1);
				relmaxpos[1] = i;
			}
		}

		/*Construct alignment up to relmaxpos*/
		if((relmaxpos[0] == (PositionA + 1)) && (relmaxpos[1] == (PositionB + 1))) { //if the next alignment is a "true" diagonal
		
			MaxAlignmentA[MaxAcounter] = strA[strAcounter];
			++strAcounter;
			++MaxAcounter;
			MaxAlignmentB[MaxBcounter] = strB[strBcounter];
			++strBcounter;
			++MaxBcounter;
			
		}

		else {
			
			if(relmaxpos[0] > (PositionA + 1)) {	//alignment skips a row so spaces must be assigned to strB

				for(i = 0; i < (relmaxpos[0] - (PositionA + 1)); ++i) {

					MaxAlignmentA[MaxAcounter] = strA[strAcounter];		//enter strA elements into maxAlignmentA
					MaxAlignmentB[MaxBcounter] = dash;					//enter dash(es) into MaxAlignmentB
					++strAcounter;										//increment all counters except strBcounter
					++MaxAcounter;										
					++MaxBcounter;
				}

				for(i = 0; i < 1; ++i) {								//make alignment at relmax position
					MaxAlignmentA[MaxAcounter] = strA[strAcounter];
					++strAcounter;
					++MaxAcounter;
					MaxAlignmentB[MaxBcounter] = strB[strBcounter];
					++strBcounter;
					++MaxBcounter;
				}
			}
			
			if(relmaxpos[1] > (PositionB + 1)) {	//alignment column a row so spaces must be assigned to strA

				for(i = 0; i < (relmaxpos[1] - (PositionB + 1)); ++i) {

					MaxAlignmentB[MaxBcounter] = strB[strBcounter];		//enter strB elements into MaxBAlignment
					MaxAlignmentA[MaxAcounter] = dash;					//enter dash(es) into MaxAlignmentA		
					++strBcounter;										//increment all counters except for strAcounter
					++MaxAcounter;
					++MaxBcounter;
				}

				for(i = 0; i < 1; ++i) {
					MaxAlignmentA[MaxAcounter] = strA[strAcounter];
					++strAcounter;
					++MaxAcounter;
					MaxAlignmentB[MaxBcounter] = strB[strBcounter];
					++strBcounter;
					++MaxBcounter;
				}
			}
		}
	PositionA = relmaxpos[0];
	PositionB = relmaxpos[1];

	}
	if(strAcounter < lenA+1) {	//letters need to be appended to MaxAlignmentA
		for(i=strAcounter; i<=lenA; ++i) {
			MaxAlignmentA[MaxAcounter] = strA[strAcounter];
			MaxAlignmentB[MaxBcounter] = dash;
			++MaxAcounter;
			++strAcounter;
			++MaxBcounter;
		}
	}
	else {

		if(strBcounter < lenB+1) {	//letters need to be appended to MaxAlignmentB
			for(i=strBcounter; i<=lenB; ++i) {
				MaxAlignmentB[MaxBcounter] = strB[strBcounter];
				MaxAlignmentA[MaxAcounter] = dash;
				++MaxBcounter;
				++strBcounter;
				++MaxAcounter;
			}
		}
	}
}
/*MAIN FUNCTIONS*/
int main()
{

	
	/*Introduction*/
	printf("Needleman-Wunsch Global Alignment Program Version 1.0\n");
	printf("Designed by: Benjamin Corum -- Marlboro College -- Spring 2001\n");
	printf("\n");
	printf("Enter:\n");

	/*Validate User Input*/
	while(intcheck == TRUE) {				
		
		/*list options*/
		printf("   1 - To align two test files\n");
		printf("   2 - To align two novel sequences\n");
		printf("   3 - To align a novel sequence with a test file\n:");

		fgets(inputC, sizeof(inputC), stdin);	//get string
		sscanf(inputC, "%i", &inputI);			//convert to integer
		
		if((inputI < 4) && (inputI > 0)){		//if its a valid user input
			intcheck = FALSE;
		}
		else {
			printf("Please enter a 1,2,or 3\n"); //if its not start over
			intcheck = TRUE;
		}
	}
	
	
	if(inputI == 1) {	//If user would like to align two test files	
		
		/*open first file*/
		ptr_file_1 = fopen("../sw_avx2_data/sw_avx2_data/ref_150.fasta", "r");
		
		/*check to see if it opened okay*/
		if(ptr_file_1 == NULL) {
			printf("Error opening 'str1.fa'\n");
			exit(8);
		}

		/*open second file*/
		ptr_file_2 = fopen("../sw_avx2_data/sw_avx2_data/reads_150.fasta", "r");
		
		/*check to see if it opened okay*/
		if(ptr_file_2 == NULL) {
			printf("Error opening 'str2.fa'\n");
			exit(8);
		}

		/*retrieve fasta information from files*/
		fgets(FASTA1, sizeof(FASTA1), ptr_file_1);
		fgets(FASTA2, sizeof(FASTA2), ptr_file_2);		

		/*measure file1 length*/
		while(holder != EOF) {
			holder = fgetc(ptr_file_1);
			if(filelen1 < 20) {
				strA[filelen1] = holder;
			}
			++filelen1;
		}

		holder = '0';

		/*measure file2 length*/
		while(holder != EOF) {
			holder = fgetc(ptr_file_2);
			if(filelen2 < 20) {
				strB[filelen2] = holder;
			}
			++filelen2;
		}

		fclose(ptr_file_1);
		fclose(ptr_file_2);								
	}
							
	if(inputI == 2) {	//two novel strings
		
		/*get string A*/
		printf("Novel strings should be no longer than 20 symbols. Program is case sensitive\n");
		printf("Please enter 1st string:");
		fgets(strA,sizeof(strA),stdin);
		printf("Please enter 2nd string:");
		fgets(strB,sizeof(strB),stdin);

		if(strA[strlen(strA)-1] = '\n') {
		strA[strlen(strA)-1] = '\0';	//edit '\n' out of string
		}
		if(strB[strlen(strB)-1] = '\n') {
		strB[strlen(strB)-1] = '\0';	//edit '\n' out of string
		}

	}

	if(inputI == 3) {	//novel string vs. test
		
		/*get string A*/
		printf("Novel strings should be no longer than 20 symbols. Program is case sensitive\n");
		printf("Please enter 1st string:");
		fgets(strA,sizeof(strA),stdin);					
		strA[strlen(strA)-1] = '\0';	//edit '\n' out of string

		/*open test file*/
		ptr_file_1 = fopen("str1.fa", "r");
		
		/*check to see if it opened okay*/
		if(ptr_file_1 == NULL) {
			printf("Error opening 'str1.fa'\n");
			exit(8);
		}
		
		/*get fasta info from test file*/
		fgets(FASTA1, sizeof(FASTA1), ptr_file_1);
	
		/*measure file1 length*/
		while(holder != EOF) {
			holder = fgetc(ptr_file_1);
			if(filelen1 < 20) {
				strB[filelen1] = holder;
			}
			++filelen1;
		}
	
	}

	/*set variables to length of strings A & B*/
	lenA = strlen(strA)-1;
	lenB = strlen(strB)-1;


	/*Fill in table*/
	for(i = 0; i <= lenA; ++i) {	//for all values of strA
		
		for(j = 0; j <= lenB; ++j) {	//& for all values of strB

			if(strA[i] == strB[j]) {	//if two characters match
				
				NWArray[i][j] = 1;		//assign a score of +1 to the pair in the NW matrix
					
			}

			else {
				NWArray[i][j] = 0;		//no match no score

			}
		}
	}

/*Rescore matrix*/	
	
	for(i = lenB; i >= 0; --i) {		//start in the bottom row moving up

		for(j = lenA; j >= 0; --j) {	//& start in the last column moving left
			
		
			MaxAdd = MaxLookUp(i,j);

			NWArray[i][j] = NWArray[i][j] + MaxAdd;

		}
	}

	/*PRINT NW matrix to screen*/
	printf(" ");	//opening space
	for(i = 0; i <= lenB; ++i) {	// print strB along the upper border of matrix
		printf("  %c",strB[i]);
	}
	printf("\n");	

	for(i = 0; i <= lenA; ++i) {
	
			printf("%c",strA[i]);	//print strA aling the left border of the matrix
	

		for(j = 0; j <= lenB; ++j) {
			printf("%3i",NWArray[i][j]);	//print NWArray[row][column]
		}

		printf("\n");
	}

	/*extract optimal global alignment from matrix*/
	PositionA = -1;
	PositionB = -1;

	MakeAlignment(PositionA,PositionB);
	printf("Here is your optimal global alignment:\n");
	printf("A)%s\nB)%s\n",MaxAlignmentA,MaxAlignmentB);

	
	return(0);

}

