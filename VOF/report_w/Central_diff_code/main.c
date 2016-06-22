/*#################################################################################
| FILE NAME:																									 |					
|																											 |																																				 			
| INTERFACE_SOLVER                                                                                                                       								 |
|																											 |
| DESCRIPTION:																								 |
| 1. Reconstructs a piecewise-linear interface from a given distribution of volume fractions using the LVIRA algorithm [1,2]			 |
| 2. Given a velocity field and an F distribution, it calculates the F distribution at the next time-step.								 |
|															              											 |
| LIMITATIONS:																								 |
| 1. Presently works only for an uniform grid.																		 |
| 2. Presently works only for 2D																					 |
| 3. Has to be modified for an open domain. See function CalculateF_Flux()													 |																										 
|                                                                                                                                         										 |		
| REFERENCES:                                                                                                                             									 |
| 1. A higher order projection method for tracking fluid interfaces in variable-density flows - JCP (130), 269-282 (1997), Puckett et al       |
| 2. Second-order volume-of-fluid algorithms for tracking material interfaces - JCP(199), 465-502 (2004), Pilliod et al                  		 |
| 3. Comparison of volume-of-fluid methods for surface-tension dominant two-phase flows - IJHMT(49), 740-754(2006), Gerlach et al       	 |
| 4. High-order surface-tension VOF model for 3D bubble flows with high density ratio - JCP (200), 153-176 (2004), Lorstad et al              |
##################################################################################


/***INCLUDE***/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <global_variables.h>
/**************/ 

/**FUNCTION DECLARATIONS*/
void        Pass_main(int, int);
void		InterfaceRecon_main();
void		Recon_Interface(int, int);
void		Write_struct_Cellattrib(int, int, long double, long double, long double, long double, long double, char[], long double, int);
void		Init_GlobalVars();
void		Calculate_InitialSlope();
void		Calculate_distance(int, int);
void		Identify_Shape(int, int);
void		Extrapolate_Line(int, int);
void		IdentifyCell_AssignArea(int, int, int, char[], long double);
long double Calculate_SumOfDiffSq(int, int);
int		Find_Quad(int, int);
void		Calculate_Theta(int, int);
void		Rotate_Line(int, int, char[]);
void		WriteToFile(int);
void		AdvanceF_field(int, int, int, int);
long double AdvanceF_WithCorrection(int, int, int, int, long double, long double);
long double AdvanceF_WithoutCorrection(int, int, int, int, long double, long double);
void		CalculateF_Flux(int);
long double CalculateX_Flux(char[], char[], long double, long double, long double, long double, long double, long double, long double, long double, long double);
long double CalculateY_Flux(char[], char[], long double, long double, long double, long double, long double, long double, long double, long double, long double);
void		AssignVal_FluxVariables(int, int, long double *, long double *, long double *, long double *, long double *, long double *, long double *, long double *, long double *, long double *, long double *);
void        BookKeeping();
/*************************/

main()
{
	   /*********DECLARATIONS***********/
	   int r, c, Step;
	   int Count;
	   int Togg;
	   int RunTimeStart, RunTimeEnd;
		RunTimeStart = time(0);
		
	   char LineFormat_debug[] = "./DebugData/VolumeFractions_%debug.dat";
  	   char LineFileName_debug[40];
	   FILE *ptrWrite_debug;

  	   char FileName_X[] = "./FError/X.dat";
	  // char FileName_Y[] = "./FError/Y.dat";
	   FILE *ptrWrite_X;
	  // FILE *ptrWrite_Y;
	   long double Tot_VolFrac;
	   /*****************/

	   /***************OPEN ALL FILE POINTERS*******/
	    ptrWrite_X =fopen(FileName_X,"w");
	   // ptrWrite_Y =fopen(FileName_Y,"w");
	   /*******************************************/
		
	   /*STEPS BEFORE CALLING RECON_INTERFACE() IN A LOOP*/	
	   Step  = 0;
	   Count = 0;
	   Togg  = 0;
	   Init_GlobalVars();	   	   
	   /**********************************************************/

 	    while(Step < 510)
 		{		
			   /*********CONSTRUCT INTERFACE*******************************/
			        //POPULATES NX AND NY INITIALLY FROM A GUESS BASED ON NEIGHBOURING CELLS VF
				Calculate_InitialSlope();
				
				for(r=0;r < NoOfRows ; r++)
				{
				   for(c=0;c < NoOfCols ; c++)
					{					  							
					  Recon_Interface(r,c);   //CALCULATE INTERFACE
					}
				}
	 			/*********************************************************/

				printf("\nRun in Progress. Step is %d\n", Step);

				//~ /***************OUTPUT FILE GENERATION**********/
				if(Step%1 == 0)
 				//if (Step>5770 && Step <=5780) 
				{	
 		       		 Tot_VolFrac = 0.0;	 	
  					 for(r=0;r < NoOfRows; r++)
					 {
					    for(c=0;c < NoOfCols ; c++)
						 {
							 Tot_VolFrac = Tot_VolFrac + (Cellattrib[r][c].VolFrac);
						 }
					 }
				
					 fprintf(ptrWrite_X,"%d\t%Lf\n", Step,Tot_VolFrac);
					fflush(ptrWrite_X);
					 //~ fprintf(ptrWrite_Y,"%Lf", Tot_VolFrac);
					 //~ fprintf(ptrWrite_Y,"\n");
 					 					
	 				 WriteToFile(Count);					// GENERATE Plot_Line{Step}.dat
					 
	// for debugging remove later				 
	 sprintf(LineFileName_debug,LineFormat_debug,Step);
	
	ptrWrite_debug =fopen(LineFileName_debug,"w");

	for(r=0;r < NoOfRows; r++)
	{
		for(c=0;c < NoOfCols; c++)
		{
			fprintf(ptrWrite_debug,"%Lf",Cellattrib[r][c].VolFrac);
			fprintf(ptrWrite_debug,"\n");
		}
	}		
	fclose(ptrWrite_debug);	
  				}
				 /***********************************************/
				
				if(Togg == 0)     //1ST HORIZONTAL AND 2ND VERTICAL PASS
				 {
					 Pass_main(1,Togg);
					InterfaceRecon_main();
					 Pass_main(2,Togg);
					 Togg = Togg + 1;
				}
				 else if (Togg == 1)   //1ST VERTICAL AND 2ND HORIZONTAL PASS
				 {
					 Pass_main(2,Togg);
					 InterfaceRecon_main();
					Pass_main(1,Togg);
					 Togg = Togg - 1;
				 }
						
				//~ /*INCREMENT COUNTERS*/
				Step++;
				 Count++;
				/******************/							

	       }//END OF WHILE

	/*********WRITE THE LAST VOLUME FRACTION DATA INTO A FILE******/
	//~ sprintf(LineFileName_debug,LineFormat_debug,Step);
	
	//~ ptrWrite_debug =fopen(LineFileName_debug,"w");

	//~ for(r=0;r < NoOfRows; r++)
	//~ {
		//~ for(c=0;c < NoOfCols; c++)
		//~ {
			//~ fprintf(ptrWrite_debug,"%d\t%d\t%Lf",r,c,Cellattrib[r][c].VolFrac);
			//~ fprintf(ptrWrite_debug,"\n");
		//~ }
	//~ }		
	//~ fclose(ptrWrite_debug);	
    
	/*#################################################################*/
	
	/***CLOSE ALL FILE POINTERS*******/
	fclose(ptrWrite_X);			
	//fclose(ptrWrite_Y);			
	/************END**************/	
	RunTimeEnd= time(0);
	RunTimeEnd = RunTimeEnd - RunTimeStart;
	printf("The total time for run is %d seconds\n",RunTimeEnd);
		//RUN IS OVER		
}

/******************************************************************************************************
NAME       : Pass_main() 
DESCRIPTION: THIS IS A DUMY FUNCTION WHICH IS USED TO REMOVE CODE FROM main()
			 TOGG IS A VARIABLE USED FOR ALTERNATING BETWEEN VERTICAL AND HORIZONTAL PASSES
			 FLAG IS A VARIABLE USED FOR DECIDING WHETHER THE CURRENT PASS IS HORIZONATAL / VERTICAL
LIMITATIONS: NONE
/******************************************************************************************************/
void Pass_main(int Flag, int Togg)
{
	int r, c;

	CalculateF_Flux(Flag);   //CALCULATE X FLUXES
	for(r=0;r < NoOfRows; r++)
	{
	   for(c=0;c < NoOfCols; c++)
		{
		  AdvanceF_field(r,c,Flag,Togg);  //ADVANCE F IN TIME BASED ON FLAG VALUE AND HORIZONTAL FLUXES
		}
	}		
	
	BookKeeping();
}

/*********************************************************************************
NAME       : InterfaceRecon_main() 
DESCRIPTION: THIS IS A DUMY FUNCTION WHICH IS USED TO REMOVE CODE FROM main()
LIMITATIONS: NONE
/********************************************************************************/
void InterfaceRecon_main()
{
	int r,c;

	Calculate_InitialSlope();
	for(r=0;r < NoOfRows ; r++)
	{
	   for(c=0;c < NoOfCols ; c++)
		{					  							
		  Recon_Interface(r,c);   
		}
	}
}

/****************************************************************************************************************
METHOD NAME: Recon_Interface
DESCRIPTION: GIVEN A CELL, THIS FUNCTION FINDS OUT ST. LINES TO REPRESENT THE INTERFACE IN THAT CELL. THE OUTPUT
			 IS WRITTEN TO A FILE NAMED PLOTLINE_{STEP}.DAT. HERE {STEP} REPRESENTS AN INTEGER, CORRESPONDING
			 TO THE ACTUAL NO. OF TIME-STEPS FOR WHICH THE CALCULATION IS PERFORMED.
LIMITATIONS:    
****************************************************************************************************************/
void Recon_Interface(int r, int c)
{
   long double Sum_OfDiff_Sq_clock, Sum_OfDiff_Sq_anticlock;
   long double Sum_OfDiff_Sq_old, Sum_OfDiff_Sq_current;   
   char sense1[20], sense2[20];
   char sense[20], Shape_current[20];
   int flag = 1;
   long double Nx_current, Ny_current, Theta_current, P_current, Nx_clock, Ny_clock, Nx_anticlock, Ny_anticlock;

	strcpy(sense1,"Clockwise");
	strcpy(sense2,"Counterclockwise");

 	if(Cellattrib[r][c].VolFrac > LOWER_LIMIT && Cellattrib[r][c].VolFrac < UPPER_LIMIT)          
	{
		 Calculate_Theta(r,c);
		
		 Nx_current = Cellattrib[r][c].Nx;   //preserve the old value
		 Ny_current = Cellattrib[r][c].Ny;   //preserve the old value  		

		 Identify_Shape(r,c);       //find the shape depending on VolFrac and Slope
		 Calculate_distance(r,c);   //find the perp. distance depending on shape
		Theta_current  = Cellattrib[r][c].Theta; //preserve the old value
		 P_current = Cellattrib[r][c].P;
		 strcpy(Shape_current,Cellattrib[r][c].Shape);
		 Extrapolate_Line(r,c);   
		 Sum_OfDiff_Sq_current = Calculate_SumOfDiffSq(r,c); 		 
		
		 /*ROTATE CLOCKWISE AND CALCULATE SUM OF SQUARES*/
		 Rotate_Line(r,c,sense1);   //change the Nx and Ny values
		 Nx_clock = Cellattrib[r][c].Nx;   //preserve the clockwise value
		 Ny_clock = Cellattrib[r][c].Ny;   //preserve the clockwise value
		 Calculate_Theta(r,c);
		 //printf("\nAfter clockwise rotation the value of Theta is %f\n", Cellattrib[r][c].Theta);
		 Identify_Shape(r,c);
		 Calculate_distance(r,c);
		 Extrapolate_Line(r,c);
		 Sum_OfDiff_Sq_clock = Calculate_SumOfDiffSq(r,c);
		 /************************************************/

		 //Write_struct_Cellattrib(r,c,-100000,-100000,Nx_current,Ny_current,-100000,"-100000",-100000,-100000);  //restore the current values
		 Write_struct_Cellattrib(r,c,-100000,-100000,Nx_current,Ny_current,P_current,Shape_current,Theta_current,-100000);  //restore the current values
		 /*ROTATE COUNTER-CLOCKWISE AND CALCULATE SUM OF SQUARES*/
		 Rotate_Line(r,c,sense2);  //change the Nx and Ny values
		 Nx_anticlock = Cellattrib[r][c].Nx;   //preserve the clockwise value
		 Ny_anticlock = Cellattrib[r][c].Ny;   //preserve the clockwise value
		 Calculate_Theta(r,c);
		 //printf("\nAfter clockwise counterclockwise rotation the value of Theta is %f\n", Cellattrib[r][c].Theta);
		 Identify_Shape(r,c);
		 Calculate_distance(r,c);
		 Extrapolate_Line(r,c);
		 Sum_OfDiff_Sq_anticlock = Calculate_SumOfDiffSq(r,c);
		 /************************************************/

		 //Write_struct_Cellattrib(r,c,-100000,-100000,Nx_current,Ny_current,-100000,"-100000",-100000,-100000);  //restore the current values
		 Write_struct_Cellattrib(r,c,-100000,-100000,Nx_current,Ny_current,P_current,Shape_current,Theta_current,-100000);  //restore the current values
		 if( Sum_OfDiff_Sq_current >= Sum_OfDiff_Sq_clock)
		   {
			 strcpy(sense,sense1);
			 Write_struct_Cellattrib(r,c,-100000,-100000,Nx_clock,Ny_clock,-100000,"-100000",-100000,-100000);  //make the clockwise value as current
			 //printf("\n THE LINE SHOULD BE ROTATED %s\n", sense);
		   }
		 else if(Sum_OfDiff_Sq_current >= Sum_OfDiff_Sq_anticlock)
		   {
			 strcpy(sense,sense2);
			 Write_struct_Cellattrib(r,c,-100000,-100000,Nx_anticlock,Ny_anticlock,-100000,"-100000",-100000,-100000);  //make the clockwise value as current
			 //printf("\n THE LINE SHOULD BE ROTATED %s\n", sense);
		   }
		 else //CURRENT VALUE IS THE MINIMUM - NO NEED FOR ITERATIONS
		   {
			 return;
		   }

		 flag = 1;		
		 while(flag)
		   {		
				Theta_current  = Cellattrib[r][c].Theta; //preserve the old value
				 P_current = Cellattrib[r][c].P;
				 strcpy(Shape_current,Cellattrib[r][c].Shape);
			   
			 Rotate_Line(r,c,sense);
			 Calculate_Theta(r,c);
			 Identify_Shape(r,c);       //Find the shape depending on VolFrac and Slope
			 Calculate_distance(r,c);   //Find the perp. distance depending on shape
			 Extrapolate_Line(r,c);

			 Sum_OfDiff_Sq_old = Sum_OfDiff_Sq_current;
			 Sum_OfDiff_Sq_current = Calculate_SumOfDiffSq(r,c);

			 if (Sum_OfDiff_Sq_old <= Sum_OfDiff_Sq_current)
			{
				 flag = 0;
				//Write_struct_Cellattrib(r,c,-100000,-100000,Nx_current,Ny_current,P_current,Shape_current,Theta_current,-100000);  //restore the current values
				 //restore the previous value - TO BE DONE
			}
		}		
	}
}

/**************************************************************************************************
METHOD NAME: Write_struct_Cellattrib
DESCRIPTION: ALL WRITE OPERATIONS TO THE STRUCTURE STRUCT_CELL_ATTRIB HAPPEN THROUGH THIS FUNCTION.
             SEND -100000 IF NOTHING IS TO BE WRITTEN FOR A GIVEN PARAMETER
LIMITATIONS: NONE  
***************************************************************************************************/
void Write_struct_Cellattrib(int r,int c,long double VolFraction, long double Area, long double Nx, long double Ny, long double P, char Shape[30], long double Theta, int Type)
  {
    if (VolFraction != -100000)
      {
       Cellattrib[r][c].VolFrac = VolFraction;
      }
    if (Area != -100000)
      {
       Cellattrib[r][c].Area = Area;
      }
    if (Nx != -100000)
      {
        Cellattrib[r][c].Nx = Nx;
      }
    if (Ny != -100000)
      {
        Cellattrib[r][c].Ny = Ny;
      }
    if (P != -100000)
      {
        Cellattrib[r][c].P = P;
      } 
    if (!strcmp(Shape,"-100000"))   //strcmp RETURNS 0 IF BOTH THE STRINGS MATCH
      {
        //DO NOTHING
      }
    else
      {
        strcpy(Cellattrib[r][c].Shape,Shape);
      }
    if(Theta != -100000)
      {
        Cellattrib[r][c].Theta = Theta;
      }
	if(Type != -100000)
	{
		Cellattrib[r][c].Type = Type;
	}
  }

/****************************************************************
METHOD NAME       : Init_GlobalVars
DESCRIPTION          : INITIALISES ALL THE GLOBAL VARIABLES 
LIMITATIONS            : NONE
****************************************************************/
void Init_GlobalVars()
  {
	/*********DECLARE LOCAL VARIABLES*************/
	int r,c;
        long double VolFrac;
        char arr[50];
	/**********************************************/

	/******************READ THE INTIAL VOLFRAC DATA FROM FILE**************/
        FILE *ptrRead;
        ptrRead =fopen("circle.dat","rt");

       for(r=0;r < NoOfRows ; r++)
      {
		for(c=0;c < NoOfCols ; c++)
	       {
		        //50 is some random no. it will read the first 50 characters of the word 
			fgets(arr, 50, ptrRead);
			sscanf (arr, "%Lf", &VolFrac);
		      
		       	Write_struct_Cellattrib(r,c,VolFrac,-100000,-100000,-100000,-100000,"-100000",-100000,-100000);
	       }
      }
     fclose(ptrRead);
     //exit(0);
     /***********************************************************************/

    /******************INITIALISE GLOBAL VARIABLES**********/
    del_sq         = pow(del,2);
    del_sqBy2      = del_sq / 2.0;
    del2           = 2.0*del;
    del3           = 3.0*del;
    delBy2         = del/2.0;
    PiBy2          = PI/2.0;
    PiBy4          = PI/4.0;
    ThreePiBy2     = 3.0*PI/2.0;
    delTBydel      = delT/del;
    delTdel        = delT*del;  	  	
     /*******************************************************/

	/********************INTIALISE THE VELOCITY/FLUX FIELD**************/
	for(r=0;r < NoOfRows ; r++)
       {
	       for(c=0;c <= NoOfCols ; c++)
	       {
			//X_Velocity[r][c] = sin(c*del*PI/30.0)*cos(r*del*PI/30.0);//(-r*del*PI/30.0) + 1.5;
			X_Velocity[r][c] = 1;
		       //X_Velocity[r][c] = sin(c*del)*cos((r+0.5)*del);
		      // X_Velocity[r][c] = -0.5*((r+0.5)*del-2);
			if(fabs(X_Velocity[r][c]) < LOWER_LIMIT)
			{
				X_Velocity[r][c] = 0.0;
			}

			X_Flux[r][c]     = 0.0;		
	       }
       }

       for(r=0;r <= NoOfRows ; r++)
       {
	       for(c=0;c < NoOfCols ; c++)
	       {
			//Y_Velocity[r][c] = -cos(c*del*PI/30.0)*sin(r*del*PI/30.0);//(c*del*PI/30.0) - 1.5;
			Y_Velocity[r][c] =0;
		       //Y_Velocity[r][c] = -cos((c+0.5)*del)*sin(r*del);
		     //  Y_Velocity[r][c] = 0.5*((c+0.5)*del-2);  
			if(fabs(Y_Velocity[r][c]) < LOWER_LIMIT)
			{
				Y_Velocity[r][c] = 0.0;
			}

			Y_Flux[r][c]     = 0.0;
		}
	}
	/**************************************************************/
  }

/******************************************************************************************************************
METHOD NAME       : Calculate_InitialSlope()
DESCRIPTION       : CALCULATES INITIAL SLOPE FROM [4,5]
KNOWN LIMITATIONS : CANNOT WORK FOR THE BOUNDARY CELLS. THIS HAS TO BE CORRECTED.      
******************************************************************************************************************/
void Calculate_InitialSlope()
{
  /*LOCAL VARIABLES**********/
  int r, c;
  long double temp1, temp2, Nx, Ny;  
  /**************************/

  for(r=0;r < NoOfRows ; r++)
   {
      for(c=0;c < NoOfCols ; c++)
       {
		 if(Cellattrib[r][c].VolFrac > 0 && Cellattrib[r][c].VolFrac < 1)
		   {
				  temp1 = Cellattrib[r+1][c+1].VolFrac + (2.0*Cellattrib[r][c+1].VolFrac) + Cellattrib[r-1][c+1].VolFrac;
				  temp2 = Cellattrib[r+1][c-1].VolFrac + (2.0*Cellattrib[r][c-1].VolFrac) + Cellattrib[r-1][c-1].VolFrac;

				  Nx =  (temp1 - temp2)/del;     

				  if(Nx != 0.0)
					{
					  Nx = -Nx; //THIS IS REQUIRED BECAUSE REF[3] FROM WHERE THE FORMULA(Nx,Ny) IS TAKEN, HAS THE DEFN OPPOSITE
					}

				  temp1 = Cellattrib[r+1][c+1].VolFrac + (2.0*Cellattrib[r+1][c].VolFrac) + Cellattrib[r+1][c-1].VolFrac;
				  temp2 = Cellattrib[r-1][c+1].VolFrac + (2.0*Cellattrib[r-1][c].VolFrac) + Cellattrib[r-1][c-1].VolFrac;
				  Ny    = (temp1 - temp2)/del;    //CHANGE POINT
              
				  if(Ny != 0.0)
					{
					  Ny    = -Ny; //CHECK THIS
					}
		
				  Write_struct_Cellattrib(r,c,-100000,-100000,Nx,Ny,-100000,"-100000",-100000,-100000);
		   }
      }
   }    
}

/**********************************************************************************************************************************************
METHOD NAME: Identify_Shape
DESCRIPTION: Takes [r,c] as a cell identifier and identifies the shape of the area based on the Slope and the VolFraction
LIMITATIONS: Check all greater than equal to signs
***********************************************************************************************************************************************/
void Identify_Shape(int r, int c) 
{
  /*LOCAL VARIABLES*/
  long double Theta;
  long double Slope;
  long double VolFrac;
  long double Area;
  char   Shape[30];
  /*****************/
       VolFrac =  Cellattrib[r][c].VolFrac;
	   Area    =  VolFrac*del_sq;
       Theta   =  Cellattrib[r][c].Theta;

        if(Theta == 0 || Theta == PiBy2 || Theta == PI || Theta == ThreePiBy2)
         {
           strcpy(Shape,"Rectangle");
           Write_struct_Cellattrib(r,c,-100000,-100000,-100000,-100000,-100000,Shape,-100000,-100000);
           return;
         }
      
        Slope   =  tan(Theta) ;   //DONT MOVE THIS UP OTHERWISE SLOPE WILL GO TO INFINITY FOR THETA = PI / 2.0
       //SHADED AREA IS A TRIANGLE
       if( (Area <= (del_sqBy2*Slope) ) && (Area <= (del_sqBy2/Slope) )  )    //both for Theta <> PI/4.0     
        {
           strcpy(Shape,"Triangle");
        } 
      //SHADED AREA IS A TRAPEZIUM WITH THETA < PI / 4
      else if ( (Area > del_sqBy2*Slope )  && Area < (del_sq - (del_sqBy2*Slope) )  )  // for Theta < PI / 4.0
        {
          strcpy(Shape,"Trapezium");
        }
      //SHADED AREA IS A TRAPEZIUM WITH THETA > PI / 4
      else if ( (Area > del_sqBy2/Slope )  && Area < (del_sq - (del_sqBy2/Slope) )  )  // for Theta > PI / 4.0
        {
          strcpy(Shape,"Trapezium");
        }
      //SHADED AREA IS A 5 SIDED FIGURE WHICH IS THE COMPLEMENT OF A TRIANGLE
      else
        {          
          strcpy(Shape,"Triangle_Complement");
        }

       //WRITE SHAPE TO THE STRUCT
       Write_struct_Cellattrib(r,c,-100000,-100000,-100000,-100000,-100000,Shape,-100000,-100000);
}

/***********************************************************************************************************************************************************
METHOD NAME: Calculate_distance
DESCRIPTION: BASED ON THE SHAPE OF THE AREA, CALCULATES THE PERPENDICULAR DISTANCE. ASSUMES EVERYTHING TO BE IN THE 1ST QUADRANT
LIMITATIONS :NO CHECKS FOR BOUNDARY VALUES AS YET / NO CHECKS FOR NX, NY = 0   
             NO CHECKS FOR SEEING IF THE STRING SHAPE IS NULL IN STRUCTCELLATTRIB
***********************************************************************************************************************************************************/
void Calculate_distance(int r, int c)
{
  /*LOCAL VARIABLES*/
  long double Theta;
  long double VolFrac;
  char   Shape[30];
  long double P;
  long double P_dash;
  /****************/

  /*INITIALIZE*/
  Theta   = 0.0;
  VolFrac = 0.0;
  strcpy(Shape,"000");
  P       = 0.0;
  P_dash  = 0.0;
  /****************************/
 
  /*FETCH THE VALUES INTO THE LOCAL VARIABLES*/
  strcpy(Shape,Cellattrib[r][c].Shape);
  Theta   =  Cellattrib[r][c].Theta;
  VolFrac =  Cellattrib[r][c].VolFrac;
  /****************************************/
  
  if(!strcmp(Shape,"Rectangle"))
    {
      P = VolFrac*del;
    }
  else if( !strcmp(Shape,"Triangle") )    //both for Theta <> PI/4.0
   {
     P = sqrt ( VolFrac*del_sq*sin(2.0*Theta) );
   }
  else if ( !strcmp(Shape,"Trapezium") && Theta < PI/4.0)  // for Theta < PI / 4.0
   {
       //THIS IS OLDER CORRECT EXPRESSION - CHANGED TO MAKE IT FASTER
       //P =  ( ( 2*VolFrac*del_sq*cos(Theta) ) + del_sq*sin(Theta) ) / del2;
       P = (VolFrac*del_sq*cos(Theta)/del) + (del*sin(Theta)*0.5);
   }
  else if ( !strcmp(Shape,"Trapezium") && Theta > PI/4.0 )  // for Theta > PI / 4.0
   {   
    //THIS IS OLDER CORRECT EXPRESSION - CHANGED TO MAKE IT FASTER	   
     //P =  ( ( 2*VolFrac*del_sq*sin(Theta) ) + del_sq*cos(Theta) ) / del2;
	P =   (VolFrac*del_sq*sin(Theta)/del ) + (del*cos(Theta)*0.5); 
   }
  else if ( !strcmp(Shape,"Triangle_Complement") )    //both for Theta <> PI/4.0
   {      
     //THIS WORKS ONLY FOR A SQUARE CELL
     P_dash = sqrt( del_sq*(1.0 - VolFrac)*sin(2.0*Theta) );
     //THIS IS OLDER CORRECT EXPRESSION - CHANGED TO MAKE IT FASTER	   
     //P = (del/cos(Theta) ) - P_dash + ( del*sin(Theta)*( 1.0-tan(Theta) ) );
	P = del*( sin(Theta) +  cos(Theta) ) - P_dash;
   }    
  else  //ERROR HANDLING
    {
      printf("\n***********************************************************\n");
      printf("\nERROR - EXECUTION TERMINATES IN SUBROUTINE CALCULATE_DISTANCE\n");
      printf("\n***********************************************************\n");
      exit(0);
    }

    //WRITE PERPENDICULAR DISTANCE TO THE STRUCT
    Write_struct_Cellattrib(r,c,-100000,-100000,-100000,-100000,P,"-100000",-100000,-100000); 
}

/*************************************************************************************************************************************************
METHOD NAME: Extrapolate_Line
DESCRIPTION: EXTRAPOLATES THE LINE IN THE "CENTRAL" CELL TO THE 3X3 BLOCK CENTRED ON IT.
             ASSIGNS AREAS BASED ON CELL NAMES EG.NORTHWEST CELL'S COORDINATES ARE [r+1,c-1] FOR 1ST QUADRANT, [r-1,c-1] FOR 2ND QUADRANT ETC.
             ROTATES THE 3X3 BLOCK CLOCKWISE TO REDUCE EVERY QUADRANT TO THE FIRST QUADRANT.
LIMITATIONS:
**************************************************************************************************************************************************
 _______________________________________________________
|                 	|                    	|       	         | 
|                 	|                    	|               	 | 
|                 	|	             	|       	       	 | 
| NorthWest	|      North     	|    NorthEast   | 
|                 	|                    	|                        |   
|                 	|                    	|                        |
|_______________	|_______________ |________________|           
|                 	|                       |                |
|                 	|                       |                | 
|                 	|                       |                |  
|    West         	|      Central     |    East        |
|                 	|                       |                |   
|                 	|                       |                |
|_______________|_______________ |________________|         
|                 |                    |                |
|                 |                    |                | 
|                 |                    |                | 
|    SouthWest    |      South         |    SouthEast   |
|                 |                    |                |  
|                 |                    |                | 
|                 |                    |                |  
|_________________|____________________|________________|          
*/

void Extrapolate_Line(int r, int c)
{
  /****LOCAL VARIABLES****/
  int    Quad;
  long double y_0, y_del, y_2del, y_3del;   //the value of y at x = 0 , x = del , x = 2*del and x = 3*del
  long double x_0, x_del, x_2del, x_3del ;  //the value of x at y = 0 , y = del , y = 2*del and y = 3*del
  long double Area;
  long double Theta; 
  long double Slope;      
  long double P;
  char   Shape[30];
  /**********************/

  strcpy(Shape,Cellattrib[r][c].Shape);
  Quad        =  Find_Quad(r,c);
  P           =  Cellattrib[r][c].P;
 
  IdentifyCell_AssignArea(r,c,Quad,"Central",Cellattrib[r][c].VolFrac);   //for the central cell, assign the same VolFrac

  //if rectangle, 
  if( !strcmp(Shape,"Rectangle") )
    {
      Area = P*del;
      if(Quad == 1)
        {
            Write_struct_Cellattrib(r+1,c+1,-100000,0.0,-100000,-100000,-100000,"-100000",-100000,-100000);
            Write_struct_Cellattrib(r  ,c+1,-100000,0.0,-100000,-100000,-100000,"-100000",-100000,-100000);
            Write_struct_Cellattrib(r-1,c+1,-100000,0.0,-100000,-100000,-100000,"-100000",-100000,-100000);
            
            Write_struct_Cellattrib(r+1,c,-100000,Area,-100000,-100000,-100000,"-100000",-100000,-100000);
            Write_struct_Cellattrib(r-1,c,-100000,Area,-100000,-100000,-100000,"-100000",-100000,-100000);

            Write_struct_Cellattrib(r+1,c-1,-100000,1.0,-100000,-100000,-100000,"-100000",-100000,-100000);
            Write_struct_Cellattrib(r  ,c-1,-100000,1.0,-100000,-100000,-100000,"-100000",-100000,-100000);
            Write_struct_Cellattrib(r-1,c-1,-100000,1.0,-100000,-100000,-100000,"-100000",-100000,-100000);
        }
      else if(Quad == 2)
        {
            Write_struct_Cellattrib(r+1,c-1,-100000,0.0,-100000,-100000,-100000,"-100000",-100000,-100000);
            Write_struct_Cellattrib(r+1,c  ,-100000,0.0,-100000,-100000,-100000,"-100000",-100000,-100000);
            Write_struct_Cellattrib(r+1,c+1,-100000,0.0,-100000,-100000,-100000,"-100000",-100000,-100000);
            
            Write_struct_Cellattrib(r,c-1,-100000,Area,-100000,-100000,-100000,"-100000",-100000,-100000);
            Write_struct_Cellattrib(r,c+1,-100000,Area,-100000,-100000,-100000,"-100000",-100000,-100000);

            Write_struct_Cellattrib(r-1,c-1,-100000,1.0,-100000,-100000,-100000,"-100000",-100000,-100000);
            Write_struct_Cellattrib(r-1,c  ,-100000,1.0,-100000,-100000,-100000,"-100000",-100000,-100000);
            Write_struct_Cellattrib(r-1,c+1,-100000,1.0,-100000,-100000,-100000,"-100000",-100000,-100000);
        }
      else if(Quad == 3)
        {
            Write_struct_Cellattrib(r+1,c-1,-100000,0.0,-100000,-100000,-100000,"-100000",-100000,-100000);
            Write_struct_Cellattrib(r  ,c-1,-100000,0.0,-100000,-100000,-100000,"-100000",-100000,-100000);
            Write_struct_Cellattrib(r-1,c-1,-100000,0.0,-100000,-100000,-100000,"-100000",-100000,-100000);
            
            Write_struct_Cellattrib(r+1,c,-100000,Area,-100000,-100000,-100000,"-100000",-100000,-100000);
            Write_struct_Cellattrib(r-1,c,-100000,Area,-100000,-100000,-100000,"-100000",-100000,-100000);

            Write_struct_Cellattrib(r+1,c+1,-100000,1.0,-100000,-100000,-100000,"-100000",-100000,-100000);
            Write_struct_Cellattrib(r  ,c+1,-100000,1.0,-100000,-100000,-100000,"-100000",-100000,-100000);
            Write_struct_Cellattrib(r-1,c+1,-100000,1.0,-100000,-100000,-100000,"-100000",-100000,-100000);
        }
      else if(Quad == 4)
        {
            Write_struct_Cellattrib(r-1,c-1,-100000,0.0,-100000,-100000,-100000,"-100000",-100000,-100000);
            Write_struct_Cellattrib(r-1,c  ,-100000,0.0,-100000,-100000,-100000,"-100000",-100000,-100000);
            Write_struct_Cellattrib(r-1,c+1,-100000,0.0,-100000,-100000,-100000,"-100000",-100000,-100000);
            
            Write_struct_Cellattrib(r,c-1,-100000,Area,-100000,-100000,-100000,"-100000",-100000,-100000);
            Write_struct_Cellattrib(r,c+1,-100000,Area,-100000,-100000,-100000,"-100000",-100000,-100000);

            Write_struct_Cellattrib(r+1,c-1,-100000,1.0,-100000,-100000,-100000,"-100000",-100000,-100000);
            Write_struct_Cellattrib(r+1,c  ,-100000,1.0,-100000,-100000,-100000,"-100000",-100000,-100000);
            Write_struct_Cellattrib(r+1,c+1,-100000,1.0,-100000,-100000,-100000,"-100000",-100000,-100000);
        } 
      return;
    }
 
         Theta       =  Cellattrib[r][c].Theta;
         Slope       =  tan(Theta);   
         /// ALL CONDITIONS ARE CHECKED USING THE FOLLOWING EQUATION
         /// y = tan(theta)[h + p*cosec(theta) - x] + h 
         /// P IS MEASURED FROM THE BOTTOM LHS CORNER

         y_0   =  Slope*( del + (P/sin(Theta))           ) + del;
         y_del =  Slope*( del + (P/sin(Theta)) - del     ) + del;
         y_2del = Slope*(       (P/sin(Theta)) - del     ) + del;
         y_3del = Slope*(       (P/sin(Theta)) - del2    ) + del;

         x_0    = del + (P/sin(Theta)) + (del/Slope);
         x_del  = del + (P/sin(Theta))              ;
         x_2del = del + (P/sin(Theta)) - (del/Slope);
         x_3del = del + (P/sin(Theta)) - (del2/Slope);

      //SHADED AREA IS A TRIANGLE

         if( !strcmp(Shape,"Triangle") )    //both for Theta <> PI/4.0     
          {
                //ASSIGN VOLFRAC TO CELLS WHICH ARE 0 OR 1
                IdentifyCell_AssignArea(r,c,Quad,"SouthWest",1.0);
                IdentifyCell_AssignArea(r,c,Quad,"North",0.0);
                IdentifyCell_AssignArea(r,c,Quad,"NorthEast",0.0);
                IdentifyCell_AssignArea(r,c,Quad,"East",0.0);
                /******************************************/

                /*****************************CALCULATE AREAS FOR NORTHWEST & WEST CELL*****************************/
                 if (y_0 <= del2)  //IF NORTHWEST IS 0
                   {
                     IdentifyCell_AssignArea(r,c,Quad,"NorthWest",0.0);

                     Area = ( (y_0 - del) + (y_del - del) )*delBy2;          //Calculate area of trapezium for West
                     IdentifyCell_AssignArea(r,c,Quad,"West",Area);
                   }
                 else if ( y_0 > del2 && y_0 <= del3 )                  //IF NORTHWEST IS A TRIANGLE
                   {
                     Area = 0.5*(y_0 - del2)*x_2del;                       //Calculate area of the triangle for NorthWest
                     IdentifyCell_AssignArea(r,c,Quad,"NorthWest",Area);

                     Area = del_sq - ( 0.5*(del - x_2del)*(del2 - y_del) ); //Calculate area of triangle and take compliment.for NorthWest
                     IdentifyCell_AssignArea(r,c,Quad,"West",Area);
                   }
                 else if( y_0 > del3 )   //IF NORTHWEST IS A TRAPEZIUM
                   {
                     Area = ( x_2del + x_3del )*delBy2;                       // Calculate the area of the trapezium for NorthWest
                     IdentifyCell_AssignArea(r,c,Quad,"NorthWest",Area);

                     Area = del_sq - ( 0.5*(del - x_2del)*(del2 - y_del) ); //Calculate area of triangle and take compliment for West
                     IdentifyCell_AssignArea(r,c,Quad,"West",Area);
                   }
                 /***********************END*****************************************************************************/

                 /********CALCULATE AREAS FOR South AND SOUTHEAST CELLS (exchange x and y in all expressions above)***********/
                  if( x_0 <= del2 )    //IF SOUTHEAST IS 0
                  {
                    IdentifyCell_AssignArea(r,c,Quad,"SouthEast",0.0);

                    Area = ( (x_0 - del) + (x_del - del) )*delBy2;       //calculate area of trapezium for South
                    IdentifyCell_AssignArea(r,c,Quad,"South",Area);
                  }
                 else if ( x_0 > del2 && x_0 <= del3 )                  //IF SOUTHEAST IS A TRIANGLE
                   {
                     Area = 0.5*(x_0 - del2)*y_2del;                       //Calculate area of the triangle for SouthEast
                     IdentifyCell_AssignArea(r,c,Quad,"SouthEast",Area);

                     Area = del_sq - ( 0.5*(del - y_2del)*(del2 - x_del) ); //Calculate area of triangle and take complement for South
                     IdentifyCell_AssignArea(r,c,Quad,"South",Area);
                   }
                 else if( x_0 > del3 )      //IF SOUTHEAST IS A TRAPEZIUM
                   {
                     Area = (y_2del + y_3del)*delBy2;                       // Calculate the area of the trapezium for SouthEast
                     IdentifyCell_AssignArea(r,c,Quad,"SouthEast",Area);

                     Area = del_sq - ( 0.5*(del - y_2del)*(del2 - x_del) ); //Calculate area of triangle and take complement for South
                     IdentifyCell_AssignArea(r,c,Quad,"South",Area);
                   }
                 /***********************END*********************************************************************************/
        }

      //SHADED AREA IS A TRAPEZIUM WITH THETA < PI / 4
         else if (  !strcmp(Shape,"Trapezium") && Theta < PI/4.0)  // for Theta < PI / 4.0
        {
                //Assign VolFRac to cells which are 0 or 1
                IdentifyCell_AssignArea(r,c,Quad,"North",0.0);
                IdentifyCell_AssignArea(r,c,Quad,"NorthEast",0.0);
                IdentifyCell_AssignArea(r,c,Quad,"SouthWest",1.0);
                IdentifyCell_AssignArea(r,c,Quad,"South",1.0);
                /******************************************/
                                
                /*****************************CALCULATE AREAS FOR NORTHWEST AND WEST CELLS*****************************/
                 if (y_0 <= del2)  //if NorthWest is 0
                   {
                     IdentifyCell_AssignArea(r,c,Quad,"NorthWest",0.0);

                     Area = ( (y_0 - del) + (y_del - del) )*delBy2;          //Calculate area of trapezium for West
                     IdentifyCell_AssignArea(r,c,Quad,"West",Area);                  
                   }
                 else if ( y_0 > del2 && y_0 <= del3 )                  //if NorthWest is a triangle
                   {
                     Area = 0.5*(y_0 - del2)*x_2del;                       //Calculate area of the triangle for NorthWest
                     IdentifyCell_AssignArea(r,c,Quad,"NorthWest",Area);                     

                     Area = del_sq - ( 0.5*(del - x_2del)*(del2 - y_del) ); //Calculate area of triangle and take compliment for West
                     IdentifyCell_AssignArea(r,c,Quad,"West",Area);                  
                   }            
                 /***********************END*****************************************************************************/

                 /*****************************CALCULATE AREAS FOR SOUTHEAST AND EAST CELLS*****************************/
                if(y_3del >= del)    //if SouthEast is 1
                  {
                     IdentifyCell_AssignArea(r,c,Quad,"SouthEast",1.0);
                     
                     Area = ( (y_2del - del) + (y_3del - del) )*delBy2;   //calculate area of the trapezium for East
                     IdentifyCell_AssignArea(r,c,Quad,"East",Area);
                  }
                else if(y_3del > 0 && y_3del < del)   //if SouthEast is a complement of a triangle
                  {
                    Area = del_sq  - ( 0.5*(del3 - x_del)*(del-y_3del) ); //calculate area of the triangle and then take complement for SouthEast
                    IdentifyCell_AssignArea(r,c,Quad,"SouthEast",Area);

                    Area = 0.5*(y_2del - del)*(x_del - del2);     //calculate area of area of the trinagle for East
                    IdentifyCell_AssignArea(r,c,Quad,"East",Area);
                  }
                 /***********************END****************************************************************************/
        }

      //SHADED AREA IS A TRAPEZIUM WITH THETA > PI / 4
         else if ( !strcmp(Shape,"Trapezium")  && Theta > PI/4.0)  // for Theta > PI / 4.0
         {
                //Assign VolFRac to cells which are 0 or 1
                IdentifyCell_AssignArea(r,c,Quad,"West",1.0);
                IdentifyCell_AssignArea(r,c,Quad,"SouthWest",1.0);
                IdentifyCell_AssignArea(r,c,Quad,"East",0.0);
                IdentifyCell_AssignArea(r,c,Quad,"NorthEast",0.0);
                /******************************************/
                
                /*****************************CALCULATE AREAS FOR NORTH & NORTHWEST*****************************/
                if(y_del >= del3)   //if NorthWest is 1
                  {
                    IdentifyCell_AssignArea(r,c,Quad,"NorthWest",1.0);
                    //printf("\nNW is 1\n");

                    Area = ( (x_2del - del) + (x_3del - del) )*delBy2;   //calculate area of the trapezium for North
                    IdentifyCell_AssignArea(r,c,Quad,"North",Area);
                    // printf("\nN %f \n",Area);
                  }
                else if(y_del < del3) //if NorthWest is complement of a triangle
                  {
                    Area = del_sq - ( 0.5*(del3 - y_del)*(del - x_3del) );   //calculate area of the triangle and take complement for NorthWest 
                    IdentifyCell_AssignArea(r,c,Quad,"NorthWest",Area);

                    Area = 0.5*(y_del - del2)*(x_2del - del);          //calculate area of the triangle for North
                    IdentifyCell_AssignArea(r,c,Quad,"North",Area);                 
                  }             
                /***********************END**************************************************************************/

                /*****************************CALCULATE AREAS FOR SOUTH & SOUTHEAST**********************************/
                if(x_0 <= del2)   //if SouthEast is 0
                  {
                    //printf("\nSE is 0\n");
                    IdentifyCell_AssignArea(r,c,Quad,"SouthEast",0.0);             

                    Area = ( (x_del - del) + (x_0 - del) )*delBy2;   //calculate area of the trapezium for South
                    IdentifyCell_AssignArea(r,c,Quad,"South",Area);                 
                    //printf("\nN %f \n",Area);
                  }
                else if(x_0 >  del2)    //if SouthEast is a triangle
                  {
                    Area = 0.5*y_2del*(x_0 - del2);    //calculate area of the triangle for SouthEast
                    IdentifyCell_AssignArea(r,c,Quad,"SouthEast",Area);            

                    Area = del_sq - (del - y_2del)*(del2 - x_del);      //calculate area of the triangle and take complement for South
                    IdentifyCell_AssignArea(r,c,Quad,"South",Area);                 
                  }             
                /***********************END**************************************************************************/
        }

      //SHADED AREA IS A 5 SIDED FIGURE WHICH IS THE COMPLEMENT OF A TRIANGLE
         else if( !strcmp(Shape,"Triangle_Complement")  )  //both for Theta <> PI/4
        {
            //Assign VolFRac to cells which are 0 or 1
            IdentifyCell_AssignArea(r,c,Quad,"West",1.0);
            IdentifyCell_AssignArea(r,c,Quad,"SouthWest",1.0);
            IdentifyCell_AssignArea(r,c,Quad,"South",1.0);
            IdentifyCell_AssignArea(r,c,Quad,"NorthEast",0.0);
            /******************************************/

          /*****************************CALCULATE AREAS FOR NORTH & NORTHWEST**********************************/
          if(x_3del >= del)   //if NorthWest is 1
           {
             //printf("NW is 1\n");
              IdentifyCell_AssignArea(r,c,Quad,"NorthWest",1.0);

              Area = ( (x_2del - del) + (x_3del - del) )*delBy2;   //calculate area of the trapezium for North
              IdentifyCell_AssignArea(r,c,Quad,"North",Area);
              //printf("\nN is %f",Area);
           }
          else if(x_3del < del  && x_3del > 0) //if NorthWest is complement of a triangle
           {
            Area = del_sq - ( 0.5*(del3 - y_del)*(del - x_3del) );   //calculate area of the triangle and take complement for NorthWest 
            IdentifyCell_AssignArea(r,c,Quad,"NorthWest",Area);
            //printf("NW is %f",Area);

            Area = 0.5*(y_del - del2)*(x_2del - del);          //calculate area of the triangle for North
            IdentifyCell_AssignArea(r,c,Quad,"North",Area);                 
           }            
          else if(x_3del  < 0)    //if NorthWest is a trapezium
           {
             Area = ( (y_0 - del2) + (y_del - del2) )*delBy2;    //calculate area of the trapezium for NorthWest
             IdentifyCell_AssignArea(r,c,Quad,"NorthWest",Area);

             Area =0.5*(y_del - del2)*(x_2del - del);   //calculate area of the triangle for North
             IdentifyCell_AssignArea(r,c,Quad,"North",Area);
           }
          /****************************************END**********************************************************/

          /*****************************CALCULATE AREAS FOR EAST & SOUTHEAST**********************************/
           if(y_3del >= del)       //if SouthEast is 1
            {
              IdentifyCell_AssignArea(r,c,Quad,"SouthEast",1.0);

              Area = ( (y_2del - del) + (y_3del - del) )*delBy2;   //calculate area of the trapezium for East
              IdentifyCell_AssignArea(r,c,Quad,"East",Area);
            }
          else if(y_3del < del && y_3del > 0)   //if SouthEast is complement of a triangle
            {
              Area = del_sq - ( 0.5*(del3 - x_del)*(del - y_3del) );   //calculate area of the triangle and take complement for SouthEast 
              IdentifyCell_AssignArea(r,c,Quad,"SouthEast",Area);

              Area = 0.5*(x_del - del2)*(y_2del - del);          //calculate area of the triangle for East
              IdentifyCell_AssignArea(r,c,Quad,"East",Area);
            }
          else if(y_3del < 0)    //if SouthEast is a trapezium
            {         
              Area = ( (x_0 - del2) + (x_del - del2) )*delBy2;    //calculate area of the trapezium for SouthEast
              IdentifyCell_AssignArea(r,c,Quad,"SouthEast",Area);

              Area =0.5*(x_del - del2)*(y_2del - del);   //calculate area of the triangle for East
              IdentifyCell_AssignArea(r,c,Quad,"East",Area);
            }
          /****************************************END**********************************************************/
        }
	 else //ERROR HANDLING
        {
	  printf("\n***********************************************************\n");
	  printf("\nERROR - ERROR IN SUBROUTINE EXTRAPOLATE_LINE\n");
	  printf("\n***********************************************************\n");
	  //exit(0);
        }
}
/************************************************************************************************************************************************************
METHOD NAME: IdentifyCell_AssignArea
DESCRIPTION: TAKES THE CELL NAME EG. NORTHWEST AND LOCATES THE POSITION OF THE CELL, BASED ON THE QUADRANT. ASSIGNS THE PASSED VALUE OF THE AREA TO THE CELL.
LIMITATIONS: NONE
*************************************************************************************************************************************************************/

void IdentifyCell_AssignArea(int r, int c,int Quad, char CellName[30],long double Area)
{
  int R, C;
  if(!strcmp(CellName,"NorthWest") )
    {
       switch (Quad)
        {
               case 1:
                 R = r+1; 
                 C = c-1;
               break;

               case 2:
                 R = r-1; 
                 C = c-1;
               break;
               
               case 3:
                 R = r-1; 
                 C = c+1;
               break;
               
               case 4:
                 R = r+1; 
                 C = c+1;
               break;
        }      
    }
  else if(!strcmp(CellName,"North") )
    {
      switch (Quad)
        {
               case 1:
                 R = r+1; 
                 C = c;
               break;

               case 2:
                 R = r; 
                 C = c-1;
               break;
               
               case 3:
                 R = r-1; 
                 C = c;
               break;
               
               case 4:
				 R =r; 
				 C =c+1;
			   break;
		}      
    }
  else if(!strcmp(CellName,"NorthEast") )
    {
        switch (Quad)
		{
			   case 1:
				 R = r+1; 
				 C = c+1;
			   break;

       		   case 2:
				R = r+1; 
				C = c-1;
			   break;
			   
       		   case 3:
				 R = r-1; 
				 C = c-1;
			   break;
			   
       		   case 4:
				 R = r-1; 
				 C = c+1;
			   break;
		} 
    }
  else if( !strcmp(CellName,"West") )
    {
		switch (Quad)
		{
			   case 1:
				R = r; 
				C = c-1;
			   break;

       		   case 2:
				R = r-1; 
				C = c;
			   break;
			   
       		   case 3:
			    R = r; 
				C = c+1;
			   break;
			   
       		   case 4:
				R = r+1; 
				C = c;
			   break;
		} 
    }   
  else if(!strcmp(CellName,"East") )
    {
		switch (Quad)
		{
			   case 1:
				R = r; 
				C = c+1;
			   break;

       		   case 2:
				 R = r+1; 
				 C = c;
			   break;
			   
       		   case 3:
				 R = r; 
				 C = c-1;
			   break;
			   
       		   case 4:
				 R = r-1; 
				 C = c;
			   break;
		} 
    }
  else if(!strcmp(CellName,"SouthWest") )
    {
		switch (Quad)
		{
			   case 1:
				R = r-1; 
				C = c-1;
			   break;

       		   case 2:
				 R = r-1; 
				 C = c+1;
			   break;
			   
       		   case 3:
				 R = r+1; 
				 C = c+1;
			   break;
			   
       		   case 4:
				 R = r+1; 
				 C = c-1;
			   break;
		} 
    }  
  else if(!strcmp(CellName,"South") )
    {
		switch (Quad)
		{
			   case 1:
				R = r-1; 
				C = c;
			   break;

       		   case 2:
				 R = r; 
				 C = c+1;
			   break;
			   
       		   case 3:
				 R = r+1; 
				 C = c;
			   break;
			   
       		   case 4:
				 R = r; 
				 C = c-1;
			   break;
		} 
    }  
  else if(!strcmp(CellName,"SouthEast") )
    {
		switch (Quad)
		{
			   case 1:
				R = r-1; 
				C = c+1;
			   break;

       		   case 2:
				 R = r+1; 
				 C = c+1;
			   break;
			   
       		   case 3:
				 R = r+1; 
				 C = c-1;
			   break;
			   
       		   case 4:
				 R = r-1; 
				 C = c-1;
			   break;
		} 
    }
  else if(!strcmp(CellName,"Central") )
    {
      R = r;
      C = c;
    }
  else //ERROR HANDLING
    {
      printf("\n\n***********************************************************\n\n");
	  printf("\nERROR IN SUBROUTINE IDENTIFYCELL_ASSIGNAREA\n");
	  printf("\n\n***********************************************************\n\n");
	  exit(0);
    }

  Write_struct_Cellattrib(R,C,-100000,Area,-100000,-100000,-100000,"-100000",-100000,-100000);
}

/*************************************************************************************************************************
NAME       : Calculate_SumOfDiffSq
DESCRIPTION: CALCULATES THE SUM OF SQUARE OF THE DIFFERENCE BETWEEN THE ACTUAL AREA AND CALCULATED AREA FOR THE 3X3 BLOCK
LIMITATIONS: NONE
**************************************************************************************************************************/
long double Calculate_SumOfDiffSq(int r, int c) 
{
  long double North_Area,South_Area,West_Area,East_Area,NorthWest_Area,NorthEast_Area,SouthWest_Area,SouthEast_Area;
  long double North_VolFrac,South_VolFrac,West_VolFrac,East_VolFrac,NorthWest_VolFrac,NorthEast_VolFrac,SouthWest_VolFrac,SouthEast_VolFrac;
  long double Square1, Square2, Square3, Square4, Square5, Square6, Square7, Square8; 
  long double Sum_Squares;

  North_Area    = Cellattrib[r+1][c].Area;
  North_VolFrac = Cellattrib[r+1][c].VolFrac*del_sq;

  South_Area    = Cellattrib[r-1][c].Area;
  South_VolFrac = Cellattrib[r-1][c].VolFrac*del_sq;

  West_Area     = Cellattrib[r][c-1].Area;
  West_VolFrac  = Cellattrib[r][c-1].VolFrac*del_sq;
    
  East_Area     = Cellattrib[r][c+1].Area;
  East_VolFrac  = Cellattrib[r][c+1].VolFrac*del_sq;

  NorthWest_Area     = Cellattrib[r+1][c-1].Area;
  NorthWest_VolFrac  = Cellattrib[r+1][c-1].VolFrac*del_sq;

  NorthEast_Area     = Cellattrib[r+1][c+1].Area;
  NorthEast_VolFrac  = Cellattrib[r+1][c+1].VolFrac*del_sq;

  SouthWest_Area     = Cellattrib[r-1][c-1].Area;
  SouthWest_VolFrac  = Cellattrib[r-1][c-1].VolFrac*del_sq;

  SouthEast_Area     = Cellattrib[r-1][c+1].Area; 
  SouthEast_VolFrac  = Cellattrib[r-1][c+1].VolFrac*del_sq;  

  Square1 = pow( North_Area      - North_VolFrac,2);
  Square2 = pow( South_Area      - South_VolFrac,2);
  Square3 = pow( West_Area       - West_VolFrac,2);
  Square4 = pow( East_Area       - East_VolFrac,2);
  Square5 = pow( NorthWest_Area  - NorthWest_VolFrac,2);
  Square6 = pow( SouthWest_Area  - SouthWest_VolFrac,2);
  Square7 = pow( NorthEast_Area  - NorthEast_VolFrac,2);
  Square8 = pow( SouthEast_Area  - SouthEast_VolFrac,2);

  Sum_Squares = Square1 + Square2 + Square3 + Square4 + Square5 + Square6 + Square7 + Square8;

  return Sum_Squares;
}

/*****************************************************************************************************************
NAME       : Find_Quad
DESCRIPTIOM: TAKES r AND c AS INPUT AND RETURNS THE QUADRANT BASED ON THE SIGNS OF Nx AND Ny
LIMITATIONS: None
*****************************************************************************************************************/
int Find_Quad(int r, int c)
{
    long double Nx, Ny;
	int Quad;

    Nx =  Cellattrib[r][c].Nx;
    Ny =  Cellattrib[r][c].Ny;

	
    if (Nx > 0 && Ny >= 0)   //FIRST QUADRANT
      {
	Quad = 1;
      }
    else if (Nx <= 0 && Ny > 0) //SECOND QUADRANT
      {
	Quad = 2;
      }
    else if (Nx < 0 && Ny <= 0) //THIRD QUADRANT
      {
	Quad = 3;
      }
    else if (Nx >= 0 && Ny < 0) //FOURTH QUADRANT
      {
	Quad = 4;
      }
     
	else //ERROR HANDLING
    {
      printf("\n\n***********************************************************\n\n");
	  printf("\nERROR IN SUBROUTINE Find_Quad\n");
	  printf("\n\n***********************************************************\n\n");
	  printf("%d\t%d\t%Lf\t%Lf",r,c,Cellattrib[r][c].Nx,Cellattrib[r][c].Ny);
	    exit(0);
    }

    return Quad;
}

/*******************************************************************************************************************
NAME       : Calculate_Theta
DESCRIPTIOM: RETURNS THE ANGLE MADE BY THE LINE WITH THE NEGATIVE DIRECTION OF THE X-AXIS, IN AN ANTICLOCKWISE SENSE
LIMITATIONS: 
********************************************************************************************************************/
void Calculate_Theta(int r, int c)
{
  /*LOCAL VARIABLES*/
  long double Nx, Ny, Theta,Quad;
  /*****************/
  Nx =  Cellattrib[r][c].Nx;
  Ny =  Cellattrib[r][c].Ny;
  Quad = Find_Quad(r,c);

  //FLAG --  IS THERE A BUG? FIND IT (corrected the theta values it looks like they are the angle for normal, but they should be the angle with interface with horizontal axis - Palas )
  if( Ny == 0.0 && Nx > 0.0 )
    {
      //Theta = 0.0; //this should be pi/2 because Theta is the  positive acute angle of the INTERFACE(not normal) made with the horizontal axis
	    Theta = PI/2;  
      Write_struct_Cellattrib(r,c,-100000,-100000,-100000,-100000,-100000,"-100000",Theta,-100000);
      return;
    }
  else if ( Nx == 0.0 && Ny > 0.0 )
    {
      Theta = PI/2.0;		// this is the second quad if you rotate the normal to quad 1, the interface will make pi/2 with x-axis
	   
      Write_struct_Cellattrib(r,c,-100000,-100000,-100000,-100000,-100000,"-100000",Theta, -100000);
      return;
    }
  else if( Ny == 0.0 && Nx < 0.0 )
    {
      //Theta = PI;   //this should be zero because Theta is the  positive acute angle of the INTERFACE(not normal) made with the horizontal axis
	  Theta = PI/2;   // third quadrant, if you rotate the normal to quad 1 angle will remain the same
      Write_struct_Cellattrib(r,c,-100000,-100000,-100000,-100000,-100000,"-100000",Theta, -100000);
      return;
    }
  else if( (Nx == 0.0 && Ny < 0.0) )
    {
      //Theta = 3.0*PiBy2;		//
	    Theta = PI/2;		//this is the fourth quad therefore it is pi/2, if you rotate the normal to 1
      Write_struct_Cellattrib(r,c,-100000,-100000,-100000,-100000,-100000,"-100000",Theta, -100000);
      return;
    }
    else 
    {
	    Theta  = (PiBy2) - fabs(atan(Ny/Nx));   
	    if(Quad == 2 || Quad == 4)
	    {
	      Theta = PiBy2 - Theta;
	    }
    }
	    
   /****/
     // Flag : I COMMENTED these because it will interfere with the rectangle cases for theta calculation, I took the non-rectangle cases in ELSE block above, the correction for QUADS 2nd and 4th are taken due consideration in rectangles case, 
    // remember that theta has to be calculated in a way that if we rotate the cell to make the normal in the first quadrant, then what angle the interface makes with the x-axis that too positive acute angle.  --------Palas
    
    
  /*fabs is required only for the 2nd and 4th Quads. atan(Nx/Ny) is positive for the 1st and 3rd Quads
  and negative for 2nd and 4th quadrants. Theta thus calculated is always a positive acute angle*/
  //~ //Theta  = (PiBy2) - fabs(atan(Ny/Nx));   // Flag : I commented this 

  //~ //FLAG - THIS IS NOT CLEAR - WHY IS THIS NECESSARY?????????????????  
  //~ if(Quad == 2 || Quad == 4)
    //~ {
      //~ Theta = PiBy2 - Theta;
    //~ }
  //~ /********/    

  Write_struct_Cellattrib(r,c,-100000,-100000,-100000,-100000,-100000,"-100000",Theta, -100000);
}

/*****************************************************************************************************************
NAME       : Rotate_Line
DESCRIPTIOM: ROTATES THE LINE CLOCKWISE / COUNTERCLOCKWISE
LIMITATIONS:
*****************************************************************************************************************/
void Rotate_Line(int r, int c, char sense[20])
{
  long double Nx, Ny;
  int Quad;

  Nx   =  Cellattrib[r][c].Nx;
  Ny   =  Cellattrib[r][c].Ny; 
  Quad =  Find_Quad(r,c);

  if(!strcmp(sense,"Clockwise"))
    {
		 if(Quad == 1)
		 {
		  Nx = Nx + Rotation_StepSize;
		  Ny = Ny - Rotation_StepSize;
		 }
		 else if(Quad == 2 )
		 {
		  Nx = Nx + Rotation_StepSize;
		  Ny = Ny + Rotation_StepSize;
		 }
		 else if(Quad == 3)
		 {
		  Nx = Nx - Rotation_StepSize;
		  Ny = Ny + Rotation_StepSize;
		 }
		 else if(Quad == 4)
		 {
		  Nx = Nx - Rotation_StepSize;
		  Ny = Ny - Rotation_StepSize;
		 }
    }
  else if(!strcmp(sense,"Counterclockwise"))
    {
		 if(Quad == 1)
		 {
		  Nx = Nx - Rotation_StepSize;
		  Ny = Ny + Rotation_StepSize;
		 }
		 else if(Quad == 2)
		 {
		  Nx = Nx - Rotation_StepSize;
		  Ny = Ny - Rotation_StepSize;
		 }
		 else if(Quad == 3)
		 {
  		  Nx = Nx + Rotation_StepSize;
		  Ny = Ny - Rotation_StepSize;
		 }
		 else if(Quad == 4)
		 {
  		  Nx = Nx + Rotation_StepSize;
		  Ny = Ny + Rotation_StepSize;
		 }
    }  
 Write_struct_Cellattrib(r,c,-100000,-100000,Nx,Ny,-100000,"-100000",-100000,-100000);  
}
/*******END************/

/*****************************************************************************************************************
NAME       : WriteToFile(int)
DESCRIPTIOM: Writes the (x1,y1) and (x2,y2) values for each interfacial cell into a file
LIMITATIONS:
*****************************************************************************************************************/
void WriteToFile(int Step)
{
  /*DECLARE LOCAL VARIABLES*/
  int r,c;
  int Quad;
  long double Theta; 
  long double P;
  char   Shape[30];
  long double X1, Y1, X2, Y2;
  FILE *ptrWrite;
  char LineFormat[] = "./Plots/PlotLine_%d.dat";
  char LineFileName[20];
  sprintf(LineFileName,LineFormat,Step);
  /*********************/

  /***********INITIALISE**************/ 
  ptrWrite =fopen(LineFileName,"w");
  /***********************************/

  for(r=0;r < NoOfRows ; r++)
   {
     for(c=0;c < NoOfCols ; c++)
       {		if(Cellattrib[r][c].VolFrac < LOWER_LIMIT && Cellattrib[r][c].VolFrac > UPPER_LIMIT)
	       { printf("%d\t%d\t%Lf\n",r,c,Cellattrib[r][c].VolFrac);}
		 if(Cellattrib[r][c].VolFrac > LOWER_LIMIT && Cellattrib[r][c].VolFrac < UPPER_LIMIT)
		   {
			 strcpy(Shape,Cellattrib[r][c].Shape);
			 Quad        =  Find_Quad(r,c);
			 P           =  Cellattrib[r][c].P;
			 Theta       =  Cellattrib[r][c].Theta;

			 if( !strcmp(Shape,"Rectangle") )
			   {
				 if(Quad == 1 || Quad == 3)
				   {
					  X1 = P + (c*del);
					  Y1 = r*del;

					  X2 = P + (c*del);
					  Y2 = (r+1)*del;
				   }
				 else if (Quad == 2 || Quad == 4)
				   {
					  X1 = (c*del);
					  Y1 = P + (r*del);

					  X2 = (c+1)*del;
					  Y2 = P + (r*del);
				   }				 
			   }
			 else if( !strcmp(Shape,"Triangle") )
			 {
				 if(Quad == 1)
				   {
					 X1 = (P/sin(Theta)) + (c*del);
					 Y1 =  r*del;

					 X2 = c*del;
					 Y2 = (P/cos(Theta)) + (r*del);
				   }
				 else if (Quad == 2)
				   {
					 X1 = del - (P/cos(Theta)) + (c*del);
					 Y1 = r*del;
					 
					 X2 = (c+1)*del;
					 Y2 = (P/sin(Theta)) + (r*del);
				   }
				 else if (Quad == 3)
				   {
					 X1 = (c+1)*del - (P/sin(Theta));
					 Y1 = (r+1)*del;

					 X2 = (c+1)*del;
					 Y2 = (r+1)*del - (P/cos(Theta));
				   }
				 else if(Quad == 4)
				   {
					 X1 = c*del;
					 Y1 = (r+1)*del - (P/sin(Theta));

					 X2 = c*del + (P/cos(Theta));
					 Y2 = (r+1)*del;
				   }
				}
				else if (  !strcmp(Shape,"Trapezium") && Theta < PiBy4)  // for Theta < PI / 4.0
				{
					 if(Quad == 1)
					   {
						 X1 = c*del;
						 Y1 = r*del + (P/cos(Theta));
						 
						 X2 = (c+1)*del;
						 Y2 = (P/cos(Theta)) - (del*tan(Theta)) + (r*del);
					   }
					 else if (Quad == 2)
					   {
						 X1 = (c+1)*del - (P/cos(Theta));
						 Y1 =  r*del;
						 
						 X2= (1+tan(Theta))*del - (P/cos(Theta)) + (c*del);
						 Y2= (r+1)*del;
					   }
					 else if (Quad == 3)
					   {
						 X1 = c*del;
						 Y1 = (1+tan(Theta))*del - (P/cos(Theta)) + (r*del);

						 X2 = (c+1)*del;
						 Y2 = (r+1)*del - (P/cos(Theta));
					   }
					 else if(Quad == 4)
					   {
						 X1 = (P/cos(Theta)) - (del*tan(Theta)) + (c*del);
						 Y1 = r*del;
						 
						 X2 = (P/cos(Theta)) + (c*del);
						 Y2 = (r+1)*del;
					   }
				   }
				 else if (  !strcmp(Shape,"Trapezium") && Theta > PiBy4)  // for Theta > PI / 4.0
				   {
					 if(Quad == 1)
					   {
						 X1 = (P/sin(Theta)) - (del/tan(Theta)) + (c*del);
						 Y1 = (r+1)*del;

						 X2 = (P/sin(Theta)) + (c*del);
						 Y2 = r*del;
					   }
					 else if (Quad == 2)
					   {
						 X1 = c*del;
						 Y1 = (P/sin(Theta)) - (del/tan(Theta)) + (r*del);
						 
						 X2 = (c+1)*del;
						 Y2 = (P/sin(Theta)) + (r*del);
					   }
					 else if (Quad == 3)
					   {
						 X1 = (c+1)*del - (P/sin(Theta));
						 Y1 = (r+1)*del;
						 
						 X2 = (1+(1/tan(Theta)))*del - (P/sin(Theta)) + (c*del);
						 Y2 = r*del;
					   }
					 else if(Quad == 4)
					   {
						 X1 = c*del;
						 Y1 = (r+1)*del - (P/sin(Theta));

						 X2 = (c+1)*del;
						 Y2 = (1+ (1/tan(Theta)))*del - (P/sin(Theta)) + (r*del);
					   }
				   }
				 else if( !strcmp(Shape,"Triangle_Complement")  )  //both for Theta <> PI/4
				   {
					 if(Quad == 1)
					   {
						 X1 = -(del/tan(Theta)) + (P/sin(Theta)) + (c*del);
						 Y1 = (r+1)*del;

						 X2 = (c+1)*del;
						 Y2 = (P/cos(Theta)) - (del*tan(Theta)) + (r*del);
					   }
					 else if (Quad == 2)
					   {
						 X1 = c*del;
						 Y1 = (P*sin(Theta)) - (del/tan(Theta)) + (P*pow(cos(Theta),2)/sin(Theta)) + (r*del);
						 
						 X2 = (1+tan(Theta))*del - (P*cos(Theta)) - (P*pow(sin(Theta),2)/cos(Theta)) + (c*del);
						 Y2 = (r+1)*del;
					   }
					 else if (Quad == 3)
					   {
						 X1 = c*del;
						 Y1 = (1+tan(Theta))*del - (P*cos(Theta)) - (P*pow(sin(Theta),2)/cos(Theta)) + (r*del);
						 
						 X2 = (1+(1/tan(Theta)))*del - (P*pow(cos(Theta),2)/sin(Theta)) - (P*sin(Theta)) + (c*del);
						 Y2 = r*del;
					   }
					 else if(Quad == 4)
					   {
						 X1 = (P/cos(Theta)) - (del*tan(Theta)) + (c*del);
						 Y1 = r*del;
						 
						 X2 = (c+1)*del;
						 Y2 = (del/tan(Theta)) - (P/sin(Theta)) + (r+1)*del;
					   }
					}

				   fprintf(ptrWrite,"%Lf \t%Lf\n", X1,X2);
				   fprintf(ptrWrite,"%Lf \t%Lf\n", Y1,Y2);
			   } 	    
       }
  }

  fclose(ptrWrite);   
}

/********************************************************************************************************************************
NAME       : AdvanceF_field
DESCRIPTION: 
LIMITATIONS:
*********************************************************************************************************************************/
void AdvanceF_field(int r,int c, int Flag, int Togg)
{
	  /*******LOCAL VARIABLES*******/
	  long double VolFrac;
	  long double Zeta, Zeta1, Zeta2;
	  long double NetXFlux, NetYFlux, NetXFlux_Corr, NetYFlux_Corr;
	  long double X_Flux_Out, X_Flux_In, X_Flux_Out_Corr; 
	  long double Y_Flux_Out, Y_Flux_In, Y_Flux_Out_Corr;
	  long double X_Flux_Right,X_Flux_Left, X_Flux_Right_Corr, X_Flux_Left_Corr;
	  long double Y_Flux_Up, Y_Flux_Down, Y_Flux_Up_Corr, Y_Flux_Down_Corr;
	  long double Ustar_Left, Ustar_Right, Vstar_Up, Vstar_Down;
	  /*****************************/
  
	  	  if(Flag == 1)    //HORIZONTAL SWEEP
		  {
			    //NetXFlux = delTBydel*(X_Flux[r][c+1] - X_Flux[r][c]);
				NetXFlux = (X_Flux[r][c+1] - X_Flux[r][c]);  // 4-6-15
			    //~ if(X_Velocity[r][c+1] >= 0.0)     //RIGHT SIDE OF THE CELL EFFLUXES
				//~ {					
					 //~ X_Flux_Right = delTBydel*X_Flux[r][c+1];
					 //~ X_Flux_Left  = delTBydel*X_Flux[r][c];
					 
					//~ //NO FLUX THROUGH THE LEFT WALL OF THE CELL - TO THE CELL OR FROM THE CELL
					//~ if(fabs(X_Flux_Left) < LOWER_LIMIT)
					//~ {
						//~ //EVERYTHING PRESENT INSIDE THE CELL LEAVES THROUGH THE RIGHT WALL
						//~ if(X_Flux_Right == Cellattrib[r][c].VolFrac) 
						//~ {	//printf("haha1");
							//~ //APPLY CORRECTION																					
							//~ Ustar_Right  = delTBydel*X_Velocity[r][c+1];

			  				//~ Zeta1 = X_Flux_Right/Ustar_Right;
							//~ Zeta2 = (Cellattrib[r][c].VolFrac - X_Flux_Right)/(1.0 - X_Flux_Right);
							
							//~ if( fabs(Zeta1 - 0.5) >= fabs(Zeta2 - 0.5) )
							//~ {
								//~ Zeta = Zeta1;
							//~ }
							//~ else if ( fabs(Zeta1 - 0.5) < fabs(Zeta2 - 0.5) )
							//~ {
								//~ Zeta = Zeta2;
							//~ }

							//~ X_Flux_Right_Corr = Cellattrib[r][c].delV*(X_Flux_Right - Ustar_Right) + (Cellattrib[r][c].delV + (1-Cellattrib[r][c].delV)*Zeta)*Ustar_Right - X_Flux_Left;
							//~ NetXFlux_Corr     = X_Flux_Right_Corr - X_Flux_Left;
					
							//~ VolFrac = AdvanceF_WithCorrection(r, c, Flag, Togg, NetXFlux_Corr, 0.0);  //Y_flux_corr is not required
						//~ }
						//~ else  //NO CORRECTION
						//~ {
							//~ VolFrac = AdvanceF_WithoutCorrection(r, c, Flag, Togg, NetXFlux, 0.0);
						//~ }
					//~ }
					//~ else  //NO CORRECTION
					//~ {
						//~ VolFrac = AdvanceF_WithoutCorrection(r, c, Flag, Togg, NetXFlux, 0.0);
					//~ }
				//~ }
				//~ else if (X_Velocity[r][c] <= 0.0)   //LEFT SIDE OF THE CELL EFFLUXES
				//~ {											
					 //~ X_Flux_Left  = delTBydel*X_Flux[r][c];
					 //~ X_Flux_Right = delTBydel*X_Flux[r][c+1];
					 
					//~ //NO FLUX THROUGH THE RIGHT WALL OF THE CELL - TO THE CELL OR FROM THE CELL
					 //~ if(fabs(X_Flux_Right) < LOWER_LIMIT)
					 //~ {
						//~ //EVERYTHING PRESENT INSIDE THE CELL LEAVES THROUGH THE LEFT WALL
						//~ if(fabs(X_Flux_Left) == Cellattrib[r][c].VolFrac) 
						//~ { //printf("haha2");
							//~ Ustar_Left   = delTBydel*X_Velocity[r][c];
							//~ //APPLY CORRECTION																					
			  				//~ Zeta1 = X_Flux_Left/Ustar_Left;
							//~ Zeta2 = (Cellattrib[r][c].VolFrac - X_Flux_Left)/(1.0 - X_Flux_Left);
							
							//~ if( fabs(Zeta1 - 0.5) >= fabs(Zeta2 - 0.5) )
							//~ {
								//~ Zeta = Zeta1;
							//~ }
							//~ else if ( fabs(Zeta1 - 0.5) < fabs(Zeta2 - 0.5) )
							//~ {
								//~ Zeta = Zeta2;
							//~ }
							
							//~ X_Flux_Left_Corr  = Cellattrib[r][c].delV*(X_Flux_Left - Ustar_Left) + (Cellattrib[r][c].delV + (1-Cellattrib[r][c].delV)*Zeta)*Ustar_Left - X_Flux_Right;
							//~ NetXFlux_Corr     = X_Flux_Right - X_Flux_Left_Corr;

							//~ VolFrac = AdvanceF_WithCorrection(r, c, Flag, Togg, NetXFlux_Corr, 0.0);  //Y_flux_corr is not required
						//~ }
						//~ else
						//~ {
							//~ VolFrac = AdvanceF_WithoutCorrection(r, c, Flag, Togg, NetXFlux, 0.0);
						//~ }
					//~ }
					//~ else 
					//~ {
						//~ VolFrac = AdvanceF_WithoutCorrection(r, c, Flag, Togg, NetXFlux, 0.0);
					//~ }

				//~ } //X_VEl < 0		
						VolFrac = AdvanceF_WithoutCorrection(r, c, Flag, Togg, NetXFlux, 0.0);
		  }
		  else if(Flag == 2)   //VERTICAL SWEEP
		  {			  
			//NetYFlux = delTBydel*(Y_Flux[r+1][c] - Y_Flux[r][c]);
			NetYFlux = (Y_Flux[r+1][c] - Y_Flux[r][c]);  // 4-6-15
				//~ if(Y_Velocity[r+1][c] >= 0.0)     //UPPER PART OF THE CELL EFFLUXES
				//~ {			
			  	    //~ Y_Flux_Up    = delTBydel*Y_Flux[r+1][c];						
					//~ Y_Flux_Down  = delTBydel*Y_Flux[r][c];

					//~ //NO FLUX THROUGH THE BOTTOM WALL OF THE CELL - TO THE CELL OR FROM THE CELL
					//~ if(fabs(Y_Flux_Down) < LOWER_LIMIT)
					//~ {
						//~ //EVERYTHING PRESENT INSIDE THE CELL LEAVES THROUGH THE UPPER WALL
						//~ if(Y_Flux_Up == Cellattrib[r][c].VolFrac) 
						//~ {	//printf("haha3");
							//~ Vstar_Up     = delTBydel*Y_Velocity[r+1][c];

							//~ Zeta1	= Y_Flux_Up/Vstar_Up;
							//~ Zeta2	= (Cellattrib[r][c].VolFrac - Y_Flux_Up)/(1.0 - Y_Flux_Up);
					
							//~ if( fabs(Zeta1 - 0.5) >= fabs(Zeta2 - 0.5) )
							//~ {
								//~ Zeta = Zeta1;
							//~ }
							//~ else if ( fabs(Zeta1 - 0.5) < fabs(Zeta2 - 0.5) )
							//~ {
								//~ Zeta = Zeta2;
							//~ }

							//~ Y_Flux_Up_Corr = Cellattrib[r][c].delV*(Y_Flux_Up - Vstar_Up) + (Cellattrib[r][c].delV + (1-Cellattrib[r][c].delV)*Zeta)*Vstar_Up - Y_Flux_Down;
							//~ NetYFlux_Corr   = Y_Flux_Up_Corr - Y_Flux_Down;

							//~ VolFrac = AdvanceF_WithCorrection(r, c, Flag, Togg, 0.0, NetYFlux_Corr);  //X_flux_corr is not required
						//~ }
						//~ else  //NO CORRECTION
						//~ {
							//~ VolFrac = AdvanceF_WithoutCorrection(r, c, Flag, Togg, 0.0, NetYFlux);  //NetXFlux is not required to be sent
						//~ }
					//~ }
					//~ else  //NO CORRECTION
					//~ {
						//~ VolFrac = AdvanceF_WithoutCorrection(r, c, Flag, Togg, 0.0, NetYFlux);
					//~ }			

				//~ }
				//~ else if (Y_Velocity[r][c] <= 0.0)   //LOWER PART OF THE CELL EFFLUXES
				//~ {
					//~ VolFrac = AdvanceF_WithoutCorrection(r, c, Flag, Togg, 0.0, NetYFlux);
				//~ }
				 VolFrac = AdvanceF_WithoutCorrection(r, c, Flag, Togg, 0.0, NetYFlux);
		  } //Vertical Sweep

		  Write_struct_Cellattrib(r,c,VolFrac,-100000,-100000,-100000,-100000,"-100000",-100000,-100000);
}

/**********************************************************************************************************
NAME       : AdvanceF_WithCorrection
DESCRIPTION: 
LIMITATIONS:
/**********************************************************************************************************/
long double AdvanceF_WithCorrection(int r, int c, int Flag, int Togg, long double NetXFlux_Corr, long double NetYFlux_Corr)
{
	/*******DECLARE******/
    long double VolFrac;
	/********************/
	printf("nonsense");
	if (Flag == 1)       // HORIZONTAL PASS
	{
		if(Togg == 0)   //HORIZONTAL PASS IS THE FIRST 
		{
			Cellattrib[r][c].delV = 1.0 - delTBydel*(X_Velocity[r][c+1] - X_Velocity[r][c]); 

			//NO CORRECTION
			Cellattrib[r][c].VolFrac = (Cellattrib[r][c].VolFrac - NetXFlux_Corr)/Cellattrib[r][c].delV;		
		}
		else if(Togg == 1)  //HORIZONTAL PASS IS THE SECOND
		{
			//NO CORRECTION
			Cellattrib[r][c].VolFrac = (Cellattrib[r][c].VolFrac*Cellattrib[r][c].delV) - NetXFlux_Corr;
			VolFrac = Cellattrib[r][c].VolFrac;
		}
	}
	else if(Flag == 2)   // VERTICAL PASS
	{
		if(Togg == 1)   //VERTICAL PASS IS THE FIRST 
		{
			Cellattrib[r][c].delV = 1.0 - delTBydel*(Y_Velocity[r+1][c] - Y_Velocity[r][c]); 

			Cellattrib[r][c].VolFrac = (Cellattrib[r][c].VolFrac - NetYFlux_Corr)/Cellattrib[r][c].delV;
			VolFrac = Cellattrib[r][c].VolFrac;
		}
		else if(Togg == 0)  //VERTICAL PASS IS THE SECOND
		{
			//WITH CORRECTION
			 Cellattrib[r][c].VolFrac = (Cellattrib[r][c].VolFrac*Cellattrib[r][c].delV) - NetYFlux_Corr;
			 VolFrac = Cellattrib[r][c].VolFrac;
		}	
	}
	else  //ERROR HANDLING
    {
      printf("\n***********************************************************\n");
      printf("\nERROR - ERROR IN SUBROUTINE AdvanceF_WithCorrection\n");
      printf("\n***********************************************************\n");
      exit(0);
    }

	return VolFrac;
}

/**********************************************************************************************************
NAME       : AdvanceF_WithoutCorrection
DESCRIPTION: 
LIMITATIONS:
/**********************************************************************************************************/
long double AdvanceF_WithoutCorrection(int r, int c, int Flag, int Togg, long double NetXFlux, long double NetYFlux)
{
	/*******DECLARE******/
    long double VolFrac;
	/********************/

	if (Flag == 1)       // HORIZONTAL PASS
	{
		if(Togg == 0)   //HORIZONTAL PASS IS THE FIRST 
		{
			//Cellattrib[r][c].delV = 1.0 - delTBydel*(X_Velocity[r][c+1] - X_Velocity[r][c]); //maximum volume change
			Cellattrib[r][c].delV = del_sq - delTdel*(X_Velocity[r][c+1] - X_Velocity[r][c]); // 4-6-15
			//NO CORRECTION
			//Cellattrib[r][c].VolFrac = (Cellattrib[r][c].VolFrac - NetXFlux)/Cellattrib[r][c].delV;
			//Cellattrib[r][c].VolFrac = (Cellattrib[r][c].VolFrac*del_sq - NetXFlux)/Cellattrib[r][c].delV; // 4-6-15   //why to update structure here ?? it is already updating in advanceF_field
			
			//VolFrac = Cellattrib[r][c].VolFrac;									
			VolFrac =  (Cellattrib[r][c].VolFrac*del_sq - NetXFlux)/Cellattrib[r][c].delV;
			
		}
		else if(Togg == 1)  //HORIZONTAL PASS IS THE SECOND
		{
			//NO CORRECTION
			//Cellattrib[r][c].VolFrac = (Cellattrib[r][c].VolFrac*Cellattrib[r][c].delV) - NetXFlux;
			//VolFrac = Cellattrib[r][c].VolFrac;
			VolFrac = (Cellattrib[r][c].VolFrac*Cellattrib[r][c].delV - NetXFlux)/del_sq;
		}

	}
	else if(Flag == 2)   // VERTICAL PASS
	{
		if(Togg == 1)   //VERTICAL PASS IS THE FIRST 
		{
			//Cellattrib[r][c].delV = 1.0 - delTBydel*(Y_Velocity[r+1][c] - Y_Velocity[r][c]); 
				Cellattrib[r][c].delV = del_sq - delTdel*(Y_Velocity[r+1][c] - Y_Velocity[r][c]); 
			//Cellattrib[r][c].VolFrac = (Cellattrib[r][c].VolFrac - NetYFlux)/Cellattrib[r][c].delV;
			//VolFrac = Cellattrib[r][c].VolFrac;
			VolFrac =  (Cellattrib[r][c].VolFrac*del_sq - NetYFlux)/Cellattrib[r][c].delV;
		}
		else if(Togg == 0)  //VERTICAL PASS IS THE SECOND
		{
			//WITH CORRECTION
			 //Cellattrib[r][c].VolFrac = (Cellattrib[r][c].VolFrac*Cellattrib[r][c].delV) - NetYFlux;
			 //VolFrac = Cellattrib[r][c].VolFrac;
			VolFrac = ((Cellattrib[r][c].VolFrac*Cellattrib[r][c].delV) - NetYFlux)/del_sq;
		}
	}
	else  //ERROR HANDLING
    {
      printf("\n***********************************************************\n");
      printf("\nERROR - ERROR IN SUBROUTINE AdvanceF_WithoutCorrection\n");
      printf("\n***********************************************************\n");
      exit(0);
    }

	return VolFrac;
}


/**********************************************************************************************************
NAME: CalculateF_Flux
DESCRIPTION: 
			1. If Flag = 1 , this subroutine does a horizontal pass and calculates all horizontal fluxes
			2. If Flag = 2 , this subroutine does a vertical pass and calculates all vertical fluxes
LIMITATIONS:
**********************************************************************************************************/
void CalculateF_Flux(int flag)
{
  /*******LOCAL VARIABLES**********/
  int		  r,c, Quad;
  long double Theta;
  char		  Shape[30];
  long double x_0,y_0,x_del,y_del,y_XVeldelT,x_YVeldelT;
  long double X_Vel, Y_Vel;
  long double XVeldelT,YVeldelT;
  /********************************/

		//HORIZONTAL PASS
		if(flag == 1)
		{
		  for(r=0;r < NoOfRows;r++)
		   {
			 for(c=0;c <= NoOfCols;c++)
			  {
				/************INITIALISE**********/
				 X_Vel = X_Velocity[r][c];
				 Y_Vel = Y_Velocity[r][c];
				 /******************************/

				 if(X_Vel > 0.0)   //THE CELL TO THE LEFT OF THE BOUNDARY DECIDES THE EFFLUX AND THAT CELL WILL FLUX FROM ITS RIGHT BOUNDARY
				  {
					if(c == 0)	
					{						
						//X_Flux[r][c] = X_Vel*delT*del/delTdel;  ////CHANGE - CURRENT ASSUMPTION IS THAT AT THE NEIGHBOURING CELL IS FULL						
						X_Flux[r][c] = 0.0;   //PRESENTLY THE ONLY THE INTERIOR PART OF THE INTERFACE CONTAINS FLUID AND HENCE NO FLUX AT THE BOUNDARY
						continue;						  	
					}
					
					if(Cellattrib[r][c-1].VolFrac > LOWER_LIMIT && Cellattrib[r][c-1].VolFrac < UPPER_LIMIT) //FOR CELLS CONTAINING INTERFACE
					{								
						//FETCH INTERFACE VALUES FROM [r][c-1] - INIT BLOCK
						Quad = Find_Quad(r,c-1);
						strcpy(Shape,Cellattrib[r][c-1].Shape);  //GET SHAPE FROM [r][c-1]			

						/**************************************/
	
						//THE PARAMETERS PASSED TO THE FUNCTION AssignVal_FluxVariables() WILL DIFFER FOR 1/3 AND 2/4 QUADRANTS. THE XVEL AND YVEL HAVE TO BE INTERCHANGED FOR 2/4 QUADRANTS
						if(Quad == 1)
						{
							AssignVal_FluxVariables(r, c-1, &Theta, &XVeldelT, &YVeldelT, &x_0, &y_0, &x_del, &y_del, &x_YVeldelT, &y_XVeldelT, &X_Vel, &Y_Vel);		 
							X_Flux[r][c] = CalculateX_Flux("Right", Shape, Theta, XVeldelT, YVeldelT, x_0, y_0, x_del, y_del, x_YVeldelT, y_XVeldelT);																				
						}
						else if(Quad == 2)
						{
							AssignVal_FluxVariables(r, c-1, &Theta, &XVeldelT, &YVeldelT, &x_0, &y_0, &x_del, &y_del, &x_YVeldelT, &y_XVeldelT, &Y_Vel, &X_Vel);		 					
							X_Flux[r][c] = CalculateY_Flux("Lower", Shape, Theta, XVeldelT, YVeldelT, x_0, y_0, x_del, y_del, x_YVeldelT, y_XVeldelT);							
						}
						else if(Quad == 3)
						{
							AssignVal_FluxVariables(r, c-1, &Theta, &XVeldelT, &YVeldelT, &x_0, &y_0, &x_del, &y_del, &x_YVeldelT, &y_XVeldelT, &X_Vel, &Y_Vel);
							X_Flux[r][c] = CalculateX_Flux("Left", Shape, Theta, XVeldelT, YVeldelT, x_0, y_0, x_del, y_del, x_YVeldelT, y_XVeldelT);
						}
						else if(Quad == 4)
						{
							AssignVal_FluxVariables(r, c-1, &Theta, &XVeldelT, &YVeldelT, &x_0, &y_0, &x_del, &y_del, &x_YVeldelT, &y_XVeldelT, &Y_Vel, &X_Vel);
							X_Flux[r][c] = CalculateY_Flux("Upper", Shape, Theta, XVeldelT, YVeldelT, x_0, y_0, x_del, y_del, x_YVeldelT, y_XVeldelT);							
						}	
					}
					else if(Cellattrib[r][c-1].VolFrac <= LOWER_LIMIT) // NO FLUX
					{
						X_Flux[r][c] = 0.0;
					}
					else if(Cellattrib[r][c-1].VolFrac >= UPPER_LIMIT) //CELL IS FULL - PARTIAL FLUX
					{
						XVeldelT     = X_Vel*delT;
						//X_Flux[r][c] = XVeldelT*del/delTdel;  // KEPT AS THE SAME FORMULA ALTHOUGH NUM. & DEN. CANCEL		
						X_Flux[r][c] = XVeldelT*del; //4-6-15
					}					
				  }
				else if (X_Vel < 0.0)   //THE CELL TO THE RIGHT OF THE BOUNDARY DECIDES THE FLUX AND THAT CELL WILL FLUX FROM ITS LEFT BOUNDARY
				  {	       			
					if(c == NoOfCols)
					{
					  //X_Flux[r][c] = X_Vel*delT*del/delTdel;  ////CHANGE - CURRENT ASSUMPTION IS THAT AT THE NEIGHBOURING CELL IS FULL						
					  X_Flux[r][c] = 0.0;   //PRESENTLY THE ONLY THE INTERIOR PART OF THE INTERFACE CONTAINS FLUID AND HENCE NO FLUX AT THE BOUNDARY
					  continue;
					}
					
					X_Vel = fabs(X_Vel);   
					Y_Vel = fabs(Y_Vel);

					if(Cellattrib[r][c].VolFrac > LOWER_LIMIT && Cellattrib[r][c].VolFrac < UPPER_LIMIT) //FOR CELLS CONTAINING INTERFACE
					{
						//FETCH INTERFACE VALUES FROM [r][c] - INIT BLOCK
						Quad = Find_Quad(r,c);							
						strcpy(Shape,Cellattrib[r][c].Shape);  //GET SHAPE FROM [r][c]
						/**************************************/			

						//THE PARAMETERS PASSED TO THE FUNCTION AssignVal_FluxVariables() WILL DIFFER FOR 1/3 AND 2/4 QUADRANTS. THE XVEL AND YVEL HAVE TO BE INTERCHANGED FOR 2/4 QUADRANTS
						if(Quad == 1)
						 {
							AssignVal_FluxVariables(r, c, &Theta, &XVeldelT, &YVeldelT, &x_0, &y_0, &x_del, &y_del, &x_YVeldelT, &y_XVeldelT, &X_Vel, &Y_Vel);
		 					X_Flux[r][c] = -CalculateX_Flux("Left",  Shape, Theta, XVeldelT, YVeldelT, x_0, y_0, x_del, y_del, x_YVeldelT, y_XVeldelT);
						 }
						 else if(Quad == 2)
						 {							 
							AssignVal_FluxVariables(r, c, &Theta, &XVeldelT, &YVeldelT, &x_0, &y_0, &x_del, &y_del, &x_YVeldelT, &y_XVeldelT, &Y_Vel, &X_Vel);
							X_Flux[r][c] = -CalculateY_Flux("Upper", Shape, Theta, XVeldelT, YVeldelT, x_0, y_0, x_del, y_del, x_YVeldelT, y_XVeldelT);						
							//if(r == 30 & c == 5)
							//	printf("\nVelocity Negative, Quad 2 %d %d %lf\n",r,c,X_Flux[r][c]);
						 }
						 else if(Quad == 3)
						 {
							AssignVal_FluxVariables(r, c, &Theta, &XVeldelT, &YVeldelT, &x_0, &y_0, &x_del, &y_del, &x_YVeldelT, &y_XVeldelT, &X_Vel, &Y_Vel);
						    X_Flux[r][c] = -CalculateX_Flux("Right", Shape, Theta, XVeldelT, YVeldelT, x_0, y_0, x_del, y_del, x_YVeldelT, y_XVeldelT);
								//if(r == 15 & c == 38)
								//printf("\nVelocity Negative Quad is 3rd %d %d %lf\n",r,c,X_Flux[r][c]);
						 }
						 else if(Quad == 4)
						 {
							AssignVal_FluxVariables(r, c, &Theta, &XVeldelT, &YVeldelT, &x_0, &y_0, &x_del, &y_del, &x_YVeldelT, &y_XVeldelT, &Y_Vel, &X_Vel);
							X_Flux[r][c] = -CalculateY_Flux("Lower", Shape, Theta, XVeldelT, YVeldelT, x_0, y_0, x_del, y_del, x_YVeldelT, y_XVeldelT);
						 }
					}
					else if(Cellattrib[r][c].VolFrac <= LOWER_LIMIT) // NO FLUX
					{
						X_Flux[r][c] = 0.0;
					}
					else if(Cellattrib[r][c].VolFrac >= UPPER_LIMIT) //CELL IS FULL - PARTIAL FLUX
					{
						XVeldelT     = X_Vel*delT;
						//X_Flux[r][c] = -XVeldelT*del/delTdel;  // KEPT AS THE SAME FORMULA ALTHOUGH NUM. & DEN. CANCEL
						X_Flux[r][c] = -XVeldelT*del;   // corrected on 4-6-15
					}
				  }
				else //NO FLUX
				  {
					X_Flux[r][c] = 0.0;					
				  }						
			  }
		   }
		}	
	  //VERTICAL PASS
		else if(flag == 2)
		{
		  for(c=0;c < NoOfCols;c++)
		  {
			  for(r=0;r <= NoOfRows;r++)
			  {				
				X_Vel = X_Velocity[r][c];
				Y_Vel = Y_Velocity[r][c];

		  		if(Y_Vel > 0.0)   //THE CELL BELOW THE BOUNDARY DECIDES THE FLUX AND THAT CELL WILL FLUX FROM ITS UPPER BOUNDARY
				  {
					if(r == 0)
					{					  
					  //Y_Flux[r][c] = Y_Vel*delT*del/delTdel;          ////CHANGE - CURRENT ASSUMPTION IS THAT AT THE BOUNDARY NEIGNBOURS ARE FULL
					  Y_Flux[r][c] = 0.0;   //PRESENTLY THE ONLY THE INTERIOR PART OF THE INTERFACE CONTAINS FLUID AND HENCE NO FLUX AT THE BOUNDARY
					  continue;	
					}
					
					if(Cellattrib[r-1][c].VolFrac > LOWER_LIMIT && Cellattrib[r-1][c].VolFrac < UPPER_LIMIT)
					{
						//FETCH INTERFACE DETAILS FROM [r-1][c] - INIT
						Quad = Find_Quad(r-1,c);									
						strcpy(Shape, Cellattrib[r-1][c].Shape);  //GET SHAPE FROM [r - 1][c]
						/**************************************/

						if(Quad == 1)
						 {		 
							AssignVal_FluxVariables(r-1, c, &Theta, &XVeldelT, &YVeldelT, &x_0, &y_0, &x_del, &y_del, &x_YVeldelT, &y_XVeldelT, &X_Vel, &Y_Vel);
							Y_Flux[r][c] = CalculateY_Flux("Upper", Shape, Theta, XVeldelT, YVeldelT, x_0, y_0, x_del, y_del, x_YVeldelT, y_XVeldelT);
						 }
						 else if(Quad == 2)
						 {
							AssignVal_FluxVariables(r-1, c, &Theta, &XVeldelT, &YVeldelT, &x_0, &y_0, &x_del, &y_del, &x_YVeldelT, &y_XVeldelT, &Y_Vel, &X_Vel);
	  					    Y_Flux[r][c] = CalculateX_Flux("Right", Shape, Theta, XVeldelT, YVeldelT, x_0, y_0, x_del, y_del, x_YVeldelT, y_XVeldelT);
						 }
						 else if(Quad == 3)
						 {
							AssignVal_FluxVariables(r-1, c, &Theta, &XVeldelT, &YVeldelT, &x_0, &y_0, &x_del, &y_del, &x_YVeldelT, &y_XVeldelT, &X_Vel, &Y_Vel);
							Y_Flux[r][c] = CalculateY_Flux("Lower", Shape, Theta, XVeldelT, YVeldelT, x_0, y_0, x_del, y_del, x_YVeldelT, y_XVeldelT);
						 }
						 else if(Quad == 4)
						 {
							AssignVal_FluxVariables(r-1, c, &Theta, &XVeldelT, &YVeldelT, &x_0, &y_0, &x_del, &y_del, &x_YVeldelT, &y_XVeldelT, &Y_Vel, &X_Vel);
							Y_Flux[r][c] = CalculateX_Flux("Left", Shape, Theta, XVeldelT, YVeldelT, x_0, y_0, x_del, y_del, x_YVeldelT, y_XVeldelT);
						 }	
					}
					else if (Cellattrib[r-1][c].VolFrac <= LOWER_LIMIT) // NO FLUX
					{
						Y_Flux[r][c] = 0.0;
					}
					else if (Cellattrib[r-1][c].VolFrac >= UPPER_LIMIT)  //CELL IS FULL - PARTIAL FLUX
					{
						YVeldelT     = Y_Vel*delT;
						//Y_Flux[r][c] = YVeldelT*del/delTdel;// 	
						Y_Flux[r][c] = YVeldelT*del;  // 4-6-15
					}
				  }
				else if (Y_Vel < 0.0)   //THE CELL ABOVE THE BOUNDARY DECIDES THE FLUX AND THAT CELL WILL FLUX FROM ITS LOWER BOUNDARY
				  {					
					if(r == NoOfRows)
					{
					  //Y_Flux[r][c] = Y_Vel*delT*del/delTdel;          ////CHANGE - CURRENT ASSUMPTION IS THAT AT THE BOUNDARY NEIGNBOURS ARE FULL
					  Y_Flux[r][c] = 0.0;   //PRESENTLY THE ONLY THE INTERIOR PART OF THE INTERFACE CONTAINS FLUID AND HENCE NO FLUX AT THE BOUNDARY
					  continue;	
					}

					X_Vel = fabs(X_Vel);   
					Y_Vel = fabs(Y_Vel);

					if(Cellattrib[r][c].VolFrac > LOWER_LIMIT && Cellattrib[r][c].VolFrac < UPPER_LIMIT)
					{
						//FETCH INTERFACE DETAILS FROM [r][c] - INIT
						Quad = Find_Quad(r,c);						
						strcpy(Shape, Cellattrib[r][c].Shape);  //GET SHAPE FROM [r][c]
						/**************************************/

						if(Quad == 1)
						 {
							AssignVal_FluxVariables(r, c, &Theta, &XVeldelT, &YVeldelT, &x_0, &y_0, &x_del, &y_del, &x_YVeldelT, &y_XVeldelT, &X_Vel, &Y_Vel);
		 				    Y_Flux[r][c] = -CalculateY_Flux("Lower", Shape, Theta, XVeldelT, YVeldelT, x_0, y_0, x_del, y_del, x_YVeldelT, y_XVeldelT);
						 }
						 else if(Quad == 2)
						 {
							AssignVal_FluxVariables(r, c, &Theta, &XVeldelT, &YVeldelT, &x_0, &y_0, &x_del, &y_del, &x_YVeldelT, &y_XVeldelT, &Y_Vel, &X_Vel);
							Y_Flux[r][c] = -CalculateX_Flux("Left", Shape, Theta, XVeldelT, YVeldelT, x_0, y_0, x_del, y_del, x_YVeldelT, y_XVeldelT);
						 }
						 else if(Quad == 3)
						 {
							AssignVal_FluxVariables(r, c, &Theta, &XVeldelT, &YVeldelT, &x_0, &y_0, &x_del, &y_del, &x_YVeldelT, &y_XVeldelT, &X_Vel, &Y_Vel);
							Y_Flux[r][c] = -CalculateY_Flux("Upper", Shape, Theta, XVeldelT, YVeldelT, x_0, y_0, x_del, y_del, x_YVeldelT, y_XVeldelT);
						 }
						 else if(Quad == 4)
						 {
							AssignVal_FluxVariables(r, c, &Theta, &XVeldelT, &YVeldelT, &x_0, &y_0, &x_del, &y_del, &x_YVeldelT, &y_XVeldelT, &Y_Vel, &X_Vel);
						    Y_Flux[r][c] = -CalculateX_Flux("Right", Shape, Theta, XVeldelT, YVeldelT, x_0, y_0, x_del, y_del, x_YVeldelT, y_XVeldelT);
						 }
					}
					else if(Cellattrib[r][c].VolFrac <= LOWER_LIMIT)  //NO FLUX
					{
						Y_Flux[r][c] = 0.0;
					}
					else if(Cellattrib[r][c].VolFrac >= UPPER_LIMIT)   //CELL IS FULL - PARTIAL FLUX
					{
						YVeldelT     = Y_Vel*delT;
						//Y_Flux[r][c] = -YVeldelT*del/delTdel;
						Y_Flux[r][c] = -YVeldelT*del; //4-6-15
					}
				  }
				else //NO FLUX
				  {
						Y_Flux[r][c] = 0.0;
				  }
			  }
		  }
		}
}

/******************************************************************************
NAME:		 CalculateX_Flux
DESCRIPTION: Calculates the fluxes for the vertical walls
LIMITATIONS:
******************************************************************************/
long double CalculateX_Flux(char Wall[20], char Shape[20], long double Theta, long double XVeldelT, long double YVeldelT, long double x_0, long double y_0, long double x_del, long double y_del, long double x_YVeldelT, long double y_XVeldelT)
{
	/***LOCAL VARIABLES***/
	long double Area_Flux, temp;
	long double tanT;
	/********************/

	/*****INITIALISE***/
	Area_Flux = 0.0;
	temp      = 0.0;
	tanT      = tan(Theta);
	/******************/

	  if(!strcmp(Wall,"Left"))      //THE WALL OF A CELL THROUGH WHICH THE FLUX LEAVES
	  {	
		  
		    if(!strcmp(Shape,"Rectangle"))
			{
				 if( (XVeldelT) >= x_0 )   //ENTIRE RECTANGLE LEAVES 
				   {
					 Area_Flux = x_0*del;
				   }
				 else                      //A TRAPEZIUM LEAVES
				   {
					 Area_Flux = (XVeldelT*del);
				   }
			}
		else  if(!strcmp(Shape,"Triangle"))
			{
				 if( (XVeldelT) >= x_0 )   //ENTIRE TRIANGLE LEAVES		 
				   {
					 Area_Flux = 0.5*x_0*y_0;
				   }
				 else                      //A TRAPEZIUM LEAVES
				   {
					 Area_Flux = ( y_0 + y_XVeldelT )*(XVeldelT*0.5);
				   }
			}
		  else if(!strcmp(Shape,"Trapezium") && Theta < PiBy4)
			{
					 Area_Flux = ( y_0 + y_XVeldelT)*(XVeldelT*0.5);  //A TRAPEZIUM LEAVES
			}
		  else  if(!strcmp(Shape,"Trapezium") && Theta > PiBy4)
			{
				 if(XVeldelT <= x_del)   //A RECTANGLE LEAVES
				   {
					 Area_Flux = XVeldelT*del;
				   }
				 else if( (XVeldelT > x_del) && (XVeldelT < x_0)  )             //A 5-SIDED FIGURE LEAVES
				   {
					 Area_Flux = ( (x_del + x_0)*delBy2 ) - ( 0.5*(x_0 - XVeldelT)*y_XVeldelT );   //TRAPEZIUM - TRIANGLE
				   }
				 else if( XVeldelT >= x_0 )    //ENTIRE TRAPEZIUM LEAVES
					{
					  Area_Flux = (x_del + x_0)*delBy2;	
					}
			}
		  else  if(!strcmp(Shape,"Triangle_Complement"))
			{
				 if(XVeldelT <= x_del)        // A RECTANGLE LEAVES
				   {
					 Area_Flux = del*XVeldelT;
				   }
				 else if(XVeldelT > x_del)   //A 5-sided figure leaves
				   {
					 Area_Flux = (x_del*del) + ( (del + y_XVeldelT)*(XVeldelT - x_del)*0.5 );  //RECTANGLE + TRAPEZIUM
				   }
			}
	  }
	  else if(!strcmp(Wall,"Right"))   //THE WALL OF A CELL THROUGH WHICH THE FLUX LEAVES
	  {						
		  if(!strcmp(Shape,"Rectangle"))
			{
				if( (del - XVeldelT) >= x_0)  // NOTHING LEAVES
				{
					Area_Flux = 0.0;
				}
				else if((del - XVeldelT) < x_0)   // A RECTANGLE LEAVES
				{
					Area_Flux = del*(x_0 - del + XVeldelT );					
				}					
			}
		  else if(!strcmp(Shape,"Triangle"))
			{
				if( (del - XVeldelT) >= x_0)  // NOTHING LEAVES
				{
					Area_Flux = 0.0;
				}
				else if((del - XVeldelT) < x_0)   // A TRIANGLE LEAVES
				{
					Area_Flux = 0.5*pow( (x_0 - del + XVeldelT), 2 )*tanT;					
				}					
			}
		  else if(!strcmp(Shape,"Trapezium") && Theta < PiBy4)
			{
					Area_Flux = ( y_0 - ((del - XVeldelT)*tanT) + y_del )*XVeldelT*0.5;  //A TRAPEZIUM LEAVES
			}
		  else  if(!strcmp(Shape,"Trapezium") && Theta > PiBy4)
			{
				if( (del - XVeldelT) < x_del )   //A TRAPEZIUM LEAVES
				{
					Area_Flux = (x_del - (2.0*(del - XVeldelT)) + x_0)*delBy2;
				}
				else if(  (del - XVeldelT) >= x_del && (del - XVeldelT) < x_0)    //A TRIANGLE LEAVES
				{
					Area_Flux = 0.5*pow( (y_0 - (del - XVeldelT)*tanT) ,2)/tanT;
				}
				else if( (del - XVeldelT) >= x_0)  //NOTHING LEAVES
				{
					Area_Flux = 0.0;
				}
			}
		  else  if(!strcmp(Shape,"Triangle_Complement"))
			{
				if( (del - XVeldelT) < x_del)   //A 5-SIDED FIGURE LEAVES
				{
					Area_Flux = (del*(x_del - del + XVeldelT)) + ((del + y_del)*(del - x_del)*0.5);  //RECTANGLE + TRAPEZIUM
				}
				else if((del - XVeldelT) >= x_del)   //A TRAPEZIUM LEAVES
				{
					Area_Flux = ( y_del + y_0 - ((del - XVeldelT)*tanT) )*XVeldelT*0.5;
				}
			}
	  }
	 //temp = Area_Flux/delTdel;

	  //return temp;
	  return Area_Flux;
}

/******************************************************************************
NAME:	     CalculateY_Flux
DESCRIPTION: Calculates the fluxes for the horizontal walls
LIMITATIONS:
******************************************************************************/
long double CalculateY_Flux(char Wall[20], char Shape[20], long double Theta, long double XVeldelT, long double YVeldelT, long double x_0, long double y_0, long double x_del, long double y_del, long double x_YVeldelT, long double y_XVeldelT)
{
	/***LOCAL VARIABLES***/
	long double Area_Flux, temp;
	long double tanT;
	/********************/

	/*****INITIALISE***/
	Area_Flux = 0.0;
	temp      = 0.0;
	tanT      = tan(Theta);
	/******************/

	  	if(!strcmp(Wall,"Lower"))    //THE WALL OF A CELL THROUGH WHICH THE FLUX LEAVES
		{	if(!strcmp(Shape,"Rectangle"))
				{
					//A rectangle always leaves
					Area_Flux = (x_0*YVeldelT);
					
				}
			else if(!strcmp(Shape,"Triangle"))
				{
					if(YVeldelT >= y_0)  //ENTIRE TRIANGLE LEAVES
					{
						Area_Flux = 0.5*x_0*y_0;
					}
					else if(YVeldelT < y_0)  //A PARALLELOGRAM LEAVES 
					{
						Area_Flux = (x_0 + x_YVeldelT)*YVeldelT*0.5;
					}
				}
			  else if(!strcmp(Shape,"Trapezium") && Theta < PiBy4)
				{
					if(YVeldelT >= y_0)   //ENTIRE TRAPEZIUM LEAVES
					{
						Area_Flux = (y_0 + y_del)*delBy2;
					}
					else if(YVeldelT > y_del && YVeldelT < y_0)  //A 5-SIDED FIGURE LEAVES
					{
						Area_Flux = (delBy2*(y_0 + y_del)) - ( (y_0 - YVeldelT)*x_YVeldelT*0.5 );  //AREA OF TRAPEZIUM - AREA OF TRIANGLE
					}
					else if(YVeldelT <= y_del) //A RECTANGLE LEAVES
					{
						Area_Flux = del*YVeldelT;
					}
				}
			  else  if(!strcmp(Shape,"Trapezium") && Theta > PiBy4)
				{
						Area_Flux = (x_0 + x_YVeldelT)*YVeldelT*0.5;	
				}
			  else  if(!strcmp(Shape,"Triangle_Complement"))
				{
					if(YVeldelT > y_del)   //A 5-SIDED FIGURE LEAVES
					{
						Area_Flux = (del*y_del) + ( (del + x_YVeldelT)*(YVeldelT - y_del)*0.5 );  //AREA OF RECTANGLE + AREA OF TRAPEZIUM
					}
					else if(YVeldelT <= y_del)  //A RECTANGLE LEAVES
					{
						Area_Flux = del*YVeldelT;
					}
				}
		}
		else if(!strcmp(Wall,"Upper"))   //THE WALL OF A CELL THROUGH WHICH THE FLUX LEAVES
		{	if(!strcmp(Shape,"Rectangle"))
				{
					//A rectangle always leaves
					Area_Flux = (x_0*YVeldelT);
					
				}
			else if(!strcmp(Shape,"Triangle"))
				{
					if(YVeldelT <= (del - y_0))   //NOTHING LEAVES
					{
						Area_Flux = 0.0;
					}
					else if(YVeldelT > (del - y_0))  //A TRIANGLE LEAVES
					{
						Area_Flux = 0.5*pow(y_0 - del + YVeldelT,2)/tanT;
					}
				}
			  else if(!strcmp(Shape,"Trapezium") && Theta < PiBy4)
				{
					if(YVeldelT <= (del - y_0))  // NOTHING LEAVES
					{
						Area_Flux = 0.0;
					}
					else if(YVeldelT > (del - y_0) && YVeldelT <= (del - y_del))  //A TRIANGLE LEAVES
					{
						Area_Flux = 0.5*pow(y_0 - del + YVeldelT,2)/tanT;
					}
					else if(YVeldelT > (del - y_del)) //A TRAPEZIUM LEAVES
					{
						Area_Flux = ((y_0 + y_del)*delBy2) - ((del - YVeldelT)*del);
					}
				}
			  else  if(!strcmp(Shape,"Trapezium") && Theta > PI/4.0)  //A TRAPEZIUM LEAVES
				{
						Area_Flux = ( x_del + ((y_0 - del + YVeldelT)/tanT) )*YVeldelT*0.5;
				}
			  else  if(!strcmp(Shape,"Triangle_Complement")) 
				{
					if(YVeldelT <= (del - y_del))    //A TAPEZIUM LEAVES
					{
						Area_Flux = ( x_del + ((y_0 - del + YVeldelT)/tanT) )*YVeldelT*0.5;
					}
					else if(YVeldelT > (del - y_del))  //A 5-SIDED FIGURE LEAVES
					{
						Area_Flux = ( (y_del - del + YVeldelT)*del) + ((x_del + del)*(del - y_del)*0.5);  //AREA OF RECTANGLE + AREA OF TRAPEZIUM
					}
				}
		}

		//temp = Area_Flux/delTdel;
		//return temp;
		return Area_Flux;
}

/******************************************************************************
NAME:		 AssignVal_FluxVariables
DESCRIPTION: Calculates the interface orientation details, for the cell whose identifier is passed.
			 The calculated values are written into the addresses passed.
LIMITATIONS:
******************************************************************************/
void AssignVal_FluxVariables(int r, int c, long double *ptr2, long double *ptr3, long double *ptr4, long double *ptr5, long double *ptr6, long double *ptr7, long double *ptr8, long double *ptr9, long double *ptr10, long double *ptr11, long double *ptr12)
{
  /*****LOOKUP FOR ARGUMENTS****
   long double *ptr2  - Theta
   long double *ptr3  - XVeldelT
   long double *ptr4  - YVeldelT
   long double *ptr5  - x_0
   long double *ptr6  - y_0
   long double *ptr7  - x_del
   long double *ptr8  - y_del
   long double *ptr9  - x_YVeldelT
   long double *ptr10 - y_XVeldelT
   long double *ptr11 - XVel
   long double *ptr12 - YVel
  /****************************/

  /*********DECLARATIONS***********/
  long double P, XVel, YVel;
  long double Theta;
  long double sinT,cosT,tanT;
  /********************************/

  /*************INITIALIZE*************/
  P				= Cellattrib[r][c].P;
  XVel          = *ptr11;                   //DO NOT RETRIEVE VELOCITY FROM CELL WALL .RETRIEVE ONLY FROM THE PASSED ADDRESS
  YVel          = *ptr12;                   //DO NOT RETRIEVE VELOCITY FROM CELL. RETRIEVE ONLY FROM THE PASSED ADDRESS SINCE CELL WALL IDENTIFIER IS DIFFRERENT FROM THE IDENTIFIER PASSED
  Theta         = Cellattrib[r][c].Theta;
  sinT			= sin(Theta);
  cosT			= cos(Theta);
  tanT			= tan(Theta);  
  /************************************/
		
		
  /*********WRITING TO THE ADDRESSES PASSED************/
  *ptr2			= Cellattrib[r][c].Theta;  
  *ptr3			= XVel*delT;
  *ptr4			= YVel*delT;
  *ptr5			= P/sinT;
  *ptr6			= P/cosT;
  *ptr7			= (*ptr6 - del)/tanT;
  *ptr8			= (*ptr6 - (del*tanT));
  *ptr9			= (*ptr6 - *ptr4)/tanT;
  *ptr10		= (*ptr6) - ((*ptr3)*tanT);
  /****************************************************/
  //FLAG: the following block is written to add the rectangle case in flux calculation
  
  if (!strcmp(Cellattrib[r][c].Shape,"Rectangle"))
  {
			 *ptr2			= Cellattrib[r][c].Theta;  
			  *ptr5			= P;
			 // *ptr6			= P/cosT;
			  *ptr7			= P;
			 // *ptr8			= (*ptr6 - (del*tanT));
			  *ptr9			= P;
			  //*ptr10		= (*ptr6) - ((*ptr3)*tanT);

			
  }
}

/**************************************************************************
NAME       : BookKeeping
DESCRIPTION: 
LIMITATIONS:
***************************************************************************/
void BookKeeping(int Flag)
{
    /*********DECLARATIONS***********/
	int r, c;
	long double VolFrac;
	/*******************************/
	
	//ASSIGN 0 AND 1 TO CELLS WHICH ARE OUT OF BOUND
	for(r=0;r < NoOfRows; r++)
	{
	   for(c=0;c < NoOfCols; c++)
		{
			if(Cellattrib[r][c].VolFrac < LOWER_LIMIT)
			  {
				//printf("\nCame here\n");
				VolFrac = 0.0;
				Write_struct_Cellattrib(r,c,VolFrac,-100000,-100000,-100000,-100000,"-100000",-100000,-100000);
			  }
			if(Cellattrib[r][c].VolFrac > UPPER_LIMIT)
			  {
				VolFrac = 1.0;
				Write_struct_Cellattrib(r,c,VolFrac,-100000,-100000,-100000,-100000,"-100000",-100000,-100000);
			  }
	   }
	}
	/***/
}
