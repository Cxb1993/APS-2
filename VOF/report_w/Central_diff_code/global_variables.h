/***********START OF USER INPUT***********************************/

/*****************GRID INPUT**********************************************/
#define NoOfRows 200
#define NoOfCols 200
/*******************************************************************/

/******************LVIRA INPUT****/
#define Rotation_StepSize 0.00031
/*********************************/

/**************************/
#define PI          3.14159265358
#define UPPER_LIMIT 0.99999
#define LOWER_LIMIT 0.00001
/**************************/
/**************NAVIER STOKES INPUT********************/
#define Rey  100
#define Fr   10
#define We   0.1 
#define del  0.02
#define delT 0.005
/*****************************************************/



/*****************END OF USER INPUT**********************************/


/*------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------*/


/*############################THE CODE BELOW DOES NOT REQUIRE ANY USER ENTRY#################################*/

typedef struct				/*2D STRUCT STORES THE ATTRIBUTES FOR A CELL*/
   {
     /***************************CELL ATTRIBUTES - LVIRA******************************************/
      long double  VolFrac;			//VOLUME FRACTION FOR EACH CELL      
      long double  Area;			//THE AREA CALCULATED FOR EACH CELL DURING INTERFACE RECONSTRUCTION
      long double  Nx;				//X-COMPONENT OF THE UNIT NORMAL TO THE LINE
      long double  Ny;				//Y-COMPONENT OF THE UNIT NORMAL TO THE LINE
      long double  P;				//PERPENDICULAR DISTANCE OF THE LINE FROM THE CENTER OF THE CELL
      char         Shape[30];       //SHAPE OF THE AREA BASED ON FITIING A STRAIGHT LINE
      long double  Theta;			//ANGLE MADE BY THE LINE WITH THE NEGATIVE DIRECTION OF X-AXIS
	  int          Type;            //TYPE OF CELL - 1. MIXED 2. ACTIVE 3. ISOLATED
	  long double  delV;            //VOLUME OF THE CELL - REQUIRED FOR BOOKKEPING DURING ADVECTION
     /********************************************************************************************/  

     /*****CELL ATTRIBUTES - MOMENTUM EQUATION***********/ 
      long double Press;
     /***************************************************/    

   } STRUCT_CELL_ATTRIB;

STRUCT_CELL_ATTRIB Cellattrib[NoOfRows][NoOfCols]; /*CREATING AN ARRAY OF STRUCTURES WHERE EACH ELEMENT OF THE ARRAY IS OF THE TYPE STRUCT_CELL_ATTRIB*/



/*X_VELOCITY[2][3] IS THE HORIZONTAL VELOCITY THROUGH THE LHS WALL OF THE CELL[2][3]
  Y_FLUX[2][3] IS THE VERTICAL FLUX THROUGH THE BOTTOM WALL OF THE CELL[2][3]

/****************WALL ATTRIBUTES***************************************************************/
long double  X_Velocity[NoOfRows][NoOfCols+1];  //X-COMPONENT OF VELOCITY AT THE LHS OF A CELL
long double  Y_Velocity[NoOfRows+1][NoOfCols];  //Y-COMPONENT OF VELOCITY AT THE BOTTOM OF A CELL
long double  X_Flux[NoOfRows][NoOfCols+1];      //THE X FLUX AT THE LHS OF A CELL
long double  Y_Flux[NoOfRows+1][NoOfCols];      //THE Y FLUX AT THE BOTTOM OF A CELL
/************************************************************************************************/

/****************GLOBAL VARIABLES****************/
long double del_sq;
long double del_sqBy2;
long double del2;
long double del3;
long double delBy2;
long double PiBy2;
long double PiBy4;
long double ThreePiBy2;
long double delTBydel;
long double delTdel;
/***********************************************/
