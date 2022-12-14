#include <stdio.h>
#include <stdlib.h>
#include<math.h>
//#include"Gene_Algorithm.h"
#include "Excel_fileopen.h"
#include "Linear_interpolation.h"
#include "Parameter.h"
#include "ParameterRg.h"

int data_length,data_width=4,TCS_width,wavelength_upper_limit,wavelength_lower_limit,wavelength_interval,Tc_lower,Tc_upper,Tc_interval,Tc_point;

int main(int argc, const char * argv[]) {
    srand(time(NULL));

    printf("Key the input of wavelength upper limit(nm)\n");
    scanf("%d",&wavelength_upper_limit);
    printf("Key the input of wavelength lower limit(nm)\n");
    scanf("%d",&wavelength_lower_limit);
    printf("Key the input of wavelength interval(nm)\n");
    scanf("%d",&wavelength_interval);
    data_length=1+(wavelength_upper_limit-wavelength_lower_limit)/wavelength_interval;
    printf("data length=%d\n",data_length);
    printf("Key the number of TCS \n");
    scanf("%d",&TCS_width);
    printf("Key the lower Tc \n");
    scanf("%d",&Tc_lower);
    printf("Key the upper Tc \n");
    scanf("%d",&Tc_upper);
    printf("Key the interval of Tc \n");
    scanf("%d",&Tc_interval);
    printf("\n\n\n");
    if(wavelength_upper_limit>830||wavelength_lower_limit<360||wavelength_interval>11||data_length<0){
        printf("error\n");
        return 0;

    }
        if(Tc_lower<2000 || Tc_upper>30000 || Tc_interval<100 || Tc_interval >1000 ){
        printf("error\n");
        return 0;

    }
    Tc_point= (Tc_upper-Tc_lower)/Tc_interval+1;


     double CIE1931_deg10[data_length][data_width],CIE1931_deg2[data_length][data_width],TCS[data_length][TCS_width],Daylight_SPD_10nm[41][4],Daylight_SPD_1nm[data_length][4];
     double XYZ[6][100];
    CIE_Excelopen_deg10(argc, argv,wavelength_lower_limit,wavelength_interval,data_length,data_width,CIE1931_deg10);
    CIE_Excelopen_deg2(argc, argv,wavelength_lower_limit,wavelength_interval,data_length,data_width,CIE1931_deg2);
    TCS_Excelopen(argc, argv,wavelength_lower_limit,wavelength_interval,data_length,TCS_width,TCS);
    Daylight_Excelopen(argc, argv,wavelength_lower_limit,1,41,4,Daylight_SPD_10nm);
    //Test_Excelopen(argc, argv,wavelength_lower_limit,1,41,4,XYZ);


    /*for(int i=0; i<401; i++){
        for(int j=0 ; j<4;j++){
           printf("%f\t",CIE1931_deg10[i][j]);
        }
        printf("\n");

    }
     printf("\n\n\n\n\n");
*/


for(int j=0;j<40;j++){
        double k=0;
    for(int  i=0;i<11;i++){
    Daylight_SPD_1nm[10*j+i][0]=Daylight_SPD_10nm[j][0]+k;
    Daylight_SPD_1nm[10*j+i][1]=Linear_interpolation(k,Daylight_SPD_10nm[j][0],Daylight_SPD_10nm[j][1],Daylight_SPD_10nm[j+1][0],Daylight_SPD_10nm[j+1][1]);
    Daylight_SPD_1nm[10*j+i][2]=Linear_interpolation(k,Daylight_SPD_10nm[j][0],Daylight_SPD_10nm[j][2],Daylight_SPD_10nm[j+1][0],Daylight_SPD_10nm[j+1][2]);
    Daylight_SPD_1nm[10*j+i][3]=Linear_interpolation(k,Daylight_SPD_10nm[j][0],Daylight_SPD_10nm[j][3],Daylight_SPD_10nm[j+1][0],Daylight_SPD_10nm[j+1][3]);
    k+=1.0;
    }
}


    Genes_Algorithm(wavelength_lower_limit,data_length,data_width,wavelength_interval,TCS_width,CIE1931_deg10,CIE1931_deg2,TCS,Daylight_SPD_1nm,Tc_lower,Tc_interval,Tc_point);



    return 0;
}
