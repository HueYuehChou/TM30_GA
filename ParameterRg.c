#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include "ParameterRg.h"

double ParameterRg(double JcaMbMT[3][100],double JcaMbMR[3][100]){


    double pi=3.14159265,A0=0.0,A1=0.0;
    double theta[99],theta_degree[99],arbrab[5][16]={0},diffar[16],meanbr[16],diffa[16],meanb[16];
    int bin[99];


    for(int i=0;i<99;i++){

            if(JcaMbMR[1][i+1]<0 && JcaMbMR[2][i+1]>0 ){

                theta[i]=atan(JcaMbMR[2][i+1]/JcaMbMR[1][i+1])+pi;
            }
            else if(JcaMbMR[1][i+1]<0 && JcaMbMR[2][i+1]<0){

                theta[i]=atan(JcaMbMR[2][i+1]/JcaMbMR[1][i+1])+pi;

            }
            else if(JcaMbMR[1][i+1]>0 && JcaMbMR[2][i+1]<0){

                theta[i]=atan(JcaMbMR[2][i+1]/JcaMbMR[1][i+1])+2*pi;

            }
            else {

                theta[i]=atan(JcaMbMR[2][i+1]/JcaMbMR[1][i+1]);

            }
        theta_degree[i]=180.0/pi* theta[i] ;
        bin[i] = floor((theta[i]/2.0)/(pi)*16.0)+1;

        //printf("%d\n",bin[i]);

        switch(bin[i]){

        case 1 :
            arbrab[0][0]+=JcaMbMR[1][i+1];
            arbrab[1][0]+=JcaMbMR[2][i+1];
            arbrab[2][0]+=JcaMbMT[1][i+1];
            arbrab[3][0]+=JcaMbMT[2][i+1];
            arbrab[4][0]+=1.0;
            break;
        case 2 :
            arbrab[0][1]+=JcaMbMR[1][i+1];
            arbrab[1][1]+=JcaMbMR[2][i+1];
            arbrab[2][1]+=JcaMbMT[1][i+1];
            arbrab[3][1]+=JcaMbMT[2][i+1];
            arbrab[4][1]+=1.0;
            break;
        case 3 :
            arbrab[0][2]+=JcaMbMR[1][i+1];
            arbrab[1][2]+=JcaMbMR[2][i+1];
            arbrab[2][2]+=JcaMbMT[1][i+1];
            arbrab[3][2]+=JcaMbMT[2][i+1];
            arbrab[4][2]+=1.0;
            break;
         case 4 :
            arbrab[0][3]+=JcaMbMR[1][i+1];
            arbrab[1][3]+=JcaMbMR[2][i+1];
            arbrab[2][3]+=JcaMbMT[1][i+1];
            arbrab[3][3]+=JcaMbMT[2][i+1];
            arbrab[4][3]+=1.0;
            break;
         case 5 :
            arbrab[0][4]+=JcaMbMR[1][i+1];
            arbrab[1][4]+=JcaMbMR[2][i+1];
            arbrab[2][4]+=JcaMbMT[1][i+1];
            arbrab[3][4]+=JcaMbMT[2][i+1];
            arbrab[4][4]+=1.0;
            break;
         case 6 :
            arbrab[0][5]+=JcaMbMR[1][i+1];
            arbrab[1][5]+=JcaMbMR[2][i+1];
            arbrab[2][5]+=JcaMbMT[1][i+1];
            arbrab[3][5]+=JcaMbMT[2][i+1];
            arbrab[4][5]+=1.0;
            break;
        case 7 :
            arbrab[0][6]+=JcaMbMR[1][i+1];
            arbrab[1][6]+=JcaMbMR[2][i+1];
            arbrab[2][6]+=JcaMbMT[1][i+1];
            arbrab[3][6]+=JcaMbMT[2][i+1];
            arbrab[4][6]+=1.0;
            break;
        case 8 :
            arbrab[0][7]+=JcaMbMR[1][i+1];
            arbrab[1][7]+=JcaMbMR[2][i+1];
            arbrab[2][7]+=JcaMbMT[1][i+1];
            arbrab[3][7]+=JcaMbMT[2][i+1];
            arbrab[4][7]+=1.0;
            break;
        case 9 :
            arbrab[0][8]+=JcaMbMR[1][i+1];
            arbrab[1][8]+=JcaMbMR[2][i+1];
            arbrab[2][8]+=JcaMbMT[1][i+1];
            arbrab[3][8]+=JcaMbMT[2][i+1];
            arbrab[4][8]+=1.0;
            break;
        case 10 :
            arbrab[0][9]+=JcaMbMR[1][i+1];
            arbrab[1][9]+=JcaMbMR[2][i+1];
            arbrab[2][9]+=JcaMbMT[1][i+1];
            arbrab[3][9]+=JcaMbMT[2][i+1];
            arbrab[4][9]+=1.0;
            break;
         case 11 :
            arbrab[0][10]+=JcaMbMR[1][i+1];
            arbrab[1][10]+=JcaMbMR[2][i+1];
            arbrab[2][10]+=JcaMbMT[1][i+1];
            arbrab[3][10]+=JcaMbMT[2][i+1];
            arbrab[4][10]+=1.0;
            break;
        case 12 :
            arbrab[0][11]+=JcaMbMR[1][i+1];
            arbrab[1][11]+=JcaMbMR[2][i+1];
            arbrab[2][11]+=JcaMbMT[1][i+1];
            arbrab[3][11]+=JcaMbMT[2][i+1];
            arbrab[4][11]+=1.0;
            break;
        case 13 :
            arbrab[0][12]+=JcaMbMR[1][i+1];
            arbrab[1][12]+=JcaMbMR[2][i+1];
            arbrab[2][12]+=JcaMbMT[1][i+1];
            arbrab[3][12]+=JcaMbMT[2][i+1];
            arbrab[4][12]+=1.0;
            break;
        case 14 :
            arbrab[0][13]+=JcaMbMR[1][i+1];
            arbrab[1][13]+=JcaMbMR[2][i+1];
            arbrab[2][13]+=JcaMbMT[1][i+1];
            arbrab[3][13]+=JcaMbMT[2][i+1];
            arbrab[4][13]+=1.0;
            break;
        case 15 :
            arbrab[0][14]+=JcaMbMR[1][i+1];
            arbrab[1][14]+=JcaMbMR[2][i+1];
            arbrab[2][14]+=JcaMbMT[1][i+1];
            arbrab[3][14]+=JcaMbMT[2][i+1];
            arbrab[4][14]+=1.0;
            break;
        case 16 :
            arbrab[0][15]+=JcaMbMR[1][i+1];
            arbrab[1][15]+=JcaMbMR[2][i+1];
            arbrab[2][15]+=JcaMbMT[1][i+1];
            arbrab[3][15]+=JcaMbMT[2][i+1];
            arbrab[4][15]+=1.0;
            break;
            }


        }
        for(int i=0 ; i<16;i++){
                for(int j=0;j<4;j++){
            arbrab[j][i]=arbrab[j][i]/arbrab[4][i];
            //printf("%f\t",arbrab[j][i]);
                }
             //   printf("\n");
        }

        for(int i=0 ;i < 16;i++){
        if(i==15){

            diffar[i]=arbrab[0][0]-arbrab[0][i];

            meanbr[i]=(arbrab[1][0]+arbrab[1][i])/2.0;

            diffa[i]=arbrab[2][0]-arbrab[2][i];

            meanb[i]=(arbrab[3][0]+arbrab[3][i])/2.0;

        }

        else{

            diffar[i]=arbrab[0][i+1]-arbrab[0][i];

            meanbr[i]=(arbrab[1][i+1]+arbrab[1][i])/2.0;

            diffa[i]=arbrab[2][i+1]-arbrab[2][i];

            meanb[i]=(arbrab[3][i+1]+arbrab[3][i])/2.0;
        }
        //printf("diffa[i]=%f\n",diffa[i]);
        //printf("meanb[i]=%f\n",meanb[i]);
        A0+=diffar[i]*meanbr[i];

        A1+=diffa[i]*meanb[i];

        }
        //printf("A1=%f\n",A1);

      double Rg= A1/A0*100.0;
      //printf("Rg=%lf\n",Rg);

      return Rg;


    }








