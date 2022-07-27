#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "Excel_fileopen.h"
#define MAX_LINE_SIZE 3500



int CIE_Excelopen_deg2(int argc, const char * argv[],int init_value ,int interval,int data_length,int data_width,double CIE1931[data_length][data_width]){


    double temp3[535][11];

    char file_name[] = "all_1nm_data.csv";
    FILE *fp;
    fp = fopen(file_name, "r");

    if (!fp) {
        fprintf(stderr, "failed to open file for reading\n");
        return 1;
    }

    char line[MAX_LINE_SIZE];
    char *result = NULL;
    int i=0;

    while(fgets(line, MAX_LINE_SIZE, fp) != NULL) {

        result = strtok(line, ",");
        int j=0;

        while( result != NULL ) {


            temp3[i][j]= atof(result);

            j++;


            result = strtok(NULL, ",");

        }

        i++;
    }
    int k=0;
   for(int i=init_value-296;i<data_length*interval+(init_value-296);i+=interval){
            CIE1931[i-(init_value-296)-(interval-1)*k][0]=temp3[i][0];
            CIE1931[i-(init_value-296)-(interval-1)*k][1]=temp3[i][5];
            CIE1931[i-(init_value-296)-(interval-1)*k][2]=temp3[i][6];
            CIE1931[i-(init_value-296)-(interval-1)*k][3]=temp3[i][7];
            k++;

    }


    fclose (fp);

}
int CIE_Excelopen_deg10(int argc, const char * argv[],int init_value ,int interval,int data_length,int data_width,double CIE1931[data_length][data_width]){


    double temp3[535][11];

    char file_name[] = "all_1nm_data.csv";
    FILE *fp;
    fp = fopen(file_name, "r");

    if (!fp) {
        fprintf(stderr, "failed to open file for reading\n");
        return 1;
    }

    char line[MAX_LINE_SIZE];
    char *result = NULL;
    int i=0;

    while(fgets(line, MAX_LINE_SIZE, fp) != NULL) {

        result = strtok(line, ",");
        int j=0;

        while( result != NULL ) {


            temp3[i][j]= atof(result);

            j++;


            result = strtok(NULL, ",");

        }

        i++;
    }
    int k=0;
   for(int i=init_value-296;i<data_length*interval+(init_value-296);i+=interval){
            CIE1931[i-(init_value-296)-(interval-1)*k][0]=temp3[i][0];
            CIE1931[i-(init_value-296)-(interval-1)*k][1]=temp3[i][8];
            CIE1931[i-(init_value-296)-(interval-1)*k][2]=temp3[i][9];
            CIE1931[i-(init_value-296)-(interval-1)*k][3]=temp3[i][10];
            k++;

    }


    fclose (fp);

}

int TCS_Excelopen(int argc, const char * argv[],int init_value ,int interval,int data_length,int TCS_width,double TCS[data_length][TCS_width]) {



    double temp1[402][100];

    char file_name1[] = "TCS_1nm_data.csv";
    FILE *fp1;
    fp1 = fopen(file_name1, "r");

    if (!fp1) {
        fprintf(stderr, "failed to open file for reading\n");
        return 1;
    }

    char line[MAX_LINE_SIZE];
    char *result = NULL;
    int i=0;

    while(fgets(line, MAX_LINE_SIZE, fp1) != NULL) {

        result = strtok(line, ",");
        int j=0;

        while( result != NULL ) {


            temp1[i][j]= atof(result);

            j++;


            result = strtok(NULL, ",");

        }

        i++;
    }
    int k=0;
   for(int i=0; i<401;i++){
    for(int j=0;j<99;j++){

        TCS[i][j]= temp1[i+1][j+1];
    }
   }

    fclose (fp1);

}
int Daylight_Excelopen(int argc, const char * argv[],int init_value ,int interval,int Daylght_length,int Daylght_width,double Daylight_SPD_10nm[Daylght_length][Daylght_width]) {



    double temp1[Daylght_length+1][Daylght_width];

    char file_name1[] = "DaylightCIECalculator.csv";
    FILE *fp1;
    fp1 = fopen(file_name1, "r");

    if (!fp1) {
        fprintf(stderr, "failed to open file for reading\n");
        return 1;
    }

    char line[MAX_LINE_SIZE];
    char *result = NULL;
    int i=0;

    while(fgets(line, MAX_LINE_SIZE, fp1) != NULL) {

        result = strtok(line, ",");
        int j=0;

        while( result != NULL ) {


            temp1[i][j]= atof(result);

            j++;


            result = strtok(NULL, ",");

        }

        i++;
    }
     int k=0;

      for(int i=init_value-379;i<Daylght_length*interval+(init_value-379);i+=interval){
        for(int j=0;j<Daylght_width;j++){

            Daylight_SPD_10nm[i-(init_value-379)-(interval-1)*k][j]=temp1[i][j];
        }


            k++;

    }


    fclose (fp1);

}
int Test_Excelopen(int argc, const char * argv[],int init_value ,int interval,int data_length,int data_width,double XYZ[6][100]) {


    double temp3[6][100];

    char file_name[] = "test.csv";
    FILE *fp;
    fp = fopen(file_name, "r");

    if (!fp) {
        fprintf(stderr, "failed to open file for reading\n");
        return 1;
    }

    char line[MAX_LINE_SIZE];
    char *result = NULL;
    int i=0;

    while(fgets(line, MAX_LINE_SIZE, fp) != NULL) {

        result = strtok(line, ",");
        int j=0;

        while( result != NULL ) {


            temp3[i][j]= atof(result);

            j++;


            result = strtok(NULL, ",");

        }

        i++;
    }
   /*for(int i=0;i<100;i++){
        for(int j=0 ; j<6; j++){
            printf("%f\t",temp3[j][i]);
        }
        printf("\n");
    }
*/

   for(int i=0;i<100;i++){
            XYZ[0][i]=temp3[0][i];
            XYZ[1][i]=temp3[1][i];
            XYZ[2][i]=temp3[2][i];
            XYZ[3][i]=temp3[3][i];
            XYZ[4][i]=temp3[4][i];
            XYZ[5][i]=temp3[5][i];


    }
  /* for(int i=0;i<100;i++){
        for(int j=0 ; j<6; j++){
            printf("%f\t",XYZ[j][i]);
        }
        printf("\n");
    }
     printf("\n");
    printf("\n");
    printf("\n");

*/
    fclose (fp);

}

