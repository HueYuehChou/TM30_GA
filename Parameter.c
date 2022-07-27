#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include "Parameter.h"


double ParameterRf(double Tristimulus_temp[3],double Tristimulus_array_TCS[3][99],double JcaMbM[3][100]){

    double LA=100.0, Yb=20.0 ,Did=1.0 , F=1.0, c=0.69, Nc=1.0, D=1.0,Aw,pi=3.14159265;
    double MCAT02[3][3]={0.7328,0.4296,-0.1624,-0.7036,1.6975,0.0061,0.003,0.0136,0.9834};
    double invCAT02[3][3]={1.09612382083551,-0.278869000218287,0.182745179382773,0.454369041975359,0.473533154307412,0.0720978037172291,-0.00962760873842936,-0.00569803121611342,1.01532563995454};
    double MHPE[3][3]= {0.38971,0.68898,-0.07868,-0.22981,1.1834,0.04641,0.0,0.0,1.0};
    //double invMHPE[3][3]={1.9102,-1.1121,0.2019,0.371,0.6291,0.0,0.0,0.0,1.0};
    double k,FL,n,Nbb,Ncb,z;
    double RGB[3][100],RwGwBw[3],RcGcBc[3][100],RcwGcwBcw[3],XcYcZc[3][100],XcwYcwZcw[3],R_G_B_prime[3][100],Rw_Gw_Bw_prime[3],RaBaGa_prime[3][100],Raw_Gaw_Baw_prime[3]
    ,a[100],b[100],h[100],A[100],J[100],Q[100],t[100],C[100],M[100],s[100],Jc[100],Mc[100],aM[100],bM[100],e[100];

    Tristimulus_temp[0]=Tristimulus_temp[0]/Tristimulus_temp[1]*100.0;
    Tristimulus_temp[2]=Tristimulus_temp[2]/Tristimulus_temp[1]*100.0;
    Tristimulus_temp[1]=100.0;

    k=1/(5*LA+1);
    //printf("k=%f\n",k);
    FL=0.2*k*k*k*k*5*LA+0.1*(1-k*k*k*k)*(1-k*k*k*k)*pow((5*LA),1.0/3.0);
   // printf("FL=%f\n",FL);
    n=Yb/Tristimulus_temp[1];
    Nbb=0.725*pow((1/n),0.2);
   // printf("Nbb=%f\n",Nbb);
    Ncb=0.725*pow((1/n),0.2);
    z=1.48+pow(n,0.5);

    for(int i=0;i<3;i++){

            RGB[i][0]=MCAT02[i][0]*Tristimulus_temp[0]+MCAT02[i][1]*Tristimulus_temp[1]+MCAT02[i][2]*Tristimulus_temp[2];
            RwGwBw[i]=MCAT02[i][0]*Tristimulus_temp[0]+MCAT02[i][1]*Tristimulus_temp[1]+MCAT02[i][2]*Tristimulus_temp[2];
            RcGcBc[i][0]=(D*(Tristimulus_temp[1]/RwGwBw[i])+1-D)*RGB[i][0];
            RcwGcwBcw[i]=(D*(Tristimulus_temp[1]/RwGwBw[i])+1-D)*RGB[i][0];
            //printf("%f\t\n\n",RwGwBw[i]);
    }
    for(int i=1;i<100;i++){
            for(int j=0; j<3;j++){

        RGB[j][i]=MCAT02[j][0]*Tristimulus_array_TCS[0][i-1]+MCAT02[j][1]*Tristimulus_array_TCS[1][i-1]+MCAT02[j][2]*Tristimulus_array_TCS[2][i-1];
        RcGcBc[j][i]=(D*(Tristimulus_temp[1]/RwGwBw[j])+1-D)*RGB[j][i];


            }
    }



    for(int i=0;i<100; i++){
        for(int j=0; j<3; j++){

            XcYcZc[j][i]=invCAT02[j][0]*RcGcBc[0][i]+invCAT02[j][1]*RcGcBc[1][i]+invCAT02[j][2]*RcGcBc[2][i];

        }
    }






    for(int i=0;i<3;i++){
        XcwYcwZcw[i]=invCAT02[i][0]*RcwGcwBcw[0]+invCAT02[i][1]*RcwGcwBcw[1]+invCAT02[i][2]*RcwGcwBcw[2];

    }




      for(int i=0;i<100; i++){
        for(int j=0; j<3; j++){

            R_G_B_prime[j][i]=MHPE[j][0]*XcYcZc[0][i]+MHPE[j][1]*XcYcZc[1][i]+MHPE[j][2]*XcYcZc[2][i];


        }

    }

    for(int i=0;i<3;i++){
        Rw_Gw_Bw_prime[i]=MHPE[i][0]*XcwYcwZcw[i]+MHPE[i][1]*XcwYcwZcw[i]+MHPE[i][2]*XcwYcwZcw[i];

    }

   for(int i=0;i<100; i++){
        for(int j=0; j<3; j++){

            RaBaGa_prime[j][i]=((400*pow((FL*R_G_B_prime[j][i]/100.0),0.42))/(pow((FL*R_G_B_prime[j][i]/100.0),0.42)+27.13))+0.1;


        }

        A[i]=(2*RaBaGa_prime[0][i]+RaBaGa_prime[1][i]+(1.0/20.0)*RaBaGa_prime[2][i]-0.305)*Nbb;

    }


   for(int i=0;i<3;i++){
        Raw_Gaw_Baw_prime[i]=((400.0*pow(FL*(Rw_Gw_Bw_prime[i]/100.0),0.42))/(pow((FL*Rw_Gw_Bw_prime[i]/100.0),0.42)+27.13))+0.1;

    }


    Aw=(2.0*Raw_Gaw_Baw_prime[0]+Raw_Gaw_Baw_prime[1]+(1.0/20.0)*Raw_Gaw_Baw_prime[2]-0.305)*Nbb;




   for(int i=0; i<100;i++){
        a[i]=RaBaGa_prime[0][i]-12.0*RaBaGa_prime[1][i]/11.0+RaBaGa_prime[2][i]/11.0;
        b[i]=(1.0/9.0)*(RaBaGa_prime[0][i]+RaBaGa_prime[1][i]-2*RaBaGa_prime[2][i]);
        if(a[i]==0 && b[i]==0){

                h[i]=0;

        }
        else if(b[i]<0){

                h[i]=360.0+(360.0/(2.0*pi))*atan(b[i]/a[i]);

        }
        else{
                h[i]=(360.0/(2.0*pi))*atan(b[i]/a[i]);

        }
        if(h[i]<0){
            h[i]=h[i]+180.0;
        }
         else if(h[i]>360){

            h[i]=h[i]-180.0;
        }
        J[i]=100.0*pow((A[i]/Aw),(c*z));
        Q[i]=(4.0/c)*(pow(J[i]/100,0.5))*(Aw+4.0)*(pow(FL,0.25));
        e[i]=((12500.0/13.0)*Nc*Ncb)*(cos((h[i]*pi/180.0)+2.0)+3.8);





    }
   // printf("\n");




    for(int i=0; i<100;i++){

        t[i]=(e[i]*pow((a[i]*a[i]+b[i]*b[i]),0.5))/(RaBaGa_prime[0][i]+RaBaGa_prime[1][i]+(21.0/20.0)*RaBaGa_prime[2][i]);

        C[i]=(pow(t[i],0.9))*(pow((J[i]/100.0),0.5))*(pow(1.64-pow(0.29,n),0.73));


        M[i]=C[i]*pow(FL,0.25);

        s[i]=100*pow((M[i]/Q[i]),0.5);

        Jc[i]=(1+100.0*0.007)*J[i]/(1+0.007*J[i]);

        Mc[i]=(1/0.0228)*log(1+0.0228*M[i]);


        aM[i]=Mc[i]*cos(h[i]*pi/180.0);

        bM[i]=Mc[i]*sin(h[i]*pi/180.0);

         //printf("%f\n",Jc[i]);

        JcaMbM[0][i]=Jc[i];
        JcaMbM[1][i]=aM[i];
        JcaMbM[2][i]=bM[i];
      /*  for(int j=0 ;j<3;j++){

            printf("%f\t",JcaMbM[j][i]);
        }
            printf("\n");*/

    }


}
