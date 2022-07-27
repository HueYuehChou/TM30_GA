#ifndef _TM30_H_
#define _TM30_H_
double TM30_calculate(double Tc,double Tc_uv_Slope_const_array[5],double LED_wavelength,
                     double LED_FWHM_L,double LED_FWHM_R,double Short_wavelength,double Short_FWHM,double Long_wavelength,
                     double Long_FWHM, double uv_interval,int TCS_width,int lower_limit,int wavelength_interval,int data_length,int data_width,
                     double CIE1931_deg10[data_length][data_width],double CIE1931_deg2[data_length][data_width],double TCS[data_length][TCS_width],double Daylight_SPD_1nm[data_length][4],double RfRg_array[2]);

double Y_normalize(int data_length,double Illuminant[data_length],double CIE1931[data_length][4]);

#endif
