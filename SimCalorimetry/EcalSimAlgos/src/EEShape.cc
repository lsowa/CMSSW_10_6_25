#include <cmath>

#include "SimCalorimetry/EcalSimAlgos/interface/EEShape.h"

EEShape::~EEShape()
{
}

EEShape::EEShape( bool   aSaveDerivative ) :
   EcalShapeBase( aSaveDerivative )
{
   buildMe() ;
}

double
EEShape::threshold() const
{
   return 0.00013 ; // magic number for barrel
}

void
EEShape::fillShape( EcalShapeBase::DVec& aVec ) const
{

   for( unsigned int i ( k1NSecBins ) ; i != k1NSecBinsTotal ; ++i )
   {
      aVec[i] = exp(2.39735 - 0.0151053* ((double)i+1.0));
//      shapeArray[i] = 10.82918*exp(-i/66.202 );
   }


   aVec[0] = 6.94068e-05 ; 
   aVec[1] = -5.03304e-05 ; 
   aVec[2] = -2.13404e-05 ; 
   aVec[3] = 6.017e-05 ; 
   aVec[4] = 2.01697e-05 ; 
   aVec[5] = 0.000114845 ; 
   aVec[6] = 2.13998e-05 ; 
   aVec[7] = 2.74476e-05 ; 
   aVec[8] = 5.2824e-05 ; 
   aVec[9] = 8.754e-05 ; 
   aVec[10] = 2.95346e-06 ; 
   aVec[11] = -7.58699e-05 ; 
   aVec[12] = -2.72224e-05 ; 
   aVec[13] = 3.10997e-06 ; 
   aVec[14] = -3.97771e-05 ; 
   aVec[15] = -1.06916e-05 ; 
   aVec[16] = -0.000113865 ; 
   aVec[17] = 6.05044e-05 ; 
   aVec[18] = -5.81202e-05 ; 
   aVec[19] = -6.58974e-06 ; 
   aVec[20] = 5.37494e-05 ; 
   aVec[21] = -0.000123729 ; 
   aVec[22] = 7.50938e-06 ; 
   aVec[23] = -1.35628e-05 ; 
   aVec[24] = 8.33725e-05 ; 
   aVec[25] = 3.19299e-05 ; 
   aVec[26] = -3.09232e-05 ; 
   aVec[27] = -7.0086e-05 ; 
   aVec[28] = 1.78937e-06 ; 
   aVec[29] = -2.20365e-05 ; 
   aVec[30] = 7.68054e-05 ; 
   aVec[31] = -2.5368e-05 ; 
   aVec[32] = 5.67291e-06 ; 
   aVec[33] = 5.87096e-05 ; 
   aVec[34] = -2.62771e-06 ; 
   aVec[35] = 4.31832e-05 ; 
   aVec[36] = 8.33616e-06 ; 
   aVec[37] = 7.27813e-05 ; 
   aVec[38] = 7.6159e-05 ; 
   aVec[39] = -1.60446e-05 ; 
   aVec[40] = -4.12127e-06 ; 
   aVec[41] = -5.93381e-05 ; 
   aVec[42] = 1.61444e-05 ; 
   aVec[43] = -5.49559e-05 ; 
   aVec[44] = 5.55254e-05 ; 
   aVec[45] = 3.32251e-05 ; 
   aVec[46] = -3.15897e-05 ; 
   aVec[47] = 7.86588e-05 ; 
   aVec[48] = -2.9704e-05 ; 
   aVec[49] = 5.66838e-05 ; 
   aVec[50] = 2.85281e-05 ; 
   aVec[51] = -3.02436e-05 ; 
   aVec[52] = -4.16265e-05 ; 
   aVec[53] = -1.63191e-05 ; 
   aVec[54] = 6.61193e-05 ; 
   aVec[55] = 9.23766e-05 ; 
   aVec[56] = 6.68903e-05 ; 
   aVec[57] = -3.20994e-05 ; 
   aVec[58] = 0.00011082 ; 
   aVec[59] = -4.07997e-05 ; 
   aVec[60] = -8.29046e-06 ; 
   aVec[61] = -7.42197e-05 ; 
   aVec[62] = -1.64386e-05 ; 
   aVec[63] = 1.02508e-05 ; 
   aVec[64] = 7.10995e-06 ; 
   aVec[65] = -5.87486e-05 ; 
   aVec[66] = -0.000101201 ; 
   aVec[67] = 1.62003e-05 ; 
   aVec[68] = -2.53093e-05 ; 
   aVec[69] = 2.65239e-05 ; 
   aVec[70] = -2.68722e-05 ; 
   aVec[71] = -4.02001e-05 ; 
   aVec[72] = 5.0674e-05 ; 
   aVec[73] = -1.75884e-05 ; 
   aVec[74] = 4.7902e-05 ; 
   aVec[75] = -1.01079e-05 ; 
   aVec[76] = 1.08427e-05 ; 
   aVec[77] = -0.000112906 ; 
   aVec[78] = 3.33076e-05 ; 
   aVec[79] = 0.000181201 ; 
   aVec[80] = 0.000426875 ; 
   aVec[81] = 0.00114222 ; 
   aVec[82] = 0.00237804 ; 
   aVec[83] = 0.00541858 ; 
   aVec[84] = 0.0089021 ; 
   aVec[85] = 0.0149157 ; 
   aVec[86] = 0.0231397 ; 
   aVec[87] = 0.0344671 ; 
   aVec[88] = 0.0471013 ; 
   aVec[89] = 0.0625517 ; 
   aVec[90] = 0.0857351 ; 
   aVec[91] = 0.108561 ; 
   aVec[92] = 0.133481 ; 
   aVec[93] = 0.163557 ; 
   aVec[94] = 0.200243 ; 
   aVec[95] = 0.225919 ; 
   aVec[96] = 0.269213 ; 
   aVec[97] = 0.302929 ; 
   aVec[98] = 0.342722 ; 
   aVec[99] = 0.378522 ; 
   aVec[100] = 0.436563 ; 
   aVec[101] = 0.467581 ; 
   aVec[102] = 0.510133 ; 
   aVec[103] = 0.550063 ; 
   aVec[104] = 0.583509 ; 
   aVec[105] = 0.619187 ; 
   aVec[106] = 0.653245 ; 
   aVec[107] = 0.686101 ; 
   aVec[108] = 0.721178 ; 
   aVec[109] = 0.745129 ; 
   aVec[110] = 0.774163 ; 
   aVec[111] = 0.799011 ; 
   aVec[112] = 0.822177 ; 
   aVec[113] = 0.838315 ; 
   aVec[114] = 0.858847 ; 
   aVec[115] = 0.875559 ; 
   aVec[116] = 0.891294 ; 
   aVec[117] = 0.90537 ; 
   aVec[118] = 0.919617 ; 
   aVec[119] = 0.930632 ; 
   aVec[120] = 0.936216 ; 
   aVec[121] = 0.947739 ; 
   aVec[122] = 0.955306 ; 
   aVec[123] = 0.961876 ; 
   aVec[124] = 0.968124 ; 
   aVec[125] = 0.97327 ; 
   aVec[126] = 0.977513 ; 
   aVec[127] = 0.984885 ; 
   aVec[128] = 0.986497 ; 
   aVec[129] = 0.990039 ; 
   aVec[130] = 0.994798 ; 
   aVec[131] = 0.994884 ; 
   aVec[132] = 0.99795 ; 
   aVec[133] = 0.99834 ; 
   aVec[134] = 0.999607 ; 
   aVec[135] = 1 ; 
   aVec[136] = 0.999047 ; 
   aVec[137] = 0.998745 ; 
   aVec[138] = 0.999219 ; 
   aVec[139] = 0.99814 ; 
   aVec[140] = 0.995082 ; 
   aVec[141] = 0.992449 ; 
   aVec[142] = 0.990418 ; 
   aVec[143] = 0.985032 ; 
   aVec[144] = 0.982308 ; 
   aVec[145] = 0.978696 ; 
   aVec[146] = 0.975656 ; 
   aVec[147] = 0.971027 ; 
   aVec[148] = 0.964811 ; 
   aVec[149] = 0.959428 ; 
   aVec[150] = 0.95096 ; 
   aVec[151] = 0.947428 ; 
   aVec[152] = 0.9419 ; 
   aVec[153] = 0.933223 ; 
   aVec[154] = 0.926482 ; 
   aVec[155] = 0.922172 ; 
   aVec[156] = 0.912777 ; 
   aVec[157] = 0.907388 ; 
   aVec[158] = 0.897289 ; 
   aVec[159] = 0.891889 ; 
   aVec[160] = 0.882056 ; 
   aVec[161] = 0.873382 ; 
   aVec[162] = 0.865442 ; 
   aVec[163] = 0.860032 ; 
   aVec[164] = 0.85202 ; 
   aVec[165] = 0.841013 ; 
   aVec[166] = 0.833802 ; 
   aVec[167] = 0.825259 ; 
   aVec[168] = 0.815013 ; 
   aVec[169] = 0.807465 ; 
   aVec[170] = 0.799428 ; 
   aVec[171] = 0.792165 ; 
   aVec[172] = 0.783088 ; 
   aVec[173] = 0.773392 ; 
   aVec[174] = 0.764982 ; 
   aVec[175] = 0.752174 ; 
   aVec[176] = 0.746487 ; 
   aVec[177] = 0.737678 ; 
   aVec[178] = 0.727396 ; 
   aVec[179] = 0.718692 ; 
   aVec[180] = 0.712737 ; 
   aVec[181] = 0.702738 ; 
   aVec[182] = 0.69559 ; 
   aVec[183] = 0.684389 ; 
   aVec[184] = 0.677989 ; 
   aVec[185] = 0.667643 ; 
   aVec[186] = 0.659009 ; 
   aVec[187] = 0.650217 ; 
   aVec[188] = 0.644479 ; 
   aVec[189] = 0.636017 ; 
   aVec[190] = 0.625257 ; 
   aVec[191] = 0.618507 ; 
   aVec[192] = 0.609798 ; 
   aVec[193] = 0.600097 ; 
   aVec[194] = 0.592788 ; 
   aVec[195] = 0.584895 ; 
   aVec[196] = 0.578228 ; 
   aVec[197] = 0.569299 ; 
   aVec[198] = 0.560576 ; 
   aVec[199] = 0.552404 ; 
   aVec[200] = 0.541405 ; 
   aVec[201] = 0.536271 ; 
   aVec[202] = 0.528734 ; 
   aVec[203] = 0.519813 ; 
   aVec[204] = 0.512264 ; 
   aVec[205] = 0.507001 ; 
   aVec[206] = 0.49828 ; 
   aVec[207] = 0.492416 ; 
   aVec[208] = 0.483181 ; 
   aVec[209] = 0.477907 ; 
   aVec[210] = 0.469623 ; 
   aVec[211] = 0.462528 ; 
   aVec[212] = 0.455099 ; 
   aVec[213] = 0.45055 ; 
   aVec[214] = 0.443576 ; 
   aVec[215] = 0.435364 ; 
   aVec[216] = 0.429789 ; 
   aVec[217] = 0.422724 ; 
   aVec[218] = 0.415621 ; 
   aVec[219] = 0.409469 ; 
   aVec[220] = 0.40401 ; 
   aVec[221] = 0.398121 ; 
   aVec[222] = 0.391079 ; 
   aVec[223] = 0.384414 ; 
   aVec[224] = 0.378214 ; 
   aVec[225] = 0.369851 ; 
   aVec[226] = 0.365966 ; 
   aVec[227] = 0.359865 ; 
   aVec[228] = 0.353505 ; 
   aVec[229] = 0.347899 ; 
   aVec[230] = 0.343829 ; 
   aVec[231] = 0.337585 ; 
   aVec[232] = 0.333089 ; 
   aVec[233] = 0.326289 ; 
   aVec[234] = 0.322249 ; 
   aVec[235] = 0.316079 ; 
   aVec[236] = 0.31061 ; 
   aVec[237] = 0.305426 ; 
   aVec[238] = 0.301885 ; 
   aVec[239] = 0.296753 ; 
   aVec[240] = 0.290931 ; 
   aVec[241] = 0.286877 ; 
   aVec[242] = 0.281831 ; 
   aVec[243] = 0.276633 ; 
   aVec[244] = 0.272283 ; 
   aVec[245] = 0.268069 ; 
   aVec[246] = 0.26399 ; 
   aVec[247] = 0.258457 ; 
   aVec[248] = 0.253549 ; 
   aVec[249] = 0.249493 ; 
}
