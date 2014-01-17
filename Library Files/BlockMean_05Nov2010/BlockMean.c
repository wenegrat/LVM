// BlockMean.c
// 2D block mean
// The mean of V*W elements along the 1st and 2nd dimensions is calculated.
// Y = BlockMean(X, V, W)
// INPUT:
//   X: UINT8 or DOUBLE array of any size.
//   V, W: Scalar numerical with integer value. Each element of the output is
//      the mean over V*W neighbouring elements of the input.
//      V and W are limited to 256 to limit memory usage.
//      A square V*V block is used, if W is omitted.
// OUTPUT:
//   Y: UINT8 or DOUBLE array, the 1st and 2nd dimensions are V and W times
//      shorter: [FLOOR(X / V) x FLOOR(Y / W) x (further dims...)]
//      If the size of the 1st or 2nd dimension is not a multiple of V and W,
//      the remaining elements at the end are skipped.
//      The empty matrix is replied for empty inputs or if the 1st or 2nd
//      dimension is shorter than V or W.
//
// NOTE: This is implemented for DOUBLE and UINT8, because I use it for
//   anti-aliasing of images stored as 3D RGB arrays.
//
// TESTING: Run uTest_BlockMean to check validity and speed.
//
// COMPILE:
//   mex -O BlockMean.c
// Linux: consider C99 comments:
//   mex -O CFLAGS="\$CFLAGS -std=c99" BlockMean.c
// Pre-compiled Mex: http://www.n-simon.de/mex
//
// Tested: Matlab 6.5, 7.7, 7.8, WinXP, 32bit
//         Compiler: LCC2.4/3.8, BCC5.5, OWC1.8, MSVC2008
// Assumed Compatibility: higher Matlab versions, Mac, Linux, 64bit
// Author: Jan Simon, Heidelberg, (C) 2009-2010 matlab.THISYEAR(a)nMINUSsimon.de

/*
% $JRev: R0D V:026 Sum:ZT2rdz6e8M5r Date:29-Oct-2010 22:18:13 $
% $License: BSD (see Docs\BSD_License.txt) $
% $UnitTest: uTest_BlockMean $
% $File: Tools\Mex\Source\BlockMean.c $
% History:
% 001: 20-Jul-2009 23:22, Generalized ChunkMeanRGB: Free trailing dimensions.
% 010: 26-Sep-2009 15:48, 40% faster with multiplicating 1/WindowSize**2.
% 013: 25-Oct-2009 21:27, Compatible with 64 bit arrays. More checks of inputs.
% 018: 12-Mar-2010 01:40, Rectangular blocks.
% 024: 29-Oct-2010 02:38, Omit access of FPU state under Linux.
%      Thanks to Dominik Brands!
*/

#include "mex.h"
#include <string.h>
#include <stdlib.h>
#include <float.h>

// Assume 32 bit addressing for Matlab 6.5:
// See MEX option "compatibleArrayDims" for MEX in Matlab >= 7.7.
#ifndef MWSIZE_MAX
#define mwSize  int32_T           // Defined in tmwtypes.h
#define mwIndex int32_T
#define MWSIZE_MAX MAX_int32_T
#endif

// Set the floating point precision from 64 to 80 bits under Windows:
// Linux runs with 80 bits as default. But even with 64 bit precision the
// results have at least the same accuracy as with Matlab's SUM. So the
// extended precision is not really necessary for BlockMean.
#if defined(__WINDOWS__) || defined(WIN32) || defined(_WIN32) || defined(_WIN64)
void C87to64bit(int Cmd);
#  define SET_PRECISION_80 C87to64bit(1);
#  define RESTORE_PRECISION C87to64bit(0);
#else  // Not windows - Linux or MacOS:
#  define SET_PRECISION_80  ;
#  define RESTORE_PRECISION ;
#endif

// Error messages do not contain the function name in Matlab 6.5! This is not
// necessary in Matlab 7, but it does not bother:
#define ERR_HEAD "*** BlockMean[mex]: "
#define ERR_ID   "JSimon:BlockMean:"

// Prototypes:
void Core_uint8(uint8_T *X, mwSize dim0, mwSize dim1,
                mwSize colSkip, mwSize rowSkip,
                mwSize V, mwSize W, mwSize numelTrail, uint8_T *Y);
void Core_double(double *X, mwSize dim0, mwSize dim1,
                 mwSize colSkip, mwSize rowSkip,
                 mwSize V, mwSize W, mwSize numelTrail, double *Y);
mxArray *EmptyAnswer(const mxArray *X);
int CheckInput(const mxArray *X, mwSize V, mwSize W);

// Main function ===============================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mwSize numelTrail, ndims, *dims, colSkip, rowSkip, V, W;
  
  // Proper inputs:
  if (nlhs > 1) {
     mexErrMsgIdAndTxt(ERR_ID   "BadNOutput",
                       ERR_HEAD "1 output allowed.");
  }
  
  if (nrhs == 2) {
    if (!mxIsNumeric(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 1) {
    mexErrMsgIdAndTxt(ERR_ID   "BadInputW",
                      ERR_HEAD "Window size W must be a scalar.");
    }
    V = (mwSize) mxGetScalar(prhs[1]);
    W = V;
    
  } else if (nrhs == 3) {
    if (!mxIsNumeric(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 1 ||
        !mxIsNumeric(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1) {
      mexErrMsgIdAndTxt(ERR_ID   "BadInputWV",
                        ERR_HEAD "Block sizes V and W must be scalar.");
    }
    V = (mwSize) mxGetScalar(prhs[1]);
    W = (mwSize) mxGetScalar(prhs[2]);
    
  } else {
    mexErrMsgIdAndTxt(ERR_ID   "BadNArgin",
                      ERR_HEAD "2 or 3 inputs required.");
  }
    
  // Reply empty DOUBLE or UINT8 matrix for empty or too small X:
  if (CheckInput(prhs[0], V, W)) {
    plhs[0] = EmptyAnswer(prhs[0]);
    return;
  } if (V * W == 1) {
    plhs[0] = mxDuplicateArray(prhs[0]);
    return;
  }
  
  // Dimensions of replied array:
  ndims = mxGetNumberOfDimensions(prhs[0]);
  if ((dims  = (mwSize *) mxMalloc(ndims * sizeof(mwSize))) == NULL) {
    mexErrMsgIdAndTxt(ERR_ID   "NoMemory",
                      ERR_HEAD "No memory for [dims].");
  }

  // Length of trailing dimensions ([M x N x ...] allowed):
  memcpy(dims, mxGetDimensions(prhs[0]), ndims * sizeof(mwSize));
  numelTrail = mxGetNumberOfElements(prhs[0]) / (dims[0] * dims[1]);
  
  // Reduced number of rows and columns:
  colSkip = dims[0] % V;   // Number of skipped final rows
  rowSkip = dims[1] % W;   // Number of skipped final columns
  dims[0] = dims[0] / V;   // Integer division with truncation
  dims[1] = dims[1] / W;

  // Create output and perform the calculations:
  if (mxIsDouble(prhs[0])) {
    plhs[0] = mxCreateNumericArray(ndims, dims, mxDOUBLE_CLASS, mxREAL);
    
    // Enable extended precision under Windows:
    SET_PRECISION_80
    
    Core_double(mxGetPr(prhs[0]), dims[0], dims[1], colSkip, rowSkip,
           V, W, numelTrail, mxGetPr(plhs[0]));
    
    RESTORE_PRECISION
   
  } else {   // UINT8:
    plhs[0] = mxCreateNumericArray(ndims, dims, mxUINT8_CLASS, mxREAL);
    Core_uint8((uint8_T *) mxGetData(prhs[0]), dims[0], dims[1],
               colSkip, rowSkip, V, W, numelTrail,
               (uint8_T *) mxGetData(plhs[0]));
  }
  
  mxFree(dims);
  
  return;
}

// =============================================================================
void Core_double(double *X, mwSize dim0, mwSize dim1,
                 mwSize colSkip, mwSize rowSkip,
                 mwSize V, mwSize W, mwSize numelTrail, double *Y)
{
  // Calculate mean over W rows and columns, DOUBLE input and output.
  double *colEnd, *rowEnd, *Xf;
  double invVW, *Yp, *Yf, *Yg;
  mwSize Block, Dim3Skip;
  
  invVW    = 1.0 / (double) (V * W);          // Q*invVW is faster than Q/VW
  Block    = W * (V * dim0 + colSkip);        // Block size with skipping
  Dim3Skip = rowSkip * (V * dim0 + colSkip);  // Skip to next 3rd dimension
  Yg       = Y + dim0 * dim1 * numelTrail;    // End of output
  
  // Loop over all trailing dimensions:
  while (Y < Yg) {
    // Loop through row blocks:
    rowEnd = X + Block * dim1;
    while (X < rowEnd) {
      // Loop over columns:
      Yf = Y + dim0;
      for (colEnd = X + Block; X < colEnd; X += colSkip) {
        // Loop over accumulator:
        for (Yp = Y; Yp < Yf; Yp++) {
          // Loop over [V] elements, sum in accumulator element:
          for (Xf = X + V; X < Xf; *Yp += *X++) ;  // empty loop
        }
      }
      
      // Normalize output by area of rectangle:
      for ( ; Y < Yf; *Y++ *= invVW) ;  // empty loop
    }
    
    X += Dim3Skip;
  }
  
  return;
}

// =============================================================================
void Core_uint8(uint8_T *X, mwSize dim0, mwSize dim1,
               mwSize colSkip, mwSize rowSkip,
               mwSize V, mwSize W, mwSize numelTrail, uint8_T *Y)
{
  // Calculate mean over W rows and columns, UINT8 input and output.
  // The accumulator vector [c] as UINT32 allows a window with side length of
  // 65536. If W would be limited to 16, a UINT16 accumulator would be enough,
  // but the acceleration is tiny.
  uint8_T  *colEnd, *rowEnd, *Xf;
  uint32_T VW, *c, *cP, *cEnd, HalfVW;
  uint32_T Block, Dim3Skip;
  mwIndex iTrail;
  
  if ((c = (uint32_T *) mxMalloc(dim0 * sizeof(uint32_T))) == NULL) {
    mexErrMsgIdAndTxt(ERR_ID   "NoMemory",
                      ERR_HEAD "No memory for accumulator.");
  }
  cEnd = c + dim0;
  
  Block    = W * (V * dim0 + colSkip);
  Dim3Skip = rowSkip * (V * dim0 + colSkip);
  
  VW     = V * W;   // No loss of bits: W <= 256
  HalfVW = VW / 2;
          
  for (iTrail = 0; iTrail < numelTrail; iTrail++) {
    rowEnd = X + Block * dim1;
    while (X < rowEnd) {
      // Initialize accumulator with V*W/2 for rounding:
      for (cP = c; cP < cEnd; *cP++ = HalfVW) ; // empty loop
        
      for (colEnd = X + Block; X < colEnd; X += colSkip) {
        for (cP = c; cP < cEnd; cP++) {
          for (Xf = X + V; X < Xf; *cP += *X++) ;  // empty loop
        }
      }
      
      // Y = uint8(floor(0.5 + sum(sum(XM, 1), 3) ./ (V * W)));
      for (cP = c; cP < cEnd; *Y++ = (uint8_T) (*cP++ / VW)) ;  // empty loop
    }
    
    X += Dim3Skip;
  }  // end: for iTrail
  
  mxFree(c);
  
  return;
}

// =============================================================================
int CheckInput(const mxArray *X, mwSize V, mwSize W)
{
  const mwSize *dims;
  
  // Check range of W (otherwise overflow of UINT32 buffer for UINT8 input):
  if (W < 1 || W > 256 || V < 1 || V > 256) {
    mexErrMsgIdAndTxt(ERR_ID   "BadWindowSize",
                      ERR_HEAD "Window size is out of range: 1 <= V <= 256.");
  }
  
  // X must be UINT8 or DOUBLE and needs 3 dimensions:
  if (!(mxIsUint8(X) || mxIsDouble(X))) {
    mexErrMsgIdAndTxt(ERR_ID   "BadClass",
                      ERR_HEAD "1st input must be UINT8 or DOUBLE array.");
  }
  
  // Reply empty matrix for any empty input:
  if (mxIsEmpty(X)) {
    return (1);
  }
  
  // Check if any dim is smaller than window size:
  dims = mxGetDimensions(X);
  if (dims[0] < V || dims[1] < W) {
    return (1);  // Or an error?
  }
  
  return (0);
}

// =============================================================================
mxArray *EmptyAnswer(const mxArray *X)
{
  // Create empty matrix with same type as the input.
  mwSize  dims[2];
  mxArray *R;
  
  if (mxIsDouble(X)) {
    R = mxCreateDoubleMatrix(0, 0, mxREAL);
  } else {
    dims[0] = 0;
    dims[1] = 0;
    R       = mxCreateNumericArray(2, dims, mxGetClassID(X), mxREAL);
  }

  return (R);
}
    
// =============================================================================
#if defined(__WINDOWS__) || defined(WIN32) || defined(_WIN32) || defined(_WIN64)
void C87to64bit(int Cmd)
{
  // Store original _control87 value and set to precision to 64 bits.
  // It seems as Matlab restores the floating point processor flags at the
  // standard exit of a Mex function. But at least Matlab 5.3 and 6.5 leave the
  // flags untouched on exits through mexErrMsgTxt.
  // INPUT: Cmd, flag, 1 means enabling, any thing else restores the former
  //        status
  
#if defined(__LCC__)           // Needs leading underscore
#define MCW_PC _MCW_PC
#define PC_64  _PC_64
#endif

  static unsigned int PC_Orig = 0;

  if (Cmd == 1) {
     PC_Orig = _control87(MCW_PC, 0);
     _control87(PC_64,   MCW_PC);
  } else {
     _control87(PC_Orig, MCW_PC);
  }
  
  return;
}
#endif
