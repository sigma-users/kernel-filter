# kernel-filter
## 1. Abstract 

These are the image data and a program that performs kernel filtering used in the following paper.  
Kohno, Y., Seki, T., Findlay, S.D. et al. Real-space visualization of intrinsic magnetic fields of an antiferromagnet. Nature 602, 234–239 (2022). https://doi.org/10.1038/s41586-021-04254-z

## 2. Image data

The image data is saved in binary format as double precision float data. Each pixel data is alined from the upper left to the lower right of the image so that the data in the same row is lined up. The number of pixels in the image is shown as a number at the end of the file name.

The following files are the unprocessed image data used in the above paper. The file with DF in the file name is the data of the ADF image. The file with Dx in the file name is the data of the deflection angle in the X (right) direction of the image. The file with Dy in the file name is the data of the deflection angle in the Y (top) direction of the image. The deflection angles are stored in rads.

### a. The following data are used in Fig. 3.

Main_DF_942_942.bin  
Main_Dx_942_942.bin  
Main_Dy_942_942.bin

### b. The following data are the room temperature data used in Fig. 4.

RT_DF_912_912.bin  
RT_Dx_912_912.bin  
RT_Dy_912_912.bin

### c. The following data are the liquid nitrogen temperature data used in Fig. 4.

N2_DF_860_860.bin  
N2_Dx_860_860.bin  
N2_Dy_860_860.bin

### 3. Processing content

In this program, the ADF image, XY direction deflection map, and B-field-filtered XY direction deflection map are calculated and saved as a bitmap image from the binary files. The binary files must be placed in the same directory as the executable file.

### 4. Prerequisites

The following development environment is required to build the program.

Visual Studio 2019  
Intel oneAPI Base Toolkit  
Intel oneAPI HPC Toolkit
