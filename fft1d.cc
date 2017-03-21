# include <iostream>
# include <fstream>
# include <sstream>
# include <string>
# include <vector>
# include <signal.h>
# include <math.h>
# include "complex.h"
# include "InputImage.h"
# define pi 3.18
using namespace std;

void Transform2d(const char *inputFN)
{
InputImage img ("inputFN");
InputImage img_op();
int c = img.GetWidth();
int r = img.GetHeight();
Complex buf1[c];
Complex buf2[c];
Complex op_img[r][c];

for(int i =0;i<r;i++)
{
buf1[].mag()= img[i][];
Transform1d(&buf1,c,&buf2);
op_img[i][].mag()= buf2[].mag();
}
img_op.SaveImageData(amit.txt,&op_img,c,r);
}

void Transform1d(Complex *h, int w, Complex *H)
{
int i,j = 0;
Complex W;
W.r=acos((2*pi)/w);
W.i=-asin((2*pi/w))); 
for(i=0;i<w;i++)
 {
  for(j=0;j<w;j++)
  {
   *H =(*h) + ((*h) * pow(W,i));
   h++;  
  
  }
   H++;
 }
}

int main (int argc, char **argv)
{
string fn("Tower.txt");
Transform2d(fn.c_str());
return 0;
}
