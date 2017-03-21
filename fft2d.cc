// Distributed two-dimensional Discrete FFT transform
// Amit Sunil Dhamne (GTID# 903223468)
// ECE8893 Project 1


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <signal.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include "Complex.h"
#include "InputImage.h"

using namespace std;
void Transform1D(Complex* h, int w, Complex* H);

void Inverse(Complex* E, int n , Complex* e);
void Transform2D(const char* inputFN,int rc)
{ // Do the 2D transform here.
  // 1) Use the InputImage object to read in the Tower.txt file and
  //    find the width/height of the input image.--- DONE---
  // 2) Use MPI to find how many CPUs in total, and which one
  //    this process is
  // 3) Allocate an array of Complex object of sufficient size to
  //    hold the 2d DFT results (size is width * height) ----DONE---
  // 4) Obtain a pointer to the Complex 1d array of input data ---DONE----
  // 5) Do the individual 1D transforms on the rows assigned to your CPU ---- DONE---
  // 6) Send the resultant transformed values to the appropriate
  //    other processors for the next phase.
  // 6a) To send and receive columns, you might need a separate
  //     Complex array of the correct size.
  // 7) Receive messages from other processes to collect your columns
  // 8) When all columns received, do the 1D transforms on the columns
  // 9) Send final answers to CPU 0 (unless you are CPU 0)
  //   9a) If you are CPU 0, collect all values from other processors
  //	   and print out with SaveImageData().
  InputImage image(inputFN);  // Create the helper object for reading the image
  // Step (1) in the comments is the line above.
  // Your code here, steps 2-9
 
  int c = image.GetWidth();
  int r = image.GetHeight();
  
  int numtasks,rank;
  int size = r*c;
  Complex* h = image.GetImageData();
  Complex buf[size] ;
  Complex* H = buf;
  int cmplex = sizeof(Complex);
  Complex buf2[size];
  Complex* v = buf2;
  Complex buf3[size];
  Complex* v1 = buf3;
  Complex buf4[size];
  Complex* v2 = buf4;
  
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int thick = c*(r/numtasks);
   printf("Number of tasks = %d My rank = %d\n",numtasks ,rank);
  if (rank>0)
    {
      for(int m =0; m<rank;m++)
         h = h +thick;
     }

  for(int sn=0;sn<16;sn++)
    {
      Transform1D(h,c,H);
      h=h+c;
      H=H+c;
     }


  if(rank>0)
  {

      H = &buf[0];
      rc = MPI_Send(H,256*16,MPI_C_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD);

      // if (rc == MPI_SUCCESS) cout<<"send success"<<endl;

      if (rc != MPI_SUCCESS)
            { cout << "Rank " << rank<< " send failed, rc " << rc << endl;
              MPI_Finalize();
              exit(1);
            }
   }

  if(rank ==0)
    {

     MPI_Status status;
     Complex buf1[r][c];
     for (int i =1;i<16;i++)
        {
         rc = MPI_Recv(H, 256*16,MPI_C_DOUBLE_COMPLEX, i,0,MPI_COMM_WORLD,&status);
         if (rc != MPI_SUCCESS)
          {cout << "Rank " << rank<< " send failed, rc " << rc << endl;}
          H = H+(256*16);
        }
       for(int i1=0; i1<r;i1++)
          {
           for(int i2 =0; i2<c; i2++)
              {
               	buf1[i1][i2] = buf[i1*256+i2];
              }
           }

	for (int i3=0; i3<c;i3++)
           {
            for(int i4=0;i4<r;i4++)
               {
                buf[i3*256+i4] = buf1[i4][i3];
               }
           }

        H=&buf[0];
         image.SaveImageData("MyAfter1D.txt",H,r,c);
        for (int a1=1;a1<16 ;a1++)
        {
            rc = MPI_Send(H,256*256,MPI_C_DOUBLE_COMPLEX,a1,0,MPI_COMM_WORLD);
            if (rc == MPI_SUCCESS) cout<<"Sent Transpose to "<<rank<<endl;
            if (rc != MPI_SUCCESS)
               {cout << "Rank " << a1<< " send failed, rc " << rc << endl; MPI_Finalize();exit(1);}
  }
	v = &buf2[0];
        for(int ns=0;ns<16;ns++)
            {
             Transform1D(H,c,v);
             v=v+c;
             H=H+c;
             cout<<"Transform sentfor 2D"<<" Rank:"<<rank<<endl;
            }
          
        for (int a3=1;a3<16;a3++)
            {
             rc = MPI_Recv(v,256*16,MPI_C_DOUBLE_COMPLEX,a3,0,MPI_COMM_WORLD,&status);
             v = v+(256*16);
            }
	v=&buf2[0];
         for(int iq1=0; iq1<r;iq1++)
          {
           for(int iq2 =0; iq2<c; iq2++)
              {
               	buf1[iq1][iq2] = buf2[iq1*256+iq2];
              }
           }

	for (int iq3=0; iq3<c;iq3++)
           {
            for(int iq4=0;iq4<r;iq4++)
               {
                buf2[iq3*256+iq4] = buf1[iq4][iq3];
               }
           }
      v =&buf2[0];
       //for(int a4=thick*16;a4<thick*16*2;a4++) buf2[a4].Print();
       image.SaveImageData("MyAfter2D.txt",v,r,c);
      
       for(int u2=1; u2<16;u2++)
         {
           rc = MPI_Send(v,256*256,MPI_C_DOUBLE_COMPLEX,u2,0,MPI_COMM_WORLD);
            if (rc != MPI_SUCCESS)
               {cout << "Rank " << u2<< " send failed, rc " << rc << endl; MPI_Finalize();exit(1);}
            if (rc==MPI_SUCCESS) {cout<<"Succesfully sent 2D for 1D to rank: "<<u2<<endl;}
         }
       v=&buf2[0];
       v1 = &buf3[0]; 
       for(int o=0;o<16;o++)
          {
           Inverse(v,c,v1);
           v=v+c;
           v1 = v1+c;
          }
       Complex  h7;
       h7.real,h7.imag = 1.0/c;
     
       for(int m1=1;m1<16;m1++)
          {
           rc = MPI_Recv(v1,256*16,MPI_C_DOUBLE_COMPLEX,m1,0,MPI_COMM_WORLD,&status);
          v1= v1+(256*16); 
          }
        for(int b1=0;b1<256*16*16;b1++) buf3[b1] =buf3[b1]*h7; 
        //for(int m9 = 0; m9<256*256;m9++) buf3[m9].Print();
        for(int t1=0; t1<r;t1++)
          {
           for(int t2 =0; t2<c; t2++)
              {
               	buf1[t1][t2] = buf3[t1*256+t2];
              }
           }

	for (int t3=0; t3<c;t3++)
           {
            for(int t4=0;t4<r;t4++)
               {
                buf3[t3*256+t4] = buf1[t4][t3];
               }
           }
       v1 = &buf3[0];
       for(int g1=1;g1<16;g1++)
          {
           rc  = MPI_Send(v1,256*256,MPI_C_DOUBLE_COMPLEX,g1,0,MPI_COMM_WORLD);
          
            if (rc != MPI_SUCCESS)
               {cout << "Rank " << rank<< " send failed, rc " << g1 << endl; MPI_Finalize();exit(1);}
            if (rc==MPI_SUCCESS) {cout<<"Succesfully sent 2D for 1D to rank: "<<g1<<endl;} 
          }
        v1 = &buf3[0];
        v2 = &buf4[0];
        for (int op1 =0; op1<16;op1++)
            {
              Inverse(v1,c,v2);
              v1 = v1+c; v2 = v2+c;
             }
         //for(int ui=0;ui<256*16;ui++) buf4[ui].Print();
        for(int m2=1;m2<16;m2++)
           {
             rc = MPI_Recv(v2, 256*16,MPI_C_DOUBLE_COMPLEX,m2,0,MPI_COMM_WORLD,&status);
             if (rc != MPI_SUCCESS)
               {cout << "Rank " << rank<< " send failed, rc " << rc << endl; MPI_Finalize();exit(1);} 
             if (rc == MPI_SUCCESS){cout<<"Recvd final Inverse from"<<m2<<endl;}
             v2 = v2+256*16;
            }
        for(int b1=0;b1<256*16*16;b1++) buf4[b1] =buf4[b1]*h7; 
      //  for (int yu = 0;yu<256*256;yu++) buf4[yu].Print();
        for(int s1=0; s1<r;s1++)
          {
           for(int s2 =0; s2<c; s2++)
              {
               	buf1[s1][s2] = buf4[s1*256+s2];
              }
           }

	for (int s3=0; s3<c;s3++)
           {
            for(int s4=0;s4<r;s4++)
               {
                buf4[s3*256+s4] = buf1[s4][s3];
               }
           }
        for (int yu = 0;yu<256*256;yu++) buf4[yu]=buf4[yu]*Complex(-1,-1);
        v2 = &buf4[0];
        image.SaveImageDataReal("MyAfterInverse.txt",&buf4[0],r,c);

     
    }//final bracket for if(rank==0)

  if (rank>0)
     {
       MPI_Status status1;
       H=&buf[0];
       rc = MPI_Recv(H, 256*256,MPI_C_DOUBLE_COMPLEX, 0,0,MPI_COMM_WORLD,&status1);
       if (rc != MPI_SUCCESS)
          {cout << "Rank " << rank<< " send failed, rc " << rc << endl; MPI_Finalize();exit(1); }
       H = H + 256*16*rank;

       for(int ns=0;ns<16;ns++)
          {
           Transform1D(H,c,v);
            v=v+c;
            H=H+c;
           cout<<"Transform senfor 2Dt"<<" Rank:"<<rank<<endl;
                // if (rank==15) c4++;
          }
           	// if (rank ==15)  cout<<"count of 5 "<<c4<<endl;//TEST CONDITION
       /*if (rank ==15)
          {
           for(int f1=0;f1<255*16;f1+=256)
              {

               for(int a5=0;a5<256;a5++)
               buf[f1+a5].Print();
               printf("\n");
             }
            }
      Test whether Transform 1D is working
      Test for negation cout<<"for loop to negate done for rank: "<<rank<<endl;*/

      v = &buf2[0];
      rc = MPI_Send(v,256*16,MPI_C_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD);
      if (rc == MPI_SUCCESS) cout<<"send success for 2D"<<endl;
      if (rc != MPI_SUCCESS)
        {
           cout << "Rank " << rank<< " send failed, rc " << rc << endl;
           MPI_Finalize();
           exit(1);
        }

      v=&buf2[0];
      rc = MPI_Recv(v,256*256,MPI_C_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD,&status1);
      if(rc == MPI_SUCCESS){cout<<"Recvd for Inverse Stage 1 Rank: "<<rank<<endl;}
      v = v+(256*16*rank);
      for(int ns1=0;ns1<16;ns1++)
         {
          Inverse(v,c,v1);
          v=v+c;
          v1= v1+c;
         }   
        v1=&buf3[0];
        rc = MPI_Send(v1,256*16,MPI_C_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD);
        if (rc == MPI_SUCCESS){cout<< "Stage 1 complete Rank: "<<rank<<endl;}
        
      v1= &buf3[0];
      rc = MPI_Recv(v1,256*256, MPI_C_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD,&status1);
      if(rc == MPI_SUCCESS){cout<<"Recvd for Inverse Stage 2 Rank: "<<rank<<endl;}
      /*if (rank ==15)
          {
           for(int f1=0;f1<255*256;f1+=256)
              {

               for(int a5=0;a5<256;a5++)
               buf3[f1+a5].Print();
               printf("\n");
             }
            }*/
      v1 = &buf3[0];
      v1 = v1+(256*16*rank);
     for(int ns3 = 0; ns3<16;ns3++)
       {
        Inverse(v1,c,v2);
       v1 = v1+c; v2 = v2=c+v2;
       }
      v2 = &buf4[0]; //if(rank==15){for(int chu=0;chu<256;chu++) buf4[chu].Print() ;}
     rc = MPI_Send(v2,256*16,MPI_C_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD);
     if(rc==MPI_SUCCESS){cout<<"Sent final 1D invrse to "<<rank<<endl;}
 


  } //Final Bracket for rank>0




                        //if (rank==0){
                        //for (int j= 0; j<((thick*16)-1);j++) buf[j].Print();}

  cout<< "Exiting Normally" << "rank"<<rank<<endl;
  MPI_Finalize();

 } //Final Bracket for function

/* Works or without MPI for(int m=0;m<r;m++)
 {

  Transform1D(h,c,H);
  h=h+256;
  H=H+256;
}
  for (int j =255; j<(256*2-1);j++) buf[j].Print(); */






 void Transform1D(Complex* h, int w, Complex* H)
 {
  // Implement a simple 1-d DFT using the double summation equation
  // given in the assignment handout.  h is the time-domain input
  // data, w is the width (N), and H is the output array.
  /* Complex twid
  twid.real = acos((2*pi/x*j*w));
  twid.imag = asin((2*pi/j*x*w)); */
  for (int x=0;x<w;x++)
  {
    H[x]=0;
    int j=0;
  for (;j<w;j++)
  {
    H[x] = H[x] + (h[j]*Complex(cos(2*M_PI*x*j/w),-sin(2*M_PI*x*j/w)));
   }
  }
}

void Inverse (Complex *E, int n, Complex* e)
{
// int m = 1/n;
// Complex  h7;
// h7.real = 1.0/n;
//h7.imag =1.0/n;
 for(int x=0;x<n;x++)
   {
     e[x]=0;
     int j=0;
    for(;j<n;j++)
     {
       e[x]=e[x]+(E[j]*Complex(cos(2*M_PI*x*j/n),sin(2*M_PI*x*j/n)));
     }
     //e[x]= e[x]*h7;
   }
}
int main(int argc, char** argv)
{
  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  // MPI initialization here
  int rc;
  rc = MPI_Init(&argc,&argv);
  if (rc!= MPI_SUCCESS){
  printf("Error starting MPI program. Terminating.\n");
  MPI_Abort(MPI_COMM_WORLD,rc);}
  else cout<<"program initialised"<<endl;
  Transform2D(fn.c_str(),rc); // Perform the transform.
  // Finalize MPI here
  //  MPI_Finalize();
}


 
       
               
            
