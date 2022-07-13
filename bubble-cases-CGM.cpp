#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath> // for fabs()
int const static n=9,mx=100,my=100; //number of latttice nodes
int c=2;//different cases
int freq=100;//output frequency
double f_r[n][mx][my],f_b[n][mx][my],f[n][mx][my],feq_r[n][mx][my],feq_b[n][mx][my],feq[n][mx][my],rho_r[mx][my],rho_b[mx][my],rho[mx][my];
double cx[n],cy[n],w[n],u[mx][my],v[mx][my],x[mx],y[my];
double CGx[mx][my],CGy[mx][my];//color-gradient
//double cosin;
double collision_2_r[n][mx][my],collision_2_b[n][mx][my];
double cosin[n];
double prodc_c_g[n][mx][my];
double c_norm[n];
double A_r=0.1,A_b=0.1;//parameter that affect the interfacial intension in second term of collision step
double beta=0.99; //parameter that affect the interface thickness in Recoloring step
double C_eq_r[n],C_eq_b[n];//aditional term in equilibrium distribution function
double a_k_r,a_k_b;//parameter to adjust density ratio in equilibrium distribution function
double rho_r0=0.5,rho_b0=1;
double radius=mx/10.0;
double a=mx/10,b=my/5;
int dx=1,dy=1; //space and time step
double const alpha=0.17;
double omega=1.0/(3.*alpha+0.5);
int mstep=5000; // The total number of time steps
void result(std::string filename,int time)
{

std::ofstream out(filename);
  out << "TITLE=\"LBM output\"" << std::endl;
  out << "VARIABLES = \"X\", \"Y\", \"U\",\"V\",\"density\",\"f0\",\"f1\",\"f2\",\"f3\",\"f4\",\"f5\",\"f6\",\"f7\",\"f8\",\"feq0\",\"feq1\",\"feq2\",\"feq3\",\"feq4\",\"feq5\",\"feq6\",\"feq7\",\"feq8\",\"CGx\",\"CGy\",\"rho_r\",\"rho_b\"" << std::endl;
  out << "ZONE T = \"fluid\", I=" << mx << ", J=" << my << ", F=POINT" << std::endl;
  //out << "SOLUTIONTIME="<< time << std::endl;
  for (unsigned j = 0; j < my; ++j)
    for (unsigned i = 0; i < mx; ++i)

    {
        out << std::scientific << std::setprecision(5) << std::setw(15) << x[i];
        out << std::scientific << std::setprecision(5) << std::setw(15) << y[j];
        out << std::scientific << std::setprecision(5) << std::setw(15) << u[i][j];
        out << std::scientific << std::setprecision(5) << std::setw(15) << v[i][j];
        out << std::scientific << std::setprecision(5) << std::setw(15) << rho[i][j];
      for(int k=0;k<9;k++)
      {
          out << std::scientific << std::setprecision(5) << std::setw(15) << f[k][i][j];
      }
      for(int k=0;k<9;k++)
      {
          out << std::scientific << std::setprecision(5) << std::setw(15) << feq[k][i][j];
      }
        out << std::scientific << std::setprecision(5) << std::setw(15) << CGx[i][j];
        out << std::scientific << std::setprecision(5) << std::setw(15) << CGy[i][j];
        out << std::scientific << std::setprecision(5) << std::setw(15) << rho_r[i][j];
        out << std::scientific << std::setprecision(5) << std::setw(15) << rho_b[i][j];
      out << std::endl;
    }
  out.close();
}

void ComputeEquilibrium()
{
    double t1,t2;
    a_k_r=0.4;
    a_k_b=(rho_b0-rho_r0+rho_r0*a_k_r)/rho_b0;
    for (int i=0;i<mx;i++)
    {
        for(int j=0;j<my;j++)
        {
                double rho_temp_r = 0.0;
                double rho_temp_b = 0.0;
                double phi=0.0;
                t1=u[i][j]*u[i][j]+v[i][j]*v[i][j];
                for(int k=0;k<9;k++)
                {
                    if(k==0)
                    {
                        C_eq_b[k]=a_k_b;
                        C_eq_r[k]=a_k_r;
                    }
                    else if(k==1||k==2||k==3||k==4)
                    {
                        C_eq_b[k]=(1-a_k_b)/5.0;
                        C_eq_r[k]=(1-a_k_r)/5.0;
                    }
                    else if(k==5||k==6||k==7||k==8)
                    {
                        C_eq_b[k]=(1-a_k_b)/20.0;
                        C_eq_r[k]=(1-a_k_r)/20.0;
                    }
                    t2=u[i][j]*cx[k]+v[i][j]*cy[k];
                    // if(k!=0)
                    // {
                    feq_r[k][i][j]=rho_r[i][j]*(C_eq_r[k]+w[k]*(3.0*t2+4.50*t2*t2-1.50*t1));
                    feq_b[k][i][j]=rho_b[i][j]*(C_eq_b[k]+w[k]*(3.0*t2+4.50*t2*t2-1.50*t1));
                    // }
                    // else
                    // {
                    //     feq_r[k][i][j]=rho_r[i][j]*(C_eq_r[k]+w[k]*(1.50*t1));
                    //     feq_b[k][i][j]=rho_b[i][j]*(C_eq_b[k]+w[k]*(1.50*t1));
                    // }
                    feq[k][i][j]=feq_b[k][i][j]+feq_r[k][i][j];
                    rho_temp_b += f_b[k][i][j];
                    rho_temp_r += f_r[k][i][j];
                }
                phi=(rho_temp_r-rho_temp_b)/(rho_temp_r+rho_temp_b);
                if (rho_temp_r < 0) std::cout << "negative density(r) at " << i << ", " << j << ", " << std::endl;
                if (rho_temp_b < 0) std::cout << "negative density(b) at " << i << ", " << j << ", " << std::endl;
                if (phi < -1 | phi>1) std::cout << "wrong phi at " << i << ", " << j << ", " << std::endl;
        }
    }
}
void initialize()
{
x[0] =0.0;
y[0] =0.0;
//coordinate
for (int i=1;i<mx;i++)
{
    x[i]=x[i-1]+dx;
}
for (int j=1;j<my;j++)
{
    y[j]=y[j-1]+dy;
}

std::cout<<"omega= "<<omega<<std::endl;
/*-------weight factor ------*/
w[0]=4./9;
    for(int i=1;i<5;i++)
    {
        w[i]=1./9;
    }
    for(int i=5;i<9;i++)
    {
        w[i]=1./36;
    }
/*---------weight factor ----------*/

/*---------streaming vector--------*/
cx[0]=0.0,cx[1]=1.0,cx[2]=0.0,cx[3]=-1.0,cx[4]=0.0,cx[5]=1.0,cx[6]=-1.0,cx[7]=-1.0,cx[8]=1.0;
cy[0]=0.0,cy[1]=0.0,cy[2]=1.0,cy[3]=0.0,cy[4]=-1.0,cy[5]=1.0,cy[6]=1.0,cy[7]=-1.0,cy[8]=-1.0;
/*---------streaming vector--------*/

/*initial condition--------------*/
    for(int i=0;i<mx;i++)
    {
        for (int j=0;j<my;j++)
        {
            switch (c)
            {
            case 1:
            {
                if((x[i]-mx/2)*(x[i]-my/2)+(y[j]-mx/2)*(y[j]-my/2)<radius*radius)
                {
                    rho_r[i][j]=rho_r0;
                    rho_b[i][j]=0.0;
                }
                else
                {
                    rho_r[i][j]=0.0;
                    rho_b[i][j]=rho_b0;
                }
                break;
            }
            case 2:
            {
                if((x[i]-mx/2+radius-1)*(x[i]-mx/2+radius-1)+(y[j]-my/2)*(y[j]-my/2)<=radius*radius
            ||(x[i]-mx/2-radius+1)*(x[i]-mx/2-radius+1)+(y[j]-my/2)*(y[j]-my/2)<radius*radius)
                {
                    rho_r[i][j]=rho_r0;
                    rho_b[i][j]=0.0;
                }
                else
                {
                    rho_r[i][j]=0.0;
                    rho_b[i][j]=rho_b0;
                }
                break;
            }
            case 3:
            {
                if((x[i]-mx/2)*(x[i]-my/2)/(a*a)+(y[j]-mx/2)*(y[j]-my/2)/(b*b)<=1)
                {
                    rho_r[i][j]=rho_r0;
                    rho_b[i][j]=0.0;
                }
                else
                {
                    rho_r[i][j]=0.0;
                    rho_b[i][j]=rho_b0;
                }
                break;
            }
            case 4:
            {
                if((x[i]-mx/2)*(x[i]-mx/2)+(y[j]-my/2+radius-1)*(y[j]-my/2+radius-1)<=radius*radius
            ||(x[i]-mx/2)*(x[i]-mx/2)+(y[j]-my/2-radius+1)*(y[j]-my/2-radius+1)<radius*radius)
                {
                    rho_r[i][j]=rho_r0;
                    rho_b[i][j]=0.0;
                }
                else
                {
                    rho_r[i][j]=0.0;
                    rho_b[i][j]=rho_b0;
                }
                break;
            }
            default:
            break;
            }
            rho[i][j]=rho_b[i][j]+rho_r[i][j];
        }
    }
//initial condition of distribution function
    for(int j=0;j<my;j++)
    {
       for(int i=0;i<mx;i++)
       {
           for(int k=0;k<9;k++)
           {
                f_r[k][i][j]=w[k]*rho_r[i][j];
                f_b[k][i][j]=w[k]*rho_b[i][j];
                f[k][i][j]=f_b[k][i][j]+f_r[k][i][j];
                u[i][j]=0.0;
                v[i][j]=0.0;
           }
       }
    }
/*initial condition--------------*/
}

void ComputeColorGradient()
{
    for (int i=0;i<mx;i++)
    {
        for(int j=0;j<my;j++)
        {

            CGx[i][j]=0.0;
            CGy[i][j]=0.0;

                int y_n = (1+j)%my;
                int x_e = (1+i)%mx;
                int y_s = my - 1 - (my- j)%my;
                int x_w = mx - 1 - (mx- i)%mx;

                CGx[i][j] = cx[1] *(rho_r[x_e][j] - rho_b[x_e][j])
                 + cx[5] *(rho_r[x_e][y_n] - rho_b[x_e][y_n])
                 + cx[8] *(rho_r[x_e][y_s] - rho_b[x_e][y_s])
                 + cx[3] *(rho_r[x_w][j] - rho_b[x_w][j])
                 + cx[6] *(rho_r[x_w][y_n] - rho_b[x_w][y_n])
                 + cx[7] *(rho_r[x_w][y_s] - rho_b[x_w][y_s]);
                CGy[i][j] = cy[2] *(rho_r[i][y_n] - rho_b[i][y_n])
                 + cy[5] *(rho_r[x_e][y_n] - rho_b[x_e][y_n])
                 + cy[6] *(rho_r[x_w][y_n] - rho_b[x_w][y_n])
                 + cy[4] *(rho_r[i][y_s] - rho_b[i][y_s])
                 + cy[7] *(rho_r[x_w][y_s] - rho_b[x_w][y_s])
                 + cy[8] *(rho_r[x_e][y_s] - rho_b[x_e][y_s]);
                // CGx[i][j]=1/6*CGx[i][j];
                // CGy[i][j]=1/6*CGy[i][j];
        //           CGx[i][j] = cx[1] *(rho_r[x_e][j] - rho_b[x_e][j])
        //           +cx[2] *(rho_r[i][y_n] - rho_b[i][y_n])
        //          + cx[5] *(rho_r[x_e][y_n] - rho_b[x_e][y_n])
        //          + cx[6] *(rho_r[x_w][y_n] - rho_b[x_w][y_n])
        //          + cx[7] *(rho_r[x_w][y_s] - rho_b[x_w][y_s])
        //          + cx[8] *(rho_r[x_e][y_s] - rho_b[x_e][y_s]);
        //         CGy[i][j] = cy[2] *(rho_r[i][y_n] - rho_b[i][y_n])
        //          + cy[5] *(rho_r[x_e][y_n] - rho_b[x_e][y_n])
        //          + cy[6] *(rho_r[x_w][y_n] - rho_b[x_w][y_n])
        //          + cy[7] *(rho_r[x_w][y_s] - rho_b[x_w][y_s])
        //          + cy[8] *(rho_r[x_e][y_s] - rho_b[x_e][y_s]);
        }
    }
}
void collision()
{
    // double B[9];
	// B[0] = -4.0/ 27.0;
	// for(int i = 1; i < 5; i++) B[i] = 2.0 / 27.0;
    // for(int i = 5; i < 9; i++) B[i] = 5.0 / 108.0;
    // double CG_norm;
    // double c_norm=0;
    // double cosin=0;
    for(int i=0;i<mx;i++)
    {
        for(int j=0;j<my;j++)
        {
            for (int k=0;k<9;k++)
            {
                // CG_norm=sqrt(CGx[i][j]*CGx[i][j]+CGy[i][j]*CGy[i][j]);
                // if(CG_norm>1e-10)
                // {
                //     double prodc_c_g=cx[k]*CGx[i][j]+cy[k]*CGy[i][j];
                //     if(k!=0)
                //     {
                //         c_norm=sqrt(cx[k]*cx[k]+cy[k]*cy[k]);
                //         cosin=prodc_c_g/(CG_norm*c_norm);
                //         collision_2_r[k][i][j]=A_r/2.0*CG_norm*(w[k]*pow((prodc_c_g/CG_norm),2)-B[k]);
                //         collision_2_b[k][i][j]=A_b/2.0*CG_norm*(w[k]*pow((prodc_c_g/CG_norm),2)-B[k]);
                //     }

                //     else
                //     {
                //         cosin=0;
                //         collision_2_r[k][i][j]=A_r/2.0*CG_norm*(-B[k]);
                //         collision_2_b[k][i][j]=A_b/2.0*CG_norm*(-B[k]);
                //     }
                // }
                // else
                // {
                //     collision_2_r[k][i][j]=0.0;
                //     collision_2_b[k][i][j]=0.0;
                //     cosin=0;
                // }
                //collision
                f_b[k][i][j] = f_b[k][i][j] * (1. - omega) + feq_b[k][i][j] * omega/*+collision_2_b[k][i][j]*/;
                f_r[k][i][j] = f_r[k][i][j] * (1. - omega) + feq_r[k][i][j] * omega/*+collision_2_r[k][i][j]*/;

                f[k][i][j]=f_b[k][i][j]+f_r[k][i][j];

                // f_b[k][i][j] = rho_b[i][j]/rho[i][j]*f[k][i][j]
                // +beta*rho_b[i][j]*rho_r[i][j]/(rho[i][j]*rho[i][j])*(C_eq_b[k]*rho_b[i][j]+C_eq_r[k]*rho_r[i][j])*cosin;

                // f_r[k][i][j] = rho_r[i][j]/rho[i][j]*f[k][i][j]
                // -beta*rho_b[i][j]*rho_r[i][j]/(rho[i][j]*rho[i][j])*(C_eq_b[k]*rho_b[i][j]+C_eq_r[k]*rho_r[i][j])*cosin;
            }
        }
    }
}

void Streaming()
{
    for (int j=0;j<my;j++)
    {
        for(int i=mx-1;i>0;i--)
        {
            f_r[1][i][j]=f_r[1][i-1][j];
            f_b[1][i][j]=f_b[1][i-1][j];
        }
        for(int i=0;i<mx-1;i++)
        {
            f_r[3][i][j]=f_r[3][i+1][j];
            f_b[3][i][j]=f_b[3][i+1][j];
        }
    }

    for(int j=my-1;j>0;j--)
    {
        for(int i=0;i<mx;i++)
        {
            f_r[2][i][j]=f_r[2][i][j-1];
            f_b[2][i][j]=f_b[2][i][j-1];
        }
        for(int i=mx-1;i>0;i--)
        {
            f_r[5][i][j]=f_r[5][i-1][j-1];
            f_b[5][i][j]=f_b[5][i-1][j-1];
        }
        for(int i=0;i<mx-1;i++)
        {
            f_r[6][i][j]=f_r[6][i+1][j-1];
            f_b[6][i][j]=f_b[6][i+1][j-1];
        }
    }

    for(int j=0;j<my-1;j++)
    {
        for(int i=0;i<mx;i++)
        {
            f_r[4][i][j]=f_r[4][i][j+1];
            f_b[4][i][j]=f_b[4][i][j+1];
        }
        for(int i=0;i<mx-1;i++)
        {
            f_r[7][i][j]=f_r[7][i+1][j+1];
            f_b[7][i][j]=f_b[7][i+1][j+1];
        }
        for(int i=mx-1;i>0;i--)
        {
            f_r[8][i][j]=f_r[8][i-1][j+1];
            f_b[8][i][j]=f_b[8][i-1][j+1];
        }
    }

//      //boundary
// 		//north boundary
// 	    int j = my-1;
//         for(int i=0;i<mx;i++)
//         {
//             f_r[2][i][j]=f_r[2][i][j-1];
//             f_b[2][i][j]=f_b[2][i][j-1];
//         }
//         for(int i=mx-1;i>0;i--)
//         {
//             f_r[5][i][j]=f_r[5][i-1][j-1];
//             f_b[5][i][j]=f_b[5][i-1][j-1];
//         }
//         for(int i=0;i<mx;i++)
//         {
//             f_r[6][i][j]=f_r[6][i+1][j-1];
//             f_b[6][i][j]=f_b[6][i+1][j-1];
//         }
//         //South boundary
// 		j = 0;
//         for(int i=0;i<mx;i++)
//         {
//             f_r[4][i][j]=f_r[4][i][j+1];
//             f_b[4][i][j]=f_b[4][i][j+1];
//         }
//         for(int i=0;i<mx;i++)
//         {
//             f_r[7][i][j]=f_r[7][i+1][j+1];
//             f_b[7][i][j]=f_b[7][i+1][j+1];
//         }
//         for(int i=mx-1;i>0;i--)
//         {
//             f_r[8][i][j]=f_r[8][i-1][j+1];
//             f_b[8][i][j]=f_b[8][i-1][j+1];
//         }
//     //east
//     int i = n-1;
//     for (int j=0;j<my;j++)
//     {
//         f_r[1][i][j]=f_r[1][i-1][j];
//         f_b[1][i][j]=f_b[1][i-1][j];
//     }

//    for(int j=my-1;j>0;j--)
//     {
//         f_r[2][i][j]=f_r[2][i][j-1];
//         f_b[2][i][j]=f_b[2][i][j-1];
//         f_r[5][i][j]=f_r[5][i-1][j-1];
//         f_b[5][i][j]=f_b[5][i-1][j-1];
//     }

//    for(int j=0;j<my;j++)
//     {
//         f_r[4][i][j]=f_r[4][i][j+1];
//         f_b[4][i][j]=f_b[4][i][j+1];
//         f_r[8][i][j]=f_r[8][i-1][j+1];
//         f_b[8][i][j]=f_b[8][i-1][j+1];
//     }


//     //west
//     i = 0;
//     for (int j=0;j<my;j++)
//     {
//         f_r[3][i][j]=f_r[3][i+1][j];
//         f_b[3][i][j]=f_b[3][i+1][j];
//     }

//     for(int j=my-1;j>0;j--)
//     {
//         f_r[2][i][j]=f_r[2][i][j-1];
//         f_b[2][i][j]=f_b[2][i][j-1];
//         f_r[6][i][j]=f_r[6][i+1][j-1];
//         f_b[6][i][j]=f_b[6][i+1][j-1];
//     }

//     for(int j=0;j<my;j++)
//     {
//         f_r[4][i][j]=f_r[4][i][j+1];
//         f_b[4][i][j]=f_b[4][i][j+1];
//         f_r[7][i][j]=f_r[7][i+1][j+1];
//         f_b[7][i][j]=f_b[7][i+1][j+1];
//     }

//     //corner
//     // north-east corner
//     i=mx-1; j=my-1;
//     f_r[1][i][j] =f_r[1][i-1][j];
// 	f_r[2][i][j] =f_r[2][i][j-1];
// 	f_r[5][i][j] =f_r[3][i-1][j-1];

//     f_b[1][i][j] =f_b[1][i-1][j];
// 	f_b[2][i][j] =f_b[2][i][j-1];
// 	f_b[5][i][j] =f_b[5][i-1][j-1];

//     //north-west corner
// 	i=0; j=my-1;
//     f_r[2][i][j] =f_r[2][i][j-1];
// 	f_r[3][i][j] =f_r[3][i+1][j];
// 	f_r[6][i][j] =f_r[6][i+1][j-1];

//     f_b[2][i][j] =f_b[2][i][j-1];
// 	f_b[3][i][j] =f_b[3][i+1][j];
// 	f_b[6][i][j] =f_b[6][i+1][j-1];

//     // south-east corner
// 	i=mx-1; j=0;
//     f_r[1][i][j] =f_r[1][i-1][j];
// 	f_r[4][i][j] =f_r[4][i][j+1];
// 	f_r[8][i][j] =f_r[8][i-1][j+1];

//     f_b[1][i][j] =f_b[1][i-1][j];
// 	f_b[4][i][j] =f_b[4][i][j+1];
// 	f_b[8][i][j] =f_b[8][i-1][j+1];

//     // south-west corner
// 	i=0; j=0;
//     f_r[3][i][j] =f_r[3][i+1][j];
// 	f_r[4][i][j] =f_r[4][i][j+1];
// 	f_r[7][i][j] =f_r[7][i+1][j+1];

//     f_b[3][i][j] =f_b[3][i+1][j];
// 	f_b[4][i][j] =f_b[4][i][j+1];
// 	f_b[7][i][j] =f_b[7][i+1][j+1];

    // for(int i=0;i<mx;i++)
    // {
    //     for(int j=0;j<my;j++)
    //     {
    //         f[1][i][j]=f_b[1][i][j]+f_r[1][i][j];
    //         f[2][i][j]=f_b[2][i][j]+f_r[2][i][j];
    //         f[3][i][j]=f_b[3][i][j]+f_r[3][i][j];
    //         f[4][i][j]=f_b[4][i][j]+f_r[4][i][j];
    //         f[5][i][j]=f_b[5][i][j]+f_r[5][i][j];
    //         f[6][i][j]=f_b[6][i][j]+f_r[6][i][j];
    //         f[7][i][j]=f_b[7][i][j]+f_r[7][i][j];
    //         f[8][i][j]=f_b[8][i][j]+f_r[8][i][j];

    //         f[1][i][j]=f_b[1][i][j]+f_r[1][i][j];
    //         f[2][i][j]=f_b[2][i][j]+f_r[2][i][j];
    //         f[3][i][j]=f_b[3][i][j]+f_r[3][i][j];
    //         f[4][i][j]=f_b[4][i][j]+f_r[4][i][j];
    //         f[5][i][j]=f_b[5][i][j]+f_r[5][i][j];
    //         f[6][i][j]=f_b[6][i][j]+f_r[6][i][j];
    //         f[7][i][j]=f_b[7][i][j]+f_r[7][i][j];
    //         f[8][i][j]=f_b[8][i][j]+f_r[8][i][j];
    //     }
    // }
}
// void Streaming()
// {
//     double f_hlp_r[9][mx][my];
//     double f_hlp_b[9][mx][my];
//     for (int i=0;i<mx;i++)
//     {
//         for(int j=0;j<my;j++)
//         {
//             int y_n = j%(my-1) + 1;
//             int x_e = i%(mx-1) + 1;
//             int y_s = my-1 - (my- j)%my;
//             int x_w = mx-1 - (mx- i)%mx;
//             f_hlp_r[1][x_e][j] = f_r[1][i][j];
//             f_hlp_r[2][i][y_n] = f_r[2][i][j];
//             f_hlp_r[3][x_w][j] = f_r[3][i][j];
//             f_hlp_r[4][i][y_s] = f_r[4][i][j];
//             f_hlp_r[5][x_e][y_n] = f_r[5][i][j];
//             f_hlp_r[6][x_w][y_n] = f_r[6][i][j];
//             f_hlp_r[7][x_w][y_s] = f_r[7][i][j];
//             f_hlp_r[8][x_e][y_s] = f_r[8][i][j];

//             f_hlp_b[1][x_e][j] = f_b[1][i][j];
//             f_hlp_b[2][i][y_n] = f_b[2][i][j];
//             f_hlp_b[3][x_w][j] = f_b[3][i][j];
//             f_hlp_b[4][i][y_s] = f_b[4][i][j];
//             f_hlp_b[5][x_e][y_n] = f_b[5][i][j];
//             f_hlp_b[6][x_w][y_n] = f_b[6][i][j];
//             f_hlp_b[7][x_w][y_s] = f_b[7][i][j];
//             f_hlp_b[8][x_e][y_s] = f_b[8][i][j];

//             for(int k=1;k<9;k++)
//             {
//                 f_r[k][i][j]=f_hlp_r[k][i][j];
//                 f_b[k][i][j]=f_hlp_b[k][i][j];
//             }
//         }
//     }
// }
void recoloring()
{
    double B[9];
	B[0] = -4.0/ 27.0;
	for(int i = 1; i < 5; i++) B[i] = 2.0 / 27.0;
    for(int i = 5; i < 9; i++) B[i] = 5.0 / 108.0;
    for(int i=0;i<mx;i++)
    {
        for(int j=0;j<my;j++)
        {
                double CG_norm=sqrt(CGx[i][j]*CGx[i][j]+CGy[i][j]*CGy[i][j]);
                //std::cout<<"CG_norm="<<CG_norm<<std::endl;
                if (CG_norm<1e-30)
                {
                    for (int k=0;k<9;k++)
                    {
                        f_b[k][i][j] = rho_b[i][j]/rho[i][j]*f[k][i][j];
                        f_r[k][i][j] = rho_r[i][j]/rho[i][j]*f[k][i][j];
                    }
                }
                else
                {
                    //second collision term
                    int k=0;
                    f_r[k][i][j] = f_r[k][i][j] + A_r/2* CG_norm *(-B[k]);
                    f_b[k][i][j] = f_b[k][i][j] + A_b/2* CG_norm *(-B[k]);
                    for(int k=1;k<9;k++)
                    {
                        prodc_c_g[k][i][j]=cx[k]*CGx[i][j]+cy[k]*CGy[i][j];
                        c_norm[k]=sqrt(cx[k]*cx[k]+cy[k]*cy[k]);
                        cosin[k]=prodc_c_g[k][i][j]/(CG_norm*c_norm[k]);
                        //std::cout<<"prodc_c_g/CG_norm"<<" "<<time<<" "<<k<<" "<<i<<" "<<j<<" "<<prodc_c_g[k][i][j]/sqrt(CGx[i][j]*CGx[i][j]+CGy[i][j]*CGy[i][j])<<std::endl;
                        f_r[k][i][j] =f_r[k][i][j]+ A_r/2*CG_norm*(w[k]*pow((prodc_c_g[k][i][j]/CG_norm),2)-B[k]);
                        f_b[k][i][j] =f_b[k][i][j]+ A_b/2*CG_norm*(w[k]*pow((prodc_c_g[k][i][j]/CG_norm),2)-B[k]);
                    }
                    //recoloring
                    int k1=0;
                    f_b[k1][i][j] = rho_b[i][j]/rho[i][j]*f[k1][i][j];
                    f_r[k1][i][j] = rho_r[i][j]/rho[i][j]*f[k1][i][j];
                    for(int k1=1;k1<9;k1++)
                    {
                        // double feq= w[k1]*rho[i][j];
                        double temp=rho_b[i][j]*rho_r[i][j]/(rho[i][j]*rho[i][j]);
                        f_r[k1][i][j] = rho_r[i][j]/rho[i][j]*f[k1][i][j]+ beta* temp * (C_eq_b[k1]*rho_b[i][j]+C_eq_r[k1]*rho_r[i][j]) *cosin[k1];
                        f_b[k1][i][j] = rho_b[i][j]/rho[i][j]*f[k1][i][j]- beta* temp * (C_eq_b[k1]*rho_b[i][j]+C_eq_r[k1]*rho_r[i][j]) *cosin[k1];
                    }
                }
        }
    }
}
void Boundarcyondition()
{
    //periodic
    // for(int j=0;j<my;j++)
    // {
    //     f_b[1][0][j]=f_b[1][mx-1][j];//left
    //     f_b[5][0][j]=f_b[5][mx-1][j];
    //     f_b[8][0][j]=f_b[8][mx-1][j];

    //     f_r[1][0][j]=f_r[1][mx-1][j];//left
    //     f_r[5][0][j]=f_r[5][mx-1][j];
    //     f_r[8][0][j]=f_r[8][mx-1][j];

    //     f_b[3][mx-1][j]=f_b[3][0][j];//right
    //     f_b[7][mx-1][j]=f_b[7][0][j];
    //     f_b[6][mx-1][j]=f_b[6][0][j];

    //     f_r[3][mx-1][j]=f_r[3][0][j];//right
    //     f_r[7][mx-1][j]=f_r[7][0][j];
    //     f_r[6][mx-1][j]=f_r[6][0][j];
    // }
    //bounce back
    for(int j=0;j<my;j++)
    {
        f_b[1][0][j]=f_b[3][mx-1][j];//left
        f_b[5][0][j]=f_b[7][mx-1][j];
        f_b[8][0][j]=f_b[6][mx-1][j];

        f_r[1][0][j]=f_r[3][mx-1][j];//left
        f_r[5][0][j]=f_r[7][mx-1][j];
        f_r[8][0][j]=f_r[6][mx-1][j];

        f_b[3][mx-1][j]=f_b[1][0][j];//right
        f_b[7][mx-1][j]=f_b[5][0][j];
        f_b[6][mx-1][j]=f_b[8][0][j];

        f_r[3][mx-1][j]=f_r[1][0][j];//right
        f_r[7][mx-1][j]=f_r[5][0][j];
        f_r[6][mx-1][j]=f_r[8][0][j];
    }
    //bounce back
    for (int i=0;i<mx;i++)
    {
        f_b[2][i][0]=f_b[4][i][0];//bottom
        f_b[5][i][0]=f_b[7][i][0];
        f_b[6][i][0]=f_b[8][i][0];

        f_r[2][i][0]=f_r[4][i][0];//bottom
        f_r[5][i][0]=f_r[7][i][0];
        f_r[6][i][0]=f_r[8][i][0];


        f_b[4][i][my-1]=f_b[2][i][my-1];//top
        f_b[7][i][my-1]=f_b[5][i][my-1];
        f_b[8][i][my-1]=f_b[6][i][my-1];

        f_r[4][i][my-1]=f_r[2][i][my-1];//top
        f_r[7][i][my-1]=f_r[5][i][my-1];
        f_r[8][i][my-1]=f_r[6][i][my-1];

    }
}

void rhouv()
{
    double ssum_r,ssum_b,usum,vsum,usum_b,usum_r,vsum_b,vsum_r;
    for(int j=0;j<my;j++)
    {
       for(int i=0;i<mx;i++)
       {
           ssum_r=0.0;
           ssum_b=0.0;
           for(int k=0;k<9;k++)
           {
                ssum_r=ssum_r+f_r[k][i][j];
                ssum_b=ssum_b+f_b[k][i][j];
           }
           rho_r[i][j]=ssum_r;
           rho_b[i][j]=ssum_b;
           rho[i][j]=rho_r[i][j]+rho_b[i][j];
       }
    }

    for(int i=0;i<mx;i++)
    {
        for(int j=0;j<my;j++)
        {
            usum=0.0;
            vsum=0.0;
            usum_b=0.0;
            usum_r=0.0;
            vsum_b=0.0;
            vsum_r=0.0;
            // TODO: check k init value
            for(int k=0;k<9;k++)
            {
                usum_b=usum_b+f_b[k][i][j]*cx[k];
                usum_r=usum_r+f_r[k][i][j]*cx[k];
                usum=usum_b+usum_r;
                vsum_b=vsum_b+f_b[k][i][j]*cy[k];
                vsum_r=vsum_r+f_r[k][i][j]*cy[k];
                vsum=vsum_b+vsum_r;
            }

            u[i][j]=usum/rho[i][j];
            v[i][j]=vsum/rho[i][j];

        }
    }
}
int main()
{
int time;
time=0;
initialize();
std::string filename = std::string("animation") +std::to_string(time)+std::string(".dat");
result(filename,time);
//main loop
for (time=1;time<mstep+1;++time)
{
// std::cout<<"a_k_b"<<a_k_b<<std::endl;
// std::cout<<"a_k_r"<<a_k_r<<std::endl;
// for(int k=0;k<9;k++)
// {
//     //std::cout<<"cosin"<<" "<<time<<" "<<k<<" "<<cosin[k]<<std::endl;
//     //std::cout<<"c_norm"<<" "<<time<<" "<<k<<" "<<c_norm[k]<<std::endl;
// }
// for(int i=0;i<mx;i++)
// {
//     for(int j=0;j<my;j++)
//     {
//         for(int k=0;k<9;k++)
//         {
//             //std::cout<<"prodc_c_g"<<" "<<time<<" "<<k<<" "<<i<<" "<<j<<" "<<prodc_c_g[k][i][j]<<std::endl;
//             //std::cout<<"prodc_c_g/CG_norm"<<" "<<time<<" "<<k<<" "<<i<<" "<<j<<" "<<prodc_c_g[k][i][j]/sqrt(CGx[i][j]*CGx[i][j]+CGy[i][j]*CGy[i][j])<<std::endl;
//         }
//     }
// }
Streaming();
rhouv();
ComputeEquilibrium();//for collision and recoloring
collision();
ComputeColorGradient();
recoloring();
Boundarcyondition();
std::string filename = std::string("animation") +std::to_string(time)+std::string(".dat");
if(time%freq==0)
{
    result(filename,time);
}

}//main loop end
return 0;
}




