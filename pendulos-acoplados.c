#include<stdio.h>
#include<math.h>

const double Pi = 3.14159265359;

double theta = 0*Pi/180.0, theta1 = 30.5*Pi/180.0, theta2 = 0, omega = 0, omega1 = 0, omega2 = 0;
double l = 0.5, L = 2, h = 0.00001, t = 0, epsilon = 0.1;
double m1 = 1, m2 = 1, m = 0.1, g = 9.8, a = 1;
double k1[6], k2[6], k3[6], k4[6];

double f(double theta, double theta1, double theta2, double omega, double omega1, double omega2, int j)
{
	if (j==0)	return omega;		//theta
	if (j==1)	return omega1;		//theta1
	if (j==2)	return omega2;		//theta2
	
	if (j==3)	return (cos(theta-theta1)*sin(theta-theta1)*l*m1*omega*omega+cos(theta-theta2)*sin(theta-theta2)*l*m2*omega*omega+m1*L*omega1*omega1*sin(theta-theta1)+m2*L*omega2*omega2*sin(theta-theta2)-cos(theta-theta1)*sin(theta1)*g*m1-cos(theta-theta2)*sin(theta2)*g*m2+g*sin(theta)*m+g*sin(theta)*m1+g*sin(theta)*m2)/(l*(cos(theta-theta1)*cos(theta-theta1)*m1+cos(theta-theta2)*cos(theta-theta2)*m2-m-m1-m2));
	if (j==4)	return -(cos(theta-theta1)*cos(theta-theta2)*sin(theta-theta2)*l*m2*omega*omega-sin(theta-theta1)*cos(theta-theta2)*cos(theta-theta2)*l*m2*omega*omega+L*cos(theta-theta1)*sin(theta-theta1)*m1*omega1*omega1+L*cos(theta-theta1)*sin(theta-theta2)*m2*omega2*omega2-cos(theta-theta1)*cos(theta-theta2)*sin(theta2)*g*m2+sin(theta-theta1)*l*m*omega*omega+sin(theta-theta1)*l*m1*omega*omega+sin(theta-theta1)*l*m2*omega*omega+sin(theta1)*cos(theta-theta2)*cos(theta-theta2)*g*m2+cos(theta-theta1)*sin(theta)*g*m+cos(theta-theta1)*sin(theta)*g*m1+cos(theta-theta1)*sin(theta)*g*m2-sin(theta1)*g*m-sin(theta1)*g*m1-sin(theta1)*g*m2)/(L*(cos(theta-theta1)*cos(theta-theta1)*m1+cos(theta-theta2)*cos(theta-theta2)*m2-m-m1-m2));
	if (j==5)	return -(-cos(theta-theta1)*cos(theta-theta1)*sin(theta-theta2)*l*m1*omega*omega+cos(theta-theta1)*sin(theta-theta1)*cos(theta-theta2)*l*m1*omega*omega+L*sin(theta-theta1)*cos(theta-theta2)*m1*omega1*omega1+L*cos(theta-theta2)*sin(theta-theta2)*m2*omega2*omega2+cos(theta-theta1)*cos(theta-theta1)*sin(theta2)*g*m1-cos(theta-theta1)*sin(theta1)*cos(theta-theta2)*g*m1+sin(theta-theta2)*l*m*omega*omega+sin(theta-theta2)*l*m1*omega*omega+sin(theta-theta2)*l*m2*omega*omega+cos(theta-theta2)*sin(theta)*g*m+cos(theta-theta2)*sin(theta)*g*m1+cos(theta-theta2)*sin(theta)*g*m2-sin(theta2)*g*m-sin(theta2)*g*m1-sin(theta2)*g*m2)/(L*(cos(theta-theta1)*cos(theta-theta1)*m1+cos(theta-theta2)*cos(theta-theta2)*m2-m-m1-m2));
}

int main()
{
	FILE *arq;
	arq = fopen("dados.txt", "w");
	fprintf(arq, "t theta theta1 theta2 omega omega1 omega2\n");
	
	for (int k=0; k<=20000000; k++)
	{
		if (k%1000 == 0)
		{
			printf(      "%f %g %g %g %g %g %g\n", t, theta*180/Pi, theta1*180/Pi, theta2*180/Pi, omega, omega1, omega2);
			fprintf(arq, "%f %g %g %g %g %g %g\n", t, theta*180/Pi, theta1*180/Pi, theta2*180/Pi, omega, omega1, omega2);
		}
		
		for (int j=0; j<=5; j++)
			k1[j] = f(theta, theta1, theta2, omega, omega1, omega2, j);
			
		for (int j=0; j<=5; j++)
			k2[j] = f(theta+0.5*h*k1[0], theta1+0.5*h*k1[1], theta2+0.5*h*k1[2], omega+0.5*h*k1[3], omega1+0.5*h*k1[4], omega2+0.5*h*k1[5], j);
			
		for (int j=0; j<=5; j++)
			k3[j] = f(theta+0.5*h*k2[0], theta1+0.5*h*k2[1], theta2+0.5*h*k2[2], omega+0.5*h*k2[3], omega1+0.5*h*k2[4], omega2+0.5*h*k2[5], j);

		for (int j=0; j<=5; j++)
			k4[j] = f(theta+h*k3[0], theta1+h*k3[1], theta2+h*k3[2], omega+h*k3[3], omega1+h*k3[4], omega2+h*k3[5], j);
			
		t = t + h;
		theta  = theta  + h/6.0 * ( k1[0] + 2*k2[0] + 2*k3[0] + k4[0] );
		theta1 = theta1 + h/6.0 * ( k1[1] + 2*k2[1] + 2*k3[1] + k4[1] );
		theta2 = theta2 + h/6.0 * ( k1[2] + 2*k2[2] + 2*k3[2] + k4[2] );
		omega =  omega  + h/6.0 * ( k1[3] + 2*k2[3] + 2*k3[3] + k4[3] );
		omega1 = omega1 + h/6.0 * ( k1[4] + 2*k2[4] + 2*k3[4] + k4[4] );
		omega2 = omega2 + h/6.0 * ( k1[5] + 2*k2[5] + 2*k3[5] + k4[5] );
	}
}