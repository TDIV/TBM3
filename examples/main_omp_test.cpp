#include <omp.h>
#include <iostream>
#include "matrix.h"
#include <unistd.h> 
#include <ctime>
#include <stdio.h>
#include <vector>
using namespace std;

// SourceTree Test
class tester{
	cmat A;
	vector<int> B;
	public:
	tester(){
		B.resize(10);
		A = cmat(2,2);

		time_t rawtime;

		time (&rawtime);
		printf ("The current local time is: %s", ctime (&rawtime));
		{
			#pragma omp parallel
			{
				#pragma omp sections
				{
					#pragma omp section
					{
						B[0] = 20;
						sleep(1);
					}
					#pragma omp section
					{
						B[1] = 20;
						sleep(1);
					}
				}
			}
		}
		time (&rawtime);
		printf ("The current local time is: %s", ctime (&rawtime));

		cout<<B[0]<<endl;
		cout<<B[1]<<endl;

	}
	void add(int i, int j, cvar val){
		A(i,j)=val;
	}
	
};

int main(){

	tester t;
	return 0;
}
