/*******************************************************
**Isabella Patnode ~ COMP233.A ~ Count Threes Experiment
*******************************************************/

#include<iostream>
#include<iomanip>
#include<time.h>
#include<omp.h>
#include<fstream>
using namespace std;

const int MAX = 1 << 26;	//64 MEG

#define STALL_CNT 8 
#define PAD 8
#define COL_WIDTH 12
#define TEST_COUNT 3
#define MAX_THREADS 5

//Function prototypes
void fillAry(int a[MAX]);
int countThreesSerial(int a[MAX]);
int countThrees1D(int a[MAX]);
int countThreesPad(int a[MAX]);
int countThreesCritical(int a[MAX]);
int countThreesAtomic(int a[MAX]);
int countThreesFor(int a[MAX]);


/********MAIN********/
int main() {
	clock_t startT, stopT; //wallclock timer
	double testTime, fastestTime; //times
	int numThrees, numThreesPar; //counts number of threes
	int* ary = new int[MAX]; //array of random MAX num ints
	int numTimes; //loop variable
	int threadNum;
	ofstream outFile;
	string fileName = "CountThreesExperiment.csv";

	cout << "Isabella Patnode ~ COMP233.A ~ Count Three's "
		<< "Experiment" << endl;
	cout << "This program times counting threes in a randomly generated"
		<< " array using a serial method and various parallel methods\n\n";

	outFile.open(fileName);
			
	/***Constructs heading of table***/

	//prints heading to console
	cout << setw(8) << "#thrds" << setw(COL_WIDTH) << "# threes" <<
		setw(COL_WIDTH) << "Method" << setw(COL_WIDTH) << "Error"
		<< setw(COL_WIDTH) << "Seconds\n";
	cout << setw(8) << "------" << setw(COL_WIDTH) << "--------" <<
		setw(COL_WIDTH) << "------" << setw(COL_WIDTH) << "-----"
		<< setw(COL_WIDTH) << "-------\n";

	//prints heading to file
	outFile << "#thrds" << "," << "# threes" <<
		"," << "Method" << "," << "Error"
		<< "," << "Seconds\n";
	outFile << "------" << "," << "--------" <<
		"," << "------" << "," << "-----"
		<< "," << "-------\n";

	fillAry(ary); //fills array with random ints between 0 and 1

	/***Serial Method***/
	//tests serial method three times
	for (numTimes = 0; numTimes < TEST_COUNT; numTimes++) {
		startT = clock();

		for (int stall = 0; stall < STALL_CNT; stall++) {
			//count's threes in array
			numThrees = countThreesSerial(ary);
		}

		stopT = clock();
		testTime = (double)(stopT - startT) / CLOCKS_PER_SEC;

		//stores fastest time
		if (numTimes == 0) {
			fastestTime = testTime;
		}
		else {
			if (testTime < fastestTime) {
				fastestTime = testTime;
			}
		}
	}

	//prints results of serial method to console 
	cout << setw(8) << 1 << setw(COL_WIDTH) << numThrees 
		<< setw(COL_WIDTH) << "Serial" << setw(COL_WIDTH) 
		<< 0 << setw(COL_WIDTH) <<fixed << setprecision(3)
		<< fastestTime << endl;

	//prints results of serial method to file
	outFile << 1 << "," << numThrees
		<< "," << "Serial" << ","
		<< 0 << "," << fixed << setprecision(3)
		<< fastestTime << endl;


	/***Parallel with 1D array***/
	//runs with number of threads of 2 - 5 inclusive
	for (threadNum = 2; threadNum <= MAX_THREADS; threadNum++) {
		//tests parallel with 1D ary 3 times for each thrd num
		for (numTimes = 0; numTimes < TEST_COUNT; numTimes++) {
			startT = clock();

			for (int stall = 0; stall < STALL_CNT; stall++) {
				//set number of threads we want to use
				omp_set_num_threads(threadNum);
				//count's threes in array
				numThreesPar = countThrees1D(ary);
			}

			stopT = clock();
			testTime = (double)(stopT - startT) / CLOCKS_PER_SEC;

			//stores fastest time
			if (numTimes == 0) {
				fastestTime = testTime;
			}
			else {
				if (testTime < fastestTime) {
					fastestTime = testTime;
				}
			}
		}

		//prints results for each number of threads to console
		cout << setw(8) << threadNum << setw(COL_WIDTH) << numThrees
			<< setw(COL_WIDTH) << "Ary1D" << setw(COL_WIDTH)
			<< numThrees - numThreesPar << setw(COL_WIDTH) << fixed
			<< setprecision(3) << fastestTime << endl;

		//prints results for each number of threads to file
		outFile << threadNum << "," << numThrees
			<< "," << "Ary1D" << ","
			<< numThrees - numThreesPar << "," << fixed
			<< setprecision(3) << fastestTime << endl;
	}


	/***Parallel with padded 2D array***/
	//runs with number of threads of 2 - 5 inclusive
	for (threadNum = 2; threadNum <= MAX_THREADS; threadNum++) {
		//tests parallel with 2D padded ary 3 times for each thrd num
		for (numTimes = 0; numTimes < TEST_COUNT; numTimes++) {
			startT = clock();

			for (int stall = 0; stall < STALL_CNT; stall++) {
				//set number of threads we want to use
				omp_set_num_threads(threadNum);
				//count's threes in array
				numThreesPar = countThreesPad(ary);
			}

			stopT = clock();
			testTime = (double)(stopT - startT) / CLOCKS_PER_SEC;

			//stores fastest time
			if (numTimes == 0) {
				fastestTime = testTime;
			}
			else {
				if (testTime < fastestTime) {
					fastestTime = testTime;
				}
			}
		}

		//prints results for each number of threads to console
		cout << setw(8) << threadNum << setw(COL_WIDTH) << numThrees 
			<< setw(COL_WIDTH)<< "Ary2DPad" << setw(COL_WIDTH) 
			<< numThrees - numThreesPar << setw(COL_WIDTH) << fixed 
			<< setprecision(3) << fastestTime << endl;

		//prints results for each number of threads to file
		outFile << threadNum << "," << numThrees
			<< "," << "Ary2DPad" << ","
			<< numThrees - numThreesPar << "," << fixed
			<< setprecision(3) << fastestTime << endl;
	}


	/***Parallel with critical section***/
	//runs with number of threads of 2 - 5 inclusive
	for (threadNum = 2; threadNum <= MAX_THREADS; threadNum++) {
		//tests parallel with critical 3 times for each thrd num
		for (numTimes = 0; numTimes < TEST_COUNT; numTimes++) {
			startT = clock();

			for (int stall = 0; stall < STALL_CNT; stall++) {
				//set number of threads we want to use
				omp_set_num_threads(threadNum);
				//count's threes in array
				numThreesPar = countThreesCritical(ary);
			}

			stopT = clock();
			testTime = (double)(stopT - startT) / CLOCKS_PER_SEC;

			//stores fastest time
			if (numTimes == 0) {
				fastestTime = testTime;
			}
			else {
				if (testTime < fastestTime) {
					fastestTime = testTime;
				}
			}
		}

		//prints results for each number of threads to console
		cout << setw(8) << threadNum << setw(COL_WIDTH) << numThrees
			<< setw(COL_WIDTH) << "Critical" << setw(COL_WIDTH)
			<< numThrees - numThreesPar << setw(COL_WIDTH) << fixed
			<< setprecision(3) << fastestTime << endl;

		//prints results for each number of threads to file
		outFile << threadNum << "," << numThrees
			<< "," << "Critical" << ","
			<< numThrees - numThreesPar << "," << fixed
			<< setprecision(3) << fastestTime << endl;
	}


	/***Parallel with atomic section***/
	//runs with number of threads of 2 - 5 inclusive
	for (threadNum = 2; threadNum <= MAX_THREADS; threadNum++) {
		//tests parallel with atomic 3 times for each thrd num
		for (numTimes = 0; numTimes < TEST_COUNT; numTimes++) {
			startT = clock();

			for (int stall = 0; stall < STALL_CNT; stall++) {
				//set number of threads we want to use
				omp_set_num_threads(threadNum);
				//count's threes in array
				numThreesPar = countThreesAtomic(ary);
			}

			stopT = clock();
			testTime = (double)(stopT - startT) / CLOCKS_PER_SEC;

			//stores fastest time
			if (numTimes == 0) {
				fastestTime = testTime;
			}
			else {
				if (testTime < fastestTime) {
					fastestTime = testTime;
				}
			}
		}

		//prints results for each number of threads to console
		cout << setw(8) << threadNum << setw(COL_WIDTH) << numThrees
			<< setw(COL_WIDTH) << "Atomic" << setw(COL_WIDTH)
			<< numThrees - numThreesPar << setw(COL_WIDTH) << fixed
			<< setprecision(3) << fastestTime << endl;

		//prints results for each number of threads to file
		outFile << threadNum << "," << numThrees
			<< "," << "Atomic" << ","
			<< numThrees - numThreesPar << "," << fixed
			<< setprecision(3) << fastestTime << endl;
	}


	/***Parallel with parallel for***/
	//runs with number of threads of 2 - 5 inclusive
	for (threadNum = 2; threadNum <= MAX_THREADS; threadNum++) {
		//tests parallel for 3 times for each thrd num
		for (numTimes = 0; numTimes < TEST_COUNT; numTimes++) {
			startT = clock();

			for (int stall = 0; stall < STALL_CNT; stall++) {
				//set number of threads we want to use
				omp_set_num_threads(threadNum);
				//count's threes in array
				numThreesPar = countThreesFor(ary);
			}

			stopT = clock();
			testTime = (double)(stopT - startT) / CLOCKS_PER_SEC;

			//stores fastest time
			if (numTimes == 0) {
				fastestTime = testTime;
			}
			else {
				if (testTime < fastestTime) {
					fastestTime = testTime;
				}
			}
		}

		//prints results for each number of threads to console
		cout << setw(8) << threadNum << setw(COL_WIDTH) << numThrees
			<< setw(COL_WIDTH) << "ParallelFor" << setw(COL_WIDTH)
			<< numThrees - numThreesPar << setw(COL_WIDTH) << fixed
			<< setprecision(3) << fastestTime << endl;

		//prints results for each number of threads to file
		outFile << threadNum << "," << numThrees
			<< "," << "ParallelFor" << ","
			<< numThrees - numThreesPar << "," << fixed
			<< setprecision(3) << fastestTime << endl;
	}

	//Finish up program
	delete[] ary;
	cout << "\n\n***Normal Termination***" << endl;
	return 0;
}

/*Load array with random ints between 0 and 100
 *Takes in array to be filled with random ints */
void fillAry(int a[MAX]) {
	int index; //loop variable
	srand((unsigned)time(0));	//seed rand()

	//fills array
	for (index = 0; index < MAX; index++) {
		a[index] = rand() % 100;
	}
}

/*count number of threes serially
 *takes in array from which threes are counted
 *returns number of threes found*/
int countThreesSerial(int a[MAX]) {
	int numThrees = 0;	//the count

	//count number of threes
	for (int index = 0; index < MAX; index++) {
		if (a[index] == 3) {
			numThrees++;
		}
	}

	return numThrees;
}

/*count number of threes in parallel using 1D array
 *takes in array from which threes are counted
 *returns number of threes found*/
int countThrees1D(int a[MAX]) {
	int numThrees[MAX_THREADS] = { 0 }; //1D ary for thrds
	int numActualThreads = 0; //num thrds we get

#pragma omp parallel 
{
		//variables for SPMD (local to each thread)
		int myID, numThreads, myStart, myStop;

		//SPMD
		myID = omp_get_thread_num();
		numThreads = omp_get_num_threads();

		myStart = (long)(myID * MAX) / numThreads;
		myStop = (long)((myID + 1) * MAX) / numThreads;

		//saves actual thread count
		if (myID == 0) {
			numActualThreads = numThreads;
		}

		//count number of threes
		for (int index = myStart; index < myStop; index++) {
			if (a[index] == 3) {
				numThrees[myID]++;
			}
		}
	} //end of parallel section

	//add each thread's count to thread 0's count
	for (int i = 1; i < numActualThreads; i++) {
		numThrees[0] += numThrees[i];
	}
	return numThrees[0];
}

/*count number of threes in parallel using 2D padded array
 *takes in array from which threes are counted
 *returns number of threes found*/
int countThreesPad(int a[MAX]) {
	int numThrees[MAX_THREADS][16] = { 0 }; //padded ary
	int numActualThreads = 0; //num thrds we get

#pragma omp parallel 
	{
		//variables for SPMD (local to each thread)
		int myID, numThreads, myStart, myStop;

		//SPMD
		myID = omp_get_thread_num();
		numThreads = omp_get_num_threads();

		myStart = (long)(myID * MAX) / numThreads;
		myStop = (long)((myID + 1) * MAX) / numThreads;

		//saves actual thread count
		if (myID == 0) {
			numActualThreads = numThreads;
		}

		//count number of threes
		for (int index = myStart; index < myStop; index++) {
			if (a[index] == 3) {
				numThrees[myID][0]++;
			}
		}
	} //end of parallel section

	//add each thread's count to thread 0's count
	for (int i = 1; i < numActualThreads; i++) {
		numThrees[0][0] += numThrees[i][0];
	}

	return numThrees[0][0];
}

/*count number of threes in parallel using critical section
 *takes in array from which threes are counted
 *returns number of threes found*/
int countThreesCritical(int a[MAX]) {
	int numThrees = 0;	//the count


	//count number of threes
#pragma omp parallel for
		for (int i = 0; i < MAX; i++) {
			if (a[i] == 3) {
#		pragma omp critical
				numThrees++;
			}
		}
	//end of parallel section

	return numThrees;
}

/*count number of threes in parallel using atomic section
 *takes in array from which threes are counted
 *returns number of threes found*/
int countThreesAtomic(int a[MAX]) {
	int numThrees = 0;	//the count


	//count number of threes
#pragma omp parallel for
		for (int i = 0; i < MAX; i++) {
			if (a[i] == 3) {
#		pragma omp atomic
				numThrees++;
			}
		}
	//end of parallel section

	return numThrees;
}

/*count number of threes in array using parallel for
 *takes in array from which threes are counted
 *returns number of threes found*/
int countThreesFor(int a[MAX]) {
	int numThrees = 0;	//the count

	
	//count number of threes
#pragma omp parallel for reduction(+: numThrees)
	for (int index = 0; index < MAX; index++) {
		if (a[index] == 3) {
			numThrees++;
		}
	}
	//end of parallel section

	return numThrees;
}
