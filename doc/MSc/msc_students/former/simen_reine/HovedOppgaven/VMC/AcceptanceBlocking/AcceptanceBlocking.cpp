#include <iostream>
#include <string>
#include <blitz/array.h>   

using namespace std;
using namespace blitz;
 
int main(int argn, char* args[]) {
 
  int numRows;
  cin >> numRows;

  int numColumns;
  cin >> numColumns;

  int blockSize;
  cin >> blockSize;


  Array<double,2> A;
  A.resize(numRows,numColumns);
  for (int i=0; i<numRows; i++)
    for (int j=0; j<numColumns; j++)
      cin >> A(i,j);

  int numBlocks = numRows/blockSize/3;
  Array<double,2> Histogram;
  Histogram.resize(blockSize+1,4);
  Histogram = 0;

  for (int column=1; column<4; column++) {
    double temp = 0;
    double temp2;
    int row = 0;
    for (int i=0; i<numBlocks; i++) {
      int BlockAccept = 0;
      for (int j=0; j<blockSize; j++) {
	temp2 = A(row, column);
	row += 3;
	if ( temp2 != temp) {
	  BlockAccept++;
	  temp = temp2;
	}
      }
      Histogram(BlockAccept, column-1) += 1;
    }
  }
  /*
  Histogram /= numBlocks;
  cerr << endl;
  for (int i=0; i<blockSize+1; i++) {
    cerr << ((double)i)/blockSize << "   ";
    for (int column=1; column<4; column++)
      cerr << Histogram(i, column-1) << " ";
    cerr << endl;
  }


  Histogram = 0;
  */


  double lowestEnergy=9999999;
  double highestEnergy=-9999999;
  for (int i=0; i<numRows; i++) {
    if (A(i, 4) < lowestEnergy) 
      lowestEnergy = A(i, 4);
    if (A(i, 4) > highestEnergy) 
       highestEnergy = A(i, 4);
    
  }

  double blockLength = (highestEnergy-lowestEnergy)/blockSize;
  for (int i=0; i<numRows; i++) {
    int histIndex = (int) ((A(i, 4)-lowestEnergy)/blockLength);
    //cerr << histIndex << " " << blockSize << endl;
    Histogram( histIndex , 3) += 1;
  }


  Range Row(0,blockSize), Column(0,2);
  Histogram(Row, Column) /= numBlocks;
  Histogram(Row, 3) /= numRows/(highestEnergy-lowestEnergy);

  cerr << endl;
  cerr << endl;
  for (int i=0; i<blockSize+1; i++) {
    cerr << ((double)i)/blockSize << "   ";
    for (int column=1; column<4; column++)
      cerr << Histogram(i, column-1) << " ";

    cerr << lowestEnergy + (0.5 + i)*blockLength <<  " "
	 << Histogram(i, 3) << endl;
  }

  return 0;
}
