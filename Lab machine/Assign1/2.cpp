#include <iostream>
#include <fstream>
using namespace std;

int main()
{
    
    // double* stars = new double[number_of_stars * 3];
    
  

  string line1[30];
  string line2[30];
  ifstream myfile("stars-9.txt");
  int a = 0;
  int b = 0;
  if(!myfile) 
  {
    cout<<"Error opening output file"<<endl;
    //system("pause");
    return -1;
  }
  while(!myfile.eof())
  {
    getline(myfile,line1[a],'\n');
    cout<<"1."<<line1[a]<<"\n";
    getline(myfile,line2[b],'\n');
    cout<<"2."<<line2[b]<<"\n";
  }

   n_file =    .size();

    for (int i=0, i < n_file, i++){

     
   }



    return 0;
}

