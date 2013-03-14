#include <iostream>
#include <limits>
#include <unistd.h>

using namespace std;
// 
// Example 1 : ./a.out 2 
//          Job : start
//          .. after two seconds
//          Job : End
//
// Example 2 :  ./a.out 2 A
//          Job A : start
//          .. after two seconds
//          Job A : End

int main(int argc, char *argv[]){

  int second = atoi(argv[1]);
  // If this job has a name, it will be printed.
  if(argc == 3){
    printf("Job %s : start\n",argv[2]);
  }else{
    printf("Job : start\n");
  }
  printf("sleep for %s second\n",argv[1]);
  sleep(second);
  if(argc == 3){
    printf("Job %s : End\n",argv[2]);
  }else{
    printf("Job : End\n");
  }
  
  return 1;
}
