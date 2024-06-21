#ifndef DECLARE_STRUCT
#define DECLARE_STRUCT

/*
CHANGE LOG
   July 3,2018: Started this code
OBJECTIVE
   Declare the Dump Structure which have 2 members:
   time: trajectory timestamp
   dumpName: Name of the trajectory file
    
*/

struct Dump{
    float time;
    char dumpName[256];
};

#endif
