#ifndef __G4SBS_TYPES_H
#define __G4SBS_TYPES_H

////////////////////////////////////////////////////////
//  Data for extracting things from GEMC
//
//  we'll hardcode them here, but it would be nice to
//  maybe get them into a database
//  I guess this could also be done through a mysql
//  interface, but I think that makes it more complicated
//  and breakable


#define NBANKS 1

// Tag numbers associated in the GEMC banks
#define __GENERATED_TAG  10
#define __CER_TAG  110

#define __GENERATED_SIZE 7

static int __g4sbs_types_datasize[NBANKS] = {21};

// FIXME:  Need to do this better,  map?
static int data_size(int tag){
    if( tag == __CER_TAG ) {
	return __g4sbs_types_datasize[0];
    }
    return 0;
}
#endif//__GEMC_TYPES_H
