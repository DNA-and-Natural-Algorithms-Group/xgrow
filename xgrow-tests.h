

#ifndef __XGROW_TESTS_H__
#define __XGROW_TESTS_H__

#include "grow.h"

typedef unsigned char ** Assembly;

void run_xgrow_tests (tube *tp,double Gmc, double Gse, int seed_i, int seed_j, int seed_n, int size);
int assembly_is_a_duplicate (void *assemblies, 
			     Assembly a, int size);
void add_assembly_to_seen (tube *tp);
void remove_assembly_from_seen (tube *tp);

void update_state_on_indicator(flake *fp, Assembly a, int size);
void update_state_off_indicator(flake *fp);
#endif
