/*
/
/ 	Interface for utilities
/ 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
/
*/
int index_two(int a, int b, int dim_a){
/*
/
/ 	Index two
/ 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
/
/ 	C++ version for pointers has -1, since the pointer refers to index 1,
/ 	which is accessed as *(h + index_two(1, 1, dim_a))
/
*/
	return -1 + dim_a*(b-1) + a;
}
