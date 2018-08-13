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
*/
	return dim_a*(b-1) + a;
}
/*
/
*/
int index_four(int a, int b, int c, int d, int dim_a, int dim_b, int dim_c){
/*
/
/ 	Index four
/ 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
/
/ 	abcd = index_two(ab,cd) = dim_ab*(cd-1)+ab
/									= dim_a*(dim_b*(dim_c*(d-1) + c - 1) + (b-1)) + a
/
*/
	return dim_a*(dim_b*(dim_c*(d-1)+c-1)+(b-1))+a;
}
