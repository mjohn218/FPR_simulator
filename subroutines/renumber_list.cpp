#include "reactions.h"
#include "utility_calls.h"

void renumber_list(int &nfree, int *freelist) {
	/*routine to sort reactants numerically to count
	 reactant species without repeats
	 */
	int i, j;
	cout << "ncount, with repeats : " << nfree << endl;
	for (i = 1; i < nfree; i++) {
		for (j = 0; j < i; j++) {
			if (freelist[j] == freelist[i]) {
				//then this interface is repeated
				nfree--;
				freelist[i] = freelist[nfree];
				j = i;
				i -= 1;
			}
		}
	}
	//now sort the interfaces numerically
	int num = nfree + 1;
	long unsigned int *index = (long unsigned int *) malloc(sizeof(long unsigned int) * (num));
	double *value = new double[num];
	for (i = 0; i < nfree; i++)
		value[i + 1] = freelist[i];
	indexx(nfree, value, index);
	int id;
	for (i = 0; i < nfree; i++) {
		id = index[i + 1] - 1; //indexing is from 1:N, skips zero
		value[i] = freelist[id];
	}

	for (i = 0; i < nfree; i++)
		freelist[i] = int(value[i]); //recopy to original array

	cout << "ncount: " << nfree << endl;
	for (i = 0; i < nfree; i++) {
		cout << freelist[i] << endl;
	}
	free(index);
	delete[] value;

}
