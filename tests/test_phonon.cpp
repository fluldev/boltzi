#include <iostream>
#include "boltzi_types.h"
#include "phonon.h"


int main() {
    using namespace std;
    using namespace boltzi;
    
    PhononData ph(
        vec3{1,2,3},
        vec3{4,5,6},
        7,
        8,
        1
    );

    ph.r.print();

    PhononList lst;
    lst.add(ph);
    lst.phonons.front().r.print();
    lst.add(ph);
    lst.remove(lst.phonons.begin());
    lst.add(ph);
    cout<<lst.phonons.size()<<endl;

    return 0;
}
