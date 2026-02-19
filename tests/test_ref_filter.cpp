#include <gatb/gatb_core.hpp>
#include <iostream>
#include <set>

using namespace gatb::core::bank;

int main() {
    std::string test_bam = "../test.bam";
    
    // First, read without filtering
    std::cout << "Reading without filtering:" << std::endl;
    IBank* bank1 = Bank::open(test_bam);
    gatb::core::tools::dp::Iterator<Sequence>* it1 = bank1->iterator();
    int count1 = 0;
    for (it1->first(); !it1->isDone(); it1->next()) {
        count1++;
        std::cout << "  Read " << count1 << ": " << it1->item().getComment() << std::endl;
    }
    delete it1;
    delete bank1;
    
    std::cout << "\nTotal reads without filtering: " << count1 << std::endl;
    
    return 0;
}
