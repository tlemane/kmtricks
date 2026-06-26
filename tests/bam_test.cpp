#include <gtest/gtest.h>
#include <gatb/gatb_core.hpp>

using namespace gatb::core::bank;

class BamTest : public ::testing::Test {
protected:
    void SetUp() override {
        test_fasta = "data/1.fasta";
        test_bam = "../test.bam";
    }

    std::string test_fasta;
    std::string test_bam;
};

TEST_F(BamTest, FormatDetection) {
    // Test that BAM and FASTA files are correctly detected
    EXPECT_EQ(Bank::getType(test_bam), "bam");
    EXPECT_EQ(Bank::getType(test_fasta), "fasta");
}

TEST_F(BamTest, SequenceEquality) {
    // Test that BAM and FASTA return identical sequences
    IBank* bam_bank = Bank::open(test_bam);
    IBank* fasta_bank = Bank::open(test_fasta);

    ASSERT_NE(bam_bank, nullptr);
    ASSERT_NE(fasta_bank, nullptr);

    gatb::core::tools::dp::Iterator<Sequence>* bam_it = bam_bank->iterator();
    gatb::core::tools::dp::Iterator<Sequence>* fasta_it = fasta_bank->iterator();

    int count = 0;
    for (bam_it->first(), fasta_it->first();
         !bam_it->isDone() && !fasta_it->isDone();
         bam_it->next(), fasta_it->next()) {
        count++;
        EXPECT_EQ(bam_it->item().toString(), fasta_it->item().toString())
            << "Sequence " << count << " should match";
    }

    EXPECT_EQ(count, 2);
    EXPECT_TRUE(bam_it->isDone() && fasta_it->isDone());

    delete bam_it;
    delete fasta_it;
    delete bam_bank;
    delete fasta_bank;
}
