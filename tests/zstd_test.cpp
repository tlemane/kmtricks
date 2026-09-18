#include <gtest/gtest.h>
#include <gatb/gatb_core.hpp>

#include <zstd.h>

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace gatb::core::bank;
using gatb::core::tools::dp::Iterator;

namespace {

std::string read_file (const std::string& path)
{
    std::ifstream ifs (path, std::ios::binary);
    std::ostringstream oss;
    oss << ifs.rdbuf();
    return oss.str();
}

/** Compress `content` into `nframes` zstd frames written back to back. Anything
 *  beyond one frame exercises the concatenated-archive path in BankZSTD::rebuffer,
 *  which is what `zstd -T` / pzstd produce. */
void write_zstd (const std::string& path, const std::string& content, int nframes = 1)
{
    std::string out;
    size_t offset = 0;

    for (int f = 0; f < nframes; f++)
    {
        size_t len = (f == nframes - 1) ? content.size() - offset : content.size() / nframes;
        size_t cap = ZSTD_compressBound (len);

        std::vector<char> buf (cap);
        size_t n = ZSTD_compress (buf.data(), cap, content.data() + offset, len, 3);
        ASSERT_FALSE (ZSTD_isError (n)) << ZSTD_getErrorName (n);

        out.append (buf.data(), n);
        offset += len;
    }

    std::ofstream ofs (path, std::ios::binary);
    ofs.write (out.data(), out.size());
    ASSERT_TRUE (ofs.good());
}

/** A record is (comment, sequence, quality). */
struct Record
{
    std::string comment, sequence, quality;

    bool operator== (const Record& o) const
    {
        return comment == o.comment && sequence == o.sequence && quality == o.quality;
    }
};

std::vector<Record> collect (IBank* bank)
{
    std::vector<Record> records;

    Iterator<Sequence>* it = bank->iterator();

    for (it->first(); !it->isDone(); it->next())
    {
        records.push_back (Record {
            it->item().getComment(), it->item().toString(), it->item().getQuality()
        });
    }

    delete it;

    return records;
}

/** A few FASTQ records, since tests/data holds no fastq fixture and we want the
 *  quality path covered. Deliberately long enough to span several zstd frames. */
std::string make_fastq (int nb_records)
{
    std::ostringstream oss;

    for (int i = 0; i < nb_records; i++)
    {
        std::string seq (100, 'A' + (i % 4 == 0 ? 0 : (i % 4 == 1 ? 2 : (i % 4 == 2 ? 6 : 19))));
        oss << "@read_" << i << " some description\n"
            << seq << "\n"
            << "+\n"
            << std::string (seq.size(), 'I') << "\n";
    }

    return oss.str();
}

} /* end of anonymous namespace */

class ZstdTest : public ::testing::Test
{
protected:
    void SetUp () override
    {
        plain_fasta = "data/1.fasta";

        fasta_zst      = "zstd_test_1.fasta.zst";
        fasta_multi_zst = "zstd_test_1.multiframe.fasta.zst";
        fastq_zst      = "zstd_test_1.fastq.zst";
        garbage_zst    = "zstd_test_garbage.zst";
        truncated_zst  = "zstd_test_truncated.zst";

        std::string fasta = read_file (plain_fasta);
        ASSERT_FALSE (fasta.empty()) << "missing fixture " << plain_fasta;

        write_zstd (fasta_zst, fasta);
        write_zstd (fasta_multi_zst, fasta, 3);

        fastq = make_fastq (50);
        write_zstd (fastq_zst, fastq);

        write_zstd (garbage_zst, "this is not a bank, just some prose.\n");

        /** Half a valid archive: the frame never completes. */
        std::string whole = read_file (fasta_zst);
        std::ofstream ofs (truncated_zst, std::ios::binary);
        ofs.write (whole.data(), whole.size() / 2);
    }

    void TearDown () override
    {
        for (const std::string& p : { fasta_zst, fasta_multi_zst, fastq_zst, garbage_zst, truncated_zst })
        {
            std::remove (p.c_str());
        }
    }

    std::string plain_fasta;
    std::string fasta_zst, fasta_multi_zst, fastq_zst, garbage_zst, truncated_zst;
    std::string fastq;
};

TEST_F (ZstdTest, FormatDetection)
{
    EXPECT_EQ (Bank::getType (fasta_zst), "zstd");
    EXPECT_EQ (Bank::getType (fastq_zst), "zstd");

    /** Registering the zstd factory must not disturb the existing formats. */
    EXPECT_EQ (Bank::getType (plain_fasta), "fasta");
}

TEST_F (ZstdTest, SequencesMatchPlainFasta)
{
    IBank* zst   = Bank::open (fasta_zst);    LOCAL (zst);
    IBank* plain = Bank::open (plain_fasta);  LOCAL (plain);

    std::vector<Record> from_zst   = collect (zst);
    std::vector<Record> from_plain = collect (plain);

    ASSERT_FALSE (from_plain.empty());
    ASSERT_EQ (from_zst.size(), from_plain.size());

    for (size_t i = 0; i < from_plain.size(); i++)
    {
        EXPECT_TRUE (from_zst[i] == from_plain[i]) << "record " << i << " differs";
    }
}

TEST_F (ZstdTest, ConcatenatedFramesReadAsOneStream)
{
    IBank* single = Bank::open (fasta_zst);        LOCAL (single);
    IBank* multi  = Bank::open (fasta_multi_zst);  LOCAL (multi);

    std::vector<Record> from_single = collect (single);
    std::vector<Record> from_multi  = collect (multi);

    ASSERT_FALSE (from_single.empty());
    ASSERT_EQ (from_multi.size(), from_single.size());

    for (size_t i = 0; i < from_single.size(); i++)
    {
        EXPECT_TRUE (from_multi[i] == from_single[i]) << "record " << i << " differs";
    }
}

TEST_F (ZstdTest, FastqQualityIsPreserved)
{
    IBank* bank = Bank::open (fastq_zst);  LOCAL (bank);

    std::vector<Record> records = collect (bank);

    ASSERT_EQ (records.size(), 50u);

    for (size_t i = 0; i < records.size(); i++)
    {
        EXPECT_EQ (records[i].comment, "read_" + std::to_string(i) + " some description");
        EXPECT_EQ (records[i].sequence.size(), 100u);
        EXPECT_EQ (records[i].quality, std::string (100, 'I'));
    }
}

TEST_F (ZstdTest, IteratingTwiceYieldsTheSameRecords)
{
    IBank* bank = Bank::open (fasta_zst);  LOCAL (bank);

    /** first() must reset the decompression stream, not resume it. */
    std::vector<Record> pass1 = collect (bank);
    std::vector<Record> pass2 = collect (bank);

    ASSERT_FALSE (pass1.empty());
    ASSERT_EQ (pass2.size(), pass1.size());

    for (size_t i = 0; i < pass1.size(); i++)
    {
        EXPECT_TRUE (pass2[i] == pass1[i]) << "record " << i << " differs";
    }
}

TEST_F (ZstdTest, NonFastaPayloadIsNotABank)
{
    /** Valid zstd, but the payload is not FASTA/FASTQ: no factory should claim it. */
    EXPECT_EQ (Bank::getType (garbage_zst), "unknown");
    EXPECT_ANY_THROW (delete Bank::open (garbage_zst));
}

TEST_F (ZstdTest, TruncatedArchiveThrows)
{
    /** Whether it fails at open() or mid-iteration depends on where the cut lands;
     *  what matters is that we never silently return a truncated bank. */
    EXPECT_ANY_THROW ({
        IBank* bank = Bank::open (truncated_zst);
        LOCAL (bank);
        collect (bank);
    });
}
