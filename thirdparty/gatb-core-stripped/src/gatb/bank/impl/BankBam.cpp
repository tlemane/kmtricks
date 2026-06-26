/*****************************************************************************
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#include <gatb/bank/impl/BankBam.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/api/StringsRepository.hpp>

#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;
using namespace gatb::core::tools::dp;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;
using namespace gatb::core::tools::misc;

#define DEBUG(a)  //printf a

/********************************************************************************/
namespace gatb {  namespace core {  namespace bank {  namespace impl {
/********************************************************************************/

/********************************************************************************/
// BGZF (Blocked GZIP Format) implementation
// BAM files use BGZF compression which is compatible with gzip but uses
// fixed-size blocks (max 64KB) for random access
/********************************************************************************/

static const int BGZF_BLOCK_SIZE = 65536;  // 64KB
static const int BGZF_MAX_BLOCK_SIZE = 65536;

/** BGZF file handle structure */
typedef struct {
    FILE* file;
    unsigned char* compressed_block;
    unsigned char* uncompressed_block;
    int block_length;
    int block_offset;
    bool eof;
    bool is_write;
} bgzf_t;

/********************************************************************************/
/** Read BGZF block header and return compressed size
 *  Returns -1 on EOF, -2 on error, or compressed block size
 */
static int bgzf_read_block_header(FILE* fp, int* block_length)
{
    unsigned char header[18];

    // Read the 18-byte BGZF header
    if (fread(header, 1, 18, fp) != 18) {
        if (feof(fp)) return -1;  // EOF
        return -2;  // Error
    }

    // Check BGZF magic bytes: 31, 139 (gzip magic)
    if (header[0] != 31 || header[1] != 139) {
        return -2;  // Not a valid gzip block
    }

    // Check extra flags
    if (header[3] != 4) {  // Should have FEXTRA flag set
        return -2;
    }

    // Get block size from BSIZE field (little-endian)
    // BSIZE is at bytes 16-17 and represents (block_size - 1)
    *block_length = ((int)header[17] << 8) | (int)header[16];
    *block_length += 1;  // BSIZE is block_size - 1

    return *block_length - 18;  // Return size of compressed data (excluding header)
}

/********************************************************************************/
/** Decompress one BGZF block */
static int bgzf_read_block(bgzf_t* fp)
{
    int block_length;
    int compressed_size = bgzf_read_block_header(fp->file, &block_length);

    if (compressed_size == -1) {
        fp->eof = true;
        return 0;
    }

    if (compressed_size < 0) {
        return -1;  // Error
    }

    // Read compressed data (excluding header which we already read)
    int remaining = block_length - 18;
    if (fread(fp->compressed_block, 1, remaining, fp->file) != (size_t)remaining) {
        return -1;
    }

    // Decompress using zlib
    z_stream zs;
    zs.zalloc = NULL;
    zs.zfree = NULL;
    zs.next_in = fp->compressed_block;
    zs.avail_in = compressed_size;
    zs.next_out = fp->uncompressed_block;
    zs.avail_out = BGZF_BLOCK_SIZE;

    // Use raw deflate (negative window bits for raw deflate without zlib header)
    if (inflateInit2(&zs, -15) != Z_OK) {
        return -1;
    }

    if (inflate(&zs, Z_FINISH) != Z_STREAM_END) {
        inflateEnd(&zs);
        return -1;
    }

    fp->block_length = zs.total_out;
    fp->block_offset = 0;

    inflateEnd(&zs);

    return fp->block_length;
}

/********************************************************************************/
/** Open BGZF file for reading */
static bgzf_t* bgzf_open(const char* path, const char* mode)
{
    bgzf_t* fp = (bgzf_t*)calloc(1, sizeof(bgzf_t));
    if (!fp) return NULL;

    fp->file = fopen(path, mode);
    if (!fp->file) {
        free(fp);
        return NULL;
    }

    fp->compressed_block = (unsigned char*)malloc(BGZF_MAX_BLOCK_SIZE);
    fp->uncompressed_block = (unsigned char*)malloc(BGZF_BLOCK_SIZE);

    if (!fp->compressed_block || !fp->uncompressed_block) {
        if (fp->compressed_block) free(fp->compressed_block);
        if (fp->uncompressed_block) free(fp->uncompressed_block);
        fclose(fp->file);
        free(fp);
        return NULL;
    }

    fp->block_length = 0;
    fp->block_offset = 0;
    fp->eof = false;
    fp->is_write = (mode[0] == 'w');

    return fp;
}

/********************************************************************************/
/** Close BGZF file */
static void bgzf_close(bgzf_t* fp)
{
    if (!fp) return;

    if (fp->file) fclose(fp->file);
    if (fp->compressed_block) free(fp->compressed_block);
    if (fp->uncompressed_block) free(fp->uncompressed_block);
    free(fp);
}

/********************************************************************************/
/** Read from BGZF file */
static int bgzf_read(bgzf_t* fp, void* data, int length)
{
    if (!fp || length < 0) return -1;
    if (length == 0) return 0;

    unsigned char* output = (unsigned char*)data;
    int bytes_read = 0;

    while (bytes_read < length) {
        // If current block is exhausted, read next block
        if (fp->block_offset >= fp->block_length) {
            if (fp->eof) break;

            if (bgzf_read_block(fp) < 0) {
                return -1;
            }

            if (fp->eof) break;
        }

        // Copy data from current block
        int available = fp->block_length - fp->block_offset;
        int to_copy = (length - bytes_read < available) ? length - bytes_read : available;

        memcpy(output + bytes_read, fp->uncompressed_block + fp->block_offset, to_copy);
        fp->block_offset += to_copy;
        bytes_read += to_copy;
    }

    return bytes_read;
}

/********************************************************************************/
// BAM format parsing
/********************************************************************************/

/** Read little-endian 32-bit integer */
static inline uint32_t read_uint32_le(const unsigned char* buf)
{
    return (uint32_t)buf[0] | ((uint32_t)buf[1] << 8) |
           ((uint32_t)buf[2] << 16) | ((uint32_t)buf[3] << 24);
}

/** Read little-endian 16-bit integer */
static inline uint16_t read_uint16_le(const unsigned char* buf)
{
    return (uint16_t)buf[0] | ((uint16_t)buf[1] << 8);
}

/** Decode BAM sequence encoding to ASCII
 *  BAM uses 4-bit encoding: =:0, A:1, C:2, M:3, G:4, R:5, S:6, V:7,
 *                            T:8, W:9, Y:10, H:11, K:12, D:13, B:14, N:15
 */
static const char BAM_NT_DECODE[16] = {
    '=', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
    'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'
};

/********************************************************************************/
// BankBam implementation
/********************************************************************************/

BankBam::BankBam (const std::string& filename)
    : _filename(filename), _filesize(0)
{
    init();
}

BankBam::~BankBam ()
{
}

void BankBam::init ()
{
    // Check that file exists and get size
    if (!System::file().doesExist(_filename)) {
        throw gatb::core::system::Exception ("BAM file does not exist: %s", _filename.c_str());
    }

    _filesize = System::file().getSize(_filename);
}

void BankBam::estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize)
{
    BankBam::Iterator it(*this);
    it.estimate(number, totalSize, maxSize);
}

void BankBam::insert (const Sequence& item)
{
    // BAM writing not implemented
    throw gatb::core::system::Exception ("BankBam::insert not implemented (read-only)");
}

void BankBam::flush ()
{
    // Not applicable for read-only
}

/********************************************************************************/
// BankBam::Iterator implementation
/********************************************************************************/

BankBam::Iterator::Iterator (BankBam& ref)
    : _ref(ref), _isDone(true), _isInitialized(false), _nIters(0), _bgzf(NULL), _index(0)
{
    _item = new Sequence(Data::ASCII);
}

BankBam::Iterator::~Iterator ()
{
    finalize();
}

void BankBam::Iterator::init ()
{
    if (_isInitialized) return;

    // Open BGZF file
    _bgzf = bgzf_open(_ref._filename.c_str(), "rb");

    if (!_bgzf) {
        throw gatb::core::system::Exception ("Failed to open BAM file: %s", _ref._filename.c_str());
    }

    // Read and validate BAM header
    unsigned char magic[4];
    if (bgzf_read((bgzf_t*)_bgzf, magic, 4) != 4) {
        throw gatb::core::system::Exception ("Failed to read BAM magic bytes");
    }

    if (memcmp(magic, "BAM\1", 4) != 0) {
        throw gatb::core::system::Exception ("Invalid BAM file format (bad magic bytes)");
    }

    // Read header text length
    unsigned char l_text_buf[4];
    if (bgzf_read((bgzf_t*)_bgzf, l_text_buf, 4) != 4) {
        throw gatb::core::system::Exception ("Failed to read BAM header length");
    }
    uint32_t l_text = read_uint32_le(l_text_buf);

    // Skip header text
    if (l_text > 0) {
        char* text = (char*)malloc(l_text);
        if (!text) {
            throw gatb::core::system::Exception ("Memory allocation failed");
        }
        if (bgzf_read((bgzf_t*)_bgzf, text, l_text) != (int)l_text) {
            free(text);
            throw gatb::core::system::Exception ("Failed to read BAM header text");
        }
        free(text);
    }

    // Read number of reference sequences
    unsigned char n_ref_buf[4];
    if (bgzf_read((bgzf_t*)_bgzf, n_ref_buf, 4) != 4) {
        throw gatb::core::system::Exception ("Failed to read number of references");
    }
    uint32_t n_ref = read_uint32_le(n_ref_buf);

    // Read reference sequence information
    for (uint32_t i = 0; i < n_ref; i++) {
        unsigned char l_name_buf[4];
        if (bgzf_read((bgzf_t*)_bgzf, l_name_buf, 4) != 4) {
            throw gatb::core::system::Exception ("Failed to read reference name length");
        }
        uint32_t l_name = read_uint32_le(l_name_buf);

        // Read and store reference name
        char* name = (char*)malloc(l_name);
        if (!name) {
            throw gatb::core::system::Exception ("Memory allocation failed");
        }
        if (bgzf_read((bgzf_t*)_bgzf, name, l_name) != (int)l_name) {
            free(name);
            throw gatb::core::system::Exception ("Failed to read reference name");
        }
        // Store reference name (null-terminated, so l_name-1, used for contig filtering)
        _ref_names.push_back(std::string(name, l_name - 1));
        free(name);

        // Skip reference length
        unsigned char l_ref_buf[4];
        if (bgzf_read((bgzf_t*)_bgzf, l_ref_buf, 4) != 4) {
            throw gatb::core::system::Exception ("Failed to read reference length");
        }
    }

    _isInitialized = true;
    _isDone = false;
}

void BankBam::Iterator::finalize ()
{
    if (_bgzf) {
        bgzf_close((bgzf_t*)_bgzf);
        _bgzf = NULL;
    }
    _isInitialized = false;
}

void BankBam::Iterator::first()
{
    init();
    _index = 0;
    _nIters = 0;
    next();
}

void BankBam::Iterator::next()
{
    if (!_isInitialized) {
        _isDone = true;
        return;
    }

    _isDone = !get_next_seq(_item->getData(), _item->_comment);

    if (!_isDone) {
        _item->setIndex(_index++);
        _nIters++;
    } else {
        finalize();
    }
}

bool BankBam::Iterator::get_next_seq (tools::misc::Vector<char>& data, std::string& comment)
{
    bgzf_t* fp = (bgzf_t*)_bgzf;

    // Read alignment block size
    unsigned char block_size_buf[4];
    int bytes_read = bgzf_read(fp, block_size_buf, 4);

    if (bytes_read == 0) {
        return false;  // EOF
    }

    if (bytes_read != 4) {
        return false;  // Error or EOF
    }

    uint32_t block_size = read_uint32_le(block_size_buf);

    if (block_size == 0) {
        return false;  // Invalid or EOF marker
    }

    // Allocate buffer for alignment block
    unsigned char* block = (unsigned char*)malloc(block_size);
    if (!block) {
        return false;
    }

    // Read alignment block
    if (bgzf_read(fp, block, block_size) != (int)block_size) {
        free(block);
        return false;
    }

    // Parse BAM alignment block
    int32_t refID = read_uint32_le(block + 0);          // Reference sequence ID
    // uint32_t pos = read_uint32_le(block + 4);        // 0-based leftmost position (skip)
    uint8_t l_read_name = block[8];                     // Length of read name
    // uint8_t mapq = block[9];                         // Mapping quality (skip)
    // uint16_t bin = read_uint16_le(block + 10);       // BAI index bin (skip)
    uint16_t n_cigar_op = read_uint16_le(block + 12);  // Number of CIGAR operations
    uint16_t flag = read_uint16_le(block + 14);         // Bitwise flags
    uint32_t l_seq = read_uint32_le(block + 16);       // Length of sequence

    // Skip secondary (0x100) and supplementary (0x800) alignments to avoid double counting
    if (flag & 0x100 || flag & 0x800) {
        free(block);
        return get_next_seq(data, comment);  // Recursively get next primary read
    }

    // Apply samtools-style flag filtering
    // -f: require flags (all specified flags must be set)
    if (_ref._require_flags != 0) {
        if ((flag & _ref._require_flags) != _ref._require_flags) {
            free(block);
            return get_next_seq(data, comment);  // Skip: required flags not all set
        }
    }
    // -F: exclude flags (none of the specified flags must be set)
    if (_ref._exclude_flags != 0) {
        if (flag & _ref._exclude_flags) {
            free(block);
            return get_next_seq(data, comment);  // Skip: excluded flag is set
        }
    }

    // Skip reads aligned to excluded reference sequences
    if (refID >= 0 && refID < (int32_t)_ref_names.size()) {
        if (_ref._excluded_refs.find(_ref_names[refID]) != _ref._excluded_refs.end()) {
            free(block);
            return get_next_seq(data, comment);
        }
    }
    // uint32_t next_refID = read_uint32_le(block + 20);// Reference ID of next read (skip)
    // uint32_t next_pos = read_uint32_le(block + 24);  // Position of next read (skip)
    // uint32_t tlen = read_uint32_le(block + 28);      // Template length (skip)

    // Extract read name (null-terminated)
    char* read_name = (char*)malloc(l_read_name + 1);
    if (!read_name) {
        free(block);
        return false;
    }
    memcpy(read_name, block + 32, l_read_name);
    read_name[l_read_name] = '\0';

    // Remove trailing null if present
    if (l_read_name > 0 && read_name[l_read_name - 1] == '\0') {
        read_name[l_read_name - 1] = '\0';
    }

    comment = read_name;
    free(read_name);

    // Calculate offset to sequence data
    uint32_t seq_offset = 32 + l_read_name + (n_cigar_op * 4);

    // Decode sequence (4-bit encoding, 2 bases per byte)
    data.resize(l_seq);
    unsigned char* seq_data = block + seq_offset;

    for (uint32_t i = 0; i < l_seq; i++) {
        uint8_t byte_val = seq_data[i / 2];
        uint8_t base_code;

        if (i % 2 == 0) {
            base_code = (byte_val >> 4) & 0x0F;  // High nibble
        } else {
            base_code = byte_val & 0x0F;         // Low nibble
        }

        data[i] = BAM_NT_DECODE[base_code];
    }

    // If read is reverse complemented (flag 0x10), reverse complement it back
    // to recover the original read sequence (important for assembly/k-mer counting)
    if (flag & 0x10) {
        // Reverse complement the sequence
        for (uint32_t i = 0; i < l_seq / 2; i++) {
            char temp = data[i];
            data[i] = data[l_seq - 1 - i];
            data[l_seq - 1 - i] = temp;
        }

        // Complement each base
        for (uint32_t i = 0; i < l_seq; i++) {
            switch (data[i]) {
                case 'A': data[i] = 'T'; break;
                case 'T': data[i] = 'A'; break;
                case 'C': data[i] = 'G'; break;
                case 'G': data[i] = 'C'; break;
                // Leave ambiguous bases (N, R, Y, etc.) unchanged
            }
        }
    }

    free(block);
    return true;
}

void BankBam::Iterator::estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize)
{
    number = 0;
    totalSize = 0;
    maxSize = 0;

    // Sample first N reads to estimate
    const int SAMPLE_SIZE = 10000;
    int count = 0;

    for (first(); !isDone() && count < SAMPLE_SIZE; next(), count++) {
        u_int64_t seqSize = item().getDataSize();
        totalSize += seqSize;
        if (seqSize > maxSize) {
            maxSize = seqSize;
        }
        number++;
    }

    if (number > 0) {
        // Extrapolate based on file size
        u_int64_t bytes_read = _nIters * 100;  // Rough estimate
        double ratio = (double)_ref._filesize / (double)bytes_read;
        number = (u_int64_t)(number * ratio);
        totalSize = (u_int64_t)(totalSize * ratio);
    }
}

/********************************************************************************/
// BankBamFactory implementation
/********************************************************************************/

IBank* BankBamFactory::createBank (const std::string& uri)
{
    // Check if this is a BAM file by reading magic bytes
    // BAM files use BGZF compression and contain "BAM\1" after decompression

    bgzf_t* fp = bgzf_open(uri.c_str(), "rb");
    if (!fp) {
        return NULL;  // Cannot open file
    }

    // Read BAM magic bytes
    unsigned char magic[4];
    int bytes_read = bgzf_read(fp, magic, 4);
    bgzf_close(fp);

    if (bytes_read == 4 && memcmp(magic, "BAM\1", 4) == 0) {
        // Valid BAM file
        try {
            return new BankBam(uri);
        } catch (...) {
            return NULL;
        }
    }

    // Not a valid BAM file
    return NULL;
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
