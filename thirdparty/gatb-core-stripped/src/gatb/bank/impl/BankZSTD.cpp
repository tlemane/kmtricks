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

#include <gatb/bank/impl/BankZSTD.hpp>

#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/api/StringsRepository.hpp>

#include <algorithm>

#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>

#include <zstd.h>

using namespace std;
using namespace gatb::core::tools::dp;
using namespace gatb::core::system;
using namespace gatb::core::system::impl;
using namespace gatb::core::tools::misc;

#define DEBUG(a)  //printf a

#define BUFFER_SIZE     (256*1024)

#define nearest_power_of_2(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

/** https://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2 */
#define is_power_of_2(v)  (((v) & ((v) - 1)) == 0)

/********************************************************************************/
namespace gatb {  namespace core {  namespace bank {  namespace impl {
/********************************************************************************/

namespace {

struct zstd_file_t
{
    FILE*           file;
    ZSTD_DStream*   dstream;

    /** Staging area for the compressed bytes read from the file. */
    unsigned char*  in_buffer;
    size_t          in_capacity;
    ZSTD_inBuffer   input;

    /** Decompressed bytes handed out to the parser. */
    unsigned char*  buffer;
    uint64_t        buffer_start, buffer_end;

    /** Total number of decompressed bytes produced so far, used by estimate(). */
    uint64_t        out_total;

    /** Return code of the last ZSTD_decompressStream call. Non-zero means the
     *  current frame is incomplete, so hitting EOF there means a truncated file. */
    size_t          last_ret;

    bool            eof;
    char            last_char;

    void reset ()
    {
        fseek (file, 0, SEEK_SET);
        ZSTD_initDStream (dstream);
        input.src    = in_buffer;
        input.size   = 0;
        input.pos    = 0;
        buffer_start = 0;
        buffer_end   = 0;
        out_total    = 0;
        last_ret     = 0;
        eof          = false;
        last_char    = 0;
    }

    /** Offset of the parser within the decompressed stream (zlib's gztell). */
    uint64_t tell () const  { return out_total - (buffer_end - buffer_start); }
};

bool has_zstd_magic (const unsigned char* bytes)
{
    uint32_t magic = (uint32_t)bytes[0]        | ((uint32_t)bytes[1] << 8)
                   | ((uint32_t)bytes[2] << 16) | ((uint32_t)bytes[3] << 24);

    if (magic == ZSTD_MAGICNUMBER)  { return true; }

    /** Skippable frames use 0x184D2A50..0x184D2A5F. */
    return (magic & ZSTD_MAGIC_SKIPPABLE_MASK) == ZSTD_MAGIC_SKIPPABLE_START;
}

zstd_file_t* zstd_open (const char* path)
{
    FILE* file = fopen (path, "rb");
    if (file == NULL)  { return NULL; }

    unsigned char magic[4];
    if (fread (magic, 1, sizeof(magic), file) != sizeof(magic) || !has_zstd_magic (magic))
    {
        fclose (file);
        return NULL;
    }

    zstd_file_t* zf = (zstd_file_t*) CALLOC (1, sizeof(zstd_file_t));

    zf->file        = file;
    zf->dstream     = ZSTD_createDStream ();
    zf->in_capacity = ZSTD_DStreamInSize ();
    zf->in_buffer   = (unsigned char*) MALLOC (zf->in_capacity);
    zf->buffer      = (unsigned char*) MALLOC (BUFFER_SIZE);

    if (zf->dstream == NULL || zf->in_buffer == NULL || zf->buffer == NULL)
    {
        if (zf->dstream != NULL)  { ZSTD_freeDStream (zf->dstream); }
        if (zf->in_buffer != NULL)  { FREE (zf->in_buffer); }
        if (zf->buffer != NULL)  { FREE (zf->buffer); }
        FREE (zf);
        fclose (file);
        return NULL;
    }

    zf->reset ();

    return zf;
}

void zstd_close (zstd_file_t* zf)
{
    if (zf == NULL)  { return; }

    if (zf->dstream != NULL)  { ZSTD_freeDStream (zf->dstream); }
    if (zf->file    != NULL)  { fclose (zf->file);              }

    FREE (zf->in_buffer);
    FREE (zf->buffer);
    FREE (zf);
}

bool rebuffer (zstd_file_t* zf)
{
    if (zf->eof)  { return false; }

    zf->buffer_start = 0;
    zf->buffer_end   = 0;

    ZSTD_outBuffer output = { zf->buffer, BUFFER_SIZE, 0 };

    while (output.pos == 0)
    {
        if (zf->input.pos >= zf->input.size)
        {
            size_t nread = fread (zf->in_buffer, 1, zf->in_capacity, zf->file);

            if (nread == 0)
            {
                if (ferror (zf->file))
                {
                    throw gatb::core::system::Exception ("BankZSTD: read error (%s)", strerror(errno));
                }

                /** last_ret != 0 means we stopped in the middle of a frame. */
                if (zf->last_ret != 0)
                {
                    throw gatb::core::system::Exception ("BankZSTD: truncated zstd stream");
                }

                zf->eof = true;
                break;
            }

            zf->input.src  = zf->in_buffer;
            zf->input.size = nread;
            zf->input.pos  = 0;
        }

        size_t ret = ZSTD_decompressStream (zf->dstream, &output, &zf->input);

        if (ZSTD_isError (ret))
        {
            throw gatb::core::system::Exception ("BankZSTD: %s", ZSTD_getErrorName (ret));
        }

        zf->last_ret = ret;
    }

    zf->buffer_end  = output.pos;
    zf->out_total  += output.pos;

    if (output.pos == 0)  { zf->eof = true;  return false; }

    return true;
}

/********************************************************************************/
struct variable_string_t
{
     variable_string_t ()  : length(0), max(0), string(0) {}
    ~variable_string_t () { if (string!=0)  { FREE(string); } }

    uint64_t length, max;
    char *string;
};

/********************************************************************************/
struct buffered_strings_t
{
     buffered_strings_t () : read(new variable_string_t), dummy(new variable_string_t), header(new variable_string_t), quality(new variable_string_t)   {}
    ~buffered_strings_t ()
    {
        delete read;
        delete dummy;
        delete header;
        delete quality;
    }

    variable_string_t *read, *dummy, *header, *quality;
};

inline signed char buffered_getc (zstd_file_t* zf)
{
    if (zf->buffer_start >= zf->buffer_end) if (!rebuffer (zf)) return -1;
    return (signed char) (zf->buffer[zf->buffer_start++]);
}

/********************************************************************************/
inline signed int buffered_gets (
    zstd_file_t*       zf,
    variable_string_t* s,
    char *dret,
    bool append,
    bool allow_spaces
)
{
    if (dret) *dret = 0;
    if (!append) s->length = 0;
    if (zf->buffer_start >= zf->buffer_end && zf->eof) return -1;
    while (1)
    {
        uint64_t i;
        if (zf->buffer_start >= zf->buffer_end) if (!rebuffer (zf)) break;
        if (allow_spaces)
        {
            for (i = zf->buffer_start; i < zf->buffer_end; i++)
                if (zf->buffer[i] == '\n') break;
        }
        else
        {
            for (i = zf->buffer_start; i < zf->buffer_end; i++)
                // isspace() answers yes for ' ', \t, \n, \v, \f, \r
                if (isspace (zf->buffer[i])) break;
        }
        if (s->max - s->length < (i - zf->buffer_start + 1))
        {
            s->max = s->length + (i - zf->buffer_start + 1);
            if (is_power_of_2(s->max))  { s->max ++; }
            nearest_power_of_2(s->max);
            s->string = (char*)  REALLOC (s->string, s->max);
        }
        memcpy (s->string + s->length, zf->buffer + zf->buffer_start, i - zf->buffer_start);
        s->length += i - zf->buffer_start;
        zf->buffer_start = i + 1;
        if (i < zf->buffer_end)
        {
            if (dret) *dret = zf->buffer[i];
            break;
        }
    }
    if (s->string == NULL)
    {
        s->max = 256;
        s->string = (char*)  CALLOC (256, 1);
    }
    else if (allow_spaces && s->length > 1 && s->string[s->length - 1] == '\r')
        s->length--;
    s->string[s->length] = '\0';
    return s->length;
}

/********************************************************************************/
/** Same heuristic as BankFastaFactory: a '>' or '@' at the start of a line, then
 *  printable characters up to a blank. Applied to the *decompressed* prefix. */
bool looks_like_fasta (const char* buffer, int size)
{
    int i = 0;
    char previous = ' ';
    bool foundStart = false;

    for (i = 0; !foundStart && i < size; i++)
    {
        if ((buffer[i]=='>' || buffer[i]=='@')  && isspace(previous))  { foundStart = true; break; }
        previous = buffer[i];
    }

    if (!foundStart)  { return false; }

    /** We look for alphanum + space. */
    bool foundSpace = false;
    for ( ; !foundSpace && i < size; i++)
    {
             if (  isspace(buffer[i]))  { foundSpace = true; }
        else if (! isprint(buffer[i]))  { break; }
    }

    return foundSpace;
}

} /* end of anonymous namespace */

/*********************************************************************
** METHOD  : BankZSTD::BankZSTD
*********************************************************************/
BankZSTD::BankZSTD (const std::string& filename)
    : _filename(filename), _filesize(0)
{
    init ();
}

/*********************************************************************
** METHOD  : BankZSTD::~BankZSTD
*********************************************************************/
BankZSTD::~BankZSTD ()
{
}

/*********************************************************************
** METHOD  : BankZSTD::init
*********************************************************************/
void BankZSTD::init ()
{
    if (!System::file().doesExist (_filename))
    {
        throw gatb::core::system::Exception ("zstd file does not exist: %s", _filename.c_str());
    }

    _filesize = System::file().getSize (_filename) * 4;
}

/*********************************************************************
** METHOD  : BankZSTD::estimate
*********************************************************************/
void BankZSTD::estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize)
{
    BankZSTD::Iterator it (*this, Iterator::NONE);

    it.estimate (number, totalSize, maxSize);
}

/*********************************************************************
** METHOD  : BankZSTD::insert
*********************************************************************/
void BankZSTD::insert (const Sequence& item)
{
    throw gatb::core::system::Exception ("BankZSTD::insert not implemented (read-only)");
}

/*********************************************************************
** METHOD  : BankZSTD::flush
*********************************************************************/
void BankZSTD::flush ()
{
    // Not applicable for read-only
}

/********************************************************************************/
// BankZSTD::Iterator implementation
/********************************************************************************/

/*********************************************************************
** METHOD  : BankZSTD::Iterator::Iterator
*********************************************************************/
BankZSTD::Iterator::Iterator (BankZSTD& ref, CommentMode_e commentMode)
    : _ref(ref), _commentsMode(commentMode), _isDone(true), _isInitialized(false), _nIters(0),
      _stream(0), _strings(0), _index(0)
{
    DEBUG (("BankZSTD::Iterator::Iterator\n"));

    /** We check that the file can be opened as a zstd stream. */
    if (zstd_file_t* zf = zstd_open (_ref._filename.c_str()))  { zstd_close (zf); }
    else  { throw gatb::core::system::ExceptionErrno (STR_BANK_unable_open_file, _ref._filename.c_str()); }
}

/*********************************************************************
** METHOD  : BankZSTD::Iterator::~Iterator
*********************************************************************/
BankZSTD::Iterator::~Iterator ()
{
    DEBUG (("BankZSTD::Iterator::~Iterator\n"));
    finalize ();
}

/*********************************************************************
** METHOD  : BankZSTD::Iterator::init
*********************************************************************/
void BankZSTD::Iterator::init ()
{
    if (_isInitialized == true)  { return; }

    zstd_file_t* zf = zstd_open (_ref._filename.c_str());

    if (zf == NULL)
    {
        fprintf (stderr, "unable to open file %s  : %s \n", _ref._filename.c_str(), strerror(errno));

        throw gatb::core::system::ExceptionErrno (STR_BANK_unable_open_file, _ref._filename.c_str());
    }

    _stream  = zf;
    _strings = new buffered_strings_t;

    _isInitialized = true;
}

/*********************************************************************
** METHOD  : BankZSTD::Iterator::finalize
*********************************************************************/
void BankZSTD::Iterator::finalize ()
{
    if (_isInitialized == false)  { return; }

    buffered_strings_t* bs = (buffered_strings_t*) _strings;
    if (bs != 0)  { delete bs;  _strings = 0; }

    zstd_close ((zstd_file_t*) _stream);
    _stream = 0;

    _isInitialized = false;
}

/*********************************************************************
** METHOD  : BankZSTD::Iterator::first
*********************************************************************/
void BankZSTD::Iterator::first()
{
    /** We may have to initialize the instance. */
    init ();

    ((zstd_file_t*) _stream)->reset ();

    _isDone = false;

    _nIters = 0;
    _index  = 0;

    next();
}

/*********************************************************************
** METHOD  : BankZSTD::Iterator::next
*********************************************************************/
void BankZSTD::Iterator::next()
{
    if (_isDone)  { return; }

    if (_commentsMode == NONE)
    {
        _isDone = get_next_seq (_item->getData()) == false;
    }
    else
    {
        _isDone = get_next_seq (_item->getData(), _item->_comment, _item->_quality, _commentsMode) == false;
    }
    _item->setIndex (_index++);

    DEBUG (("BankZSTD::Iterator::next  _isDone=%d\n", _isDone));
}

/*********************************************************************
** METHOD  : BankZSTD::Iterator::get_next_seq
*********************************************************************/
bool BankZSTD::Iterator::get_next_seq (Vector<char>& data)
{
    string dummy;
    return get_next_seq (data, dummy, dummy, NONE);
}

/*********************************************************************
** METHOD  : BankZSTD::Iterator::get_next_seq
*********************************************************************/
bool BankZSTD::Iterator::get_next_seq (Vector<char>& data, string& comment, string& quality, CommentMode_e mode)
{
    buffered_strings_t* bs = (buffered_strings_t*) _strings;

    signed char c;
    zstd_file_t* zf = (zstd_file_t*) _stream;

    if (zf->last_char == 0)
    {
        while ((c = buffered_getc (zf)) != -1 && c != '>' && c != '@')
            ; // go to next header
        if (c == -1) return false; // eof
        zf->last_char = c;
    }
    bs->quality->length = bs->read->length = bs->dummy->length = 0;

    if (buffered_gets (zf, bs->header, (char *) &c, false, false) < 0)
        return false; // eof

    if (c != '\n')
    {
        if (mode == IDONLY)
        {
            variable_string_t dummy;
            buffered_gets (zf, &dummy, NULL, true, true); // read header
        }
        else
        {
            // We add the last read character (likely a space)
            bs->header->string[bs->header->length] = c;
            bs->header->length++;
            bs->header->string[bs->header->length] = 0;

            buffered_gets (zf, bs->header, NULL, true, true); // read header
        }
    }

    if (bs->read->string == NULL)
    {
        bs->read->max = 256;
        bs->read->string = (char*)  MALLOC (bs->read->max);
    }
    while ((c = buffered_getc (zf)) != -1 && c != '>' && c != '+' && c != '@')
    {
        if (c == '\n') continue; // empty line
        bs->read->string[bs->read->length++] = c;
        buffered_gets (zf, bs->read, NULL, true, true);
    }
    if (c == '>' || c == '@') zf->last_char = c;
    if (bs->read->length + 1 >= bs->read->max)
    {
        bs->read->max = bs->read->length + 2;
        nearest_power_of_2(bs->read->max);
        bs->read->string = (char*)  REALLOC (bs->read->string, bs->read->max);
    }
    bs->read->string[bs->read->length] = '\0';
    if (c == '+') // fastq
    {
        if (bs->quality->max < bs->read->max) // resize quality to match read length
        {
            bs->quality->max = bs->read->max;
            bs->quality->string = (char*)  REALLOC (bs->quality->string, bs->quality->max);
        }
        while ((c = buffered_getc (zf)) != -1 && c != '\n')
            ; // read rest of quality comment
        while (buffered_gets (zf, bs->quality, NULL, true, true) >= 0 && bs->quality->length < bs->read->length)
            ; // read rest of quality
        zf->last_char = 0;

        quality.assign (bs->quality->string, bs->quality->length);
    }

    /** We update the data of the sequence. */
    data.set (bs->read->string, bs->read->length);

    comment.assign (bs->header->string, bs->header->length);

    return true;
}

/*********************************************************************
** METHOD  : BankZSTD::Iterator::estimate
*********************************************************************/
void BankZSTD::Iterator::estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize)
{
    /** We may have to initialize the instance. */
    init ();

    Vector<char> data;

    /** We rewind the stream. */
    ((zstd_file_t*) _stream)->reset ();

    /** We initialize the provided arguments. */
    number    = 0;
    totalSize = 0;
    maxSize   = 0;

    while (get_next_seq (data)  &&  number <= _ref.getEstimateThreshold())
    {
        number ++;
        if (data.size() > maxSize)  { maxSize = data.size(); }
        totalSize += data.size ();
    }

    u_int64_t actualPosition = ((zstd_file_t*) _stream)->tell ();

    if (actualPosition > 0)
    {
        // linear extrapolation
        number    = (number    * _ref.getSize()) / actualPosition;
        totalSize = totalSize * ( (float)_ref.getSize() / (float)actualPosition) ;
    }
}

/********************************************************************************/
// BankZSTDFactory implementation
/********************************************************************************/

IBank* BankZSTDFactory::createBank (const std::string& uri)
{
    zstd_file_t* zf = zstd_open (uri.c_str());
    if (zf == NULL)  { return NULL; }

    bool isFasta = false;

    try
    {
        if (rebuffer (zf))
        {
            int size = (int) min ((uint64_t) 256, zf->buffer_end);
            isFasta = looks_like_fasta ((const char*) zf->buffer, size);
        }
    }
    catch (...)
    {
        zstd_close (zf);
        return NULL;
    }

    zstd_close (zf);

    if (!isFasta)  { return NULL; }

    try
    {
        return new BankZSTD (uri);
    }
    catch (...)
    {
        return NULL;
    }
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
