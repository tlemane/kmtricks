/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
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

#include <gatb/bank/impl/BankFasta.hpp>
#include <gatb/bank/impl/BankComposite.hpp>

#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/api/StringsRepository.hpp>
#include <gatb/tools/misc/impl/Tokenizer.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>

#include <algorithm>
#include <string.h>
#include <errno.h>
#include <zlib.h> // Added by Pierre Peterlongo on 02/08/2012.

using namespace std;
using namespace gatb::core::tools::dp;
using namespace gatb::core::tools::dp::impl;
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

size_t BankFasta::_dataLineSize = 70;

/********************************************************************************/
// heavily inspired by kseq.h from Heng Li (https://github.com/attractivechaos/klib)
typedef struct
{
    gzFile stream;
    unsigned char *buffer;
    uint64_t buffer_start, buffer_end;
    bool eof;
    char last_char;

    void rewind ()
    {
        gzrewind (stream);
        last_char    = 0;
        eof          = 0;
        buffer_start = 0;
        buffer_end   = 0;
    }

} buffered_file_t;

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

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankFasta::BankFasta (const std::string& filename, bool output_fastq, bool output_gz)
    : filesizes(0), nb_files(0), _insertHandle(0), _gz_insertHandle(0)
{
    _output_fastq = output_fastq;
    _output_gz= output_gz;
    _filenames.push_back (filename);
    init ();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankFasta::~BankFasta ()
{
    if (_insertHandle    != 0)  { fclose  (_insertHandle);    }
    if (_gz_insertHandle != 0)  { gzclose (_gz_insertHandle); }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankFasta::finalize ()
{
    if (_insertHandle    != 0)  { fclose  (_insertHandle);    _insertHandle    = 0; }
    if (_gz_insertHandle != 0)  { gzclose (_gz_insertHandle); _gz_insertHandle = 0; }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankFasta::init ()
{
    /** We check that we don't exceed the number of allowed files. */
    if (_filenames.empty() || _filenames.size() > getMaxNbFiles())
    {
        /** We send an exception. */
	    fprintf(stderr,"unable to open file(s)"); // normally you get an information exception, but not always in some test programs (e.g. GATB training day), so i'm printing an explicit message
        throw gatb::core::system::Exception (STR_BANK_bad_file_number, _filenames.size(), getMaxNbFiles());
    }

    nb_files  = _filenames.size();
    filesizes = 0;

    // estimate total size of files
    for (size_t i=0; i<nb_files; i++)
    {
        /** Shortcut. */
        const char* fname = _filenames[i].c_str();

        bool compressed = false;
        u_int64_t estimated_filesize;

        if (strstr (fname, "gz") == (fname + strlen (fname) - 2))
            compressed = true;

        if (compressed)
            // crude hack, based on Quip paper reporting compression ratio (~0.3).
            // gzseek(SEEK_END) isn't supported. need to read whole file otherwise :/

            estimated_filesize = System::file().getSize (fname) * 4;
        else
            estimated_filesize = System::file().getSize (fname);

        filesizes += estimated_filesize;
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankFasta::estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize)
{
    /** We create an iterator for the bank. */
    BankFasta::Iterator it (*this, Iterator::NONE);

    /** We delegate the computation to the iterator. */
    it.estimate (number, totalSize, maxSize);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankFasta::flush ()
{
    if (_insertHandle    != 0)  { fflush  (_insertHandle);              }
    if (_gz_insertHandle != 0)  { gzflush (_gz_insertHandle, Z_FINISH); }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankFasta::insert (const Sequence& item)
{
    /** We open the last file if needed. */
    if (_insertHandle == 0  &&  _filenames.empty()==false)
    {
        _insertHandle = fopen (_filenames[_filenames.size()-1].c_str(), "w");
    }

    if (_output_gz && _gz_insertHandle == 0  &&  _filenames.empty()==false)
    {
        _gz_insertHandle =  gzopen (_filenames[_filenames.size()-1].c_str(), "w");
    }
    
    if (_insertHandle != 0)
    {
        if(_output_fastq)
        {
            if(_output_gz)
            {
                gzprintf (_gz_insertHandle, "@%s\n", item.getComment().c_str());
                gzprintf (_gz_insertHandle, "%.*s\n",(int)item.getDataSize(),  item.getDataBuffer());
                gzprintf (_gz_insertHandle, "+\n");
                gzprintf (_gz_insertHandle, "%s\n", item.getQuality().c_str());
            }
            else
            {
                fprintf (_insertHandle, "@%s\n", item.getComment().c_str());
                fprintf (_insertHandle, "%.*s\n",(int)item.getDataSize(),  item.getDataBuffer());
                fprintf (_insertHandle, "+\n");
                fprintf (_insertHandle, "%s\n", item.getQuality().c_str());
            }

        }
        else
        {
        /** We add the sequence into the bag. */
        fprintf (_insertHandle, ">%s\n", item.getComment().c_str());

#if 1
        fprintf (_insertHandle, "%.*s\n",(int)item.getDataSize(),  item.getDataBuffer());
#else
        // We dump the data with fixed sized columns
        char line[4*1024];
        char* loop = line;
        size_t actualDataLineSize = _dataLineSize > 0 ? _dataLineSize : (size_t)(~0);

        size_t len = item.getDataSize();

        for (size_t i=0; i<len; )
        {
            size_t j=0;
            for (j=0; j<actualDataLineSize && i<len; j++, i++)
            {
                *(loop++) = item.getDataBuffer() [i];
                if (j >= sizeof(line)-2)
                {
                    *(loop++) = 0;   fprintf (_insertHandle, "%s", line);     loop = line;
                }
            }
            *(loop++) = 0;    fprintf (_insertHandle, "%s\n", line);    loop = line;
        }
#endif
        }
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankFasta::Iterator::Iterator (BankFasta& ref, CommentMode_e commentMode)
    : _ref(ref), _commentsMode(commentMode), _isDone(true), _isInitialized(false), _nIters(0),
      index_file(0), buffered_file(0), buffered_strings(0), _index(0)
{
    DEBUG (("Bank::Iterator::Iterator\n"));

    /** We check that the file can be opened. */
    if (gzFile stream = gzopen (_ref._filenames[0].c_str(), "r"))  {  gzclose (stream);  }
    else  {  
        throw gatb::core::system::ExceptionErrno (STR_BANK_unable_open_file, _ref._filenames[0].c_str());  }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
BankFasta::Iterator::~Iterator ()
{
    DEBUG (("Bank::Iterator::~Iterator\n"));
    finalize ();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankFasta::Iterator::first()
{
    /** We may have to initialize the instance. */
    init  ();

    for (u_int64_t i = 0; i < _ref.nb_files; i++)
    {
        buffered_file_t* bf = (buffered_file_t *) buffered_file[i];
        if (bf != 0)  { bf->rewind(); }
    }

    index_file = 0;
    _isDone = false;
    
    _nIters = 0;
    _index  = 0;
    
    next();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankFasta::Iterator::next()
{
    if (_isDone)  { return; }

    if (_commentsMode == NONE)
    {
        _isDone = get_next_seq (_item->getData()) == false;
    }
    else
    {
        _isDone = get_next_seq (_item->getData(), _item->_comment,_item->_quality, _commentsMode) == false;
    }
    _item->setIndex (_index++);
    DEBUG (("Bank::Iterator::next  _isDone=%d\n", _isDone));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
// the following functions are adapted from kseq.h by Heng Li (https://github.com/attractivechaos/klib)
inline bool rebuffer (buffered_file_t *bf)
{
    if (bf->eof) return false;
    bf->buffer_start = 0;
    bf->buffer_end = gzread (bf->stream, bf->buffer, BUFFER_SIZE);
    if (bf->buffer_end < BUFFER_SIZE) bf->eof = 1;
    if (bf->buffer_end == 0) return false;
    return true;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
inline signed char buffered_getc (buffered_file_t *bf)
{
    if (bf->buffer_start >= bf->buffer_end) if (!rebuffer (bf)) return -1;
    return (signed char) (bf->buffer[bf->buffer_start++]);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
inline signed int buffered_gets (
    buffered_file_t*   bf,
    variable_string_t* s,
    char *dret,
    bool append,
    bool allow_spaces
)
{
    if (dret) *dret = 0;
    if (!append) s->length = 0;
    if (bf->buffer_start >= bf->buffer_end && bf->eof) return -1;
    while (1)
    {
        uint64_t i;
        if (bf->buffer_start >= bf->buffer_end) if (!rebuffer (bf)) break;
        if (allow_spaces)
        {
            for (i = bf->buffer_start; i < bf->buffer_end; i++)
                if (bf->buffer[i] == '\n') break;
        }
        else
        {
            for (i = bf->buffer_start; i < bf->buffer_end; i++)
                // isspace() answers yes for ' ', \t, \n, \v, \f, \r
                if (isspace (bf->buffer[i])) break;
        }
        if (s->max - s->length < (i - bf->buffer_start + 1))
        {
            s->max = s->length + (i - bf->buffer_start + 1);
            if (is_power_of_2(s->max))  { s->max ++; }
            nearest_power_of_2(s->max);
            s->string = (char*)  REALLOC (s->string, s->max);
            //std::cout << "realloc of size " << s->max << " res: " << (uint64_t)(s->string) << std::endl;
        }
         memcpy (s->string + s->length, bf->buffer + bf->buffer_start, i - bf->buffer_start);
        s->length += i - bf->buffer_start;
        bf->buffer_start = i + 1;
        if (i < bf->buffer_end)
        {
            if (dret) *dret = bf->buffer[i];
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

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool BankFasta::Iterator::get_next_seq_from_file (Vector<char>& data, string& comment, string& quality, int file_id, CommentMode_e mode)
{
    
    buffered_strings_t* bs = (buffered_strings_t*) buffered_strings;
   // printf("%i -\n",bs->header->length);

    signed char c;
    buffered_file_t *bf = (buffered_file_t *) buffered_file[file_id];
    if (bf->last_char == 0)
    {
        while ((c = buffered_getc (bf)) != -1 && c != '>' && c != '@')
            ; // go to next header
        if (c == -1) return false; // eof
        bf->last_char = c;
    }
    bs->quality->length = bs->read->length = bs->dummy->length = 0;

    if (buffered_gets (bf, bs->header, (char *) &c, false, false) < 0) //ici
        return false; // eof

    if (c != '\n')
    {
        if (mode == IDONLY)
        {
            variable_string_t dummy;
            buffered_gets (bf, &dummy, NULL, true, true); // read header
        }
        else
        {
            // We add the last read character (likely a space)
            bs->header->string[bs->header->length] = c;
            bs->header->length++;
            bs->header->string[bs->header->length] = 0;

            buffered_gets (bf, bs->header, NULL, true, true); // read header
        }
	}

    if (bs->read->string == NULL)
    {
        bs->read->max = 256;
        bs->read->string = (char*)  MALLOC (bs->read->max);
    }
    while ((c = buffered_getc (bf)) != -1 && c != '>' && c != '+' && c != '@')
    {
        if (c == '\n') continue; // empty line
        bs->read->string[bs->read->length++] = c;
        buffered_gets (bf, bs->read, NULL, true, true);
    }
    if (c == '>' || c == '@') bf->last_char = c;
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
        while ((c = buffered_getc (bf)) != -1 && c != '\n')
            ; // read rest of quality comment
        while (buffered_gets (bf, bs->quality, NULL, true, true) >= 0 && bs->quality->length < bs->read->length)
            ; // read rest of quality
        bf->last_char = 0;
        
        quality.assign (bs->quality->string, bs->quality->length);
       // printf("%i  %i\n",bs->quality->length,bs->header->length);
    }

    /** We update the data of the sequence. */
#if 1
    data.set (bs->read->string, bs->read->length);
#else
    data.setRef (bs->read->string, bs->read->length);
#endif

    //if (comment.empty() == false)
    {
        comment.assign (bs->header->string, bs->header->length);
    }

    return true;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool BankFasta::Iterator::get_next_seq_from_file (Vector<char>& data, int file_id)
{
    string dummy;
    return get_next_seq_from_file (data, dummy,dummy, file_id, NONE);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool BankFasta::Iterator::get_next_seq (Vector<char>& data, string& comment,string& quality, CommentMode_e mode)

{
    bool success = get_next_seq_from_file (data, comment,quality, index_file, mode);
    if (success) return true;

    // cycle to next file if possible
    if ((u_int64_t)index_file < _ref.nb_files - 1)
    {
        index_file++;
        return get_next_seq (data, comment,quality, mode);
    }
    return false;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool BankFasta::Iterator::get_next_seq (Vector<char>& data)
{
    string  dummy;
    return get_next_seq (data, dummy,dummy, NONE);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankFasta::Iterator::init ()
{

    if (_isInitialized == true)  { return ;}

    /** We initialize the array of files. */
    buffered_file = (void**) CALLOC (getMaxNbFiles(), sizeof(void*));

    /** Shortcut. */
    vector<string>& fnames = _ref._filenames;

    // open each file for reading
    for (size_t i=0; i<_ref.nb_files; i++)
    {
        /** Shortcut. */
        const char* fname = fnames[i].c_str();

        buffered_file_t** bf = (buffered_file_t **) buffered_file + i;
        *bf = (buffered_file_t *)  CALLOC (1, sizeof(buffered_file_t));
        (*bf)->buffer = (unsigned char*)  MALLOC (BUFFER_SIZE);
        (*bf)->stream = gzopen (fname, "r");
        gzbuffer((*bf)->stream,2*1024*1024);
		
        /** We check that we can open the file. */
        if ((*bf)->stream == NULL)
        {
            // there used to be some cleanup here but what's the point, we're going to throw an exception anyway
        
			//GR : dunno why this exception does not show up, adding a message here
            //RC: the exceptino doesn't even trigger, or there is an exception but the message doesn't show?
			fprintf(stderr,"unable to open file %s  : %s \n",fname,strerror(errno));
			
            /** We launch an exception. */
            throw gatb::core::system::ExceptionErrno (STR_BANK_unable_open_file, fname);

        }
    }

    index_file = 0; // initialize the get_next_seq iterator to the first file

    // init read and dummy (for readname and quality)
    buffered_strings_t* bs = new buffered_strings_t;
    buffered_strings = bs;

    _isInitialized = true;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankFasta::Iterator::finalize ()
{
    if (_isInitialized == false)  { return; }

    buffered_strings_t* bs = (buffered_strings_t*) buffered_strings;

    if (bs != 0)  { delete bs; }

    for (u_int64_t i = 0; i < _ref.nb_files; i++)
    {
        buffered_file_t* bf = (buffered_file_t *) buffered_file[i];

        if (bf != 0)
        {
            /** We close the handle of the file. */
            if (bf->stream != NULL)  {  gzclose (bf->stream);  bf->stream = 0; }

            /** We delete the buffer. */
            FREE (bf->buffer);

            /** We delete the buffered file itself. */
            FREE (bf);
        }
    }

    /** We release the array of files. */
    FREE (buffered_file);

    /** We reset the initialization flag. */
    _isInitialized = false;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void BankFasta::Iterator::estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize)
{
    /** We may have to initialize the instance. */
    init  ();

    Vector<char> data;

    /** We rewind the files. */
    for (u_int64_t i = 0; i < _ref.nb_files; i++)
    {
        buffered_file_t* bf = (buffered_file_t *) buffered_file[i];
        if (bf != 0)  { bf->rewind(); }
    }

    /** We initialize the provided arguments. */
    number    = 0;
    totalSize = 0;
    maxSize   = 0;

    number = 0;
    while (get_next_seq (data)  &&  number <= _ref.getEstimateThreshold())
    {
        number ++;
        if (data.size() > maxSize)  { maxSize = data.size(); }
        totalSize += data.size ();
    }

    u_int64_t actualPosition = 0;

    /** We compute the aggregated size from the files having been read until we
     * reached our limit number of sequences. */
    for (int i=0; i<=index_file; i++)
    {
        buffered_file_t* current = (buffered_file_t *) buffered_file[i];

        actualPosition += gztell (current->stream);
    }

    if (actualPosition > 0)
    {
        // linear extrapolation
        number    = (number    * _ref.getSize()) / actualPosition;
        totalSize = totalSize * ( (float)_ref.getSize() / (float)actualPosition) ;
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IBank* BankFastaFactory::createBank (const std::string& uri)
{
    bool isFASTA = false;

    /** We check whether the uri looks like a FASTA bank. */
    gzFile file = gzopen (uri.c_str(), "r");
    if (file != 0)
    {
        char buffer[256];
        int res = gzread (file, buffer, sizeof(buffer));
        if (res > 0)
        {
            int i=0;
            char previous=' ';
            bool foundStart = false;
            for (i=0; !foundStart && i<res; i++)
            {
                if ((buffer[i]=='>' || buffer[i]=='@')  && isspace(previous))    {  foundStart = true;  break; }
                previous=buffer[i];
            }
            if (foundStart)
            {
                /** We look for alphanum + space. */
                bool foundSpace = false;
                for ( ; !foundSpace && i<res; i++)
                {
                         if (  isspace(buffer[i]))  { foundSpace=true; }
                    else if (! isprint(buffer[i]))  { break; }
                }
                isFASTA = foundSpace;
            }
        }

        gzclose (file);
    }

    return (isFASTA ? new BankFasta (uri) : NULL);
}

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
