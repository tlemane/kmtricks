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

/** \file IteratorFile.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Iterator implementation for file
 */

#ifndef _GATB_CORE_TOOLS_COLLECTIONS_IMPL_ITERATOR_FILE_HPP_
#define _GATB_CORE_TOOLS_COLLECTIONS_IMPL_ITERATOR_FILE_HPP_

/********************************************************************************/

#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <gatb/tools/collections/api/Iterable.hpp>
#include <gatb/system/impl/System.hpp>

#include <string>
#include <vector>
#include <zlib.h>

/********************************************************************************/
namespace gatb          {
namespace core          {
namespace tools         {
namespace collections   {
/** Implementation of API from collections */
namespace impl          {
/********************************************************************************/

#define DEBUG_ITERATORFILE(x) {}//{x}

#define BUFFER_SIZE (128*1024)

/** \brief Iterator implementation for file
 */
template <class Item> class IteratorFile : public dp::Iterator<Item>
{
public:

    /** Constructor. */
    IteratorFile () : _file(0), _buffer(0), _cpt_buffer(0), _idx(0), _cacheItemsNb(0), _isDone(true) {}

    IteratorFile (const IteratorFile& it):
        _filename(it._filename), _file(0),  _buffer(0), _cpt_buffer(0), _idx(0), _cacheItemsNb(it._cacheItemsNb), _isDone(true)
    {
        _file    = system::impl::System::file().newFile (_filename, "rb");
        _buffer  = (Item*) MALLOC (sizeof(Item) * _cacheItemsNb);
        DEBUG_ITERATORFILE(std::cout << "iteratorfile constructed" << std::endl;)
    }

    /** Constructor. */
    IteratorFile (const std::string& filename, size_t cacheItemsNb=10000) :
        _filename(filename), _file(0),  _buffer(0), _cpt_buffer(0), _idx(0), _cacheItemsNb(cacheItemsNb), _isDone(true)

    {
        _file    = system::impl::System::file().newFile (filename, "rb");
        _buffer  = (Item*) MALLOC (sizeof(Item) * _cacheItemsNb);
    }

    /** Destructor. */
    ~IteratorFile ()
    {
        if (_file)  { delete _file;  }
        if (_buffer) { FREE (_buffer); }
    }

    /** Affectation. */
    IteratorFile& operator= (const IteratorFile& it)
    {
        if (this != &it)
        {
            if (_file)    { delete _file; }
            if (_buffer)  { FREE(_buffer); }

            _filename     = it._filename;
            _cpt_buffer   = it._cpt_buffer;
            _idx          = it._idx;
            _cacheItemsNb = it._cacheItemsNb;
            _isDone       = it._isDone;

            _file    = system::impl::System::file().newFile (it._filename, "rb");
            _buffer  = (Item*) MALLOC (sizeof(Item) * it._cacheItemsNb);
        }
        return *this;
    }

    /** \copydoc dp::Iterator::first */
    void first()
    {
        _file->seeko (0, SEEK_SET);
        _cpt_buffer = 0;
        _idx        = 0;
        _isDone     = false;
        next ();
    }

    /** \copydoc dp::Iterator::next */
    void next()
    {
        if (_cpt_buffer==0)
        {
            _idx = 0;
            DEBUG_ITERATORFILE(std::cout << "(next) doing a fread of " << _cacheItemsNb << " items at position " << _file->tell() << " file size " << _file->getSize() << std::endl;)
            _cpt_buffer = _file->fread (_buffer, sizeof(Item), _cacheItemsNb);
            if (_cpt_buffer==0)  { _isDone = true;  return; }
        }

        *(this->_item) =  _buffer[_idx];
        _cpt_buffer --;
        _idx ++;
    }

    /** \copydoc dp::Iterator::isDone */
    bool isDone()  { return _isDone; }

    /** \copydoc dp::Iterator::item */
    Item& item ()  { return *(this->_item); }

    /** */
    size_t fill (std::vector<Item>& vec, size_t len=0)
    {
        if (len==0)  { len = vec.size(); }
        DEBUG_ITERATORFILE(std::cout << "(fill) doing a fread of " << len << " items at position " << _file->tell() << " file size " << _file->getSize() << std::endl;)
        return _file->fread (vec.data(), sizeof(Item), len);
    }

private:
    std::string     _filename;
    system::IFile*  _file;
    Item*           _buffer;
    int             _cpt_buffer;
    int             _idx;
    size_t          _cacheItemsNb;
    bool            _isDone;
};

/********************************************************************************/

/** \brief Implementation of the Iterable interface as a file
 *
 * This implementation uses a file as the source of the items that can be iterated.
 *
 */
template <class Item> class IterableFile : public tools::collections::Iterable<Item>, public virtual system::SmartPointer
{
public:

    /** Constructor
     * \param[in] filename : name of the file to be iterated.
     * \param[in] cacheItemsNb : number of items in the cache memory
     */
    IterableFile (const std::string& filename, size_t cacheItemsNb=10000)
        :   _filename(filename), _cacheItemsNb (cacheItemsNb), 
        _file(0)  // hacking my own iterator, for getItems, separate from IteratorFile. dirty, but nothing used to work at all. _file is used in getItems() only
    {
        // if the file doesn't exist (meaning that BagFile hasn't created it yet), let's create it just for the sake of it. but then we'll open it just for reading
        if (!system::impl::System::file().doesExist(filename))
            {
                auto   _file2 = system::impl::System::file().newFile (filename, "wb");
                delete _file2;
            }
        /* _file should be initialized here but actually, the iterator() method will also create its own file.
         * so, instead of opening _file here, let's wait until getItems() is actually called (sometimes it won't).
         */
    }

    /** Destructor. */
    ~IterableFile () {
        if (_file)  { delete _file;  }
    }

    /** \copydoc Iterable::iterator */
    dp::Iterator<Item>* iterator ()  { return new IteratorFile<Item> (_filename, _cacheItemsNb); }

    /** \copydoc Iterable::getNbItems */
    int64_t getNbItems ()   {  
        DEBUG_ITERATORFILE(std::cout << "IteratorFile::getNbItems called (file size: "<< system::impl::System::file().getSize(_filename) << "), returning " << system::impl::System::file().getSize(_filename) / sizeof(Item) << std::endl;)
        return system::impl::System::file().getSize(_filename) / sizeof(Item); }

    /** \copydoc Iterable::estimateNbItems */
    int64_t estimateNbItems ()   {  return getNbItems(); }
    
    Item* getItems (Item*& buffer)
    {
        std::cout << "IteratorFile::getItems(buffer) not implemented" << std::endl; exit(1);
    }
    
    /* from ../src/gatb/tools/collections/api/Iterable.hpp:
       Return a buffer of items.
        * \param[out] buffer : the buffer
        * \param[in] start : index where to start in the buffer --> NOTE: it is ignored here because we're not buffering anything (this ties in with a different behavior in istream of Storage with CurrentIdx, depending on STORAGE_FILE, currentIdx should be set to 0 when calling getItems, which in turns calls this function)
        * \param[in] nb : number of items to be retrieved
        * \return the number of items retrieved 
    */
    size_t getItems (Item*& buffer, size_t start, size_t nb)
    {
        if (_file == 0) 
            _file = system::impl::System::file().newFile (_filename, "rb"); 
        DEBUG_ITERATORFILE(std::cout << "want to read " << nb << " elements of size " << sizeof(Item) << " at position " << _file->tell() << " file size " << _file->getSize() /*<< " then write them to buffer at position " << (sizeof(Item) * start) << std::endl*/;)
        size_t n = _file->fread (buffer /*+ sizeof(Item) * start*/, sizeof(Item), nb);
        DEBUG_ITERATORFILE(std::cout << "read " << n << " elements" << std::endl;)
        return n;
    }

private:
    std::string     _filename;
    size_t          _cacheItemsNb;
    system::IFile*  _file;
};
    
/********************************************************************************/
/* EXPERIMENTAL (not documented). */
template <class Item> class IteratorGzFile : public dp::Iterator<Item>
{
public:
    
    /** Constructor. */
    IteratorGzFile () : _gzfile(0), _buffer(0), _cpt_buffer(0), _idx(0), _cacheItemsNb(0), _isDone(true) {}
    
    IteratorGzFile (const IteratorGzFile& it):
    _filename(it._filename), _gzfile(0),  _buffer(0), _cpt_buffer(0), _idx(0), _cacheItemsNb(it._cacheItemsNb), _isDone(true)
    {
        _gzfile =   gzopen(_filename.c_str(),"rb");
        gzbuffer(_gzfile,2*1024*1024);
        _buffer  = (Item*) MALLOC (sizeof(Item) * _cacheItemsNb);
    }
    
    /** Constructor. */
    IteratorGzFile (const std::string& filename, size_t cacheItemsNb=10000) :
    _filename(filename), _gzfile(0),  _buffer(0), _cpt_buffer(0), _idx(0), _cacheItemsNb(cacheItemsNb), _isDone(true)
    
    {
        _gzfile =   gzopen(_filename.c_str(),"rb");
        gzbuffer(_gzfile,2*1024*1024);
        _buffer  = (Item*) MALLOC (sizeof(Item) * _cacheItemsNb);
    }
    
    /** Destructor. */
    ~IteratorGzFile ()
    {
        if (_gzfile)  { gzclose(_gzfile);   }
        if (_buffer) { FREE (_buffer); }
    }
    
    /** Affectation. */
    IteratorGzFile& operator= (const IteratorGzFile& it)
    {
        if (this != &it)
        {
            if (_gzfile)    {  gzclose(_gzfile);  }
            if (_buffer)  { FREE(_buffer); }
            
            _filename     = it._filename;
            _cpt_buffer   = it._cpt_buffer;
            _idx          = it._idx;
            _cacheItemsNb = it._cacheItemsNb;
            _isDone       = it._isDone;
            
            _gzfile =   gzopen(_filename.c_str(),"r");
            gzbuffer(_gzfile,2*1024*1024);
            _buffer  = (Item*) MALLOC (sizeof(Item) * it._cacheItemsNb);
        }
        return *this;
    }
    
    /** \copydoc dp::Iterator::first */
    void first()
    {
        gzseek(_gzfile,0,SEEK_SET);
        _cpt_buffer = 0;
        _idx        = 0;
        _isDone     = false;
        next ();
    }
    
    /** \copydoc dp::Iterator::next */
    void next()
    {
        if (_cpt_buffer==0)
        {
            _idx = 0;
            _cpt_buffer = (gzread(_gzfile, _buffer, sizeof(Item)*_cacheItemsNb)  /  sizeof(Item))  ;  // gzread returns number of bytes uncompressed returned in _buffer
            //printf("refreshing buffer from gzread(), cptbuffer now: %d\n",_cpt_buffer);
            if (_cpt_buffer < 0)
            { // On other errors, gzread() shall return a value less than 0 and and applications may examine the cause using gzerror().
                // FIXME: more tests are needed but on my system (R: linux64), gzread returns the fixed number of items, then a lower number of items (as it reaches the end), then it returns -1 instead of 0 and prints a "data error" but at this point it's fine to just return
                int err;
                fprintf(stderr, "gzread error: %s\n", gzerror(_gzfile, &err));
            }
            if (_cpt_buffer<=0)  { _isDone = true;  return; }
        }
        *(this->_item) =  _buffer[_idx];
        _cpt_buffer --;
        _idx ++;
    }
    
    /** \copydoc dp::Iterator::isDone */
    bool isDone()  { return _isDone; }
    
    /** \copydoc dp::Iterator::item */
    Item& item ()  { return *(this->_item); }
    
    /** */
    size_t fill (std::vector<Item>& vec, size_t len=0)
    {
        if (len==0)  { len = vec.size(); }
        
        return (gzread(_gzfile,vec.data(), sizeof(Item)*len)  /  sizeof(Item));
 
    }
    
private:
    std::string     _filename;
    gzFile  _gzfile;
    Item*           _buffer;
    int             _cpt_buffer;
    int             _idx;
    size_t          _cacheItemsNb;
    bool            _isDone;
};

/********************************************************************************/
/* EXPERIMENTAL (not documented).
 * reading from a  sorted compressed file, read is buffered
 */
template <class Item> class IteratorCountCompressedFile : public dp::Iterator<Item>
{
public:
    
    /** Constructor. */
    IteratorCountCompressedFile () : _file(0), _buffer(0), _cpt_buffer(0), _idx(0), _cacheItemsNb(0), _isDone(true),_abundance(0) {}
    
    IteratorCountCompressedFile (const IteratorCountCompressedFile& it):
    _filename(it._filename), _file(0),  _buffer(0), _cpt_buffer(0), _idx(0), _cacheItemsNb(it._cacheItemsNb), _isDone(true),_abundance(0)
    {
        _file    = system::impl::System::file().newFile (_filename, "rb");
        _buffer  = (u_int8_t*) MALLOC (sizeof(u_int8_t) * _cacheItemsNb);
    }
    
    /** Constructor. */
    IteratorCountCompressedFile (const std::string& filename, size_t cacheItemsNb=10000) :
    _filename(filename), _file(0),  _buffer(0), _cpt_buffer(0), _idx(0), _cacheItemsNb(cacheItemsNb), _isDone(true),_abundance(0)
    
    {
        _file    = system::impl::System::file().newFile (filename, "rb");
        _buffer  = (u_int8_t*) MALLOC (sizeof(u_int8_t) * _cacheItemsNb);
    }
    
    /** Destructor. */
    ~IteratorCountCompressedFile ()
    {
        if (_file)  { delete _file;  }
        if (_buffer) { FREE (_buffer); }
    }
    
    /** Affectation. */
    IteratorCountCompressedFile& operator= (const IteratorCountCompressedFile& it)
    {
        if (this != &it)
        {
            if (_file)    { delete _file; }
            if (_buffer)  { FREE(_buffer); }
            
            _filename     = it._filename;
            _cpt_buffer   = it._cpt_buffer;
            _idx          = it._idx;
            _cacheItemsNb = it._cacheItemsNb;
            _isDone       = it._isDone;
            _abundance    = it._abundance;
            _file    = system::impl::System::file().newFile (it._filename, "rb");
            _buffer  = (u_int8_t*) MALLOC (sizeof(u_int8_t) * it._cacheItemsNb);
        }
        return *this;
    }
    
    /** \copydoc dp::Iterator::first */
    void first()
    {
        _file->seeko (0, SEEK_SET);
        _cpt_buffer = 0;
        _idx        = 0;
        _abundance  = 0;
        _isDone     = false;
        next ();
    }
    
    /** \copydoc dp::Iterator::next */
    void next()
    {
        if(_abundance)
        {
            *(this->_item) =  _previous;
            _abundance--;
        }
        else
        {
            if (!readChunkIfNeeded (1)) return;

            //read a couple (byte, Item)
            _abundance =  _buffer[_idx];
            _cpt_buffer --;  _idx ++;

            if (!readChunkIfNeeded (sizeof(Item))) return; // this one should succeed (ie, in the file always a couple  (byte, elem))
            _previous =  *((Item *) (_buffer + _idx ));
            _cpt_buffer -= sizeof(Item); _idx += sizeof(Item);
            

            *(this->_item) =  _previous;
            _abundance--;
        }
    }
    
    /** \copydoc dp::Iterator::isDone */
    bool isDone()  { return _isDone; }
    
    /** \copydoc dp::Iterator::item */
    Item& item ()  { return *(this->_item); }
    
    /** */
    size_t fill (std::vector<Item>& vec, size_t len=0)
    {
        printf("Not yet implemented \n");
        return 0;
    }
    
private:
    std::string     _filename;
    system::IFile*  _file;
    u_int8_t*       _buffer;
    int             _cpt_buffer; // how many unread bytes are remaining in the buffer
    int             _idx; // where we should read the next elem in the buffer
    size_t          _cacheItemsNb; //in bytes for this  compressed type file
    bool            _isDone;
    u_int8_t        _abundance ;
    Item            _previous;

    bool readChunkIfNeeded (size_t needNBytes)
    {
        if (_cpt_buffer < (int)needNBytes)
        {
            // printf("Read new chunk  _cacheItemsNb %zu B nedd %zu have %i  \n",_cacheItemsNb,needNBytes,_cpt_buffer);
            // _idx = 0;
            // printf(" A new pos in file %p  %llu  \n",_file, _file->tell());

            //deplacer ce quil reste au debut avant
            memcpy (_buffer,_buffer + _idx,_cpt_buffer ); _idx = 0;
            int remaining = _cpt_buffer;
            //std::cout << "(readChunkIfNeeded) doing a fread of " << (_cacheItemsNb - remaining) << " items at position " << _file->tell() << std::endl;
            _cpt_buffer += _file->fread (_buffer+ remaining , sizeof(u_int8_t), _cacheItemsNb - remaining  );
            // printf(" B new pos in file %p  %llu  \n",_file, _file->tell());

            if (_cpt_buffer==0)  {
                _isDone = true;
                //printf("should end iterating file \n");
            }
            // printf("end read have %i \n",_cpt_buffer);
        }
        return !_isDone;
    }
};

/********************************************************************************/
/* EXPERIMENTAL (not documented). */
template <class Item> class IterableGzFile : public tools::collections::Iterable<Item>, public virtual system::SmartPointer
{
public:
    
    /** */
    IterableGzFile (const std::string& filename, size_t cacheItemsNb=10000)
    :   _filename(filename), _cacheItemsNb (cacheItemsNb)  {}
    
    /** */
    ~IterableGzFile () {}
    
    /** */
    dp::Iterator<Item>* iterator ()  { return new IteratorGzFile<Item> (_filename, _cacheItemsNb); }
    
    /** */
    int64_t getNbItems ()   {  return -1; } // does not know value
    
    int64_t estimateNbItems ()   {  return 3* (system::impl::System::file().getSize(_filename) / sizeof(Item)); }

private:
    std::string     _filename;
    size_t          _cacheItemsNb;
};

/********************************************************************************/
/* EXPERIMENTAL (not documented). */
template <class Item> class IterableCountCompressedFile : public tools::collections::Iterable<Item>, public virtual system::SmartPointer
{
public:
    
    /** */
    IterableCountCompressedFile (const std::string& filename, size_t cacheItemsNb=10000)
    :   _filename(filename), _cacheItemsNb (cacheItemsNb)  {}
    
    /** */
    ~IterableCountCompressedFile () {}
    
    /** */
    dp::Iterator<Item>* iterator ()  { return new IteratorCountCompressedFile<Item> (_filename, _cacheItemsNb); }
    
    /** */
    int64_t getNbItems ()   {  return -1; } // does not know value
    
    int64_t estimateNbItems ()   {  return 2* (system::impl::System::file().getSize(_filename) / sizeof(Item)); }
    
private:
    std::string     _filename;
    size_t          _cacheItemsNb;
};

    
/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_COLLECTIONS_IMPL_ITERATOR_FILE_HPP_ */
