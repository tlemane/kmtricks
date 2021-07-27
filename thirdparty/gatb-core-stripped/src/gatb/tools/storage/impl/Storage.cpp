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

/********************************************************************************/

#include <gatb/tools/storage/impl/Storage.hpp>

#include <iostream>
#include <fstream>

/********************************************************************************/
namespace gatb { namespace core {  namespace tools {  namespace storage {  namespace impl {
/********************************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Storage::Storage (StorageMode_e mode, const std::string& name, bool autoRemove)
    : Cell(0, ""), _name(name), _factory(0), _root(0), _autoRemove(autoRemove)
{
    setFactory (new StorageFactory (mode));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Storage::~Storage ()
{
    setRoot    (0);
    setFactory (0);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Group* Storage::getRoot ()
{
    if (_root == 0)  { setRoot    (_factory->createGroup (this, ""));   _root->setCompressLevel (this->getCompressLevel()); }
    return _root;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Group& Storage::operator() (const std::string name)
{
    if (name.empty())  { return *getRoot(); }
    else               { return getRoot ()->getGroup (name);  }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Storage::remove ()
{
    getRoot()->remove();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void Storage::setFactory (StorageFactory* factory)
{
    SP_SETATTR(factory);
}

/********************************************************************************
             #####   #######  ######   #######     #     #     #
            #     #     #     #     #  #          # #    ##   ##
            #           #     #     #  #         #   #   # # # #
             #####      #     ######   #####    #     #  #  #  #
                  #     #     #   #    #        #######  #     #
            #     #     #     #    #   #        #     #  #     #
             #####      #     #     #  #######  #     #  #     #
********************************************************************************/

/* WRAPPERS BETWEEN STORAGE AND C++ STREAMS
 *
 * Got inspiration from :
 *      http://savingyoutime.wordpress.com/2009/04/21/using-c-stl-streambufostream-to-create-time-stamped-logging-class/
 *      http://www.mr-edd.co.uk/blog/beginners_guide_streambuf
 */

/* Output stream buffer implementation. */
class Storage_ostreambuf : public std::streambuf
{
protected:

    static const int bufferSize = 4*1024;   // size of data buffer
    char buffer[bufferSize];                // data buffer

public:
    Storage_ostreambuf (Group& group, const std::string& name) : _nbWritten(0)
    {
        setp (buffer, buffer+(bufferSize-1));
        _collection = & group.getCollection<math::NativeInt8> (name);
    }

    virtual ~Storage_ostreambuf() { 
        sync(); 
        //std::cout << "ostream destructor" << std::endl; 
        }

protected:

    collections::Collection<math::NativeInt8>* _collection;

    // flush the characters in the buffer
    int flushBuffer ()
    {
        int num = pptr()-pbase();
        _collection->insert ((math::NativeInt8*)buffer, num);
        _collection->flush(); // important!
        _nbWritten += num;
        pbump(-num); // reset put pointer accordingly
        return num;
    }

    virtual int overflow ( int c = EOF )
    {
        if (c != EOF) {
            *pptr() = c;    // insert character into the buffer
            pbump(1);
        }
        if (flushBuffer() == EOF)
            return EOF;
        return c;
    }

    virtual int sync()
    {
        if (flushBuffer() == EOF) {  return -1; }  // ERROR
        return 0;
    }

    virtual pos_type  seekoff (off_type off, std::ios_base::seekdir dir,  std::ios_base::openmode mode)
    {
        sync ();  // We may have to flush the current buffer first
        return _nbWritten;
    }

    pos_type _nbWritten;
};

/*********************************************************************
*********************************************************************/

/* Input stream buffer implementation. */
class Storage_istreambuf : public std::streambuf
{
    public:
    Storage_istreambuf (Group& group, const std::string& name, std::size_t buff_sz = 1024, std::size_t put_back = 64) :
            put_back_(std::max(put_back, size_t(1))),
            buffer_(std::max(buff_sz, put_back_) + put_back_), currentIdx(0)
        {
            char *end = &buffer_.front() + buffer_.size();
            setg(end, end, end);
            _collection = & group.getCollection<math::NativeInt8> (name);
        }

    private:
        // overrides base class underflow()
        int_type underflow()
        {
            if (gptr() < egptr()) // buffer not exhausted
                return traits_type::to_int_type(*gptr());

            char *base = &buffer_.front();
            char *start = base;

            if (eback() == base) // true when this isn't the first fill
            {
                // Make arrangements for putback characters
                std::memmove(base, egptr() - put_back_, put_back_);
                start += put_back_;
            }
            // start is now the start of the buffer, proper.
            // Read from fptr_ in to the provided buffer
            math::NativeInt8* start2 = (math::NativeInt8*) start;
            //std::cout << "(storage) istream calling getItems with params: start2=" << start2 << " currentIdx=" << currentIdx << " buffer.size()-start+base=" << (buffer_.size() - (start - base)) << ", total buffer size " << buffer_.size()  << std::endl;
            if ( ((int64_t)base > (int64_t)start) || (((int64_t)start - (int64_t)base) > (int64_t)buffer_.size()))
            {
                std::cout << "Error: trying to read more elements " << (start - base) << " = (" << start << " - " << base << ") than the buffer size" << std::endl; exit(1);
            }
            size_t offset = currentIdx ; // in h d f 5: needs to be currentIdx; in file: needs to be 0 (but it will be fixed in IteratorFile, ok, not here)
            size_t n = _collection->getItems (start2, offset, buffer_.size() - (start - base));
            currentIdx += n;

            if (n == 0)  {   return traits_type::eof();  }

            // Set buffer pointers
            setg(base, start, start + n);

            return traits_type::to_int_type(*gptr());
        }

        // copy ctor and assignment not implemented;
        // copying not allowed
        Storage_istreambuf(const Storage_istreambuf &);
        Storage_istreambuf &operator= (const Storage_istreambuf &);

    private:
        collections::Collection<math::NativeInt8>* _collection;

        const std::size_t put_back_;
        std::vector<char> buffer_;
        size_t currentIdx;
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Storage::ostream::ostream (Group& group, const std::string& name)
    : std::ios(0), std::ostream(new Storage_ostreambuf(group,name))
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Storage::ostream::~ostream()
{
    delete rdbuf();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Storage::istream::istream (Group& group, const std::string& name)
    : std::ios(0), std::istream(new Storage_istreambuf(group,name))
{
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
Storage::istream::~istream ()
{
    delete rdbuf();
}
	
	
///////////////////////////////////////
////////// SuperKmerBinFiles //////////
///////////////////////////////////////
	
SuperKmerBinFiles::SuperKmerBinFiles(const std::string& path,const std::string& name, size_t nb_files, bool lz4) : _basefilename(name), _path(path),_nb_files(nb_files) 
{
	_nbKmerperFile.resize(_nb_files,0);
	_FileSize.resize(_nb_files,0);
	
	openFiles("wb"); //at construction will open file for writing
	// then use close() and openFiles() to open for reading
	
}

// Rayan: added for dsk separation into phase1 and phase2
SuperKmerBinFiles::SuperKmerBinFiles(const std::string& prefix, bool lz4)
{
    std::ifstream myfile (prefix+"/SuperKmerBinInfoFile");
    std::string line;
    if (myfile.is_open())
    {
        getline (myfile,line); _basefilename = line; 
        getline (myfile,line); _path = line;
        getline (myfile,line); _nb_files = atoi(line.c_str());
	
        _nbKmerperFile.resize(_nb_files,0);
    	_FileSize.resize(_nb_files,0);
	    _files.resize(_nb_files,0);
 	   	_synchros.resize(_nb_files,0);
        for(unsigned int ii=0;ii<_files.size();ii++)
        {
            getline (myfile,line); _nbKmerperFile[ii] = atol(line.c_str());
            getline (myfile,line); _FileSize[ii] = atol(line.c_str());
        }
        myfile.close();
    }
    else
    {
        std::cout << "error reading infoFile at " << prefix << std::endl; exit(1);
    }
}

void SuperKmerBinFiles::saveInfoFile(const std::string& prefix)
{
  // noob C++ style :) -Rayan
  std::ofstream myfile;
  myfile.open (prefix + "/SuperKmerBinInfoFile");
  myfile << _basefilename << "\n";
  myfile << _path << "\n";
  myfile << _nb_files << "\n";
  for(unsigned int ii=0;ii<_files.size();ii++)
  {
	myfile << _nbKmerperFile[ii] << "\n";
	myfile << _FileSize[ii] << "\n";
  }
  myfile.close();
}

void SuperKmerBinFiles::openFile( const char* mode, int fileId)
{
	std::stringstream ss;
	ss << _basefilename << "." << fileId;
	
	_files[fileId] = system::impl::System::file().newFile (_path, ss.str(), mode);
	_synchros[fileId] = system::impl::System::thread().newSynchronizer();
	_synchros[fileId]->use();
}
	
void SuperKmerBinFiles::openFiles( const char* mode)
{
	_files.resize(_nb_files,0);
	_synchros.resize(_nb_files,0);
	
	system::impl::System::file().mkdir(_path, 0755);

	for(unsigned int ii=0;ii<_files.size();ii++)
	{
		std::stringstream ss;
		ss << _basefilename << "." << ii;
		_files[ii] = system::impl::System::file().newFile (_path, ss.str(), mode);
		_synchros[ii] = system::impl::System::thread().newSynchronizer();
		_synchros[ii]->use();

	}
}

	
std::string SuperKmerBinFiles::getFileName(int fileId)
{
	
	std::stringstream ss;
	ss << _path << "/" <<_basefilename << "." << fileId;
	
	return ss.str();
}

	
int SuperKmerBinFiles::readBlock(unsigned char ** block, unsigned int* max_block_size, unsigned int* nb_bytes_read, int file_id)
{
	_synchros[file_id]->lock();
	
	int nbr = _files[file_id]->fread(nb_bytes_read, sizeof(*max_block_size),1);
	std::cerr << "part " << file_id << std::endl;
	std::cerr << "nbr " << nbr << std::endl;
	std::cerr << "max " << *max_block_size << std::endl;
	std::cerr << "read " << *nb_bytes_read << std::endl;
	if(nbr == 0)
	{
		//printf("__ end of file %i __\n",file_id);
		_synchros[file_id]->unlock();
		return 0;
	}
	
	if(*nb_bytes_read > *max_block_size)
	{
		*block = (unsigned char *) realloc(*block, *nb_bytes_read);
		*max_block_size = *nb_bytes_read;
	}
	
	//block
	_files[file_id]->fread(*block, sizeof(unsigned char),*nb_bytes_read);
	
	_synchros[file_id]->unlock();
	
	return *nb_bytes_read;
}

int SuperKmerBinFiles::getNbItems(int fileId)
{
	return _nbKmerperFile[fileId];
}
	
	
u_int64_t SuperKmerBinFiles::getFileSize(int fileId)
{

	return _FileSize[fileId];
}

void SuperKmerBinFiles::getFilesStats(u_int64_t & total, u_int64_t & biggest, u_int64_t & smallest, float & mean)
{
	total =0;
	smallest = ~0;
	biggest  = 0;
	mean=0;
	for(unsigned int ii=0;ii<_FileSize.size();ii++)
	{
		smallest = std::min (smallest, _FileSize[ii]);
		biggest  = std::max (biggest,  _FileSize[ii]);
		total+=_FileSize[ii];
	}
	if(_FileSize.size()!=0)
		mean= total/_FileSize.size();
	
}
	
	
void SuperKmerBinFiles::writeBlock(unsigned char * block, unsigned int block_size, int file_id, int nbkmers)
{
	_synchros[file_id]->lock();
	
	_nbKmerperFile[file_id]+=nbkmers;
	_FileSize[file_id] += block_size+sizeof(block_size);
	//block header
	_files[file_id]->fwrite(&block_size, sizeof(block_size),1);

	//block
	_files[file_id]->fwrite(block, sizeof(unsigned char),block_size);
	std::cerr << "P " << file_id << " block size " << block_size << std::endl;
	_synchros[file_id]->unlock();
}
	
void SuperKmerBinFiles::flushFiles()
{
	for(unsigned int ii=0;ii<_files.size();ii++)
	{
		_synchros[ii]->lock();

		if(_files[ii]!=0 || _files[ii]!=0)
		{
			_files[ii]->flush();
		}
		
		_synchros[ii]->unlock();
	}
}

void SuperKmerBinFiles::eraseFiles()
{
	for(unsigned int ii=0;ii<_files.size();ii++)
	{
		std::stringstream ss;
		ss << _path << "/" <<_basefilename << "." << ii;
		system::impl::System::file().remove(ss.str());
	}
	system::impl::System::file().rmdir(_path);
}

void SuperKmerBinFiles::eraseFile(int fileId)
{
	std::stringstream ss;
	ss << _path << "/" << _basefilename << "." << fileId;
	system::impl::System::file().remove(ss.str());
}

void SuperKmerBinFiles::closeFile(  int fileId)
{
	if(_files[fileId]!=0)
	{
		delete _files[fileId];
		_files[fileId] = 0;
		_synchros[fileId]->forget();
	}
}

	
void SuperKmerBinFiles::closeFiles()
{
	for(unsigned int ii=0;ii<_files.size();ii++)
	{
		if(_files[ii]!=0)
		{
			delete _files[ii];
			_files[ii] = 0;
			_synchros[ii]->forget();
		}
	}
}
	
SuperKmerBinFiles::~SuperKmerBinFiles()
{
	this->closeFiles();
	//this->eraseFiles();
}
	
int SuperKmerBinFiles::nbFiles()
{
	return _files.size();
}

////////////////////////////////////////////
//////////  CacheSuperKmerBinFiles /////////
////////////////////////////////////////////



	
CacheSuperKmerBinFiles::CacheSuperKmerBinFiles(SuperKmerBinFiles * ref, size_t buffsize )
{
	_ref = ref;

	_nb_files = _ref->nbFiles();
	_nbKmerperFile.resize(_nb_files,0);

	_buffer_max_capacity = buffsize; // this is per file, per thread
	//printf("buffsize %i per file per thread \n",_buffer_max_capacity);

	_max_superksize= 255; // this is extra size from regular kmer; ie total max superksize is kmersize +  _max_superksize
	
	_buffers.resize(_nb_files);
	_buffers_idx.resize(_nb_files,0);
	
	for(unsigned int ii=0; ii<_buffers.size();ii++ )
	{
		_buffers[ii] = (u_int8_t*) MALLOC (sizeof(u_int8_t) * _buffer_max_capacity);
	}
	
}
	
//copy construc : alloc own buffer for new object
CacheSuperKmerBinFiles::CacheSuperKmerBinFiles (const CacheSuperKmerBinFiles& p)
{
	_ref = p._ref;
	_nb_files= p._nb_files;
	_buffer_max_capacity= p._buffer_max_capacity;
	_max_superksize= p._max_superksize;
	_nbKmerperFile.resize(_nb_files,0);

	_buffers.resize(_nb_files);
	_buffers_idx.resize(_nb_files,0);
	
	for(unsigned int ii=0; ii<_buffers.size();ii++ )
	{
		_buffers[ii] = (u_int8_t*) MALLOC (sizeof(u_int8_t) * _buffer_max_capacity);
	}
}
	
void CacheSuperKmerBinFiles::flushAll()
{
	//printf("flush all buffers\n");
	for(unsigned int ii=0; ii<_buffers.size();ii++ )
	{
		flush(ii);
	}
}
	
	
void CacheSuperKmerBinFiles::flush(int file_id)
{
	if(_buffers_idx[file_id]!=0)
	{
		_ref->writeBlock(_buffers[file_id],_buffers_idx[file_id],file_id,_nbKmerperFile[file_id]);
		
		_buffers_idx[file_id]=0;
		_nbKmerperFile[file_id] = 0;
	}
}
	
	
void CacheSuperKmerBinFiles::insertSuperkmer(u_int8_t* superk, int nb_bytes, u_int8_t nbk, int file_id)
{
	std::cerr << "part " << file_id << std::endl;
	std::cerr << "bytes" << nb_bytes << std::endl;
	if( (_buffers_idx[file_id]+nb_bytes+1) > _buffer_max_capacity)
	{
		flush(file_id);
	}
	
	_buffers[file_id][_buffers_idx[file_id]++] = nbk;
	
	memcpy(_buffers[file_id] + _buffers_idx[file_id]  , superk,nb_bytes);
	_buffers_idx[file_id] += nb_bytes;
	_nbKmerperFile[file_id]+=nbk;
	
}
	
CacheSuperKmerBinFiles::~CacheSuperKmerBinFiles()
{
	this->flushAll();
	for(unsigned int ii=0; ii<_buffers.size();ii++ )
	{
		FREE (_buffers[ii]);
	}
}
/********************************************************************************/
} } } } } /* end of namespaces. */
/********************************************************************************/

