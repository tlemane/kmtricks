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

/** \file IFileSystem.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Operating System abstraction for file system management.
 */

#ifndef _GATB_CORE_SYSTEM_IFILE_HPP_
#define _GATB_CORE_SYSTEM_IFILE_HPP_

/********************************************************************************/

#include <gatb/system/api/types.hpp>
#include <string>

/********************************************************************************/
namespace gatb      {
namespace core      {
 namespace system   {
/********************************************************************************/

/** \brief Abstraction of what we need about file
 *
 *  We define here a few methods we need for handling files.
 *
 *  NOTE: deleting a IFile instance won't delete the file on the file system; one
 *  should see IFile as a logical handle on a physical file.
 */
class IFile
{
public:

    /** Tells whether the file is opened or not.
     * \return open status.
     */
    virtual bool isOpen () = 0;

    /** Tells whether or not we are at the end of a file.
     *  \return true if we are at the end of the file, false otherwise.
     */
    virtual bool isEOF () = 0;

    /** Locates the cursor into the file (similar to 'fseeko' functions)
     * \param[in] offset : bytes number we want to go, relatively from the 'whence' paramater
     * \param[in] whence : SEEK_SET, SEEK_END, or SEEK_CUR
     * \return 0 if successful, -1 otherwise
     */
    virtual int seeko (u_int64_t offset, int whence) = 0;

    /** Get the position in a file.
     * \return the current position
     */
    virtual u_int64_t tell () = 0;

    /** Get the currently pointed character in the file.
     * \return the current character.
     */
    virtual int get () = 0;

    /** Unget the currently pointed character in the file.
     * \param[in] c : the current character.
     * \return c on success, EOF on error
     */
    virtual int unget (int c) = 0;

    /** Reads a line in the file.
     *  \param[in] s : buffer where to put the read line.
     *  \param[in] size : maximum number of characters to be read
     *  \return the actual size of the read buffer (0 if nothing read)
     */
    virtual int gets (char* s, int size) = 0;

    /** Writes a buffer (0 terminated) into the file with a new line.
     *  \param[in] format : format of the data to be dumped into the file.
     *  \param[in] ...    : arguments (as an ellipsis) to b dumped
     */
    virtual void print (const char* format, ...) = 0;

    /** Reads a buffer from the file.
     * \param[in] ptr : the buffer to be read
     * \param[in] size : size of the buffer
     * \param[in] nmemb : number of elements to be written
     * \return number of  items successfully written */
    virtual size_t fread (void* ptr, size_t size, size_t nmemb) = 0;

    /** Writes a buffer into the file.
     * \param[in] ptr : the buffer to be written
     * \param[in] size : size of the buffer
     * \param[in] nmemb : number of elements to be written
     * \return number of  items successfully written */
    virtual size_t fwrite (const void* ptr, size_t size, size_t nmemb) = 0;

    /** Flush the file.
     */
    virtual void flush () = 0;

    /** Get the size of a file.
     * \return the size of the file.
     */
    virtual u_int64_t getSize () = 0;

    /** Get the file URI
     * \return the URI of the file */
    virtual const std::string& getPath () const = 0;

    /** Destructor. */
    virtual ~IFile () {}
};

/********************************************************************************/

/** \brief interface for some operations at file system level.
 *
 * This interface define a few operations like creating/deleting directories.
 *
 * It can provide information like the max number of files that can be used at the
 * same time, or the available space at some location in the FFS.
 *
 * It also acts as a factory that creates IFile instances.
 *
 * \see IFile
 */
 class IFileSystem
 {
 public:

     /** Alias type for a file system path. */
     typedef std::string Path;

     /** Return the maximum number of files that can be used at the same time.
      * \return the max number of files. */
     virtual size_t getMaxFilesNumber () = 0;

     /** Return the available space at the location given by the provided path.
      * \param[in] path : the location from where the space size is computed.
      * \return the available size (in KBytes). */
     virtual u_int64_t  getAvailableSpace (const Path& path) = 0;

     /** Retrieve the current directory absolute path.
      * \return the current directory to be retrieved.
      */
     virtual Path getCurrentDirectory () = 0;

     /** Retrieve the directory of the provided name.
      * \return the directory
      */
     virtual Path getDirectory (const Path& path) = 0;

     /** Retrieve the default temporary directory absolute path.
      * \return the temporary directory to be retrieved.
      */
     virtual Path getTemporaryDirectory () = 0;

     /** Return the base of the URI
      * \param[in] path : uri from which we want to extract the base name.
      * \param[in] cutToFirstDot : by default, getBaseName("a.b.c.fq") = "a.b.c". If cutToFirstDot is true, will return just "a"
      * \return the base name of the uri
      */
     virtual Path getBaseName (const Path& path, bool cutToFirstDot = false) = 0;

     /** Return the canonical path to the given file, ie replace symbolic links or relative path.
      * \param[in] file : the file we want the canonical path.
      * \return the real path. */
     virtual Path getRealPath (const Path& file) = 0;

     /** Return the extension to the given file,
      * \param[in] file : file
      * \return extension */
     virtual std::string getExtension (const Path& file) = 0;
 
     /** Get a temporary file name. One may provide an argument; in such a case some prefix/suffix will
      * be appended to this name in order to make it unique.
      * \param[in] filename : file name (may be empty)
      * \return a unique file name. */
     virtual std::string getTemporaryFilename (const std::string& filename="") = 0;

     /** Tells whether a file exists or not.
      * \return true if file exists, false otherwise
      */
     virtual bool doesExist (const Path& path) = 0;

     /** Tells whether a folder exists or not.
      * \return true if folder exists, false otherwise
      */
     virtual bool doesExistDirectory (const Path& path) = 0;
 
     /** Tells whether path is a folder that ends with a certain string
      * \return true if folder ends with string, false otherwise
      */
     virtual bool isFolderEndingWith (const Path& path, const std::string &ending) = 0;

     /** Return the size of the file given by the provided path.
      * \param[in] path : the location from where the space size is computed.
      * \return the available size (in Bytes). */
     virtual u_int64_t  getSize (const Path& path) = 0;

     /** Clear the file system cache. */
     virtual int clearCache () = 0;

     /** Create a directory for the provided path and mode.
      * \param[in] path : path of the directory to be created.
      * \param[in] mode : mode of creation (see 'mkdir' system function documentation). */
     virtual int mkdir (const Path& path, u_int64_t mode) = 0;

     /** Delete the entry given its path.
      * \param[in] path : path of the entry to be removed from the file system. */
     virtual int rmdir (const Path& path) = 0;

     /** Delete the file given its path.
      * \param[in] path : path of the file to be removed from the file system. */
     virtual int remove (const Path& path) = 0;

     /** Rename the provided uri
      * \param[in] from : initial name.
      * \param[in] to : final name. */
     virtual int rename (const Path& from, const Path& to) = 0;

     /** Iterates entries of a directory.
      * \param[in] path : path of the directory to be iterated
      * \param[in] callback : callback called for each found entry
      * \param[in] data : private data given to the callback. */
     virtual void iterate (const Path& path, void (*callback) (const Path& entry, void* data), void* data) = 0;

     /** Get entries of a directory.
      * \param[in] path : path of the explored directory */
     virtual std::vector<Path> listdir(const Path& path) = 0;

     /** Creates a new IFile instance (equivalent to 'fopen' function)
      * \param[in] path : uri of the file to be opened.
      * \param[in] mode : mode of the file (like fopen)
      * \return instance of IFile, 0 otherwise.
      */
     virtual IFile* newFile (const Path& path, const char* mode) = 0;

     /** Creates a new IFile instance (equivalent to 'fopen' function)
      * \param[in] dirpath  : uri of the directory where the file is meant to be.
      * \param[in] filename : name of the file
      * \param[in] mode : mode of the file (like fopen)
      * \return instance of IFile, 0 otherwise.
      */
     virtual IFile* newFile (const Path& dirpath, const Path& filename, const char* mode) = 0;

     /** Get metadata associated with the file for a given key.
      * \param[in] filename : name of the file.
      * \param[in] key : key for which we want the value
      * \param[in] value : value associated to the key
      * \return -1 if KO, 0 otherwise, otherwise length of the retrieved value. */
     virtual ssize_t getAttribute (const Path& filename, const char* key, std::string& value) = 0;

     /** Set metadata associated with the file for a given key.
      * \param[in] filename : name of the file.
      * \param[in] key : key for which we want the value
      * \param[in] fmt : format of the string (in the 'printf' way)
      * \param[in] ... : paramters for the format parameter
      * \return -1 if error, 0 otherwise. */
     virtual ssize_t setAttribute (const Path& filename, const char* key, const char* fmt, ...) = 0;

     /** Destructor. */
     virtual ~IFileSystem () {}
};

/********************************************************************************/
} } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_SYSTEM_IFILE_HPP_ */
