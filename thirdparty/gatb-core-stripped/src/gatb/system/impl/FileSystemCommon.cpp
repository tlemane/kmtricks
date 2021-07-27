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

#include <gatb/system/impl/FileSystemCommon.hpp>
#include <gatb/system/impl/System.hpp>

#include <sys/resource.h>
#include <sys/statvfs.h>
#include <stdlib.h>
#include <dirent.h>
#include <libgen.h>

#include <string>
#include <sstream>

using namespace std;

/********************************************************************************/
namespace gatb { namespace core { namespace system { namespace impl {
/********************************************************************************/

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
u_int64_t  FileSystemCommon::getAvailableSpace (const Path& path)
{
    struct statvfs buffer;

    statvfs (path.c_str(), &buffer);

    u_int64_t available = (buffer.f_bavail * buffer.f_bsize) / 1024;

    return available;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IFileSystem::Path FileSystemCommon::getCurrentDirectory ()
{
    char path[1000];
    char* buffer = getcwd (path, sizeof(path));

    if (buffer == 0)  {  throw ExceptionErrno ("unable to get current directory");  }

    return path;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IFileSystem::Path FileSystemCommon::getDirectory (const Path& path)
{
     size_t pos = path.find_last_of("\\/");
     return (std::string::npos == pos)  ? "."  : path.substr(0, pos);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IFileSystem::Path FileSystemCommon::getTemporaryDirectory ()
{
    const char* dir = 0;

         if ( (dir = getenv ("TMPDIR"))  != 0)  {  return dir;    }
    else if ( (dir = getenv ("TMP"))     != 0)  {  return dir;    }
    else if ( (dir = getenv ("TEMPDIR")) != 0)  {  return dir;    }
    else                                        {  return "/tmp"; }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS : Warning! this method isn't exactly basename() as you'd expect in C++. It returns the base name but cuts everything after the last dot.
*********************************************************************/
IFileSystem::Path FileSystemCommon::getBaseName (const Path& path, bool cutToFirstDot)
{
    /** We duplicate the provided path. */
    char* reads_path = strdup (path.c_str());

    /** We build the basename; it still may have a suffix. */
    std::string reads_name (basename(reads_path)); // posix basename() may alter reads_path

    /** We release the duplicated path. */
    free (reads_path);

	//string prefix = System::file().getBaseName(_inputFilename);;
	while (reads_name.find('.') != string::npos){ // make sure there is a dot in the file, else the basename is the file itself

	    /** We look for the beginnin of the suffix. */
		int lastindex = reads_name.find_last_of(".");

	    /** We build the result. */
		reads_name = reads_name.substr(0, lastindex);

        if (cutToFirstDot == false)
            break;
	}

    //int lastindex = reads_name.find_last_of (".");
    //Path result = reads_name.substr(0, lastindex);

    /** We return the base name, without suffix. */
    return reads_name;
}


std::string FileSystemCommon::getExtension(const Path &path)
{
    std::string spath = string(path.c_str());
    return spath.substr(spath.find_last_of(".") + 1);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IFileSystem::Path FileSystemCommon::getRealPath (const Path& file)
{
    /** We must use a big buffer. Previously was 1024 but we got crashes on
     * the CI server with this value... */
    char buf [4*1024];

    if (realpath (file.c_str(), buf) != 0)
    {
        return buf;
    }
    throw Exception ("Unable to get the real path for '%s'", file.c_str());
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
string FileSystemCommon::getTemporaryFilename (const std::string& filename)
{
    stringstream ss;
    ss << tmp_prefix() << "_" << System::thread().getProcess();
    if (filename.empty()==false)  { ss << "_" << filename; }
    return ss.str();
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
bool FileSystemCommon::doesExist (const Path& path)
{
    FILE* fp = fopen (path.c_str(), "rb");

   if (fp != NULL)
   {
       fclose (fp);
       return true;
   }
   else
   {
       return false;
   }
}

bool FileSystemCommon::doesExistDirectory (const Path& path)
{
   DIR* dir = opendir(path.c_str());
   if (dir)
   {
           /* Directory exists. */
           closedir(dir);
           return true;
   }
   return false;
}


bool FileSystemCommon::isFolderEndingWith (const Path& path, const std::string &ending)
{
    if (!doesExistDirectory(path))
        return false;

    // strip '/' at end of path
    std::string path_without_slashes = path;
   if (path_without_slashes.length() > 0)
   {
       std::string::iterator it = path_without_slashes.end() - 1;
       if (*it == '/')
           path_without_slashes.erase(it);
   }

    // http://stackoverflow.com/questions/874134/find-if-string-ends-with-another-string-in-c
    if (path_without_slashes.length() >= ending.length()) {
        if (0 == path_without_slashes.compare (path_without_slashes.length() - ending.length(), ending.length(), ending))
            return true;
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
u_int64_t  FileSystemCommon::getSize (const Path& path)
{
    struct stat st;

    if (stat (path.c_str(), &st) == 0) return st.st_size;

    return 0;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void FileSystemCommon::iterate (const Path& path, void (*callback) (const Path&, void* data), void* data)
{
    DIR* dp = opendir (path.c_str());

    if (dp)
    {
        struct dirent* dirp = 0;

        while ( (dirp = readdir(dp)) != 0)
        {
            callback (dirp->d_name, data);
        }

        closedir (dp);
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
std::vector<std::string> FileSystemCommon::listdir (const Path& path)
{
	std::vector<std::string> filenames;
    DIR* dp = opendir (path.c_str());

    if (dp)
    {
        struct dirent* dirp = 0;

        while ( (dirp = readdir(dp)) != 0)
        {
        	filenames.push_back(std::string(dirp->d_name));
        }

        closedir (dp);
    }

    return filenames;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
#if 0
IFile* FileSystemCommon::newFile (const Path& dirpath, const Path& filename, const char* mode)
{
    /** We build the full file path. */
    stringstream ss;
    ss << dirpath << "/" << filename;

    /** We create the file handle. */
    return newFile (ss.str(), mode);
}
#endif

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
