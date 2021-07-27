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

#ifdef __linux__

#include <gatb/system/impl/FileSystemLinux.hpp>

#include <sys/resource.h>
#include <sys/statvfs.h>
#include <stdlib.h>
#include <dirent.h>

#include <sys/xattr.h>

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

class FileLinux : public CommonFile
{
public:

    /** Constructor. */
    FileLinux (const char* path, const char* mode) : CommonFile(path,mode)  { }

    /** \copydoc IFile::tell */
    u_int64_t tell ()  { return (isOpen() ? ftello64 (_handle) : 0); }

	/** \copydoc IFile::seeko */
    int seeko (u_int64_t offset, int whence)  {  return (isOpen() ?  fseek /* cygwin doesnt like fseeko and fseek/fseeko seems similar */(_handle, offset, whence) : -1);  }
};

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
size_t FileSystemLinux::getMaxFilesNumber ()
{
    size_t result = 0;

    struct rlimit64 lim;

    if (getrlimit64 (RLIMIT_NOFILE, &lim) == 0)   {  result = lim.rlim_cur;  }

    return result;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IFile* FileSystemLinux::newFile (const Path& path, const char* mode)
{
    return new FileLinux (path.c_str(), mode);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
IFile* FileSystemLinux::newFile (const Path& dirpath, const Path& filename, const char* mode)
{
    /** We build the full file path. */
    stringstream ss;
    ss << dirpath << "/" << filename;

    /** We create the file handle. */
    return newFile (ss.str(), mode);
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
int FileSystemLinux::clearCache ()
{
    int result = EXIT_FAILURE;

    //result = ::system ("echo 3 > /proc/sys/vm/drop_caches");

    return result;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
ssize_t FileSystemLinux::getAttribute (const Path& filename, const char* key, string& value)
{
    char buffer[4*1024];

    value.clear();

    ssize_t res = ::getxattr (filename.c_str(), (string("user.") + key).c_str(), buffer, sizeof(buffer));

    if (res >= 0)   { value.assign (buffer, res); }

    return res;
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
ssize_t FileSystemLinux::setAttribute (const Path& filename, const char* key, const char* fmt, ...)
{
    char buffer[4*1024];

    va_list ap;
    va_start (ap, fmt);
    vsnprintf (buffer, sizeof(buffer), fmt, ap);
    va_end (ap);

    return ::setxattr (filename.c_str(), (string("user.") + key).c_str(), buffer, strlen(buffer), XATTR_CREATE);
}


/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* __LINUX__ */
