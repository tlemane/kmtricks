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

/** \file FileSystemMacos.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Implementation for MacOs.
 */

#ifndef _GATB_CORE_SYSTEM_IMPL_FILE_SYSTEM_MACOS_HPP_
#define _GATB_CORE_SYSTEM_IMPL_FILE_SYSTEM_MACOS_HPP_

/********************************************************************************/

#include <gatb/system/impl/FileSystemCommon.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace system    {
/** \brief Implementation of Operating System abstraction layer */
namespace impl      {
/********************************************************************************/

/** \brief default implementation
 */
class FileSystemMacos : public FileSystemCommon
{
public:

	/** \copydoc IFileSystem::getMaxFilesNumber */
    size_t getMaxFilesNumber ();

    /** \copydoc IFileSystem::getAttribute */
    ssize_t getAttribute (const Path& filename, const char* key, std::string& value);

    /** \copydoc IFileSystem::setAttribute */
    ssize_t setAttribute (const Path& filename, const char* key, const char* fmt, ...);

    /** \copydoc IFileSystem::newFile */
    IFile* newFile (const Path& path, const char* mode);

    /** \copydoc IFileSystem::newFile(const Path&, const Path&, const char*) */
    IFile* newFile (const Path& dirpath, const Path& filename, const char* mode);

    /** \copydoc IFileSystem::clearCache */
    int clearCache ();
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_SYSTEM_IMPL_FILE_SYSTEM_MACOS_HPP_ */
