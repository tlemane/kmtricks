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

/** \file StringsRepository.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Pool of strings providing information to end users
 */

#ifndef _GATB_CORE_TOOLS_MISC_STRINGS_REPOSITORY_HPP_
#define _GATB_CORE_TOOLS_MISC_STRINGS_REPOSITORY_HPP_

/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Tools package */
namespace tools     {
/** \brief Misc interfaces */
namespace misc      {
/********************************************************************************/

/* \brief Pool of strings
 *
 * This class provides constant strings used throughout the code. It may be interesting
 * in order to have a central point for constant strings management:
 *      - ease translation in different languages
 *      - entry point for strings obfuscation if needed
 *
 * // rayan's remark: I respectfully disagree, this is not useful: we won't translate and we won't obfuscate
 *
 * It could also be possible to read the strings from a configuration file.
 *
 * The class defines one (static) method per constant string to be used. Note that we
 * also (see below) define a macro definition that eases the use of such a facility.
 */

class MessageRepository
{
public:

    /** \brief Singleton method.
     *
     * This method could return different types for the string repository, for translation
     * for instance.
     *
     * \return the singleton instance.
     */
    static MessageRepository& singleton()  { static MessageRepository instance; return instance; }

    const char* BANK_bad_file_number    () { return "bank files number is %d but should be in [1..%d]"; }
    const char* BANK_bad_file_path      () { return "unable to find file '%s'"; }
    const char* BANK_unable_open_file   () { return "error opening file: %s"; }
    const char* BANK_unable_write_file  () { return "unable to write into file"; }
};

/********************************************************************************/

/* \brief Pool of strings
 */
class StringRepository
{
public:
    static StringRepository& singleton()  { static StringRepository instance; return instance; }

    const char* db             ()  { return "-db";             }
    const char* file           ()  { return "-file";           }
    const char* graph          ()  { return "-graph";          }
    const char* kmer_size      ()  { return "-kmer-size";      }
    const char* minimizer_size ()  { return "-minimizer-size"; }
    const char* edge_km_representation ()  { return "-edge-km"; }
    const char* all_abundance_counts   ()  { return "-all-abundance-counts"; }
    const char* kmer_abundance ()  { return "-abundance"; }
    const char* kmer_abundance_min ()  { return "-abundance-min"; }
    const char* kmer_abundance_min_threshold ()  { return "-abundance-min-threshold"; }
    const char* kmer_abundance_max ()  { return "-abundance-max"; }
    const char* max_memory     ()  { return "-max-memory";     }
    const char* max_disk       ()  { return "-max-disk";       }
    const char* kmer_solid     ()  { return "-kmer-solid";     }
    const char* kmer_cFP       ()  { return "-kmer-cFP";       }
    const char* prefix         ()  { return "-prefix";         }
    const char* progress_bar   ()  { return "-bargraph";       }
    const char* nb_cores       ()  { return "-nb-cores";       }
    const char* partition_type ()  { return "-partition-type"; }
    const char* histogram_max  ()  { return "-histo-max";      }
    const char* uri_debloom    ()  { return "-debloom";        }
    const char* uri_input      ()  { return "-in";             }
    const char* uri_output     ()  { return "-out";            }
    const char* uri_output_dir ()  { return "-out-dir";        }
    const char* uri_output_tmp ()  { return "-out-tmp";        }
    const char* verbose        ()  { return "-verbose";        }
    const char* help           ()  { return "-help";           }
	const char* help_short     ()  { return "-h";              }
    const char* version        ()  { return "-version";        }
    const char* bloom_type     ()  { return "-bloom";          }
    const char* debloom_type   ()  { return "-debloom";        }
    const char* debloom_impl   ()  { return "-debloom-impl";   }
    const char* branching_type ()  { return "-branching-nodes";}
    const char* topology_stats ()  { return "-topology-stats";}
    const char* uri_solid_kmers()  { return "-solid-kmers-out";    }
    const char* bank_convert_type ()  { return "-bank-convert";   }
    const char* integer_precision ()  { return "-integer-precision";}
    const char* solidity_kind  ()  { return "-solidity-kind"; }
	const char* solidity_custom  ()  { return "-solidity-custom"; }
	const char* histo2D  ()  { return "-histo2D"; }
	const char* histo  ()  { return "-histo"; }
    const char* minimizer_type ()  { return "-minimizer-type"; }
    const char* repartition_type() { return "-repartition-type"; }
    const char* compress_level()   { return "-out-compress"; }
    const char* config_only()      { return "-config-only"; }
    const char* storage_type()     { return "-storage-type"; }

    const char* attr_uri_input      ()  { return "input";           }
    const char* attr_kmer_size      ()  { return "kmer_size";       }
    const char* attr_kmer_abundance ()  { return "abundance";       }
    const char* attr_bloom_type     ()  { return "bloom_kind";      }
    const char* attr_debloom_type   ()  { return "debloom_kind";    }
};

/********************************************************************************/

/** Shortcuts. */
#define STR_URI_DB              gatb::core::tools::misc::StringRepository::singleton().db ()
#define STR_URI_FILE            gatb::core::tools::misc::StringRepository::singleton().file ()
#define STR_URI_GRAPH           gatb::core::tools::misc::StringRepository::singleton().graph ()
#define STR_KMER_SIZE           gatb::core::tools::misc::StringRepository::singleton().kmer_size ()
#define STR_MINIMIZER_SIZE      gatb::core::tools::misc::StringRepository::singleton().minimizer_size ()
#define STR_EDGE_KM_REPRESENTATION      gatb::core::tools::misc::StringRepository::singleton().edge_km_representation ()
#define STR_ALL_ABUNDANCE_COUNTS        gatb::core::tools::misc::StringRepository::singleton().all_abundance_counts ()
#define STR_INTEGER_PRECISION   gatb::core::tools::misc::StringRepository::singleton().integer_precision ()
#define STR_KMER_ABUNDANCE      gatb::core::tools::misc::StringRepository::singleton().kmer_abundance ()
#define STR_KMER_ABUNDANCE_MIN  gatb::core::tools::misc::StringRepository::singleton().kmer_abundance_min ()
#define STR_KMER_ABUNDANCE_MIN_THRESHOLD  gatb::core::tools::misc::StringRepository::singleton().kmer_abundance_min_threshold ()
#define STR_KMER_ABUNDANCE_MAX  gatb::core::tools::misc::StringRepository::singleton().kmer_abundance_max ()
#define STR_MAX_MEMORY          gatb::core::tools::misc::StringRepository::singleton().max_memory ()
#define STR_MAX_DISK            gatb::core::tools::misc::StringRepository::singleton().max_disk ()
#define STR_KMER_SOLID          gatb::core::tools::misc::StringRepository::singleton().kmer_solid ()
#define STR_KMER_CFP            gatb::core::tools::misc::StringRepository::singleton().kmer_cFP ()
#define STR_PREFIX              gatb::core::tools::misc::StringRepository::singleton().prefix ()
#define STR_PROGRESS_BAR        gatb::core::tools::misc::StringRepository::singleton().progress_bar ()
#define STR_NB_CORES            gatb::core::tools::misc::StringRepository::singleton().nb_cores ()
#define STR_PARTITION_TYPE      gatb::core::tools::misc::StringRepository::singleton().partition_type ()
#define STR_HISTOGRAM_MAX       gatb::core::tools::misc::StringRepository::singleton().histogram_max ()
#define STR_URI_DEBLOOM         gatb::core::tools::misc::StringRepository::singleton().uri_debloom ()
#define STR_URI_INPUT           gatb::core::tools::misc::StringRepository::singleton().uri_input ()
#define STR_URI_OUTPUT          gatb::core::tools::misc::StringRepository::singleton().uri_output ()
#define STR_URI_OUTPUT_DIR      gatb::core::tools::misc::StringRepository::singleton().uri_output_dir ()
#define STR_URI_OUTPUT_TMP      gatb::core::tools::misc::StringRepository::singleton().uri_output_tmp ()
#define STR_VERBOSE             gatb::core::tools::misc::StringRepository::singleton().verbose ()
#define STR_HELP                gatb::core::tools::misc::StringRepository::singleton().help ()
#define STR_HELP_SHORT          gatb::core::tools::misc::StringRepository::singleton().help_short ()
#define STR_VERSION             gatb::core::tools::misc::StringRepository::singleton().version ()
#define STR_BLOOM_TYPE          gatb::core::tools::misc::StringRepository::singleton().bloom_type()
#define STR_DEBLOOM_TYPE        gatb::core::tools::misc::StringRepository::singleton().debloom_type()
#define STR_DEBLOOM_IMPL        gatb::core::tools::misc::StringRepository::singleton().debloom_impl()
#define STR_BRANCHING_TYPE      gatb::core::tools::misc::StringRepository::singleton().branching_type()
#define STR_TOPOLOGY_STATS      gatb::core::tools::misc::StringRepository::singleton().topology_stats()
#define STR_URI_SOLID_KMERS     gatb::core::tools::misc::StringRepository::singleton().uri_solid_kmers()
#define STR_BANK_CONVERT_TYPE   gatb::core::tools::misc::StringRepository::singleton().bank_convert_type()
#define STR_SOLIDITY_KIND       gatb::core::tools::misc::StringRepository::singleton().solidity_kind()
#define STR_SOLIDITY_CUSTOM     gatb::core::tools::misc::StringRepository::singleton().solidity_custom()
#define STR_HISTO2D             gatb::core::tools::misc::StringRepository::singleton().histo2D()
#define STR_HISTO               gatb::core::tools::misc::StringRepository::singleton().histo()
#define STR_MINIMIZER_TYPE      gatb::core::tools::misc::StringRepository::singleton().minimizer_type()
#define STR_REPARTITION_TYPE    gatb::core::tools::misc::StringRepository::singleton().repartition_type()
#define STR_COMPRESS_LEVEL      gatb::core::tools::misc::StringRepository::singleton().compress_level()
#define STR_CONFIG_ONLY         gatb::core::tools::misc::StringRepository::singleton().config_only()
#define STR_STORAGE_TYPE        gatb::core::tools::misc::StringRepository::singleton().storage_type ()

/********************************************************************************/

#define ATTR_URI_INPUT          gatb::core::tools::misc::StringRepository::singleton().attr_uri_input()
#define ATTR_KMER_SIZE          gatb::core::tools::misc::StringRepository::singleton().attr_kmer_size()
#define ATTR_KMER_ABUNDANCE     gatb::core::tools::misc::StringRepository::singleton().attr_kmer_abundance()
#define ATTR_BLOOM_TYPE         gatb::core::tools::misc::StringRepository::singleton().attr_bloom_type()
#define ATTR_DEBLOOM_TYPE       gatb::core::tools::misc::StringRepository::singleton().attr_debloom_type()

/********************************************************************************/

/** Shortcuts. */
#define STR_BANK_bad_file_number    gatb::core::tools::misc::MessageRepository::singleton().BANK_bad_file_number ()
#define STR_BANK_bad_file_path      gatb::core::tools::misc::MessageRepository::singleton().BANK_bad_file_path ()
#define STR_BANK_unable_open_file   gatb::core::tools::misc::MessageRepository::singleton().BANK_unable_open_file ()
#define STR_BANK_unable_write_file  gatb::core::tools::misc::MessageRepository::singleton().BANK_unable_write_file ()

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MISC_STRINGS_REPOSITORY_HPP_ */
