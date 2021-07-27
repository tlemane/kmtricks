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

#ifndef _CONFIGURATION_ALGORITHM_HPP_
#define _CONFIGURATION_ALGORITHM_HPP_

/********************************************************************************/

#include <gatb/tools/misc/impl/Algorithm.hpp>
#include <gatb/tools/misc/impl/Property.hpp>

#include <gatb/kmer/impl/Model.hpp>
#include <gatb/kmer/impl/Configuration.hpp>

#include <gatb/bank/api/IBank.hpp>

/********************************************************************************/
namespace gatb      {
namespace core      {
namespace kmer      {
namespace impl      {
/********************************************************************************/

template<size_t span>
class ConfigurationAlgorithm : public gatb::core::tools::misc::impl::Algorithm
{
public:

    /** */
    ConfigurationAlgorithm (bank::IBank* bank, tools::misc::IProperties* input);

    /** */
    ~ConfigurationAlgorithm ();

    /** */
    void execute ();

    /** */
    const Configuration&  getConfiguration() const { return _config; }

private:
    /** */
    static std::vector<tools::misc::CountRange> getSolidityThresholds (tools::misc::IProperties* params);

	static std::vector<bool> getSolidityCustomVector (tools::misc::IProperties* params);

    /** Shortcut. */
    typedef typename Kmer<span>::Type Type;

    Configuration _config;

    bank::IBank* _bank;
    void setBank (bank::IBank* bank) { SP_SETATTR(bank); }

    tools::misc::IProperties* _input;
    void setInput (tools::misc::IProperties* input)  { SP_SETATTR(input); }
};

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _CONFIGURATION_ALGORITHM_HPP_ */
