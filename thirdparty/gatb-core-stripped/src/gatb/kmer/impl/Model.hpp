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

/** \file Model.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Kmer management
 */

#ifndef _GATB_CORE_KMER_IMPL_MODEL_HPP_
#define _GATB_CORE_KMER_IMPL_MODEL_HPP_

/********************************************************************************/

#include <gatb/system/api/Exception.hpp>

#include <gatb/kmer/api/IModel.hpp>
#include <gatb/tools/collections/api/Bag.hpp>

#include <gatb/tools/designpattern/api/Iterator.hpp>
#include <gatb/tools/designpattern/impl/IteratorHelpers.hpp>
#include <gatb/tools/misc/api/Data.hpp>
#include <gatb/tools/misc/api/Abundance.hpp>

#include <gatb/tools/math/Integer.hpp>

#include <gatb/tools/storage/impl/Storage.hpp>

#include <vector>
#include <algorithm>
#include <iostream>
#include <bitset>

extern const char bin2NT[] ;
extern const char binrev[] ;
extern const unsigned char revcomp_4NT[];
extern const unsigned char comp_NT    [];

/********************************************************************************/
namespace gatb      {
namespace core      {
/** \brief Package for genomic databases management. */
namespace kmer      {
/** \brief Implementation for genomic databases management. */
namespace impl      {
/********************************************************************************/

/** We get the nth defined kmer size. */
#define KMER_SPAN(n)  (boost::mpl::at<gatb::core::tools::math::IntegerList, boost::mpl::int_<n> >::type::value)

/** We get the first value in the kmer size list. */
#define KMER_DEFAULT_SPAN  KMER_SPAN(0)

/********************************************************************************/

/** \brief Entry point for kmer management.
 *
 * This structure is only a container for other types defined inside. The specificity is
 * that this structure is templated by a 'span' integer that represents the maximal kmer
 * size supported (actually, the max value is 'span-1').
 *
 * Inside this structure, we have the following main elements:
 *      - 'Type'  : this is the integer type representing kmer values
 *      - 'Model' : provides many services for managing kmers
 *      - 'Count' : roughly speaking, this a kmer value with an associated abundance
 *
 * This structure must be used only with for 4 values (32,64,96,128 for instance, see Model.cpp), otherwise
 * a compilation error occurs (more values could be added in the future).
 *
 * A default value of 32 is defined for the template parameter, so writing 'Kmer<>::Model'
 * represents a model that supports kmers of size up to 31 (included).
 */
template <size_t span=KMER_DEFAULT_SPAN>
struct Kmer
{
    /************************************************************/
    /***********************     TYPE     ***********************/
    /************************************************************/

    /** Alias type for the integer value of a kmer. We use the LargeInt class for supporting big integers.
     * Note that the template parameter 'span' represents the maximal kmer size supported by the Kmer class.
     * A conversion to the template parameter of LargeInt is done.
     */
    typedef tools::math::LargeInt<(span+31)/32> Type;


    /************************************************************/
    /***********************     MODEL    ***********************/
    /************************************************************/

    /** Forward declarations. */
    class ModelDirect;
    class ModelCanonical;
    template<class Model, class Comparator> class ModelMinimizer;

    /** Now, we need to define what is a kmer for each kind of model.
     *
     * The simple case is KmerDirect, where only the value of the kmer is available
     * as a method 'value' returning a Type object.
     *
     * The second case is KmerCanonical, which is the same as KmerDirect, but with
     * two other methods 'forward' and 'revcomp'
     *
     * The third case is KmerMinimizer<Model> which allows to handle minimizers associated
     * to a kmer. This class inherits from the Model::Kmer type and adds methods specific
     * to minimizers, such as 'minimizer' itself (ie the Model::Kmer object holding the
     * minimizer), 'position' giving the position of the minimizer whithin the kmer and
     * 'hasChanged' telling whether a minimizer has changed during iteration of kmers from
     * some data source (a sequence data for instance).
     */

    /** \brief Kmer type for the ModelDirect class.
     *
     * This class represent direct kmers, ie. a mere Type value.
     *
     *  NOTE: this class is not intended to be used directly by end users. Instead, the typedef definition
     *  \ref ModelDirect::Kmer should be preferred.
     *
     * Example of use:
     * \snippet kmer2.cpp  snippet1_direct
     */
    class KmerDirect
    {
    public:
        /** Returns the value of the kmer.
         * \return the kmer value as a Type object. */
        const Type& value  () const { return _value;   }

        /** This is a dummy function that always returns the value of the kmer in the forward direction, even when "which" is 1.
         * it's provided for API compatibility with KmerCanonical
         * We could compute revcomp(_value) when which=1, but KmerDirect doesn't know about its k-mer size.
         * \param[in] which: dummy parameter
         * \return the kmer value as a Type object. */
        const Type& value  (int which) const { if (which==1){ std::cout << "unsupported call to value(which) for KmerDirect" << std::endl; exit(1); } 
                                               return _value;   }

        /*_ Comparison operator between two instances.
         * \param[in] t : object to be compared to
         * \return true if the values are the same, false otherwise. */
        bool operator< (const KmerDirect& t) const  { return this->_value < t._value; };

        /** Set the value of the kmer
         * \param[in] val : value to be set. */
        void set (const Type& val) { _value=val; }

        /** Tells whether the kmer is valid or not. It may be invalid if some unwanted
         * nucleotides characters (like N) have been used to build it.
         * \return true if valid, false otherwise. */
        bool isValid () const { return _isValid; }

        /* compatibility with KmerCanonical API */
        bool which () const { return true; }

        /* compatibility with KmerCanonical API */
        Strand strand() const { return STRAND_FORWARD;  }

        /* compatibility with KmerCanonical API */
        const Type& forward() const { return value(); }

        /** Returns the reverse complement value of this canonical kmer.
         * \return the reverse complement value */
//        const Type& revcomp() const { return  }


    protected:
        Type _value;
        bool _isValid;
        friend class ModelDirect;

        /** Extract a mmer from a kmer. This is done by using a mask on the kmer.
         * \param[in] mask : mask to be applied to the current kmer
         * \param[in] size : shift size (needed for some kmer classes but not all)
         * \param[in] mmer_lut : lookup table of minimizers
         * \return the extracted kmer.
         */

        KmerDirect extract      (const Type& mask, size_t size, Type * mmer_lut)  {  KmerDirect output;  output.set (mmer_lut[(this->value() & mask).getVal()]);  return output;  }
        KmerDirect extractShift (const Type& mask, size_t size, Type * mmer_lut)  {  KmerDirect output = extract(mask,size,mmer_lut);  _value = _value >> 2;  return output;  }
    };

    /** \brief Kmer type for the ModelCanonical class.
     *
     * This class represent canonical kmers, ie. a value that is the minimal value of the forward kmer
     * and its reverse complement.
     *
     * The implementation maintains a table of two Type objects, the first one for the forward kmer and
     * the second one for the reverse complement.
     *
     * We can know which object is the canonical one (ie. the minimum) by using the
     *  method \ref value.*  One can also retrieve the strand used for the canonical form with the method
     *  \ref strand.
     *
     * It is still possible to get the forward kmer with \ref forward and the reverse complement with \ref
     *  revcomp.
     *
     *  NOTE: this class is not intended to be used directly by end users. Instead, the typedef definition
     *  \ref ModelCanonical::Kmer should be preferred.
     *
     * Example of use:
     * \snippet kmer2.cpp  snippet1_canonical
     */
    class KmerCanonical
    {
    public:

        /** Returns the value of the kmer.
         * \return the kmer value as a Type object. */
        const Type& value  () const { return table[(int)choice];   }

        /** Returns the value of the kmer.
         * \param[in] which: forward or reverse strand
         * \return the kmer value as a Type object. */
        const Type& value  (int which) const { return table[which];   }

        /** Comparison operator between two instances.
         * \param[in] t : object to be compared to
         * \return true if the values are the same, false otherwise. */
        bool operator< (const KmerDirect& t) const  { return this->value() < t.value(); };

        /** Set the value of the kmer. IMPORTANT: Not really a forward/revcomp couple,
         * but may be useful for the minimizer default value.
         * \param[in] val : value to be set (set to both forward and reverse complement). */
        void set (const Type& val)
        {
            table[0]=val;
            table[1]=val;
            choice = 0;
        }

        void set (const u_int64_t& val)
        {
            table[0].setVal(val);
            table[1].setVal(val);
            choice = 0;
        }


        /** Set the forward/revcomp attributes. The canonical form is computed here.
         * \param[in] forward : forward value
         * \param[in] revcomp : reverse complement value.
         */
        void set (const Type& forward, const Type& revcomp)
        {
            table[0]=forward;
            table[1]=revcomp;
            updateChoice ();
        }
		
        /** Tells whether the kmer is valid or not. It may be invalid if some unwanted
         * nucleotides characters (like N) have been used to build it.
         * \return true if valid, false otherwise. */
        bool isValid () const { return _isValid; }

        /** Returns the forward value of this canonical kmer.
         * \return the forward value */
        const Type& forward() const { return table[0]; }

        /** Returns the reverse complement value of this canonical kmer.
         * \return the reverse complement value */
        const Type& revcomp() const { return table[1]; }

        /** Tells which strand is used for the kmer.
         * \return true if the kmer value is the forward value, false if it is the reverse complement value
         */
        bool which () const { return choice==0 ? true : false; }

        /** Tells which strand is used.
         * \return the used strand. */
        Strand strand() const { return which() ? STRAND_FORWARD : STRAND_REVCOMP; }

        /* tells whether a kmer and its revcomp are identical */
        bool isPalindrome () const { return table[0] == table[1]; }

    protected:
        Type table[2];  char choice;
		
        bool _isValid;
        void updateChoice () { choice = (table[0] < table[1]) ? 0 : 1; }
        friend class ModelCanonical;

        /** Extract a mmer from a kmer. This is done by using a mask on the kmer.
         * \param[in] mask : mask to be applied to the current kmer
         * \param[in] size : shift size (needed for some kmer classes but not all)
         * \param[in] mmer_lut : lookup table of minimizers
         * \return the extracted kmer.
         */
        KmerCanonical extract (const Type& mask, size_t size, Type* mmer_lut)
        {

            KmerCanonical output;
			
			output.set(mmer_lut[(this->table[0] & mask).getVal()]); //no need to recomp updateChoice with this
			//mmer_lut takes care of revcomp and forbidden mmers
			//output.set (this->table[0] & mask, (this->table[1] >> size) & mask);
            //output.updateChoice();
            return output;
        }


        KmerCanonical extractShift (const Type& mask, size_t size, Type * mmer_lut)
        {
            KmerCanonical output = extract (mask, size,mmer_lut);
            table[0] = table[0] >> 2;   table[1] = table[1] << 2;  updateChoice();
            return output;
        }
       
     };

    /** \brief Kmer type for the ModelMinimizer class.
     *
     * This class associates a kmer and its minimizer. It inherits from the Model::Kmer type
     * and adds methods specific to minimizers, such as \ref minimizer itself
     * (ie the Model::Kmer object holding the minimizer),
     * \ref position giving the position of the minimizer whithin the kmer and
     * \ref hasChanged telling whether a minimizer has changed during iteration of kmers from
     * some data source (a sequence data for instance).
     *
     * NOTE: this class is not intended to be used directly by end users. Instead, the typedef definition
     * \ref ModelMinimizer::Kmer should be preferred.
     *
     * Example of use:
     * \snippet kmer2.cpp  snippet1_minimizer
     */
    template<class Model, class Comparator>
    class KmerMinimizer : public Model::Kmer
    {
    public:

        /** Returns the minimizer of the current kmer as a Model::Kmer object
         * \return the Model::Kmer instance */
        const typename Model::Kmer& minimizer() const  {  return _minimizer; }

        /** Returns the position of the minimizer within the kmer. By convention,
         * a negative value means that there is no minimizer inside the kmer.
         * \return the position of the minimizer. */
        int position () const  {  return _position;  }

        /** Tells whether the minimizer has changed; useful while iterating kmers
         * \return true if changed, false otherwise */
        bool hasChanged () const  {  return _changed;  }

    protected:

        typename Model::Kmer _minimizer;
        int16_t              _position;
        bool                 _changed;
        friend class ModelMinimizer<Model,Comparator>;
    };

    /** Abstract class that provides kmer management.
     *
     * This class is the base class for kmer management. It provides several services on this purpose
     * like getting kmer information from some nucleotides sequence, or iterate kmers through such
     * a sequence.
     *
     * This class has two templates types :
     *
     *      1) ModelImpl : ModelAbstract is design for static polymorphism and ModelImpl is the implementation
     *                     that must be provided to it
     *
     *      2) T : type of kmers handled by the class (ie KmerDirect, KmerCanonical...); I was not successful
     *             in trying to hide KmerXXX classes in the dedicated ModelXXX classes because of mutual
     *             dependencies while template specializations (maybe a solution one day)
     *
     * End user will be given instances of Kmer class, delivering more or less information according to the
     * specific type of ModelImpl
     */
    template <class ModelImpl, typename T>
    class ModelAbstract : public system::SmartPointer
    {
    public:

        /** Type of kmers provided by the class. It can be KmerDirect, KmerCanonical or KmerMinimizer.
         *
         * The simple way to get the value of the kmer is done with the 'value' method.
         *
         * Note that, according to the true type of T, this Kmer typedef may have more or less methods. */
        typedef T Kmer;

        /** (default) Constructor. The provided (runtime) kmer size must be coherent with the span (static) value.
         * \param[in] sizeKmer : size of kmers handled by the instance.*/
        ModelAbstract (size_t sizeKmer=span-1) : _kmerSize(sizeKmer)
        {
            /** We check that the Type precision is enough for the required kmers span. */
            if (sizeKmer >= span)
            {
                throw system::Exception ("Type '%s' has too low precision (%d bits) for the required %d kmer size",
                    Type().getName(), Type().getSize(), sizeKmer
                );
            }

            /** We compute the mask of the kmer. Useful for computing kmers in a recursive way. */
            Type un;
            un.setVal(1);
            _kmerMask = (un << (_kmerSize*2)) - un;

            size_t shift = 2*(_kmerSize-1);

            /** The _revcompTable is a shortcut used while computing revcomp recursively. */
            /** Important: don't forget the Type cast, otherwise the result in only on 32 bits. */
            for (size_t i=0; i<4; i++)   {  Type tmp; tmp.setVal(comp_NT[i]);  _revcompTable[i] = tmp << shift;  }
        }

        /** Returns the span of the model
         * \return the model span. */
        size_t getSpan () const { return span; }

        /** Get the memory size (in bytes) of a Kmer<span>::Type object.
         * \return the memory size of a kmer. */
        size_t getMemorySize ()  const  { return sizeof (Type); }

        /** Gives the kmer size for this model.
         * \return the kmer size. */
        size_t getKmerSize () const { return _kmerSize; }

        /** Gives the maximum value of a kmer for the instance.
         * \return the maximum kmer value. */
        const Type& getKmerMax () const { return _kmerMask; }

        /** Returns an ascii representation of the kmer value.
         * \param[in] kmer : the kmer we want an ascii representation for
         * \return a string instance holding the ascii representation. */
        std::string toString (const Type& kmer) const  {  return kmer.toString(_kmerSize);  }
        std::string toString (u_int64_t kmer) const  {  tools::math::LargeInt<1> km; km.setVal(kmer); return km.toString(_kmerSize);  }

        /** Compute the reverse complement of a kmer.
         * \param[in] kmer : the kmer to be reverse-completed.
         * \return the reverse complement. */
        Type reverse (const Type& kmer)  const  { return revcomp (kmer, this->_kmerSize); }

        /** Build a kmer from a Data object (ie a sequence of nucleotides), starting at an index in the nucleotides sequence.
         * The result is a pair holding the built kmer and a boolean set to yes if the built kmer has to be understood in
         * the forward sense, false otherwise.
         * \param[in] data : the data from which we extract a kmer
         * \param[in] startIndex : start index in the data object (default to 0)
         * \return a pair with the built kmer and a boolean set to yes if the kmer is understood in the forward strand
         */
        Kmer getKmer (const tools::misc::Data& data, size_t startIndex=0)  const
        {
            return codeSeed (data.getBuffer(), data.getEncoding(), startIndex);
        }

        /** Iteration of the kmers from a data object through a functor (so lambda expressions can be used).
         * Note that the functor takes the currently iterated Kmer object and its index during the iteration.
         *
         *  Example of use:
         * \snippet kmer8.cpp  snippet1_iterate
         * \param[in] data : the sequence of nucleotides as a Data object.
         * \param[in] callback  : functor that handles one kmer */
        template<typename Callback>
        bool iterate (tools::misc::Data& data, Callback callback) const
        {
            return execute <Functor_iterate<Callback> > (data.getEncoding(), Functor_iterate<Callback>(data,callback));
        }

        /** Compute the kmer given some nucleotide data.
         *  Note that we don't check if we have enough nucleotides in the provided data.
         * \param[in] seq : the sequence
         * \param[in] encoding : encoding mode of the sequence
         * \param[in] startIndex : index of the first nucleotide in the seq buffer
         * \return the kmer for the given nucleotides. */
        Kmer codeSeed (const char* seq, tools::misc::Data::Encoding_e encoding, size_t startIndex=0) const
        {
            return execute<Functor_codeSeed> (encoding, Functor_codeSeed(seq, startIndex));
        }

        /** Compute the next right kmer given a current kmer and a nucleotide.
         * \param[in] kmer : the current kmer as a starting point
         * \param[in] nucl : the next nucleotide
         * \param[in] encoding : encoding mode of the sequence
         * \return the kmer on the right of the given kmer. */
        Kmer codeSeedRight (const Kmer& kmer, char nucl, tools::misc::Data::Encoding_e encoding)  const
        {
            return execute<Functor_codeSeedRight> (encoding, Functor_codeSeedRight(kmer,nucl));
        }

        /** Build a vector of successive kmers from a given sequence of nucleotides provided as a Data object.
         * \param[in] data : the sequence of nucleotides.
         * \param[out] kmersBuffer : the successive kmers built from the data object.
         * \return true if kmers have been extracted, false otherwise. */
		//GR : est ce quon pourrait passer un pointeur (de taille suffisante) au lieu  vector, pour pas avoir a faire resize dessus dans tas
		// si taille pas suffisante,
        bool build (tools::misc::Data& data, std::vector<Kmer>& kmersBuffer)  const
        {
            /** We compute the number of kmers for the provided data. Note that we have to check that we have
             * enough nucleotides according to the current kmer size. */
            int32_t nbKmers = data.size() - this->getKmerSize() + 1;
            if (nbKmers <= 0)  { return false; }

            /** We resize the resulting kmers vector. */
            kmersBuffer.resize (nbKmers);

            /** We fill the vector through a functor. */
            this->iterate (data, BuildFunctor<Kmer>(kmersBuffer));

            return true;
        }

        /** Iterate the neighbors of a given kmer; these neighbors are:
         *  - 4 outgoing neighbors (with nt A,C,T,G)
         *  - 4 incoming neighbors (with nt A,C,T,G)
         *  This method uses a functor that will be called for each possible neighbor of the source kmer.
         *  \param[in] source : the kmer from which we want neighbors.
         *  \param[in] fct : a functor called for each neighbor.
         *  \param[in] mask : holds 8 bits for each possible neighbor (1 means that the neighbor is computed)
         */
        template<typename Functor>
        void iterateNeighbors (const Type& source, const Functor& fct, const std::bitset<8>& mask = 0xFF)  const
        {
            // hacky to cast Functor& instead of const Functor&, but don't wanna break API yet want non-const functor
            iterateOutgoingNeighbors(source, (Functor&) fct, std::bitset<4> ( (mask.to_ulong() >> 0) & 15));
            iterateIncomingNeighbors(source, (Functor&) fct, std::bitset<4> ( (mask.to_ulong() >> 4) & 15));
        }

        /** Iterate the neighbors of a given kmer; these neighbors are:
         *  - 4 outcoming neighbors
         *  This method uses a functor that will be called for each possible neighbor of the source kmer.
         *  \param[in] source : the kmer from which we want neighbors.
         *  \param[in] fct : a functor called for each neighbor.
         *  \param[in] mask : mask of the neighbors to be used
         */
        template<typename Functor>
        void iterateOutgoingNeighbors (const Type& source, Functor& fct, const std::bitset<4>& mask = 0x0F)  const
        {
            /** We compute the 4 possible neighbors. */
            for (size_t nt=0; nt<4; nt++)
            {
                if (mask[nt] == true)
                {
                    Type next1 = (((source) * 4 )  + nt) & getKmerMax();
                    Type next2 = revcomp (next1, getKmerSize());
                    fct (std::min (next1, next2));
                }
            }
        }

        /** Iterate the neighbors of a given kmer; these neighbors are:
         *  - 4 incoming neighbors
         *  This method uses a functor that will be called for each possible neighbor of the source kmer.
         *  \param[in] source : the kmer from which we want neighbors.
         *  \param[in] fct : a functor called for each neighbor.
         *  \param[in] mask : mask of the neighbors to be used
         */
        template<typename Functor>
        void iterateIncomingNeighbors (const Type& source, Functor& fct, const std::bitset<4>& mask = 0x0F)  const
        {
            Type rev = core::tools::math::revcomp (source, getKmerSize());

            /** We compute the 4 possible neighbors. */
            for (size_t nt=0; nt<4; nt++)
            {
                /** Here, we use the complement of the current nucleotide 'nt', the idea is to have the same
                 * nucleotide iteration than the iterateOutgoingNeighbors method.
                 * Remember : A=0, C=1, T=2, G=3  (each coded on 2 bits)
                 * => we can get the complement by negating the most significant bit (ie "nt^2") */

                if (mask[nt] == true)
                {
                    Type next1 = (((rev) * 4 )  + (nt^2)) & getKmerMax();
                    Type next2 = revcomp (next1, getKmerSize());
                    fct (std::min (next1, next2));
                }
            }
        }

        /************************************************************/
        /* \brief Iterator on successive kmers
         *
         * This class will iterate successive kmers extracted from a Data object.
         * It is similar to the Model::build, except that here we don't have a container
         * holding all the successive kmers (ie. we have here only sequential access and
         * not direct access).
         *
         * To be used, such an iterator must be initialized with some sequence of nucleotides,
         * which is done with the 'setData' method.
         */
        class Iterator : public tools::dp::impl::VectorIterator<Kmer>
        {
        public:
            /** Constructor.
             * \param[in] ref : the associated model instance.
             */
            Iterator (ModelAbstract& ref)  : _ref(ref)   {}

            /** Set the data to be iterated.
             * \param[in] d : the data as information source for the iterator
             */
            void setData (tools::misc::Data& d) // TODO: should this be const? I feel like it should
            {
                /** We fill the vector with the items to be iterated. */
                _ref.build (d, this->_items);

                /** We set the vector size. */
                this->_nb = this->_items.size();
            }

        private:
            /** Reference on the underlying model; called for its 'build' method. */
            ModelAbstract& _ref;
        };

    protected:

        /* Shortcuts. */
        typedef tools::misc::Data::ConvertChar      ConvertChar;
        typedef tools::misc::Data::ConvertASCII     ConvertASCII;
        typedef tools::misc::Data::ConvertInteger   ConvertInteger;
        typedef tools::misc::Data::ConvertBinary    ConvertBinary;

        /* Size of a kmer for this model. */
        size_t  _kmerSize;

        /* Mask for the kmer. Used for computing recursively kmers. */
        Type  _kmerMask;

        /* Shortcut for easing/speeding up the recursive revcomp computation. */
        Type _revcompTable[4];

        /* \return -1 if valid, otherwise index of the last found bad character. */
        template<class Convert>
        int polynom (const char* seq, Type& kmer, size_t startIndex)  const
        {
            ConvertChar c;
            int badIndex = -1;

            /** We iterate 'kmersize" nucleotide to build the first kmer as a polynomial evaluation. */
            kmer.setVal(0);
            for (size_t i=0; i<_kmerSize; ++i)
            {
                /** We get the current nucleotide (and its invalid status). */
                c = Convert::get(seq,i+startIndex);

                /** We update the polynome value. */
                kmer = (kmer<<2) + c.first;

                /** We update the 'invalid' status: a single bad character makes the result invalid. */
                if (c.second)  { badIndex = i; }
            }

            return badIndex;
        }

        /* Generic function that switches to the correct implementation according to the encoding scheme
         * of the provided Data parameter; the provided functor class is specialized with the correct data conversion type
         * and the called.
         */
        template<class Functor>
        typename Functor::Result execute (tools::misc::Data::Encoding_e encoding, Functor action) const
        {
            switch (encoding)
            {
                case  tools::misc::Data::ASCII:    return action.template operator()<ConvertASCII  > (this);
                case  tools::misc::Data::INTEGER:  return action.template operator()<ConvertInteger> (this);
                case  tools::misc::Data::BINARY:   return action.template operator()<ConvertBinary>  (this);
                default:  throw system::Exception ("BAD FORMAT IN 'execute'");
            }
        }

        /* Adaptor between the 'execute' method and the 'codeSeed' method. */
        struct Functor_codeSeed
        {
            typedef typename ModelImpl::Kmer Result;
            const char* buffer;
            size_t startIndex;
            Functor_codeSeed (const char* buffer, size_t startIndex) : buffer(buffer), startIndex(startIndex) {}
            template<class Convert>  Result operator() (const ModelAbstract* model)
            {
                Result result;
                static_cast<const ModelImpl*>(model)->template first <Convert> (buffer, result, startIndex);
                return result;
            }
        };

        /* Adaptor between the 'execute' method and the 'codeSeedRight' method. */
        struct Functor_codeSeedRight
        {
            typedef typename ModelImpl::Kmer Result;
            const Kmer& kmer; char nucl;
            Functor_codeSeedRight (const Kmer& kmer, char nucl) : kmer(kmer), nucl(nucl) {}
            template<class Convert>  Result operator() (const ModelAbstract* model)
            {
                ConvertChar c = Convert::get(&nucl,0);
                Result result=kmer;
                static_cast<const ModelImpl*>(model)->template next <Convert> (c.first, result, c.second==0);
                return result;
            }
        };

        /* Adaptor between the 'execute' method and the 'iterate' method. */
        template<class Callback>
        struct Functor_iterate
        {
            typedef bool Result;
            tools::misc::Data& data; Callback callback;
            Functor_iterate (tools::misc::Data& data, Callback callback) : data(data), callback(callback) {}
            template<class Convert>  Result operator() (const ModelAbstract* model)
            {
                return static_cast<const ModelImpl*>(model)->template iterate<Callback, Convert> (data.getBuffer(), data.size(), callback);
            }
        };

        /** Template method that iterates the kmer of a given sequence (provided as a buffer and its length).
         *  Note : we use static polymorphism here (http://en.wikipedia.org/wiki/Template_metaprogramming)
         *  \param[in] seq : the sequence to be iterated
         *  \param[in] length : length of the sequence
         *  \param[in] callback : functor called on each found kmer in the sequence
         *  \return true if kmers have been found, false otherwise.
         */
        template<typename Callback, typename Convert>
        bool iterate (const char* seq, size_t length, Callback callback) const
        {
            /** We compute the number of kmers for the provided data. Note that we have to check that we have
             * enough nucleotides according to the current kmer size. */
            int32_t nbKmers = length - _kmerSize + 1;
            if (nbKmers <= 0)  { return false; }

            /** We create a result instance. */
            typename ModelImpl::Kmer result;

            /** We compute the initial seed from the provided buffer. */
            int indexBadChar = static_cast<const ModelImpl*>(this)->template first<Convert> (seq, result, 0);

            /** We need to keep track of the computed kmers. */
            size_t idxComputed = 0;

            /** We notify the result. */
            this->notification<Callback> (result, idxComputed, callback);

            /** We compute the following kmers from the first one.
             * We have consumed 'kmerSize' nucleotides so far for computing the first kmer,
             * so we start the loop with idx=_kmerSize.
             */
            for (size_t idx=_kmerSize; idx<length; idx++)
            {
                /** We get the current nucleotide. */
                ConvertChar c = Convert::get (seq, idx);

                if (c.second)  { indexBadChar = _kmerSize-1; }
                else           { indexBadChar--;     }

                /** We compute the next kmer from the previous one. */
                static_cast<const ModelImpl*>(this)->template next<Convert> (c.first, result, indexBadChar<0);

                /** We notify the result. */
                this->notification<Callback> (result, ++idxComputed, callback);
            }

            return true;
        }

        template <class Callcack>
        void  notification (const Kmer& value, size_t idx, Callcack callback) const {  callback (value, idx);  }

        /** */
        template<typename Type>
        struct BuildFunctor
        {
            std::vector<Type>& kmersBuffer;
            BuildFunctor (std::vector<Type>& kmersBuffer) : kmersBuffer(kmersBuffer) {}
            void operator() (const Type& kmer, size_t idx)  {  kmersBuffer[idx] = kmer;  }
        };
		
    };

    /********************************************************************************/

    /** \brief Model that handles "direct" kmers, ie sequences of nucleotides.
     * The associated value of such a kmer is computed as a polynom P(X) with X=4
     * and where the coefficients are in [0..3].
     * By convention, we use A=0, C=1, T=2 and G=3
     *
     * Example of use:
     * \snippet kmer2.cpp  snippet1_direct
     */
    class ModelDirect :  public ModelAbstract<ModelDirect, Kmer<span>::KmerDirect>
    {
    public:

        /** Kmer type for this kind of model.  */
        typedef Kmer<span>::KmerDirect Kmer;

        /** Constructor.
         * \param[in] kmerSize : size of the kmers handled by the model. */
        ModelDirect (size_t kmerSize=span-1) : ModelAbstract<ModelDirect, Kmer> (kmerSize) {}

        /** Computes a kmer from a buffer holding nucleotides encoded in some format.
         * The way to interpret the buffer is done through the provided Convert template class.
         * \param[in] buffer : holds the nucleotides sequence from which the kmer has to be computed
         * \param[out] value : kmer as a result
         * \param[in] startIndex : index of the first nucleotide of the kmer to retrieve in the buffer
         */
        template <class Convert>
        int first (const char* buffer, Kmer& value, size_t startIndex)   const
        {
           int result = this->template polynom<Convert> (buffer, value._value, startIndex);
            value._isValid = result < 0;
            return result;
        }

        /** Computes a kmer in a recursive way, ie. from a kmer and the next
         * nucleotide. The next nucleotide is computed from a buffer and the
         * index of the nucleotide within the buffer.
         * The way to interpret the buffer is done through the provided Convert template class.
         * \param[in] c : next nucleotide
         * \param[out] value : kmer to be updated with the provided next nucleotide
         * \param[in] isValid : tells whether the updated kmer is valid or not
         */
        template <class Convert>
        void  next (char c, Kmer& value, bool isValid)   const
        {
            value._value   = ( (value._value << 2) +  c) & this->_kmerMask;
            value._isValid = isValid;
        }
	};

    /********************************************************************************/

    /** \brief Model that handles "canonical" kmers, ie the minimum value of the
     * direct kmer and its reverse complement.
     *
     * Example of use:
     * \snippet kmer2.cpp  snippet1_canonical
     */
    class ModelCanonical :  public ModelAbstract<ModelCanonical, Kmer<span>::KmerCanonical>
    {
    public:

        /** Kmer type for this kind of model.  */
        typedef Kmer<span>::KmerCanonical Kmer;

        /** Constructor.
         * \param[in] kmerSize : size of the kmers handled by the model. */
        ModelCanonical (size_t kmerSize=span-1) : ModelAbstract<ModelCanonical, Kmer> (kmerSize) {}

        /** Computes a kmer from a buffer holding nucleotides encoded in some format.
         * The way to interpret the buffer is done through the provided Convert template class.
         * \param[in] seq : holds the nucleotides sequence from which the kmer has to be computed
         * \param[out] value : kmer as a result
         * \param[in] startIndex : index of the first nucleotide of the kmer to retrieve in the buffer
         */
        template <class Convert>
        int first (const char* seq, Kmer& value, size_t startIndex)   const
        {

            int result = this->template polynom<Convert> (seq, value.table[0], startIndex);
            value._isValid = result < 0;
            value.table[1] = this->reverse (value.table[0]);
            value.updateChoice();
            return result;
        }

        /** Computes a kmer in a recursive way, ie. from a kmer and the next
         * nucleotide. The next nucleotide is computed from a buffer and the
         * index of the nucleotide within the buffer.
         * The way to interpret the buffer is done through the provided Convert template class.
         * \param[in] c : next nucleotide
         * \param[out] value : kmer to be updated with the provided next nucleotide
         * \param[in] isValid : tells whether the updated kmer is valid or not
         */
        template <class Convert>
        void  next (char c, Kmer& value, bool isValid)   const
        {
            value.table[0] = ( (value.table[0] << 2) +  c                          ) & this->_kmerMask;
            value.table[1] = ( (value.table[1] >> 2) +  this->_revcompTable[(int)c]) & this->_kmerMask;
            value._isValid = isValid;

            value.updateChoice();
        }

        /* compute kmer hash */
        uint64_t getHash(const Type &k) const
        {
            return hash1(k, 0);
        }

        /* for profiling only, compute kmer hash using hash2 */
        void getHash2(const Type &k) const
        {
            hash2(k, 1LL);
        }

    };

    /********************************************************************************/

    struct ComparatorMinimizer
    {
        template<class Model>  void init (const Model& model, Type& best) const { best = model.getKmerMax(); }
        bool operator() (const Type& current, const Type& best) const { return current < best; }
    };

    /* compare the minimizers by frequency, if information is available, else lexicographical */
    /* maybe this code can be factorized with ComparatorMinimizer or also one can say that it subsumes it 
     * (at the cost of accessing has_frequency for each comparison) */
    struct ComparatorMinimizerFrequencyOrLex
    {
        template<class Model>  void init (const Model& model, Type& best) 
        {   
            best = model.getKmerMax(); 
            has_frequency = false;
        }

        void include_frequency (uint32_t *freq_order)
        { 
            _freq_order = freq_order;
            has_frequency = true;
        }

		template<class Model>   Type computeLargest (const Model& model,int mmersize)
		{
			Type largest;
			if(has_frequency)
			{
				u_int64_t nbminims_total = ((u_int64_t)1 << (2*mmersize));

				Type  mmer_max;
				mmer_max.setVal(0);
				uint32_t _freq_max = _freq_order[mmer_max.getVal()];
				for(uint32_t ii=0; ii< nbminims_total; ii++)
				{
					Type Tii;
					Tii.setVal(ii);
					
					if( ! (*this)(Tii,mmer_max ))
					{
						mmer_max.setVal(ii);
						_freq_max = _freq_order[ii];
					}
				}
				printf("largest freq is %i   for %s\n",_freq_max,mmer_max.toString(mmersize).c_str());
				largest = mmer_max;
			}
			else
			{
				largest = model.getKmerMax();
			}
			
			return largest;
		}
		
        bool operator() (const Type& a_t, const Type& b_t) const {
            u_int64_t a = a_t.getVal();
            u_int64_t b = b_t.getVal();

            if (has_frequency)
            {
                //printf("testing freq order of %d %d: %d %d, min is gonna be: %d\n",a,b,_freq_order[a], _freq_order[b], (_freq_order[a] < _freq_order[b]) ? a : b);
				//printf("freq order %llu %llu        %i %i  %i\n",a,b,_freq_order[a],_freq_order[b], _freq_order[a] > _freq_order[b]);
                if (_freq_order[a] == _freq_order[b])
                    return a < b;
                return _freq_order[a] < _freq_order[b];
            }
            else
            {
                return a < b; 
            }
        }

        private:
        uint32_t* _freq_order;
        bool has_frequency;
    };


    /** \brief Model that handles kmers of the Model type + a minimizer
     *
     * This model supports the concept of minimizer. It acts as a Model instance (given as a
     * template class) and add minimizer information to the Kmer type.
     *
     * Example of use:
     * \snippet kmer2.cpp  snippet1_minimizer
     */
    template<class ModelType, class Comparator=Kmer<span>::ComparatorMinimizerFrequencyOrLex> 
    class ModelMinimizer :  public ModelAbstract <ModelMinimizer<ModelType,Comparator>, KmerMinimizer<ModelType,Comparator> >
    {
    public:

        /** Type of the model for kmer and mmers.  */
        typedef ModelType Model;

        /** Kmer type for this kind of model.  */
        typedef KmerMinimizer<ModelType,Comparator> Kmer;

        /** Return a reference on the model used for managing mmers.
         * \return the minimizer model. */
        const ModelType& getMmersModel() const { return _miniModel; }

        /** Constructor.
         * \param[in] kmerSize      : size of the kmers handled by the model.
         * \param[in] minimizerSize : size of the mmers handled by the model.
         * \param[in] cmp : functor that compares two minizers
         * \param[in] freq_order : a 4^m table containing the frequency of each minimizer 
         */
        ModelMinimizer (size_t kmerSize, size_t minimizerSize, Comparator cmp=Comparator(), uint32_t *freq_order=NULL)
            : ModelAbstract <ModelMinimizer<ModelType,Comparator>, Kmer > (kmerSize),
              _kmerModel(kmerSize), _miniModel(minimizerSize), _cmp(cmp), _freq_order(freq_order)
        {
            if (kmerSize < minimizerSize)  { throw system::Exception ("Bad values for kmer %d and minimizer %d", kmerSize, minimizerSize); }

            _minimizerSize = minimizerSize;
			
            /** We compute the number of mmers found in a kmer. */
            _nbMinimizers = _kmerModel.getKmerSize() - minimizerSize + 1;

            /** We need a mask to extract a mmer from a kmer. */

            _mask.setVal(((u_int64_t)1 << (2*_minimizerSize)) - 1);
            _shift = 2*(_nbMinimizers-1);

            /** We initialize the default value of the minimizer.
             * The value is actually set by the Comparator instance provided as a template of the class. */
            Type tmp;
            _cmp.template init<ModelType> (getMmersModel(), tmp);
            _minimizerDefault.set (tmp); //////////max value of minim
			
			u_int64_t nbminims_total = ((u_int64_t)1 << (2*_minimizerSize));
			_mmer_lut = (Type *) MALLOC(sizeof(Type) * nbminims_total ); //free that in destructor
                
            /* if it's ModelDirect, don't do a revcomp; also, use slow method */
            ModelCanonical* isModelCanonical_p = dynamic_cast<ModelCanonical*>(&_kmerModel);
            bool isModelCanonical = isModelCanonical_p != NULL;
            _defaultFast = isModelCanonical;

			for(u_int64_t ii=0; ii< nbminims_total; ii++)
			{
				Type mmer;
				mmer.setVal(ii);
				
				// if(!is_allowed(mmer.getVal(),minimizerSize)) mmer = _mask;
				// if(!is_allowed(rev_mmer.getVal(),minimizerSize)) rev_mmer = _mask;
				
				if (isModelCanonical)
				{
					Type rev_mmer = revcomp(mmer, minimizerSize);
					if(rev_mmer < mmer) mmer = rev_mmer;
					
					//may be cleaner with this
					//if (_cmp (rev_mmer, mmer ) == true)
					//	mmer = rev_mmer;
				}
				
				//std:: cout << "ii " << ii << " is allowed " << is_allowed(mmer.getVal(),minimizerSize) <<  " mmer getval " << mmer.getVal() << " is model canonical " << isModelCanonical << std::endl;
				
				if (!is_allowed(mmer.getVal(),minimizerSize))
					mmer = _mask;
				
				_mmer_lut[ii] = mmer;
			}

            if (freq_order)
                setMinimizersFrequency(freq_order);

			
        }



        /** Destructor */
        ~ModelMinimizer ()
        {
            if (_mmer_lut != 0)  { FREE (_mmer_lut); }
        }

        /** Computes a kmer from a buffer holding nucleotides encoded in some format.
         * The way to interpret the buffer is done through the provided Convert template class.
         * \param[in] seq : holds the nucleotides sequence from which the kmer has to be computed
         * \param[out] kmer : kmer as a result
         * \param[in] startIndex : index of the first nucleotide of the kmer to retrieve in the buffer
         */
        template <class Convert>
        int first (const char* seq, Kmer& kmer, size_t startIndex) const 
        {
            /** We compute the first kmer. */
            int result = _kmerModel.template first<Convert> (seq, kmer, startIndex);

            /** We compute the minimizer of the kmer. */
            computeNewMinimizer (kmer, _defaultFast);

            return result;
        }

        /** Computes a kmer in a recursive way, ie. from a kmer and the next
         * nucleotide. The next nucleotide is computed from a buffer and the
         * index of the nucleotide within the buffer.
         * The way to interpret the buffer is done through the provided Convert template class.
         * \param[in] c : next nucleotide
         * \param[out] kmer : kmer to be updated with the provided next nucleotide
         * \param[in] isValid : tells whether the updated kmer is valid or not
         */
        template <class Convert>
        void  next (char c, Kmer& kmer, bool isValid) const
        {
            /** We compute the next kmer. */
            _kmerModel.template next<Convert> (c, kmer, isValid);

            /** We set the valid status according to the Convert result. */
            kmer._isValid = isValid;

            /** We extract the new mmer from the kmer. also applies the mmer_lut */
            typename ModelType::Kmer mmer = kmer.extract (this->_mask, this->_shift,_mmer_lut);

            /** We update the position of the previous minimizer. */
            kmer._position--;

            /** By default, we consider that the minimizer is still the same. */
            kmer._changed  = false;
                
            /** We have to update the minimizer in the following case:
             *      1) the new mmer is the new minimizer
             *      2) the previous minimizer is invalid or out from the new kmer window.
             */
            if (_cmp (mmer.value(), kmer._minimizer.value()) == true) // .value()
            {
                kmer._minimizer = mmer; // extract() above has already done the job of querying mmer_lut for revcomp / forbidden kmers
                kmer._position  = _nbMinimizers - 1;
                kmer._changed   = true;
            }

            else if (kmer._position < 0)
            {
                computeNewMinimizer (kmer, _defaultFast);
            }
        }

        /** Get the minimizer value of the provided kmer. Note that minimizers are supposed to be
         * of small sizes, so their values can fit a u_int64_t type.
         * \return the miminizer value as an integer. */
        u_int64_t getMinimizerValue (const Type& k, bool fastMethod = true) const
        {
            Kmer km; km.set(k);  this->computeNewMinimizer (km, fastMethod);
            return km.minimizer().value().getVal();
        }


        /* for profiling purpose only */ 
        u_int64_t getMinimizerValueDummy (const Type& k) 
        {
            Kmer km; km.set(k);  /* don't execute anything, this function is here to get a baseline time */
            km._minimizer = this->_minimizerDefault;
            km._position = 0;
            return km.minimizer().value().getVal();
        }
        
        /** Get the minimizer string of the provided kmer. (used for debugging purposes only)
            \return the miminizer as a nucleotide string. */
        std::string getMinimizerString (const Type& k, bool fastMethod = true) const
        {
            Kmer km; km.set(k);  this->computeNewMinimizer (km, fastMethod);
            return _miniModel.toString(km.minimizer().value());
        }

         /** Get the minimizer position of the provided kmer. (used for debugging purposes only)
            \return the miminizer position */
        int getMinimizerPosition (const Type& k, bool fastMethod = true) const
        {
            Kmer km; km.set(k);  this->computeNewMinimizer(km, fastMethod);
            return km.position();
        }

        /* for profiling only, sweep the kmer to count the number of AA's */
        void sweepForAA(const Type &k) const
        {
            unsigned int dummy;
            Kmer km; km.set(k);
            justSweepForAA(km.value(0), _nbMinimizers, dummy);
        }

		//return value of larger mmer in the freq order
        void setMinimizersFrequency (uint32_t *freq_order)
        {
            _cmp.include_frequency(freq_order);
			
			//Type tmp =_cmp.template computeLargest<ModelType>(getMmersModel(),_minimizerSize);
			//_minimizerDefault.set (tmp);
        }

        // hack to access compare int's, for bcalm, needs to be made cleaner later
        bool compareIntMinimizers( size_t a, size_t b)
        {
            if (!_freq_order) return a <= b;
            if (_freq_order[a] == _freq_order[b])
                return a <= b;
            return _freq_order[a] <= _freq_order[b];
        }

    private:
        ModelType  _kmerModel;
        ModelType  _miniModel;
        size_t     _minimizerSize;
        Comparator _cmp;
        size_t     _nbMinimizers;
        Type       _mask;

		Type * _mmer_lut;
        size_t     _shift;
        typename ModelType::Kmer _minimizerDefault;
        bool       _defaultFast;

        uint32_t *_freq_order;
		

        /** Tells whether a minimizer is valid or not, in order to skip minimizers
         *  that are too frequent. */
        bool is_allowed (uint32_t mmer, uint32_t len)
		{
			if (_freq_order) return true; // every minimizer is allowed in freq order
			
			u_int64_t  _mmask_m1  ;
			u_int64_t  _mask_0101 ;
			u_int64_t  _mask_ma1 ;
			
			//code to ban mmer with AA inside except if at the beginnning
			// A C T G        00   01   10   11
			_mmask_m1  = (1 << ((len-2)*2)) -1 ; //vire 2 premieres lettres m = 8  donne    00 00 11 11 11 11 11 11
			_mask_0101 = 0x5555555555555555  ; //         01 01 01 01 01 01 01 01
			_mask_ma1  = _mask_0101 & _mmask_m1;//        00 00 01 01 01 01 01 01
			
			u_int64_t a1 = mmer; //
			a1 =   ~(( a1 )   | (  a1 >>2 ));  //
			a1 =((a1 >>1) & a1) & _mask_ma1 ;  //
			
			if(a1 != 0) return false;
			
			// if ((mmer & 0x3f) == 0x2a)   return false;   // TTT suffix
			// if ((mmer & 0x3f) == 0x2e)   return false;   // TGT suffix
			// if ((mmer & 0x3c) == 0x28)   return false;   // TT* suffix
			// for (uint32_t j = 0; j < len - 3; ++j)       // AA inside
			//      if ((mmer & 0xf) == 0)  return false;
			//      else                    mmer >>= 2;
			// if (mmer == 0)               return false;   // AAA prefix
			// if (mmer == 0x04)            return false;   // ACA prefix
			// if ((mmer & 0xf) == 0)   return false;       // *AA prefix
			
			return true;
		}
		
        /** Returns the minimizer of the provided vector of mmers. */
        void computeNewMinimizerOriginal(Kmer& kmer) const
        {
            /** We update the attributes of the provided kmer. Note that an invalid minimizer is
             * memorized by convention by a negative minimizer position. */
            kmer._minimizer = this->_minimizerDefault;
            kmer._position  = -1;
            kmer._changed   = true;

            typename ModelType::Kmer mmer;

            /** We compute each mmer and memorize the minimizer among them. */

            Type kmer_minimizer_value = kmer._minimizer.value();
            Type val = kmer.value(0);

            for (int16_t idx=_nbMinimizers-1; idx>=0; idx--)
            {

                /** We extract the most left mmer in the kmer. */
                Type candidate_minim = _mmer_lut[(val & _mask).getVal()];
				

                /** We check whether this mmer is the new minimizer. */
                if (_cmp (candidate_minim, kmer_minimizer_value ) == true)  
                {
                    mmer.set(candidate_minim);
                    kmer._minimizer = mmer;   
                    kmer._position = idx; 
                    kmer_minimizer_value = candidate_minim; 
                }
            
                val >>= 2;    
            }
        }
   
        /** Returns the minimizer of the provided vector of mmers, fast method (may fallback to normal method)
         * Note: only used for KmerCanonicals */
        void computeNewMinimizer(KmerMinimizer<ModelCanonical, Comparator>& kmer, bool fastMethod = true) const 
        {
            //if (!fastMethod || _freq_order) // fast method doesn't work with frequency order
			//temporarily desactivate fastmode, seem to returns wrong minimizer and leads to duplicated kmer output
            {
                computeNewMinimizerOriginal(kmer);
                return;
            }

            bool validResult;
            size_t position = -1;
            u_int32_t minim;
            fastLexiMinimizer(kmer.value(0), _nbMinimizers, _minimizerSize, minim, position, validResult);

            if (!validResult)
            {
                computeNewMinimizerOriginal(kmer);
                //_invalidMinimizersCounter++; // have to get rid of this metric else that function isn't "const" anymore, and it's a cascade (next() uses it for instance)
            }
            else
            {
                kmer._minimizer.set(minim);
                kmer._changed = true;
                kmer._position = position; 
                // _minimizersCounter++; // have to get rid of this metric else that function isn't "const" anymore, and it's a cascade
            }
        }

        /* direct model cannot use the speed-up version because it handles reverse complements in a hardcoded way */
        void computeNewMinimizer(KmerMinimizer<ModelDirect, Comparator>  &kmer, bool fastMethod = true) const
        {
            computeNewMinimizerOriginal(kmer);
            //std::cout << "specialized kmer minimizer " << _miniModel.toString(kmer.minimizer().value()) << std::endl;
        }



    };


    /************************************************************/
    /*********************  SUPER KMER    ***********************/
    /************************************************************/
	
	//now with  vector containing the overlapping kmers of the superkmers (instead of reference to large external vector buffer)
    class SuperKmer
    {
    public:

        //typedef Type SType[2];

#ifdef NONCANONICAL
        typedef ModelMinimizer<ModelDirect> Model;
#else
        typedef ModelMinimizer<ModelCanonical> Model;
#endif
        typedef typename Model::Kmer           Kmer;

        static const u_int64_t DEFAULT_MINIMIZER = 1000000000 ;

        SuperKmer (size_t kmerSize, size_t miniSize)
            : minimizer(DEFAULT_MINIMIZER), kmerSize(kmerSize), miniSize(miniSize)
        {
			_max_size_sk = 1000;
			kmers.clear();
			_sk_buffer = (u_int8_t *) malloc(_max_size_sk);
			_sk_buffer_idx=0;
          //  if (kmers.empty())  { kmers.resize(kmerSize); range.second = kmers.size()-1; }
        }

        u_int64_t                minimizer;

        Kmer& operator[] (size_t idx)  {
			return kmers[idx];
		}

		size_t size() const {
			return kmers.size();
		}

        bool isValid() const { return minimizer != DEFAULT_MINIMIZER; }

		void addKmer(Kmer newkmer)
		{
			kmers.push_back(newkmer);
		}
		
		void reset()
		{
			kmers.clear();
			//binrep.clear();
			_sk_buffer_idx =0;
		}
		
		//save superkmer to CacheSuperKmerBinFiles

        // for kmtricks
        template<typename Storage>
        void save(int file_id, Storage* storage)
        {
			size_t superKmerLen = size();

			int required_bytes = (superKmerLen + kmerSize +3) /4 ;
			if(required_bytes > _max_size_sk)
			{
				_sk_buffer = (u_int8_t *) realloc(_sk_buffer, _max_size_sk);
				_max_size_sk = required_bytes;
			}
			_sk_buffer_idx =0;
			Type basekmer = (*this)[0].forward();
			int rem_size = kmerSize;
			u_int8_t newbyte=0;
			u_int64_t mask4nt  = 255;
			u_int64_t mask1nt  = 3;
			while(rem_size>=4)
			{
				newbyte = basekmer.getVal() & mask4nt ; // ptet un getVal et cast to u_int8_t
				rem_size -= 4;
				basekmer = basekmer >> 8;
				_sk_buffer[_sk_buffer_idx++]= newbyte;
			}
			//reste du kmer
			newbyte = basekmer.getVal() & mask4nt;
			int uid = rem_size; //uid = nb nt used in this newbyte

			uint skid =1;
			while(true)
			{
				while(uid<4 && skid < superKmerLen)
				{
					u_int8_t newnt = ((*this)[skid].forward()).getVal() & mask1nt ;
					newbyte |=  newnt << (uid*2);
					uid++; skid++;
				}
				if(uid > 0)
					_sk_buffer[_sk_buffer_idx++]= newbyte;
				if(skid >= superKmerLen) break;

				newbyte=0; uid=0;
			}
			//storage->writeBlock(_sk_buffer, _sk_buffer_idx, file_id, kmers.size());
			storage->insertSuperkmer(_sk_buffer, _sk_buffer_idx, kmers.size(), file_id);
        }

		void save(tools::storage::impl::CacheSuperKmerBinFiles  & cacheSuperkFile, int file_id)
		{
//			printf("saving superk to file %i \n",file_id);
//			//debug
//			for (size_t ii=0 ; ii < kmers.size(); ii++)
//			{
//			
//				printf("%s\n",	(((*this)[ii].forward()).toString(kmerSize)).c_str());
//
//			}
//			//
			size_t superKmerLen = size();

			
			int required_bytes = (superKmerLen + kmerSize +3) /4 ;
			if(required_bytes > _max_size_sk)
			{
				_sk_buffer = (u_int8_t *) realloc(_sk_buffer, _max_size_sk);
				_max_size_sk = required_bytes;
			}
			
			
			//binrep.clear();
			_sk_buffer_idx =0;
			
			Type basekmer = (*this)[0].forward();
			
			int rem_size = kmerSize;
			u_int8_t newbyte=0;
			u_int64_t mask4nt  = 255;
			u_int64_t mask1nt  = 3;
			
			while(rem_size>=4)
			{
				newbyte = basekmer.getVal() & mask4nt ; // ptet un getVal et cast to u_int8_t
				rem_size -= 4;
				basekmer = basekmer >> 8;
				_sk_buffer[_sk_buffer_idx++]= newbyte;
				//binrep.push_back(newbyte);
//				//debug pushing
//				Type dd; dd.setVal(newbyte);
//				printf("pushing %s\n",	(dd.toString(4)).c_str());
//				//
			}
			
			//reste du kmer
			newbyte = basekmer.getVal() & mask4nt;
			int uid = rem_size; //uid = nb nt used in this newbyte

			//reste du newbyte avec le superk

			uint skid =1;
			
			while(true)
			{

				while(uid<4 && skid < superKmerLen)
				{
					
					u_int8_t newnt = ((*this)[skid].forward()).getVal() & mask1nt ;
					
					newbyte |=  newnt << (uid*2);
					uid++; skid++;
				}
				
				if(uid > 0)
					_sk_buffer[_sk_buffer_idx++]= newbyte;
				
				//binrep.push_back(newbyte);
//				//debug pushing
//				Type dd; dd.setVal(newbyte);
//				printf("pushing %s\n",	(dd.toString(4)).c_str());
//				//
				
				if(skid >= superKmerLen) break;

				newbyte=0; uid=0;
			}
			
			
			//printf("insert superK %i  _sk_buffer_idx %i \n",kmers.size(),_sk_buffer_idx);
			
		//	cacheSuperkFile.insertSuperkmer(binrep.data(), binrep.size(), kmers.size(),  file_id);
			cacheSuperkFile.insertSuperkmer(_sk_buffer, _sk_buffer_idx, kmers.size(),  file_id);
			
		}
		
		
        /** */
        void save (tools::collections::Bag<Type>& bag)
        {
            size_t superKmerLen = size();

            int64_t zero = 0;
            Type masknt;
            masknt.setVal((int64_t) 3);
            Type nbK;
            nbK.setVal((int64_t) size());
            Type compactedK;
            compactedK.setVal(zero);

            for (size_t ii=1 ; ii < superKmerLen; ii++)
            {
                compactedK = compactedK << 2  ;
                compactedK = compactedK | ( ((*this)[ii].forward()) & masknt) ;
            }

            int maxs = (compactedK.getSize() - 8 ) ;

            compactedK = compactedK | (  nbK << maxs ) ;

            bag.insert (compactedK);
            bag.insert ((*this)[0].forward());
        }

        /** NOT USED YET. */
#if 0
        void load (tools::dp::Iterator<Type>& iter)
        {
            Type superk = iter.item(); iter.next();
            Type seedk = iter.item();

            u_int8_t        nbK, rem ;
            Type compactedK;
            int ks = kmerSize;
            Type un ; un.setVal( 1);
            size_t _shift_val = Type::getSize() -8;
            Type kmerMask = (un << (ks*2)) - un;
            size_t shift = 2*(ks-1);

            compactedK =  superk;
            nbK = (compactedK >> _shift_val).getVal() & 255; // 8 bits poids fort = cpt //todo for large k values
            rem = nbK;

            Type temp = seedk;
            Type rev_temp = revcomp(temp,ks);
            Type newnt ;
            Type mink;

            /** We loop over each kmer of the current superkmer. */
            for (int ii=0; ii<nbK; ii++,rem--)
            {
                mink = std::min (rev_temp, temp);

                /** We set the current (canonical) kmer. */
                kmers[ii].set (rev_temp, temp);

                if(rem < 2) break;
                newnt =  ( superk >> ( 2*(rem-2)) ) & 3 ;

                temp = ((temp << 2 ) |  newnt   ) & kmerMask;
                newnt .setVal(comp_NT[newnt.getVal()]) ;
                rev_temp = ((rev_temp >> 2 ) |  (newnt << shift) ) & kmerMask;
            }
        }
#endif

		~SuperKmer()
		{
			free(_sk_buffer);
		}
		
    private:
        size_t              kmerSize;
        size_t              miniSize;
        std::vector<Kmer>  kmers;
		//std::vector<u_int8_t> binrep;
		u_int8_t * _sk_buffer;
		int _sk_buffer_idx;
		int _max_size_sk;
    };

    /************************************************************/
    /***********************     COUNT    ***********************/
    /************************************************************/
    /** \brief Structure associating a kmer value with an abundance value.
     *
     * This structure is useful for methods that counts kmer, such as the SortingCount algorithm.
     * It is also interesting to save [kmer,abundance] in a h d f 5 format.
     *
     * By default, the abundance value is coded on 32 bits, so abundance up to 1<<32 can be used.
     */
    struct Count : tools::misc::Abundance<Type,CountNumber>
    {
        /** Constructor.
         * \param[in] val : integer value of the kmer
         * \param[in] abund : abundance for the kmer */
        Count(const Type& val, const CountNumber& abund) : tools::misc::Abundance<Type,CountNumber>(val, abund) {}

        /** Default constructor. */
        Count() : tools::misc::Abundance<Type,CountNumber>() {}

        /** Copy constructor. */
        Count(const Count& val) : tools::misc::Abundance<Type,CountNumber>(val.value, val.abundance) {}

        /** Comparison operator
         * \param[in] other : object to be compared to
         * \return true if the provided kmer value is greater than the current one. */
        bool operator< (const Count& other) const {  return this->value < other.value; }
        
        /** Equal operator
         * \param[in] other : object to be compared to
         * \return true if the provided kmer value is greater than the current one. */
        bool operator== (const Count& other) const {  return (this->value == other.value && this->abundance == other.abundance); }
    };

};  // struct Kmer

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/

#endif /* _GATB_CORE_KMER_IMPL_MODEL_HPP_ */
