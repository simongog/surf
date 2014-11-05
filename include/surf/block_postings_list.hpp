#ifndef SURF_BLOCK_POSTINGS_LIST_H
#define SURF_BLOCK_POSTINGS_LIST_H

#include <limits>
#include <stdexcept>

#include "util.h"
#include "memutil.h"
#include "codecs.h"
#include "codecfactory.h"
#include "bitpacking.h"
#include "simdfastpfor.h"
#include "deltautil.h"

#include "sdsl/int_vector.hpp"

using namespace sdsl;

namespace surf {

struct vbyte_coder {
    static size_t encode_num(uint32_t num,uint8_t* out) {
        size_t written_bytes = 0;
        uint8_t w = num & 0x7F;
        num >>= 7;
        while (num > 0) {
            w |= 0x80; // mark overflow bit
            *out = w;
            ++out;
            w = num & 0x7F;
            num >>= 7;
            written_bytes++;
        }
        *out = w;
        ++out;
        written_bytes++;
        return written_bytes;
    }
    static uint32_t decode_num(const uint8_t*& in) {
        uint32_t num = 0;
        uint8_t w=0;
        uint32_t shift=0;
        do {
            w = *in;
            in++;    
            num |= (((uint32_t)(w&0x7F))<<shift);
            shift += 7;
        } while ((w&0x80) > 0);
        return num;
    }
    static void encode(const uint32_t* A,size_t n,uint32_t* out,size_t& written_u32s) {
        uint8_t* out_bytes = (uint8_t*) out;
        size_t written_bytes = 0;
        for(size_t i=0;i<n;i++) {
            size_t written = encode_num(A[i],out_bytes);
            out_bytes += written;
            written_bytes += written;
        }
        written_u32s = written_bytes/4;
        if(written_bytes%4 != 0) written_u32s++;
    }
    static void decode(const uint32_t* in,size_t n,uint32_t* out) {
        const uint8_t* in_bytes = (const uint8_t*) in;
        for(size_t i=0;i<n;i++) {
            *out = decode_num(in_bytes);
            out++;
        }
    }
};


template<uint64_t t_block_size>
class block_postings_list;

template<uint64_t t_block_size>
class plist_iterator
{
    public:
        typedef block_postings_list<t_block_size>                         list_type;
        typedef typename list_type::size_type                             size_type;
        typedef uint64_t                                                 value_type;
    public: // default implementation used. not necessary to list here
        plist_iterator() = default;
        plist_iterator(const plist_iterator& pi) = default;
        plist_iterator(plist_iterator&& pi) = default;
        plist_iterator& operator=(const plist_iterator& pi) = default;
        plist_iterator& operator=(plist_iterator&& pi) = default;
    public:
        plist_iterator(const list_type& l,size_t pos);
        plist_iterator& operator++();
        bool operator ==(const plist_iterator& b) const;
        bool operator !=(const plist_iterator& b) const;
        uint64_t docid() const;
        uint64_t freq() const;
        void skip_to_id(uint64_t id);
        void skip_to_block_with_id(uint64_t id);
        uint64_t block_rep() const { return m_plist_ptr->block_rep(m_cur_block_id); }
        size_t size() const { return m_plist_ptr->size(); }
        size_t remaining() const { return size() - m_cur_pos; }
        size_t offset() const { return m_cur_pos; }
    private:
        void access_and_decode_cur_pos() const;
    private:
        size_type m_cur_pos = std::numeric_limits<uint64_t>::max();
        mutable size_type m_cur_block_id = std::numeric_limits<uint64_t>::max();
        mutable size_type m_last_accessed_block = std::numeric_limits<uint64_t>::max()-1;
        mutable size_type m_last_accessed_id = std::numeric_limits<uint64_t>::max()-1;
        mutable value_type m_cur_docid = 0;
        mutable value_type m_cur_freq = 0;
        const list_type* m_plist_ptr = nullptr;
        mutable std::vector<uint32_t, FastPForLib::cacheallocator> m_decoded_ids;
        mutable std::vector<uint32_t, FastPForLib::cacheallocator> m_decoded_freqs;
};


template<uint64_t t_block_size=128>
class block_postings_list {
	static_assert(t_block_size % 32 == 0,"blocksize must be multiple of 32.");
public: // types
	friend class plist_iterator<t_block_size>;
	using comp_codec = FastPForLib::OPTPFor<t_block_size/32,FastPForLib::Simple16<false>>;
	using size_type = sdsl::int_vector<>::size_type;
	using const_iterator = plist_iterator<t_block_size>;
	using pfor_data_type = std::vector<uint32_t, FastPForLib::cacheallocator>;
	#pragma pack(push, 1)
	struct block_data {
		uint32_t max_block_id = 0;
		uint32_t id_offset = 0;
		uint32_t freq_offset = 0;
	};
	#pragma pack(pop)
private: // actual data
	uint32_t m_size = 0;
	double m_list_maximuim = std::numeric_limits<double>::lowest();
	double m_max_doc_weight = std::numeric_limits<double>::lowest();
	std::vector<block_data> m_block_data;
    pfor_data_type m_docid_data;
    pfor_data_type m_freq_data;
public: // default 
    block_postings_list() {
    	m_block_data.resize(1);
    }
    block_postings_list(const block_postings_list& pl) = default;
    block_postings_list(block_postings_list&& pl) = default;
    block_postings_list& operator=(const block_postings_list& pi) = default;
    block_postings_list& operator=(block_postings_list&& pi) = default;
    double list_max_score() const { return m_list_maximuim; };
    double max_doc_weight() const { return m_max_doc_weight; };
public: // constructors
    block_postings_list(std::istream& in) {
        load(in);
    }
    template<class t_rank> 
    block_postings_list(const t_rank& ranker,sdsl::int_vector<>& D,size_t sp,size_t ep) {
	    if (ep<sp) {
	        std::cerr << "ERROR: trying to create empty postings list.\n";
	        throw std::logic_error("trying to create empty postings list.");
	    }

	    std::sort(D.begin()+sp,D.begin()+ep+1);

	    // count uniq docs
	    size_t unique = 1;
	    for(size_t i=sp+1;i<=ep;i++) {
	        if(D[i] != D[i-1]) unique++;
	    }

	    // extract doc_ids and freqs
    	m_size = unique;
	    sdsl::int_vector<32> tmp_data(unique);
	    sdsl::int_vector<32> tmp_freq(unique);
	    size_type j = 0; size_type freq = 1;
	    for (size_type i=sp+1; i<=ep; i++) {
	        if (D[i] != D[i-1]) {
	            tmp_data[j] = D[i-1];
	            tmp_freq[j] = freq;
	            j++;
	            freq = 0;
	        }
	        freq++;
	    }
	    tmp_data[j] = D[ep];
	    tmp_freq[j] = freq;


	    // create block max structure first
	    create_block_support(tmp_data);

	    // create rank support structure
	    create_rank_support(tmp_data,tmp_freq,ranker);

	    // compress postings
	    compress_postings_data(tmp_data,tmp_freq);
    }
    template<class t_rank> 
    block_postings_list(const t_rank& ranker,
    			  std::vector<std::pair<uint64_t,uint64_t>>& pre_sorted_data) 
    	: block_postings_list(pre_sorted_data)
    {
    	m_size = pre_sorted_data.size();

	    // extract doc_ids and freqs
	    sdsl::int_vector<32> tmp_data(pre_sorted_data.size());
	    sdsl::int_vector<32> tmp_freq(pre_sorted_data.size());
	    for (size_type i=0; i<pre_sorted_data.size(); i++) {
	        tmp_data[i] = pre_sorted_data[i].first;
	        tmp_freq[i] = pre_sorted_data[i].second;
	    }

	    // create block max structure first
	    create_block_support(tmp_data);

	    // create rank support structure
	    create_rank_support(tmp_data,tmp_freq,ranker);

	    // compress postings
	    compress_postings_data(tmp_data,tmp_freq);
    }
    block_postings_list(std::vector<std::pair<uint64_t,uint64_t>>& pre_sorted_data) {
    	m_size = pre_sorted_data.size();

	    // extract doc_ids and freqs
	    sdsl::int_vector<32> tmp_data(pre_sorted_data.size());
	    sdsl::int_vector<32> tmp_freq(pre_sorted_data.size());
	    for (size_type i=0; i<pre_sorted_data.size(); i++) {
	        tmp_data[i] = pre_sorted_data[i].first;
	        tmp_freq[i] = pre_sorted_data[i].second;
	    }

	    // create block max structure first
	    create_block_support(tmp_data);

	    // compress postings
	    compress_postings_data(tmp_data,tmp_freq);
    }
private: // functions used during construction
	void create_block_support(const sdsl::int_vector<32>& ids)
	{
	    size_t num_blocks = ids.size() / t_block_size;
	    if (ids.size() % t_block_size != 0) num_blocks++;
	    m_block_data.resize(num_blocks);
	    size_t j = 0;
	    for (size_t i=t_block_size-1; i<ids.size(); i+=t_block_size) {
	        m_block_data[j++].max_block_id = ids[i];
	    }
	    if (ids.size() % t_block_size != 0) m_block_data[j].max_block_id = ids[ids.size()-1];
	}
	template<class t_rank>
	void create_rank_support(const sdsl::int_vector<32>& ids,
							 const sdsl::int_vector<32>& freqs,
							 const t_rank& ranker)
	{
		auto F_t = std::accumulate(freqs.begin(),freqs.end(),0);
		auto f_t = ids.size();
	    for (size_t l=0; l<ids.size(); l++) {
	        auto id = ids[l];
	        auto f_dt = freqs[l];
	        double W_d = ranker.doc_length(id);
	        double doc_weight = ranker.calc_doc_weight(W_d);
	        double score = ranker.calculate_docscore(1.0f,f_dt,f_t,F_t,W_d,true);
	        m_list_maximuim = std::max(m_list_maximuim,score);
	        m_max_doc_weight = std::max(m_max_doc_weight,doc_weight);
	    }
	}
	void compress_postings_data(const sdsl::int_vector<32>& ids,
						        sdsl::int_vector<32>& freqs)
	{
		// delta compress ids first
		uint32_t* id_input = (uint32_t*) ids.data();
		FastPForLib::Delta::fastDelta(id_input,ids.size());

        // substract one from all freqs
        for(size_t i=0;i<freqs.size();i++) freqs[i]--;

	    // encode ids and freqs using pfor
	    static comp_codec c;
	    m_docid_data.resize(2 * ids.size() + 1024);
	    uint32_t* id_out = m_docid_data.data();
	    m_freq_data.resize(2 * freqs.size() + 1024);
	    uint32_t* freq_out = m_freq_data.data();
	    uint32_t* freq_input = (uint32_t*) freqs.data();

	    size_type cur_block = 0;
	    uint64_t id_offset = 0;
	    uint64_t freq_offset = 0;
	    size_t encoded_id_size = 0;
	    size_t encoded_freq_size = 0;
	    for (size_t i=0; i<ids.size(); i+=t_block_size) {
	    	if(i+t_block_size > ids.size()) break;
	    	m_block_data[cur_block].id_offset = id_offset;
	    	m_block_data[cur_block].freq_offset = freq_offset;
	    	c.encodeBlock(&id_input[i],&id_out[id_offset],encoded_id_size);
	    	c.encodeBlock(&freq_input[i],&freq_out[freq_offset],encoded_freq_size);
	    	id_offset += encoded_id_size;
	    	freq_offset += encoded_freq_size;
	    	cur_block++;
	    }

	    // any non-full blocks?
	    size_type n = ids.size() % t_block_size;
	    if( n != 0 ) {
	    	m_block_data[cur_block].id_offset = id_offset;
	    	m_block_data[cur_block].freq_offset = freq_offset;
	    	size_type i = ids.size() - n;
	    	vbyte_coder::encode(&id_input[i],n,&id_out[id_offset],encoded_id_size);
	    	vbyte_coder::encode(&freq_input[i],n,&freq_out[freq_offset],
	    						encoded_freq_size);
	    	id_offset += encoded_id_size;
	    	freq_offset += encoded_freq_size;
	    }
	    m_docid_data.resize(id_offset);
	    m_docid_data.shrink_to_fit();
	    m_freq_data.resize(freq_offset);
	    m_freq_data.shrink_to_fit();
	}
public: // functions used during processing
	void decompress_block(size_t block_id,
						  pfor_data_type& id_data,
						  pfor_data_type& freq_data) const
	{
		uint32_t delta_offset = 0;
		if(block_id != 0) {
			delta_offset = m_block_data[block_id-1].max_block_id;
		}

		const uint32_t* id_start = m_docid_data.data() + 
							m_block_data[block_id].id_offset;
		const uint32_t* freq_start = m_freq_data.data() + 
							m_block_data[block_id].freq_offset;
		auto block_size = postings_in_block(block_id);

	    if (id_data.size() != block_size) { // did we allocate space already?
	        id_data.resize(block_size);
	        freq_data.resize(block_size);
	    }

	    size_t rec_ids;
	    size_t rec_freqs;
		if(block_size == t_block_size) { // PFor
			static comp_codec c;
			c.decodeBlock(id_start,id_data.data(),rec_ids);
			c.decodeBlock(freq_start,freq_data.data(),rec_freqs);
		} else { // vbyte
			vbyte_coder::decode(id_start,block_size,id_data.data());
			vbyte_coder::decode(freq_start,block_size,freq_data.data());
			rec_ids = rec_freqs = block_size;
		}

		// undo delta compression
		id_data[0] += delta_offset;
        freq_data[0]++;
		for(size_t i=1;i<block_size;i++) {
			id_data[i] += id_data[i-1];
            freq_data[i]++;
		}

		if( rec_ids != rec_freqs ) {
	        std::cerr << "ERROR: number of decoded ids and freqs is not equal. "
	                  << rec_ids << " != " << rec_freqs << "\n";
	        throw std::logic_error("number of decoded ids and freqs is not equal.");
		}
	}

	size_type find_block_with_id(uint64_t id,size_t start_block) const {
	    size_t block_id = start_block;
	    size_t nblocks = m_block_data.size();
	    while (block_id < nblocks && m_block_data[block_id].max_block_id < id) {
	        block_id++;
	    }
	    return block_id;
	}
	size_type size() const {
		return m_size;
	}
	uint32_t block_rep(size_t bid) const {
		return m_block_data[bid].max_block_id;
	}
	size_type num_blocks() const {
		return m_block_data.size();
	}
	size_type postings_in_block(size_type block_id) const {
		size_type block_size = t_block_size;
		size_type mod = m_size % t_block_size;
		if(block_id == m_block_data.size()-1 && mod != 0) {
			block_size = mod;
		}
		return block_size;
	}
    const_iterator begin() const {
        return const_iterator(*this,0);
    }
    const_iterator end() const {
        return const_iterator(*this,m_size);
    }
    auto serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, 
    			   std::string name = "") const -> size_type 
	{
	    size_type written_bytes = 0;

		sdsl::structure_tree_node* child;
	    if(m_size <= t_block_size) { // only one block
	    	child = sdsl::structure_tree::add_child(v,"single block list",sdsl::util::class_name(*this));
	    } else {
	    	child = sdsl::structure_tree::add_child(v, "multi block list",sdsl::util::class_name(*this));
	    }

	    written_bytes += sdsl::write_member(m_size,out,child,"size");

	    if(m_size <= t_block_size) { // only one block
	    	written_bytes += sdsl::write_member(m_block_data[0].max_block_id,out,child,"max block id");
	    } else {
	    	auto* blockdata = sdsl::structure_tree::add_child(child, "block data","block data");
	    	out.write((const char*)m_block_data.data(), m_block_data.size()*sizeof(block_data));
	    	written_bytes += m_block_data.size()*sizeof(block_data);
	    	sdsl::structure_tree::add_size(blockdata, m_block_data.size()*sizeof(block_data));
	    }

        uint32_t docidu32 = m_docid_data.size();
        uint32_t frequ32 = m_freq_data.size();
        written_bytes += sdsl::write_member(docidu32,out,child,"docid u32s");
        written_bytes += sdsl::write_member(frequ32,out,child,"freq u32s");

    	auto* idchild = sdsl::structure_tree::add_child(child, "id data","delta compressed");
        out.write((const char*)m_docid_data.data(), m_docid_data.size()*sizeof(uint32_t));
        sdsl::structure_tree::add_size(idchild, m_docid_data.size()*sizeof(uint32_t));
        written_bytes +=  m_docid_data.size()*sizeof(uint32_t);

        auto* fchild = sdsl::structure_tree::add_child(child, "freq data", "compressed");
        out.write((const char*)m_freq_data.data(), m_freq_data.size()*sizeof(uint32_t));
        written_bytes +=  m_freq_data.size()*sizeof(uint32_t);
    	sdsl::structure_tree::add_size(fchild, m_freq_data.size()*sizeof(uint32_t));

	    written_bytes += sdsl::write_member(m_list_maximuim,out,child,"list max score");
	    written_bytes += sdsl::write_member(m_max_doc_weight,out,child,"max doc weight");

	    sdsl::structure_tree::add_size(child, written_bytes);
	    return written_bytes;
	}
	void load(std::istream& in) {
		read_member(m_size,in);
		if(m_size <= t_block_size) { // only one block
			uint32_t max_block_id;
			read_member(max_block_id,in);
			m_block_data.resize(1);
			m_block_data[0].max_block_id = max_block_id;
		} else {
			uint64_t num_blocks = m_size / t_block_size;
			if(m_size % t_block_size != 0) num_blocks++;
			m_block_data.resize(num_blocks);
			in.read((char*)m_block_data.data(),num_blocks*sizeof(block_data));
		}

		// load compressed data
        uint32_t docidu32;
        uint32_t frequ32;
        read_member(docidu32,in);
        read_member(frequ32,in);
        m_docid_data.resize(docidu32);
        m_freq_data.resize(frequ32);
        in.read((char*)m_docid_data.data(),docidu32*sizeof(uint32_t));
        in.read((char*)m_freq_data.data(),frequ32*sizeof(uint32_t));

	    read_member(m_list_maximuim,in);
	    read_member(m_max_doc_weight,in);
	}
};


template<uint64_t t_bs>
plist_iterator<t_bs>::plist_iterator(const list_type& l,size_t pos) : plist_iterator()
{
    m_cur_pos = pos;
    m_plist_ptr = &l;
}

template<uint64_t t_bs>
plist_iterator<t_bs>& plist_iterator<t_bs>::operator++()
{
    if (m_cur_pos != size()) { // end?
        (*this).m_cur_pos++;
    } else {
        std::cerr << "ERROR: trying to advance plist iterator beyond list end.\n";
        throw std::out_of_range("trying to advance plist iterator beyond list end");
    }
    return (*this);
}

template<uint64_t t_bs>
bool plist_iterator<t_bs>::operator ==(const plist_iterator& b) const
{
    return ((*this).m_cur_pos == b.m_cur_pos) && ((*this).m_plist_ptr == b.m_plist_ptr);
}

template<uint64_t t_bs>
bool plist_iterator<t_bs>::operator !=(const plist_iterator& b) const
{
    return !((*this)==b);
}

template<uint64_t t_bs>
typename plist_iterator<t_bs>::value_type plist_iterator<t_bs>::docid() const
{
    if (m_cur_pos == m_plist_ptr->size()) { // end?
        std::cerr << "ERROR: plist iterator dereferenced at list end.\n";
        throw std::out_of_range("plist iterator dereferenced at list end");
    }
    if (m_cur_pos == m_last_accessed_id) {
        return m_cur_docid;
    }
    access_and_decode_cur_pos();
    return m_cur_docid;
}

template<uint64_t t_bs>
typename plist_iterator<t_bs>::value_type plist_iterator<t_bs>::freq() const
{
    if (m_cur_pos == m_plist_ptr->size()) { // end?
        std::cerr << "ERROR: plist iterator dereferenced at list end.\n";
        throw std::out_of_range("plist iterator dereferenced at list end");
    }
    if (m_cur_pos == m_last_accessed_id) {
        return m_cur_freq;
    }
    access_and_decode_cur_pos();
    return m_cur_freq;
}

template<uint64_t t_bs>
void plist_iterator<t_bs>::access_and_decode_cur_pos() const
{
    m_cur_block_id = m_cur_pos / t_bs;
    if (m_cur_block_id != m_last_accessed_block) {  // decompress block
        m_last_accessed_block = m_cur_block_id;
        m_plist_ptr->decompress_block(m_cur_block_id,m_decoded_ids,m_decoded_freqs);
    }
    size_t in_block_offset = m_cur_pos % t_bs;
    m_cur_docid = m_decoded_ids[in_block_offset];
    m_cur_freq = m_decoded_freqs[in_block_offset];
    m_last_accessed_id = m_cur_pos;
}

template<uint64_t t_bs>
void plist_iterator<t_bs>::skip_to_block_with_id(uint64_t id)
{
    size_t old_block = m_cur_block_id;
    m_cur_block_id = m_plist_ptr->find_block_with_id(id,m_cur_block_id);

    // we now go to the first id in the new block!
    if (old_block != m_cur_block_id) {
        m_cur_pos = m_cur_block_id*t_bs;
        if (m_cur_pos > m_plist_ptr->size()) { // don't go past the end!
            m_cur_pos = m_plist_ptr->size();
        }
    }
}

template<uint64_t t_bs>
void plist_iterator<t_bs>::skip_to_id(uint64_t id)
{
    if(id == m_cur_docid) {
        return;
    }

    skip_to_block_with_id(id);
    // check if we reached list end!
    if (m_cur_block_id >= m_plist_ptr->num_blocks()) {
        m_cur_pos = m_plist_ptr->size();
        return;
    }
    if (m_last_accessed_block != m_cur_block_id) {
        m_last_accessed_block = m_cur_block_id;
        //std::cout << "skip decompress" << std::endl;
        m_plist_ptr->decompress_block(m_cur_block_id,m_decoded_ids,m_decoded_freqs);
        //std::cout << "size = " << m_decoded_ids.size() << std::endl;
        // new block -> find from the beginning
        auto block_itr = std::lower_bound(m_decoded_ids.begin(),m_decoded_ids.end(),id);
        m_cur_pos = (t_bs*m_cur_block_id) + std::distance(m_decoded_ids.begin(),block_itr);
    } else {
        size_t in_block_offset = m_cur_pos % t_bs;
        auto block_itr = std::lower_bound(m_decoded_ids.begin()+in_block_offset,m_decoded_ids.end(),id);
        m_cur_pos = (t_bs*m_cur_block_id) + std::distance(m_decoded_ids.begin(),block_itr);
    }
    size_t inblock_offset = m_cur_pos % t_bs;
    m_cur_docid = m_decoded_ids[inblock_offset];
    m_cur_freq = m_decoded_freqs[inblock_offset];
    m_last_accessed_id = m_cur_pos;
}


}

#endif
