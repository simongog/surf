#ifndef SURF_INVFILE_POSTINGS_LIST_HPP
#define SURF_INVFILE_POSTINGS_LIST_HPP

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

enum class compression_codec : std::uint8_t {
    fastbinarypacking8 = 0,
    fastbinarypacking16 = 1,
    fastbinarypacking32 = 2,
    BP32 = 3,
    vsencoding = 4,
    fastpfor = 5,
    simdfastpfor = 6,
    simplepfor = 7,
    pfor = 8,
    pfor2008 = 9,
    newpfor = 10,
    optpfor = 11,
    vbyte = 12,
    simple8b = 13,
    varintg8iu = 14,
    simdbinarypacking = 15
};

static std::map<compression_codec,std::string> compression_codec_names {
    {compression_codec::fastbinarypacking8,"fastbinarypacking8"},
    {compression_codec::fastbinarypacking16,"fastbinarypacking16"},
    {compression_codec::fastbinarypacking32,"fastbinarypacking32"},
    {compression_codec::BP32,"BP32"},
    {compression_codec::vsencoding,"vsencoding"},
    {compression_codec::fastpfor,"fastpfor"},
    {compression_codec::simdfastpfor,"simdfastpfor"},
    {compression_codec::simplepfor,"simplepfor"},
    {compression_codec::pfor,"pfor"},
    {compression_codec::pfor2008,"pfor2008"},
    {compression_codec::newpfor,"newpfor"},
    {compression_codec::optpfor,"optpfor"},
    {compression_codec::vbyte,"vbyte"},
    {compression_codec::simple8b,"simple8b"},
    {compression_codec::varintg8iu,"varintg8iu"},
    {compression_codec::simdbinarypacking,"simdbinarypacking"}
};

template<compression_codec t_codec,uint64_t t_block_size>
class postings_list;

template<compression_codec t_codec,uint64_t t_block_size>
class plist_iterator
{
    public:
        typedef postings_list<t_codec,t_block_size>                       list_type;
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
        double block_max() const {
            return m_plist_ptr->block_max(m_cur_block_id); 
        }
        double block_max_doc_weight() const {
            return m_plist_ptr->block_max_doc_weight(m_cur_block_id); 
        }
        uint64_t block_rep() const {
            return m_plist_ptr->block_rep(m_cur_block_id); 
        }
        size_t size() const {
            return m_plist_ptr->size();
        }
        size_t remaining() const {
            return size() - m_cur_pos;
        }
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
        mutable std::vector<uint32_t, FastPFor::cacheallocator> m_decoded_ids;
        mutable std::vector<uint32_t, FastPFor::cacheallocator> m_decoded_freqs;
};


template<compression_codec t_codec=compression_codec::optpfor,uint64_t t_block_size=128>
class postings_list
{
    public:
        friend class plist_iterator<t_codec,t_block_size>;
        typedef sdsl::int_vector<>::size_type                             size_type;
        typedef std::pair<uint64_t,uint64_t>                             value_type;
        typedef plist_iterator<t_codec,t_block_size>                 const_iterator;
        typedef std::vector<uint32_t, FastPFor::cacheallocator>      pfor_data_type;
    public: // not necessary. just listed to be verbose
        postings_list() = default;
        postings_list(const postings_list& pl) = default;
        postings_list(postings_list&& pl) = default;
        postings_list& operator=(const postings_list& pi) = default;
        postings_list& operator=(postings_list&& pi) = default;
    public:
        postings_list(std::istream& in) {
            load(in);
        }
        template<class t_rank> 
        postings_list(const t_rank& ranker,sdsl::int_vector<>& data,size_t sp,size_t ep);
        template<class t_rank> 
        postings_list(const t_rank& ranker,std::vector<std::pair<uint64_t,uint64_t>>& pre_sorted_data);
        postings_list(std::vector<std::pair<uint64_t,uint64_t>>& pre_sorted_data);
        const_iterator begin() const {
            return const_iterator(*this,0);
        }
        const_iterator end() const {
            return const_iterator(*this,m_size);
        }
        size_t num_blocks() const {
            return m_block_representatives.size();
        }
        size_t size() const {
            return m_size;
        }
        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="") const;
        void load(std::istream& in);
        double max_doc_weight() const {
            return m_max_doc_weight;
        }
        double block_max(uint64_t bid) const {
            return m_block_maximums[bid];
        }
        double block_max_doc_weight(uint64_t bid) const {
            return m_block_max_doc_weights[bid];
        }
        uint64_t block_rep(uint64_t bid) const {
            return m_block_representatives[bid];
        }
        double list_max_score() const {
            return m_list_maximuim;
        }
    private:
        size_t m_size = 0;
        pfor_data_type m_docid_data;
        pfor_data_type m_freq_data;
        sdsl::int_vector<> m_id_block_ptr;
        sdsl::int_vector<> m_freq_block_ptr;
        sdsl::int_vector<> m_block_representatives;
    private: // rank function dependent
        std::vector<double> m_block_maximums;
        std::vector<double> m_block_max_doc_weights;
        double m_list_maximuim = std::numeric_limits<double>::lowest();
        double m_max_doc_weight = std::numeric_limits<double>::lowest();
    private:
        void decompress_block(size_t bid,pfor_data_type& id_data,pfor_data_type& freq_data) const;
        void create_block_support(const sdsl::int_vector<32>& ids,const sdsl::int_vector<32>& freqs);
        template<class t_rank> void create_rank_support(const sdsl::int_vector<32>& ids,const sdsl::int_vector<32>& freqs,double F_t,const t_rank& ranker);
        void compress_postings_data(const sdsl::int_vector<32>& ids,const sdsl::int_vector<32>& freqs);
        size_t find_block_with_id(uint64_t id,size_t start_block = 0) const;
};

template<compression_codec t_codec,uint64_t t_bs>
template<class t_rank>
postings_list<t_codec,t_bs>::postings_list(const t_rank& ranker,std::vector<std::pair<uint64_t,uint64_t>>& pre_sorted_data)
{
    // extract doc_ids and freqs
    sdsl::int_vector<32> tmp_data(pre_sorted_data.size());
    sdsl::int_vector<32> tmp_freq(pre_sorted_data.size());
    for (size_type i=0; i<pre_sorted_data.size(); i++) {
        tmp_data[i] = pre_sorted_data[i].first;
        tmp_freq[i] = pre_sorted_data[i].second;
    }

    // create block max structure first
    create_block_support(tmp_data,tmp_freq);

    // create rank support structure
    create_rank_support(tmp_data,tmp_freq,pre_sorted_data.size(),ranker);

    // compress postings
    compress_postings_data(tmp_data,tmp_freq);
}

template<compression_codec t_codec,uint64_t t_bs>
postings_list<t_codec,t_bs>::postings_list(std::vector<std::pair<uint64_t,uint64_t>>& pre_sorted_data)
{
    // extract doc_ids and freqs
    sdsl::int_vector<32> tmp_data(pre_sorted_data.size());
    sdsl::int_vector<32> tmp_freq(pre_sorted_data.size());
    for (size_type i=0; i<pre_sorted_data.size(); i++) {
        tmp_data[i] = pre_sorted_data[i].first;
        tmp_freq[i] = pre_sorted_data[i].second;
    }

    // create block max structure first
    create_block_support(tmp_data,tmp_freq);

    // compress postings
    compress_postings_data(tmp_data,tmp_freq);
}



template<compression_codec t_codec,uint64_t t_bs>
template<class t_rank>
postings_list<t_codec,t_bs>::postings_list(const t_rank& ranker,
                                           sdsl::int_vector<>& D,
                                           size_t sp,size_t ep)
{
    if (ep<sp) {
        std::cerr << "ERROR: trying to create empty postings list.\n";
        throw std::logic_error("trying to create empty postings list.");
    }

    std::sort(D.begin()+sp,D.begin()+ep);

    // count uniq docs
    size_t unique = 1;
    for(size_t i=sp+1;i<=ep;i++) {
        if(D[i] != D[i-1]) unique++;
    }

    // extract doc_ids and freqs
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
    create_block_support(tmp_data,tmp_freq);

    // create rank support structure
    create_rank_support(tmp_data,tmp_freq,ep-sp+1,ranker);

    // compress postings
    compress_postings_data(tmp_data,tmp_freq);
}

template<compression_codec t_codec,uint64_t t_bs>
auto
postings_list<t_codec,t_bs>::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const -> size_type 
{
    sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written_bytes = 0;

    // id and freq stuff
    written_bytes += sdsl::write_member(m_size,out,child,"size");
    written_bytes += m_id_block_ptr.serialize(out,child,"id ptrs");
    written_bytes += m_freq_block_ptr.serialize(out,child,"freq ptrs");
    written_bytes += sdsl::write_member(m_docid_data.size(),out,child,"docid u32s");
    written_bytes += sdsl::write_member(m_freq_data.size(),out,child,"freq u32s");
    sdsl::structure_tree_node* idchild = sdsl::structure_tree::add_child(child, "id data", compression_codec_names[t_codec]);
    out.write((const char*)m_docid_data.data(), m_docid_data.size()*sizeof(uint32_t));
    sdsl::structure_tree::add_size(idchild, m_docid_data.size()*sizeof(uint32_t));
    written_bytes +=  m_docid_data.size()*sizeof(uint32_t);
    out.write((const char*)m_freq_data.data(), m_freq_data.size()*sizeof(uint32_t));
    written_bytes +=  m_freq_data.size()*sizeof(uint32_t);
    sdsl::structure_tree_node* fchild = sdsl::structure_tree::add_child(child, "freq data", compression_codec_names[t_codec]);
    sdsl::structure_tree::add_size(fchild, m_freq_data.size()*sizeof(uint32_t));

    // block skip stuff
    written_bytes += m_block_representatives.serialize(out,child,"block skip");

    // rank dependent stuff
    sdsl::structure_tree_node* bmchild = sdsl::structure_tree::add_child(child, "blockmax", "blockmax");
    sdsl::structure_tree_node* bmmchild = sdsl::structure_tree::add_child(bmchild, "block maximums", "double");
    size_t bm_written_bytes = sdsl::write_member(m_block_maximums.size(),out,bmchild,"num max_scores");
    out.write((const char*)m_block_maximums.data(), m_block_maximums.size()*sizeof(double));
    out.write((const char*)m_block_max_doc_weights.data(), m_block_max_doc_weights.size()*sizeof(double));
    bm_written_bytes += 2*m_block_maximums.size()*sizeof(double);
    sdsl::structure_tree::add_size(bmmchild,2*m_block_maximums.size()*sizeof(double));
    sdsl::structure_tree::add_size(bmchild,bm_written_bytes);
    written_bytes += bm_written_bytes;
    written_bytes += sdsl::write_member(m_list_maximuim,out,child,"list max score");
    written_bytes += sdsl::write_member(m_max_doc_weight,out,child,"max doc weight");

    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<compression_codec t_codec,uint64_t t_bs>
void postings_list<t_codec,t_bs>::load(std::istream& in)
{
    read_member(m_size,in);
    m_id_block_ptr.load(in);
    m_freq_block_ptr.load(in);
    size_t docidu32;
    size_t frequ32;
    read_member(docidu32,in);
    read_member(frequ32,in);
    m_docid_data.resize(docidu32);
    m_freq_data.resize(frequ32);
    in.read((char*)m_docid_data.data(),docidu32*sizeof(uint32_t));
    in.read((char*)m_freq_data.data(),frequ32*sizeof(uint32_t));
    m_block_representatives.load(in);
    size_t num_block_max_scores;
    read_member(num_block_max_scores,in);
    m_block_maximums.resize(num_block_max_scores);
    in.read((char*)m_block_maximums.data(),num_block_max_scores*sizeof(double));
    m_block_max_doc_weights.resize(num_block_max_scores);
    in.read((char*)m_block_max_doc_weights.data(),num_block_max_scores*sizeof(double));
    read_member(m_list_maximuim,in);
    read_member(m_max_doc_weight,in);
}


template<compression_codec t_codec,uint64_t t_bs>
template<class t_rank>
void postings_list<t_codec,t_bs>::create_rank_support(const sdsl::int_vector<32>& ids,const sdsl::int_vector<32>& freqs,double F_t,const t_rank& ranker)
{
    size_t num_blocks = ids.size() / t_bs;
    if (ids.size() % t_bs != 0) num_blocks++;

    m_block_maximums.resize(num_blocks);
    m_block_max_doc_weights.resize(num_blocks);
    double max_score = 0.0f;
    double max_doc_weight = std::numeric_limits<double>::lowest();
    double f_t = ids.size();
    size_t i = 1;
    size_t j = 0;
    m_max_doc_weight = std::numeric_limits<double>::lowest();
    m_list_maximuim = std::numeric_limits<double>::lowest();
    for (size_t l=0; l<ids.size(); l++) {
        auto id = ids[l];
        auto f_dt = freqs[l];
        double W_d = ranker.doc_length(id);
        double doc_weight = ranker.calc_doc_weight(W_d);
        double score = ranker.calculate_docscore(1.0f,f_dt,f_t,F_t,W_d);
        max_score = std::max(max_score,score);
        max_doc_weight = std::max(max_doc_weight,doc_weight);
        if (i % t_bs == 0) {
            m_block_maximums[j] = max_score;
            m_block_max_doc_weights[j] = max_doc_weight;
            m_list_maximuim = std::max(m_list_maximuim,max_score);
            m_max_doc_weight = std::max(m_max_doc_weight,max_doc_weight);
            i = 0;
            max_score = 0.0f;
            max_doc_weight = std::numeric_limits<double>::lowest();
            j++;
        }
        i++;
    }
    if (ids.size() % t_bs != 0) {
        m_block_maximums[num_blocks-1] = max_score;
        m_block_max_doc_weights[num_blocks-1] = max_doc_weight;
    }
    m_list_maximuim = std::max(m_list_maximuim,max_score);
    m_max_doc_weight = std::max(m_max_doc_weight,max_doc_weight);
}


template<compression_codec t_codec,uint64_t t_bs>
void postings_list<t_codec,t_bs>::create_block_support(const sdsl::int_vector<32>& ids,const sdsl::int_vector<32>& freqs)
{
    size_t num_blocks = ids.size() / t_bs;
    if (ids.size() % t_bs != 0) num_blocks++;

    m_block_representatives.resize(num_blocks);
    size_t j = 0;
    for (size_t i=t_bs-1; i<ids.size(); i+=t_bs) {
        m_block_representatives[j++] = ids[i];
    }
    m_block_representatives[j] = ids[ids.size()-1];
    sdsl::util::bit_compress(m_block_representatives);
}

template<compression_codec t_codec,uint64_t t_bs>
void postings_list<t_codec,t_bs>::compress_postings_data(const sdsl::int_vector<32>& ids,const sdsl::int_vector<32>& freqs)
{
    // encode in blocks
    size_t num_blocks = ids.size() / t_bs;
    if (ids.size() % t_bs != 0) num_blocks++;
    m_id_block_ptr.resize(num_blocks+1);
    m_freq_block_ptr.resize(num_blocks+1);
    m_size = ids.size();

    // encode using pfor
    static std::shared_ptr<FastPFor::IntegerCODEC> codec =
        FastPFor::CODECFactory::getFromName(compression_codec_names[t_codec]);
    FastPFor::IntegerCODEC& c = *codec;

    // encode ids
    size_t id_start = 0;
    m_docid_data.resize(2 * ids.size() + 1024);

    // delta compress ids first
    uint32_t* input = (uint32_t*) ids.data();
    for (size_t i=0; i<ids.size(); i+=t_bs) {
        size_t n = t_bs;
        if (i+t_bs >= ids.size()) n = ids.size()-i;
        FastPFor::Delta::fastDelta(&input[i],n);
    }

    // pfor compress ids
    uint32_t* id_data = m_docid_data.data();
    size_t j = 0;
    for (size_t i=0; i<ids.size(); i+=t_bs) {
        m_id_block_ptr[j] = id_start;
        size_t n = t_bs;
        if (i+t_bs >= ids.size()) n = ids.size()-i;
        size_t encoded_size = m_docid_data.size() - id_start; // how much space we have left
        c.encodeArray(&input[i],n,&id_data[id_start],encoded_size);
        id_start += encoded_size;
        j++;
    }
    m_id_block_ptr[j] = id_start;
    m_docid_data.resize(id_start);

    // encode freqs
    size_t freq_start = 0;
    m_freq_data.resize(2 * freqs.size() + 1024);
    uint32_t* freq_data = m_freq_data.data();
    uint32_t* finput = (uint32_t*) freqs.data();
    j = 0;
    for (size_t i=0; i<freqs.size(); i+=t_bs) {
        m_freq_block_ptr[j] = freq_start;
        size_t n = t_bs;
        if (i+t_bs >= freqs.size()) n = freqs.size()-i;
        size_t encoded_size =  m_freq_data.size() - id_start; // how much space we have left
        c.encodeArray(&finput[i],n,&freq_data[freq_start],encoded_size);
        freq_start += encoded_size;
        j++;
    }
    m_freq_block_ptr[j] = freq_start;
    m_freq_data.resize(freq_start);

    sdsl::util::bit_compress(m_id_block_ptr);
    sdsl::util::bit_compress(m_freq_block_ptr);
}

template<compression_codec t_codec,uint64_t t_bs>
void postings_list<t_codec,t_bs>::decompress_block(size_t bid,pfor_data_type& id_data,pfor_data_type& freq_data) const
{
    if (id_data.size() != t_bs) { // did we allocate space already?
        id_data.resize(t_bs);
        freq_data.resize(t_bs);
    }
    static std::shared_ptr<FastPFor::IntegerCODEC> codec =
        FastPFor::CODECFactory::getFromName(compression_codec_names[t_codec]);
    FastPFor::IntegerCODEC& c = *codec;
    size_t block_start = m_id_block_ptr[bid];
    size_t block_stop = m_id_block_ptr[bid+1];
    size_t compressed_block_size = block_stop-block_start;

    // decode doc ids
    uint32_t* aligned_docids = (uint32_t*)(m_docid_data.data()+block_start);
    size_t num_recovered_ids = id_data.size();
    uint32_t* result_ids = id_data.data();
    c.decodeArray(aligned_docids,compressed_block_size,result_ids,num_recovered_ids);
    FastPFor::Delta::fastinverseDelta(result_ids,num_recovered_ids);

    // decode freqs (no d-gap as non-increasing)
    block_start = m_freq_block_ptr[bid];
    block_stop = m_freq_block_ptr[bid+1];
    compressed_block_size = block_stop-block_start;
    size_t num_recovered_freqs = freq_data.size();
    uint32_t* aligned_freqs = (uint32_t*)(m_freq_data.data()+block_start);
    uint32_t* result_freqs = freq_data.data();
    c.decodeArray(aligned_freqs , compressed_block_size ,result_freqs,num_recovered_freqs);

    if (num_recovered_ids != num_recovered_freqs) {
        std::cerr << "ERROR: number of decoded ids and freqs is not equal. "
                  << num_recovered_ids << " != " << num_recovered_freqs << "\n";
        throw std::logic_error("number of decoded ids and freqs is not equal.");
    }
    if (num_recovered_ids != t_bs) {
        freq_data.resize(num_recovered_ids);
        id_data.resize(num_recovered_ids);
    }
}

template<compression_codec t_codec,uint64_t t_bs>
size_t postings_list<t_codec,t_bs>::find_block_with_id(uint64_t id,size_t start_block) const
{
    size_t block_id = start_block;
    size_t nblocks = m_block_representatives.size();
    while (block_id < nblocks && m_block_representatives[block_id] < id) {
        block_id++;
    }
    return block_id;
}

template<compression_codec t_codec,uint64_t t_bs>
plist_iterator<t_codec,t_bs>::plist_iterator(const list_type& l,size_t pos) : plist_iterator()
{
    m_cur_pos = pos;
    m_plist_ptr = &l;
}

template<compression_codec t_codec,uint64_t t_bs>
plist_iterator<t_codec,t_bs>& plist_iterator<t_codec,t_bs>::operator++()
{
    if (m_cur_pos != size()) { // end?
        (*this).m_cur_pos++;
    } else {
        std::cerr << "ERROR: trying to advance plist iterator beyond list end.\n";
        throw std::out_of_range("trying to advance plist iterator beyond list end");
    }
    return (*this);
}

template<compression_codec t_codec,uint64_t t_bs>
bool plist_iterator<t_codec,t_bs>::operator ==(const plist_iterator& b) const
{
    return ((*this).m_cur_pos == b.m_cur_pos) && ((*this).m_plist_ptr == b.m_plist_ptr);
}

template<compression_codec t_codec,uint64_t t_bs>
bool plist_iterator<t_codec,t_bs>::operator !=(const plist_iterator& b) const
{
    return !((*this)==b);
}

template<compression_codec t_codec,uint64_t t_bs>
typename plist_iterator<t_codec,t_bs>::value_type plist_iterator<t_codec,t_bs>::docid() const
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

template<compression_codec t_codec,uint64_t t_bs>
typename plist_iterator<t_codec,t_bs>::value_type plist_iterator<t_codec,t_bs>::freq() const
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

template<compression_codec t_codec,uint64_t t_bs>
void plist_iterator<t_codec,t_bs>::access_and_decode_cur_pos() const
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

template<compression_codec t_codec,uint64_t t_bs>
void plist_iterator<t_codec,t_bs>::skip_to_block_with_id(uint64_t id)
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

template<compression_codec t_codec,uint64_t t_bs>
void plist_iterator<t_codec,t_bs>::skip_to_id(uint64_t id)
{
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

} // end surf namespace

#endif
