#include "provided.h"
#include "Trie.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
using namespace std;


class GenomeMatcherImpl
{
public:
    GenomeMatcherImpl(int minSearchLength);
    void addGenome(const Genome& genome);
    int minimumSearchLength() const;
    bool findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const;
    bool findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const;
private:
    struct DNA_position
    {
        DNA_position(string name, int pos, int index) :m_genome_name(name), m_pos(pos), m_index(index) {}
        string m_genome_name;
        int m_pos;
        int m_index;
    };
    static bool percentComp(GenomeMatch i, GenomeMatch j)      //used to sort the results in findRelatedGenomes
        {
            if(i.percentMatch > j.percentMatch)
                return true;
            else if(j.percentMatch > i.percentMatch)
                return false;
            else      //the percentages are equal, so return true if i comes before j alphabetically
                return i.genomeName < j.genomeName;
        }
    int m_minSearchLength;
    vector<Genome> m_genomes;
    Trie<DNA_position> m_DNA_positions;
};

GenomeMatcherImpl::GenomeMatcherImpl(int minSearchLength)
{
    m_minSearchLength = minSearchLength;
}

void GenomeMatcherImpl::addGenome(const Genome& genome)
{
    m_genomes.push_back(genome);
    for(int i = 0; i <= genome.length()-minimumSearchLength(); i++)
    {
        string DNA;
        genome.extract(i, minimumSearchLength(), DNA);
        m_DNA_positions.insert(DNA, DNA_position(genome.name(), i, static_cast<int>(m_genomes.size())-1));       //this genome is in the back of the vector so its index is size-1
    }
}

int GenomeMatcherImpl::minimumSearchLength() const
{
    return m_minSearchLength;
}

bool GenomeMatcherImpl::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    //minimumLength is the minimum length for each match, not the search
    if(fragment.length() < minimumLength || minimumLength < minimumSearchLength())
        return false;
    vector<DNA_position> hits = m_DNA_positions.find(fragment.substr(0, minimumSearchLength()), exactMatchOnly);      //find all hits
    vector<int> correspondingIndex (m_genomes.size(), -1);     //Indices line up w/ m_genome. Will store the location of each genome in the matches vector
                                                                //avoids having to loop through the matches vector to search for the correpsonding match to a genome
    vector<DNAMatch> m;       //will store DNAMatches
    for(int i = 0; i < hits.size(); i++)
    {
        int matchLen = 0;      //count the length of current match corresponding with the ith hit
        int misMatches = 0;    //count the number of mismatches corresponding with the ith hit and check it w exactMatchOnly
        int k = 0;      //counter for fragment position
        string genomeBase;      //will store base of genome to compare to fragment
        while(k < fragment.length())      //iterate through the characters of fragment
        {
            m_genomes[hits[i].m_index].extract(hits[i].m_pos + k, 1, genomeBase);
            if(genomeBase == fragment.substr(k, 1))
                matchLen++;
            else
            {
                misMatches++;
                if( (exactMatchOnly && misMatches > 0) || (!exactMatchOnly && misMatches >1) )
                    break;
                else      //have not exceded the max number of misMatches, the length of the match increases
                    matchLen++;
            }
            k++;
        }
        
        if(matchLen < minimumLength)     //match is not long enough, do not need to process this match
            continue;
        
        int indexOfMatch = correspondingIndex[hits[i].m_index];       //corresponding index of the for this genome in the match vector
        if(indexOfMatch == -1)      // no match has been found for this genome yet, so push this match
        {
            m.push_back(DNAMatch());
            m.back().genomeName = hits[i].m_genome_name;
            m.back().position = hits[i].m_pos;
            m.back().length = matchLen;
            correspondingIndex[hits[i].m_index] = static_cast<int>(m.size())-1;       //update where the corresponding index of this genome is in the match vector
        }
        else        //there is an existing match for this genome, we want to compare them and update accordingly
        {
            if(matchLen > m[indexOfMatch].length)
            {
                m[indexOfMatch].length = matchLen;
                m[indexOfMatch].position = hits[i].m_pos;
            }
        }
    }
    if(m.empty())
        return false;
    matches = m;
    return true;
}

bool GenomeMatcherImpl::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    if(fragmentMatchLength < minimumSearchLength())
        return false;
    int nFragments = query.length()/fragmentMatchLength;
    string queryFrag;
    vector<DNAMatch> DNAmatches;
    vector<GenomeMatch> genomeMatches;
    vector<int> nMatches (m_genomes.size(), 0);      //indexes line up with m_genomes
    for(int i = 0; i < nFragments; i++)
    {
        query.extract(i*fragmentMatchLength, fragmentMatchLength, queryFrag);
        findGenomesWithThisDNA(queryFrag, fragmentMatchLength, exactMatchOnly, DNAmatches);
        for(int k = 0; k < DNAmatches.size(); k++)       //iterates through matches
        {
            for(int j = 0; j < m_genomes.size(); j++)      //iterates through genomes to find index corresponding genome to match
            {
                if(DNAmatches[k].genomeName == m_genomes[j].name())
                {
                    nMatches[j]++;      //increment the number of matches that that genome has
                    break;
                }
            }
        }
    }
    
    for(int i = 0; i < nMatches.size(); i++)
    {
        double matchPct = (static_cast<double>(nMatches[i])/nFragments) * static_cast<double>(100);
        if(matchPct > matchPercentThreshold)
        {
            genomeMatches.push_back(GenomeMatch());
            genomeMatches.back().genomeName = m_genomes[i].name();
            genomeMatches.back().percentMatch = matchPct;
        }
    }
    sort(genomeMatches.begin(), genomeMatches.end(), percentComp);
    if (genomeMatches.empty())
        return false;
    results = genomeMatches;
    return true;
}

//******************** GenomeMatcher functions ********************************

// These functions simply delegate to GenomeMatcherImpl's functions.
// You probably don't want to change any of this code.

GenomeMatcher::GenomeMatcher(int minSearchLength)
{
    m_impl = new GenomeMatcherImpl(minSearchLength);
}

GenomeMatcher::~GenomeMatcher()
{
    delete m_impl;
}

void GenomeMatcher::addGenome(const Genome& genome)
{
    m_impl->addGenome(genome);
}

int GenomeMatcher::minimumSearchLength() const
{
    return m_impl->minimumSearchLength();
}

bool GenomeMatcher::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    return m_impl->findGenomesWithThisDNA(fragment, minimumLength, exactMatchOnly, matches);
}

bool GenomeMatcher::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    return m_impl->findRelatedGenomes(query, fragmentMatchLength, exactMatchOnly, matchPercentThreshold, results);
}

