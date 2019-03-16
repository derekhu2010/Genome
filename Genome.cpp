#include "provided.h"
#include <string>
#include <vector>
#include <iostream>
#include <istream>
#include <fstream>
#include <cctype>
using namespace std;

class GenomeImpl
{
public:
    GenomeImpl(const string& nm, const string& sequence);
    static bool load(istream& genomeSource, vector<Genome>& genomes);
    int length() const;
    string name() const;
    bool extract(int position, int length, string& fragment) const;
private:
    string m_name;
    string m_sequence;
};

GenomeImpl::GenomeImpl(const string& nm, const string& sequence)
{
    m_name = nm;
    m_sequence = sequence;
}

bool GenomeImpl::load(istream& genomeSource, vector<Genome>& genomes)
{
    string s;       //will store a line
    string sequence = "";        //will store a line if it is a sequence
    string name = "";       //will store a line if it is a name;
    bool hasBaseLine = false;
    while(getline(genomeSource, s))
    {
        if(s.empty())
            return false;
        if(s[0] == '>')      //the line is a name line
        {
            if(name != "" && hasBaseLine)       //insert previous genome if there is a valid one
            {
                genomes.push_back(Genome(name, sequence));
                sequence = "";
                hasBaseLine = false;
            }
            else if(name != "" && !hasBaseLine)       //2 name lines in a row
                return false;
            if(s.size() <= 1)      //only contains >
                return false;
            name = s.substr(1);
        }
        else     //the line must be a subsequence of DNA sequence
        {
            for (int i = 0; i < s.size(); i++)
            {
                s[i] = toupper(s[i]);
                switch (s[i])
                {
                    case 'A':
                    case 'C':
                    case 'T':
                    case 'G':
                    case 'N':
                        sequence += s[i];
                        break;
                    default:      //invalid DNA base return false
                        return false;
                }
            }
            if(name == "")     //a name was not inputted, thus the first line read was a DNA sequence
                return false;
            hasBaseLine = true;
        }
    }
    if(name != "" && hasBaseLine)       //inserts the last genome into the vector
    {
        genomes.push_back(Genome(name, sequence));
        return true;
    }
    return false;      //last line is a name followed by no baseline
}

int GenomeImpl::length() const
{
    return static_cast<int>(m_sequence.length());
}

string GenomeImpl::name() const
{
    return m_name;  // This compiles, but may not be correct
}

bool GenomeImpl::extract(int position, int length, string& fragment) const
{
    int l = this->length();
    if((position + length) > l)
        return false;
    fragment = m_sequence.substr(position, length);
    return true;
}

//******************** Genome functions ************************************

// These functions simply delegate to GenomeImpl's functions.
// You probably don't want to change any of this code.

Genome::Genome(const string& nm, const string& sequence)
{
    m_impl = new GenomeImpl(nm, sequence);
}

Genome::~Genome()
{
    delete m_impl;
}

Genome::Genome(const Genome& other)
{
    m_impl = new GenomeImpl(*other.m_impl);
}

Genome& Genome::operator=(const Genome& rhs)
{
    GenomeImpl* newImpl = new GenomeImpl(*rhs.m_impl);
    delete m_impl;
    m_impl = newImpl;
    return *this;
}

bool Genome::load(istream& genomeSource, vector<Genome>& genomes)
{
    return GenomeImpl::load(genomeSource, genomes);
}

int Genome::length() const
{
    return m_impl->length();
}

string Genome::name() const
{
    return m_impl->name();
}

bool Genome::extract(int position, int length, string& fragment) const
{
    return m_impl->extract(position, length, fragment);
}
