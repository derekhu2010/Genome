#ifndef TRIE_INCLUDED
#define TRIE_INCLUDED

#include <string>
#include <vector>
#include <cctype>
#include <iostream>
using namespace std;

template<typename ValueType>
class Trie
{
public:
    Trie();
    ~Trie();
    void reset();
    void insert(const std::string& key, const ValueType& value);
    std::vector<ValueType> find(const std::string& key, bool exactMatchOnly) const;
    
    // C++11 syntax for preventing copying and assignment
    Trie(const Trie&) = delete;
    Trie& operator=(const Trie&) = delete;
    
private:
    struct TrieNode
    {
    public:
        vector<ValueType> m_values;
        vector<char> m_labels;
        vector<TrieNode*> m_children;
    };
    
    TrieNode* m_root;
    void resetHelper(TrieNode* start);
    std::vector<ValueType> findHelper(const std::string& key, bool exactMatchOnly, TrieNode* start, int nMismatches) const;
};

template<typename ValueType>
Trie<ValueType>::Trie()
{
    m_root = new TrieNode;
}

template<typename ValueType>
Trie<ValueType>::~Trie()
{
    reset();
    delete m_root;
}

template<typename ValueType>
void Trie<ValueType>::reset()
{
    resetHelper(m_root);
    m_root = new TrieNode;
}

template<typename ValueType>
void Trie<ValueType>::insert(const std::string &key, const ValueType &value)
{
    TrieNode* cur = m_root;
    for (int i = 0; i < key.size(); i++)
    {
        bool charfound = false;
        for(int k = 0; k < cur->m_labels.size(); k++)
        {
            if(key[i] == cur->m_labels[k])
            {
                charfound = true;
                cur = cur->m_children[k];
                break;
            }
        }
        if(charfound)      //if char is found, move on to the next character
            continue;
        else
        {
            cur->m_labels.push_back(key[i]);
            cur->m_children.push_back(new TrieNode);
            cur = cur->m_children.back();
        }
    }
    cur->m_values.push_back(value);
}

template<typename ValueType>
void Trie<ValueType>::resetHelper(TrieNode *start)
{
    if(start->m_children.size() == 0)
    {
        delete start;
        return;
    }
    while( !start->m_children.empty())
    {
        resetHelper(start->m_children.back());
        start->m_children.pop_back();
    }
    delete start;
}

template<typename ValueType>
std::vector<ValueType> Trie<ValueType>::find(const std::string& key, bool exactMatchOnly) const
{
    TrieNode* start = nullptr;
    vector<ValueType> t;
    for (int i = 0; i < m_root->m_labels.size(); i++)      //this loop finds the child pointer that corresponds with the first letter
    {
        if(m_root->m_labels[i] == key[0])
            start = m_root->m_children[i];
    }
    if (start == nullptr)      //first letter is not found in the root node, return an empty vector
        return t;
    t = findHelper(key.substr(1), exactMatchOnly, start, 0);
    return t;
}

template<typename ValueType>
std::vector<ValueType> Trie<ValueType>::findHelper(const std::string& key, bool exactMatchOnly, TrieNode* start, int nMismatches) const
{
    vector<ValueType> t, temp;
    int startIndex;
    if((exactMatchOnly && nMismatches > 0) || (!exactMatchOnly && nMismatches > 1) )       //break the recursion, return an empty vector, do not process anymore nodes if this is met.
        return t;
    if(key.size() == 0)       //traversed through all the letters, get the values and return
    {
        for(int i = 0; i < start->m_values.size(); i++)
            t.push_back(start->m_values[i]);
        return t;
    }
    for(startIndex = 0; startIndex < start->m_children.size(); startIndex++)
    {
        if(start->m_labels[startIndex] == key[0])      //no matter what, if the letter matches a child, we want to add the values to our final answer
            temp = findHelper(key.substr(1), exactMatchOnly, start->m_children[startIndex], nMismatches);
        else if(!exactMatchOnly)       //we only want to add mismatches if exactMatchOnly is false
            temp = findHelper(key.substr(1), exactMatchOnly, start->m_children[startIndex], nMismatches+1);      //increment nMismatches by 1; only allowed 1 mismatch
        
        for(int i = 0; i < temp.size(); i++)
            t.push_back(temp[i]);
        temp.clear();
    }
    return t;
}
#endif // TRIE_INCLUDED
