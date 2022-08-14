// UMBC - CMSC 341 - Spring 2022 - Proj2
// Author: Sriram Vema
// Date: 5/10/2022
// Section:5
// File: dnadb.cpp
// Description: dnadb file that holds dnadb class functions
#include "dnadb.h"
DnaDb::DnaDb(int size, hash_fn hash){
    if (size < MINPRIME) {
        size = MINPRIME;                   // Makes sure the size is within bounds
    } else if (size > MAXPRIME) {
        size = MAXPRIME;
    } else {
        if (!isPrime(size)) {
            size = findNextPrime(size);
        }
    }
    m_hash = hash;
    m_currentCap = size; // Sets all private members
    m_currentSize = 0;
    m_currentTable = new DNA[size]; // sets the current table
    m_currNumDeleted = 0;
    m_oldTable = nullptr;
    m_oldCap = 0;
    m_oldNumDeleted = 0;
    m_oldSize = 0;
    for (int i = 0; i < m_currentCap; i++) { // sets current table to empty objects
        m_currentTable[i] = EMPTY;
    }
}

DnaDb::~DnaDb(){
    delete [] m_currentTable;
    if (m_oldTable != nullptr) {  // deleted both tables
        delete [] m_oldTable;
    }
    m_currentTable = nullptr; // sets private members to base values
    m_oldTable = nullptr;
    m_oldCap = 0;
    m_oldNumDeleted = 0;
    m_oldSize = 0;
    m_currNumDeleted = 0;
    m_currentSize = 0;
    m_currentCap = 0;
}

bool DnaDb::checkExists(string sequence, int location) {
    for (int i = 0; i <  m_currentCap; i++) {   // checks to see if any of the dna objects are duplicates
        if (m_currentTable[i].getSequence() == sequence && m_currentTable[i].getLocId() == location) {
            return true;
        }
    }
    return false;
}

bool DnaDb::insert(DNA dna){
    bool result = false;
    if (checkExists(dna.m_sequence, dna.m_location)) {  // doesnt insert if it is a duplicate
        return false;
    }
    int index = m_hash(dna.m_sequence) % m_currentCap;  // calculates the first index
    if (m_currentSize == 0) {   // adds the object straight away if the table is empty
        m_currentTable[index] = dna;
        m_currentSize ++;
        rehash();
        result = true;
    } else {
        if (m_currentTable[index].getSequence() != "") {        // if the index is occupied
            int i = 0;
            while(m_currentTable[index].getSequence() != "") {          // while loop keeps running quadratic probing until it finds an index that isn't occupied
                index = ((m_hash(dna.m_sequence) % m_currentCap) + (i*i)) % m_currentCap;
                i++;
            }
        }
        m_currentTable[index] = dna; // adds the the object to the final index
        m_currentSize ++;
        rehash();
        result = true;
    } 
    rehash();      
    return result;
}

void DnaDb::rehash(){
    if (lambda() > 0.5 || deletedRatio() > 0.8) { // if either of these conditions are met and oldtable is empty, it copies the values over
        if (m_oldTable == nullptr) {
            m_oldCap = m_currentCap;
            m_oldSize = m_currentSize;
            m_oldNumDeleted = m_currNumDeleted;
            m_oldTable = new DNA[m_oldCap];
            for (int i = 0; i < m_oldCap; i++) {
                m_oldTable[i] = m_currentTable[i];
            }
            m_currentCap = findNextPrime((m_oldCap * 4));
            m_currentSize = 0;
            m_currNumDeleted = 0;
            delete [] m_currentTable;
            m_currentTable = new DNA[m_currentCap];
        }
    }
    if(m_oldTable != nullptr){       // if oldtable is filled it moves 25% of data over
        int hashedInsert = 0;
        int i = 0;
        while (hashedInsert <= (m_oldSize/4) && i < m_oldCap) {
            if(m_oldTable[i].getSequence() != "" && m_oldTable[i].getSequence() != DELETEDKEY) {
                currentInsert(m_oldTable[i]);
                m_oldTable[i] = DELETED;
                m_oldNumDeleted++;
                hashedInsert++;
            }
            i++;
        }
    }
    bool filled = false;
    for (int i = 0; i < m_oldCap; i++) {  // checks to see if old table has any non deleted strands
        if (m_oldTable[i].getSequence() != "" && m_oldTable[i].getSequence() != "DELETED") {
            filled = true;
        }
    }
    if (m_oldNumDeleted == m_oldSize) {  // makes oldtable nullptr again if everything has been rehashed
        delete [] m_oldTable;
        m_oldCap = 0;
        m_oldNumDeleted = 0;
        m_oldSize = 0;
        m_oldTable = nullptr;
    }
}

bool DnaDb::remove(DNA dna) {                          
    int index = m_hash(dna.m_sequence) % m_currentCap;            // first checks the hash index of the current and old tables
    if (m_currentTable[index].m_location == dna.m_location) {
        m_currentTable[index] = DELETED;
        m_currNumDeleted ++;
        rehash();
        return true;
    } else if (m_oldTable != nullptr) {
        if (m_oldTable[index].m_location == dna.m_location) {
            m_oldTable[index] = DELETED;
            m_oldNumDeleted++;
            rehash();
            return true;            
        }
    } else {                                        //if the object isnt found it goes through the array to find it
        for (int i = 0; i < m_currentCap; i++) {
            if (m_currentTable[i].m_location == dna.m_location) {
                m_currentTable[i] = DELETED;    // sets the object to deleted
                m_currNumDeleted++;
                rehash();    // rehashes after every remove
                return true;
            }
        }
        if (m_oldTable != nullptr) {
            for (int n = 0; n < m_oldCap; n++) {
                if (m_oldTable[n].m_location == dna.m_location) {
                    m_oldTable[n] = DELETED;
                    m_oldNumDeleted++;
                    rehash();
                    return true;
                }
            }            
        }

    }
    rehash();
    return false;
}

DNA DnaDb::getDNA(string sequence, int location){
    if (m_currentSize > 0) {                      // iterates through the arrays to find the dna strand
        for (int i = 0; i < m_currentCap; i++) {
            if (m_currentTable[i].m_location == location) {
                return m_currentTable[i];
            }
        }
    }
    if (m_oldTable != nullptr) {
        for (int i = 0; i < m_oldCap; i++) {
            if (m_oldTable[i].m_location == location) {
                return m_oldTable[i];
            }
        }
    }
    return EMPTY;
}

float DnaDb::lambda() const {
    return (float) m_currentSize/m_currentCap;
}

float DnaDb::deletedRatio() const {
    return (float) m_currNumDeleted/m_currentSize;
}

void DnaDb::dump() const {
    cout << "Dump for current table: " << endl;
    if (m_currentTable != nullptr)
        for (int i = 0; i < m_currentCap; i++) {
            cout << "[" << i << "] : " << m_currentTable[i] << endl;
        }
    cout << "Dump for old table: " << endl;
    if (m_oldTable != nullptr)
        for (int i = 0; i < m_oldCap; i++) {
            cout << "[" << i << "] : " << m_oldTable[i] << endl;
        }
}

bool DnaDb::isPrime(int number){
    bool result = true;
    for (int i = 2; i <= number / 2; ++i) {
        if (number % i == 0) {
            result = false;
            break;
        }
    }
    return result;
}

int DnaDb::findNextPrime(int current){
    //we always stay within the range [MINPRIME-MAXPRIME]
    //the smallest prime starts at MINPRIME
    if (current < MINPRIME) current = MINPRIME-1;
    for (int i=current; i<MAXPRIME; i++) { 
        for (int j=2; j*j<=i; j++) {
            if (i % j == 0) 
                break;
            else if (j+1 > sqrt(i) && i != current) {
                return i;
            }
        }
    }
    //if a user tries to go over MAXPRIME
    return MAXPRIME;
}

DNA::DNA(string sequence, int location) {
    if ((location >= MINLOCID && location <= MAXLOCID) ||
        (location == 0 && sequence == "DELETED")){
        // this is a normal or a DELETED object
        m_sequence = sequence;
        m_location = location;
    }
    else{
        // this is the empty object
        m_sequence = "";
        m_location = 0;
    }
}

string DNA::getSequence() const {
    return m_sequence;
}

int DNA::getLocId() const {
    return m_location;
}

// Overloaded assignment operator
const DNA& DNA::operator=(const DNA& rhs){
    if (this != &rhs){
        m_sequence = rhs.m_sequence;
        m_location = rhs.m_location;
    }
    return *this;
}

// Overloaded insertion operator.  Prints DNA's sequence (key),
// and the location ID. This is a friend function in DNA class.
ostream& operator<<(ostream& sout, const DNA &dna ) {
    if (!dna.m_sequence.empty())
        sout << dna.m_sequence << " (Location ID " << dna.m_location << ")";
    else
        sout << "";
  return sout;
}

// Overloaded equality operator. This is a friend function in DNA class.
// To test inequality we may negate the results of this operator.
bool operator==(const DNA& lhs, const DNA& rhs){
    return ((lhs.m_sequence == rhs.m_sequence) && (lhs.m_location == rhs.m_location));
}



void DnaDb::currentInsert(DNA dna) {        // inserts into the current array without rehashing
    bool result = false;                   // only called in rehash
    int index = m_hash(dna.m_sequence) % m_currentCap;
    if (m_currentSize == 0) {
        m_currentTable[index] = dna;
        m_currentSize ++;
        result = true;
    } else {
        if (m_currentTable[index].getSequence() != "") {
            int i = 0;
            while(m_currentTable[index].getSequence() != "") {
                index = ((m_hash(dna.m_sequence) % m_currentCap) + (i*i)) % m_currentCap;
                i++;
            }
        }
        m_currentTable[index] = dna;
        m_currentSize ++;
        result = true;
    } 
}

