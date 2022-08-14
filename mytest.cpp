// UMBC - CMSC 341 - Spring 2022 - Proj2
// Author: Sriram Vema
// Date: 5/10/2022
// Section:5
// File: mytest.cpp
// Description: mytest file that holds Tester class functions
#include "dnadb.h"
#include <random>
#include <vector>
enum RANDOM {UNIFORMINT, UNIFORMREAL, NORMAL};
class Random {
public:
    Random(int min, int max, RANDOM type=UNIFORMINT, int mean=50, int stdev=20) : m_min(min), m_max(max), m_type(type)
    {
        if (type == NORMAL){
            //the case of NORMAL to generate integer numbers with normal distribution
            m_generator = std::mt19937(m_device());
            //the data set will have the mean of 50 (default) and standard deviation of 20 (default)
            //the mean and standard deviation can change by passing new values to constructor 
            m_normdist = std::normal_distribution<>(mean,stdev);
        }
        else if (type == UNIFORMINT) {
            //the case of UNIFORMINT to generate integer numbers
            // Using a fixed seed value generates always the same sequence
            // of pseudorandom numbers, e.g. reproducing scientific experiments
            // here it helps us with testing since the same sequence repeats
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_unidist = std::uniform_int_distribution<>(min,max);
        }
        else{ //the case of UNIFORMREAL to generate real numbers
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_uniReal = std::uniform_real_distribution<double>((double)min,(double)max);
        }
    }
    void setSeed(int seedNum){
        // we have set a default value for seed in constructor
        // we can change the seed by calling this function after constructor call
        // this gives us more randomness
        m_generator = std::mt19937(seedNum);
    }

    int getRandNum(){
        // this function returns integer numbers
        // the object must have been initialized to generate integers
        int result = 0;
        if(m_type == NORMAL){
            //returns a random number in a set with normal distribution
            //we limit random numbers by the min and max values
            result = m_min - 1;
            while(result < m_min || result > m_max)
                result = m_normdist(m_generator);
        }
        else if (m_type == UNIFORMINT){
            //this will generate a random number between min and max values
            result = m_unidist(m_generator);
        }
        return result;
    }

    double getRealRandNum(){
        // this function returns real numbers
        // the object must have been initialized to generate real numbers
        double result = m_uniReal(m_generator);
        // a trick to return numbers only with two deciaml points
        // for example if result is 15.0378, function returns 15.03
        // to round up we can use ceil function instead of floor
        result = std::floor(result*100.0)/100.0;
        return result;
    }
    
    private:
    int m_min;
    int m_max;
    RANDOM m_type;
    std::random_device m_device;
    std::mt19937 m_generator;
    std::normal_distribution<> m_normdist;//normal distribution
    std::uniform_int_distribution<> m_unidist;//integer uniform distribution
    std::uniform_real_distribution<double> m_uniReal;//real uniform distribution

};
class Tester{
    public:
    bool insertNormal();   // makes sure everything is inserted in the proper index
    bool getDNAerror();    // makes sure getDNA returns false if the dna doesn't exist
    bool getDNAedge();     // sees if getDNA works without needing to rehash
    bool getDNAnormal();   // sees if getDNA works after rehashing is called
    bool removeNormal();   // sees if remove works properly
    bool removeEdge();     // sees if remove works without rehashing
    bool removeError();    // sees if remove returns false if the strand doesn't exist
    bool rehashComplete();  // tests if rehashing is done completely and old table is cleared
    bool rehashRemove();   // tests that remove calls rehash
    bool rehashDeletedRatio();  // tests that rehashing is called when deleted ratio is exceeded
};

unsigned int hashCode(const string str);
string sequencer(int size, int seedNum);

int main() {
    Tester tester;
    cout << "Testing normal case for insert:" << endl;
    if (tester.insertNormal()) {
        cout << "\tNormal case for insert passed!" << endl;
    } else {
        cout << "\tNormal case for insert failed!" << endl;
    }
    cout << endl;
    cout << "Testing error case for getDNA (dna does not exist):" << endl;
    if (tester.getDNAerror()) {
        cout << "\tError case for getDNA passed!" << endl;
    } else {
        cout << "\tError case for getDNA failed!" << endl;
    }
    cout << endl;
    cout << "Testing edge case for getDNA:" << endl;
    if (tester.getDNAedge()) {
        cout << "\tEdge case for getDNA passed!" << endl;
    } else {
        cout << "\tEdge case for getDNA failed!" << endl;
    }
    cout << endl;
    cout << "Testing normal case for getDNA:" << endl;
    if (tester.getDNAnormal()) {
        cout << "\tNormal case for getDNA passed!" << endl;
    } else {
        cout << "\tNormal case for getDNA failed!" << endl;
    }
    cout << endl;
    cout << "Testing normal case for remove:" << endl;
    if (tester.removeNormal()) {
        cout << "\tNormal case for remove passed!" << endl;
    } else {
        cout << "\tNormal case for remove failed!" << endl;
    }
    cout << endl;
    cout << "Testing edge case for remove:" << endl;
    if (tester.removeEdge()) {
        cout << "\tEdge case for remove passed!" << endl;
    } else {
        cout << "\tEdge case for remove failed!" << endl;
    }
    cout << endl;
    cout << "Testing error case for remove (removing something that doesn't exist):" << endl;
    if (tester.removeError()) {
        cout << "\tError case for remove passed!" << endl;
    } else {
        cout << "\tError case for remove failed!" << endl;
    }
    cout << endl;
    cout << "Testing complete case for rehash (every data point gets rehashed and old table is emptied):" << endl;
    if (tester.rehashComplete()) {
        cout << "\tComplete case for rehash passed!" << endl;
    } else {
        cout << "\tComplete case for rehash failed!" << endl;
    }
    cout << endl;
    cout << "Testing rehash for deleted ratio:" << endl;
    if (tester.rehashDeletedRatio()) {
        cout << "\tNormal case for deletedRatio rehash passed!" << endl;
    } else {
        cout << "\tNormal case for deletedRatio rehash failed!" << endl;
    }
    cout << endl;
    cout << "Testing normal case for removal rehash:" << endl;
    if (tester.rehashRemove()) {
        cout << "\tNormal case for removal rehash passed!" << endl;
    } else {
        cout << "\tNormal case for removal rehash failed!" << endl;
    }
}


bool Tester::insertNormal() {
    DnaDb dnadb(MINPRIME, hashCode);
    DNA dataObj = DNA("AAAAA", 1001);  // insert a very small data set
    if (!dnadb.insert(dataObj)) {
        return false;
    }
    dataObj = DNA("ATTTT", 1002);
    if (!dnadb.insert(dataObj)) {
        return false;
    }
    dataObj = DNA("AGGGG", 1003);
    if (!dnadb.insert(dataObj)) {
        return false;
    }
    dataObj = DNA("ACCCC", 1004);
    if (!dnadb.insert(dataObj)) {
        return false;
    }
    dataObj = DNA("GGGGG", 1005);
    if (!dnadb.insert(dataObj)) {
        return false;
    }
    for (int i = 0; i < dnadb.m_currentCap; i ++) {
        if (dnadb.m_currentTable[i].m_sequence != "" && dnadb.m_currentTable[i].m_sequence != DELETEDKEY) {
            if (i != hashCode(dnadb.m_currentTable[i].m_sequence) % dnadb.m_currentCap) {  // makes sure every index abides by the proper index
                return false;
            }
        }
    }
    return true;
}

bool Tester::getDNAerror() {
    Random RndLocation(MINLOCID,MAXLOCID);
    vector<DNA> dataList;
    DnaDb dnadb(MINPRIME, hashCode);
    bool result = true;
    for (int i=0;i<49;i++){
        // generating random data
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());   // insert objects without needing to rehash
        // saving data for later use
        dataList.push_back(dataObj);
        // inserting data in to the DnaDb object
        dnadb.insert(dataObj);
    }
    DNA nonExist = DNA(sequencer(5, 3), 0);  // try to remove a nonexistent strand
    if (dnadb.getDNA(nonExist.getSequence(), nonExist.m_location).m_sequence != "") {   
        return false;
    }
    return true;
}

bool Tester::getDNAedge() {
    Random RndLocation(MINLOCID,MAXLOCID);
    vector<DNA> dataList;
    DnaDb dnadb(MINPRIME, hashCode);
    bool result = true;
    for (int i=0;i<5;i++){
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());  // add noncolliding data points
        dataList.push_back(dataObj);
        dnadb.insert(dataObj);
    }
    DNA dna = dataList.at(2);  // find strand with noncolliding points
    if (dnadb.getDNA(dna.getSequence(), dna.m_location).m_sequence != dna.m_sequence && dnadb.getDNA(dna.getSequence(), dna.m_location).m_location != dna.m_location) {   
        return false;
    }
    return true;
}

bool Tester::getDNAnormal() {
    Random RndLocation(MINLOCID,MAXLOCID);
    vector<DNA> dataList;
    DnaDb dnadb(MINPRIME, hashCode);
    bool result = true;
    for (int i=0;i<49;i++){
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());   //insert 49 so that rehashing is not needed
        dataList.push_back(dataObj);
        dnadb.insert(dataObj);
    }
    DNA dna = dataList.at(2);   // try to find a normal data point
    if (dnadb.getDNA(dna.getSequence(), dna.m_location).m_sequence != dna.m_sequence && dnadb.getDNA(dna.getSequence(), dna.m_location).m_location != dna.m_location) {   
        return false;
    }
    return true;
}

bool Tester::removeNormal() {
    DnaDb dnadb(MINPRIME, hashCode);
    DNA dataObj = DNA("AAAAA", 1001);
    dnadb.insert(dataObj);
    dataObj = DNA("ATTTT", 1002);
    dnadb.insert(dataObj);
    dataObj = DNA("AGGGG", 1003);
    dnadb.insert(dataObj);
    dataObj = DNA("ACCCC", 1004);
    dnadb.insert(dataObj);
    dataObj = DNA("GGGGG", 1005);
    dnadb.insert(dataObj);
    return dnadb.remove(dataObj);  // remove with noncolliding data points
}

bool Tester::removeEdge() {
    Random RndLocation(MINLOCID,MAXLOCID);
    vector<DNA> dataList;
    DnaDb dnadb(MINPRIME, hashCode);
    DNA dataObj;
    for (int i=0;i<49;i++){
        dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());  // insert 49 so that rehashing is not needed
        dataList.push_back(dataObj);
        dnadb.insert(dataObj);
    }
    return dnadb.remove(dataObj);  // remove  without calling rehash
}

bool Tester::removeError() {
    Random RndLocation(MINLOCID,MAXLOCID);
    vector<DNA> dataList;
    DnaDb dnadb(MINPRIME, hashCode);
    DNA dataObj;
    for (int i=0;i<48;i++){
        dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        dataList.push_back(dataObj);
        dnadb.insert(dataObj);
    }
    dataObj = DNA(sequencer(5, 48), RndLocation.getRandNum());  // try to remove a nonexistent object from dnadb
    if (dnadb.remove(dataObj)) {
        return false;
    }
    return true;
}

bool Tester::rehashComplete() {
    Random RndLocation(MINLOCID,MAXLOCID);
    DnaDb dnadb(MINPRIME, hashCode);
    for (int i=0;i<56;i++){
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        dnadb.insert(dataObj);
    }
    if (dnadb.m_currentSize == 56 && dnadb.m_oldTable == nullptr) {  // make sure all data points have been moved over and old table is cleared
        return true;
    }
    return false;
}

bool Tester::rehashRemove() {
    Random RndLocation(MINLOCID,MAXLOCID);
    vector<DNA> dataList;
    DnaDb dnadb(MINPRIME, hashCode);         
    bool result = true;
    for (int i=0;i<40;i++){
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        dataList.push_back(dataObj);
        dnadb.insert(dataObj);
    } 
    for (int i = 0; i < 33; i++) {   // remove 33 so that 0.8 is exceeded with deleted ratio
        dnadb.remove(dataList.at(i));
    }
    
    if (dnadb.m_currentCap != MINPRIME) {  // make sure current table is resized which would mean that rehash has been called
        return true;
    }
    return false;
}

bool Tester::rehashDeletedRatio() {
    Random RndLocation(MINLOCID,MAXLOCID);
    vector<DNA> dataList;
    DnaDb dnadb(MINPRIME, hashCode);
    bool result = true;
    for (int i=0;i<40;i++){
        DNA dataObj = DNA(sequencer(5, i), RndLocation.getRandNum());
        dataList.push_back(dataObj);
        dnadb.insert(dataObj);
    } 
    for (int i=0;i<38;i++){
        DNA dataObj = dataList.at(i);  // remove 38/40 so that tables should be completely rehashed
        dnadb.remove(dataObj);
    }
    if (dnadb.m_oldTable == nullptr) {  // makes sure it is completely rehashed and old table is cleared
        return true;
    }
    return false;
}

unsigned int hashCode(const string str) {
   unsigned int val = 0 ;
   const unsigned int thirtyThree = 33 ;  // magic number from textbook
   for ( int i = 0 ; i < str.length(); i++)
      val = val * thirtyThree + str[i] ;
   return val ;
}
string sequencer(int size, int seedNum){
    //this function returns a random DNA sequence
    string sequence = "";
    Random rndObject(0,3);
    rndObject.setSeed(seedNum);
    for (int i=0;i<size;i++){
        sequence = sequence + ALPHA[rndObject.getRandNum()];
    }
    return sequence;
}