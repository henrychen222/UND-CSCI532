/**
 *  To compile on the linux lab machines:
 *  g++ -std=c++11 unordered_map_test.cxx -o unordered_map_test
 *
 *  and run:
 *  ./unordered_map_tes
 */


#include<iostream>
using std::cout;
using std::endl;

#include<map>
using std::map;

#include<unordered_map>
using std::unordered_map;

#include <string>
using std::string;


int main(int argc, char **argv) {

    map<string, int> my_map;
    //unordered_map<string, int> my_map;

    cout << "my_map[\"ACGTA\"] is: " << my_map["ACGTA"] << endl;

    my_map["ACGTA"]++;

    cout << "my_map[\"ACGTA\"] is: " << my_map["ACGTA"] << endl;

    my_map["ACGTA"] = 100;

    cout << "my_map[\"ACGTA\"] is: " << my_map["ACGTA"] << endl;

    my_map["GGGTA"] = 25;
    my_map["CCCCCATAGA"] = 3;

    if (my_map.count("ACGTA") > 0) {
        cout << "my_map contained \"ACGTA\"" << endl;
    } else {
        cout << "my_map DID NOT CONTAIN \"ACGTA\"" << endl;
    }

    if (my_map.count("AGGGA") > 0) {
        cout << "my_map contained \"AGGGA\"" << endl;
    } else {
        cout << "my_map DID NOT CONTAIN \"AGGGA\"" << endl;
    }

    for (auto i = my_map.begin(); i != my_map.end(); i++) {
        cout << "i->first: " << i->first << ", i->second: " << i->second << endl;
    }
}